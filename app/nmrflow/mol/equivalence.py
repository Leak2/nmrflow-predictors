"""等价性判定：化学等价、非对映异位、磁等价。"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

from .prepare import PreparedMol, ring_sizes


@dataclass
class SpinGroup:
    """一组化学等价的核。"""

    group_id: int
    nucleus: str
    atom_indices: list[int]        # 核本身的原子索引（H 为氢原子；13C 为碳）
    heavy_indices: list[int]       # 对应的重原子索引（H 时为其所连重原子）
    n_nuclei: int
    diastereotopic_pair: int = -1  # >=0 表示与另一个组是同一碳上的非对映异位对


def symmetry_classes(mol: Chem.Mol) -> list[int]:
    """拓扑对称类（同类的原子化学等价）。includeChirality=True 以区分手性环境。"""
    return list(Chem.CanonicalRankAtoms(mol, breakTies=False, includeChirality=True))


def has_stereochemistry(mol: Chem.Mol) -> bool:
    """是否存在已指定的立体化学（手性中心或 E/Z 双键）。"""
    try:
        rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:
        pass
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            return True
    for bond in mol.GetBonds():
        if bond.GetStereo() not in (Chem.BondStereo.STEREONONE,):
            return True
    return False


def potential_stereocenters(mol: Chem.Mol) -> set[int]:
    """可能的手性中心，**包含未指定构型的**。

    ``has_stereochemistry`` 只认已经画出来的 ``@``；但 CF₂Cl–CFClI 这类分子，
    样品里那个碳就是手性中心，只是 SMILES 里没写 @。判断"CF₂ 上两个 F 是不是
    非对映异位"必须用这一版，否则未标注的谱子会被当成两个 F 等价。
    """
    centers: set[int] = set()
    try:
        for info in Chem.FindPotentialStereo(mol, cleanIt=False, flagPossible=True):
            if info.type == Chem.StereoType.Atom_Tetrahedral:
                centers.add(info.centeredOn)
    except Exception:  # 老版本 RDKit 没有该 API：当作无手性中心，退回旧行为
        return centers
    return centers


def _hydrogens_of(mol: Chem.Mol, heavy_idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(heavy_idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]


def _is_rigid_ring(idx: int, rings: dict[int, list[int]]) -> bool:
    return any(size <= 4 for size in rings.get(idx, []))


def _ch2_is_diastereotopic(
    mol: Chem.Mol, heavy_idx: int, ranks: list[int], chiral: bool,
    rings: dict[int, list[int]],
) -> bool:
    """判断 sp3-CH2 的两个 H 是否非对映异位（不等价）。"""
    atom = mol.GetAtomWithIdx(heavy_idx)

    # sp2 上的 =CH2：两个 H 分别为 cis / trans，恒不等价
    if atom.GetHybridization() == Chem.HybridizationType.SP2 and atom.GetIsAromatic() is False:
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                return True

    neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
    if len(neighbors) < 2:
        return False
    if ranks[neighbors[0]] == ranks[neighbors[1]]:
        return False           # 两个取代基相同 -> 同topic，等价
    # 取代基不同：只有存在立体元素或处于刚性小环时才表现为非对映异位
    return chiral or _is_rigid_ring(heavy_idx, rings)


def hydrogen_groups(prep: PreparedMol, start: int = 0) -> list[SpinGroup]:
    """把分子中所有氢划分成化学等价组。"""
    mol = prep.mol                    # 含显式氢
    heavy = prep.heavy                # 不含显式氢
    ranks = symmetry_classes(heavy)
    chiral = has_stereochemistry(heavy)
    rings = ring_sizes(heavy)

    buckets: dict[tuple[int, int], list[int]] = {}
    order: list[tuple[int, int]] = []

    for heavy_idx in range(heavy.GetNumAtoms()):
        hs = _hydrogens_of(mol, heavy_idx)
        if not hs:
            continue
        cls = ranks[heavy_idx]
        if len(hs) == 1:
            keys = [(cls, 0)] * 1
            for h, k in zip(hs, keys):
                buckets.setdefault(k, []).append(h)
                if k not in order:
                    order.append(k)
            continue
        split = _ch2_is_diastereotopic(mol, heavy_idx, ranks, chiral, rings) if len(hs) == 2 else False
        if split:
            for sub, h in enumerate(hs):
                k = (cls, sub)
                buckets.setdefault(k, []).append(h)
                if k not in order:
                    order.append(k)
        else:
            k = (cls, 0)
            buckets.setdefault(k, []).extend(hs)
            if k not in order:
                order.append(k)

    groups: list[SpinGroup] = []
    for offset, key in enumerate(order):
        gid = start + offset
        hs = buckets[key]
        heavy_of = [next(n.GetIdx() for n in mol.GetAtomWithIdx(h).GetNeighbors()) for h in hs]
        groups.append(
            SpinGroup(gid, "1H", hs, heavy_of, len(hs))
        )
    return groups


def carbon_groups(prep: PreparedMol, start: int = 0) -> list[SpinGroup]:
    """按对称类划分碳原子（13C 用）。"""
    heavy = prep.heavy
    ranks = symmetry_classes(heavy)
    buckets: dict[int, list[int]] = {}
    for idx in range(heavy.GetNumAtoms()):
        if heavy.GetAtomWithIdx(idx).GetAtomicNum() == 6:
            buckets.setdefault(ranks[idx], []).append(idx)
    return [
        SpinGroup(start + gid, "13C", idxs, idxs, len(idxs))
        for gid, idxs in enumerate(buckets.values())
    ]


def hetero_groups(prep: PreparedMol, element: str,
                  start: int = 0) -> list[SpinGroup]:
    """19F / 31P 等异核：每个原子一组（再用对称类合并）。"""
    heavy = prep.heavy
    ranks = symmetry_classes(heavy)
    buckets: dict[int, list[int]] = {}
    atomic = {"F": 9, "P": 15}[element]
    for idx in range(heavy.GetNumAtoms()):
        if heavy.GetAtomWithIdx(idx).GetAtomicNum() == atomic:
            buckets.setdefault(ranks[idx], []).append(idx)
    return [
        SpinGroup(start + gid, {"F": "19F", "P": "31P"}[element],
                  idxs, idxs, len(idxs))
        for gid, idxs in enumerate(buckets.values())
    ]
