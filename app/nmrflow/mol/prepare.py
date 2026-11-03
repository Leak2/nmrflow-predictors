"""SMILES 预处理：解析、去盐、加氢、立体与芳香性感知。"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


@dataclass
class FragmentInfo:
    """一个组分的独立信息（盐的阳离子/阴离子各算一个）。"""

    smiles: str
    counts: dict[str, int]
    charge: int


@dataclass
class PreparedMol:
    """经过标准化处理的分子。"""

    smiles_in: str
    mol: Chem.Mol            # 已去盐、已加显式氢、含构象信息占位
    heavy: Chem.Mol          # 同一分子去掉显式氢（用于对称性与 SMARTS）
    smiles_out: str
    formula: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    fragments: list[FragmentInfo] = field(default_factory=list)

    @property
    def num_atoms(self) -> int:
        return self.mol.GetNumAtoms()


def _largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    """取最大的有机碎片（去盐 / 去溶剂 / 去对离子）。"""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol

    def score(m: Chem.Mol) -> tuple[int, int]:
        organic = sum(
            1 for a in m.GetAtoms()
            if a.GetAtomicNum() in (1, 6, 7, 8, 9, 15, 16, 17, 35, 53)
        )
        return (organic, m.GetNumAtoms())

    return max(frags, key=score)


def _organic_fragments(mol: Chem.Mol) -> Chem.Mol:
    """保留**所有**含碳碎片，只丢掉真正的无机对离子。

    NMR 测的是整份样品，光致产酸剂这类盐的阳离子（三苯基锍 15H）和
    阴离子（全氟链 12F）都会出峰，只取"最大碎片"会把一半信号弄丢。
    相反 Na+ / Cl- / BF4- / PF6- / 硫酸根这些不含碳，应当剔除。
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1:
        return mol
    kept = [
        f for f in frags
        if any(a.GetAtomicNum() == 6 for a in f.GetAtoms())
    ]
    if not kept:                # 全是无机物，退化成取最大碎片
        return _largest_fragment(mol)
    if len(kept) == 1:
        return kept[0]

    combined = kept[0]
    for frag in kept[1:]:
        combined = Chem.CombineMols(combined, frag)
    # CombineMols 不会重建环信息 / 芳香性缓存。
    # 不做这一步，合并后所有环查询都会失手 —— 表现为芳香氢全部退回
    # 苯的基准值 7.26、13C 全部 128.5（看起来像"没有取代基"）。
    rw = Chem.RWMol(combined)
    Chem.SanitizeMol(rw)
    return rw.GetMol()


def prepare(smiles: str, *, add_hs: bool = True,
            all_fragments: bool = True) -> PreparedMol:
    """解析并标准化一个 SMILES。

    处理：多组分去盐、芳香性感知、立体信息保留、加显式氢。
    """
    warnings: list[str] = []
    raw = Chem.MolFromSmiles(smiles)
    if raw is None:
        raise ValueError(f"无法解析 SMILES：{smiles!r}")

    try:
        Chem.SanitizeMol(raw)
    except Exception as exc:  # pragma: no cover - 极端输入
        raise ValueError(f"SMILES 结构不合理（{exc}）：{smiles!r}") from exc

    n_frags = len(Chem.GetMolFrags(raw))
    if all_fragments:
        mol = _organic_fragments(raw)
    else:
        mol = _largest_fragment(raw)
    if n_frags > 1:
        n_org = sum(
            1 for f in Chem.GetMolFrags(raw, asMols=True, sanitizeFrags=True)
            if any(a.GetAtomicNum() == 6 for a in f.GetAtoms())
        )
        if all_fragments and n_org > 1:
            warnings.append(
                f"输入含 {n_frags} 个组分，其中 {n_org} 个含碳，"
                f"已合并全部有机碎片（无机对离子已剔除）"
            )
        else:
            warnings.append(
                f"输入含 {n_frags} 个组分，已取最大有机碎片（去盐/去对离子）"
            )

    n_rad = sum(1 for a in mol.GetAtoms() if a.GetNumRadicalElectrons() > 0)
    if n_rad:
        warnings.append(f"分子含 {n_rad} 个自由基电子，谱图预测可能不可靠")

    charge = Chem.GetFormalCharge(mol)
    if charge != 0:
        warnings.append(f"母体净电荷 {charge:+d}，位移会显著偏离中性参考值")

    heavy = Chem.Mol(mol)
    with_h = Chem.Mol(mol)
    if add_hs:
        with_h = Chem.AddHs(with_h)

    formula = {}
    for atom in with_h.GetAtoms():
        formula[atom.GetSymbol()] = formula.get(atom.GetSymbol(), 0) + 1

    # 各组分单独的元素计数与电荷：质谱对盐必须分别报阳离子/阴离子，
    # 把 C18H15S+ 和 C14H7F12O9S4- 相加成 C32H22F12O9S5 再算 [M+H]+ 没有化学意义。
    fragments: list[FragmentInfo] = []
    if n_frags > 1:
        for frag in Chem.GetMolFrags(raw, asMols=True, sanitizeFrags=True):
            frag_h = Chem.AddHs(frag)
            counts: dict[str, int] = {}
            for atom in frag_h.GetAtoms():
                counts[atom.GetSymbol()] = counts.get(atom.GetSymbol(), 0) + 1
            fragments.append(
                FragmentInfo(
                    smiles=Chem.MolToSmiles(frag),
                    counts=counts,
                    charge=Chem.GetFormalCharge(frag),
                )
            )

    return PreparedMol(
        smiles_in=smiles,
        mol=with_h,
        heavy=heavy,
        smiles_out=Chem.MolToSmiles(mol),
        formula=formula,
        warnings=warnings,
        fragments=fragments,
    )


def ring_sizes(mol: Chem.Mol) -> dict[int, list[int]]:
    """原子索引 -> 其所在环的大小列表。"""
    info = mol.GetRingInfo()
    out: dict[int, list[int]] = {}
    for ring in info.AtomRings():
        for idx in ring:
            out.setdefault(idx, []).append(len(ring))
    return out


def mol_weight(mol: Chem.Mol) -> float:
    return rdMolDescriptors.CalcExactMolWt(mol)
