"""1H-1H 耦合常数预测。"""

from __future__ import annotations

from rdkit import Chem

from ..core.types import Coupling, Signal
from ..data.loader import load_toml
from ..mol.envs import aromatic_ring_of, bond_paths, ring_distance
from ..mol.prepare import PreparedMol, ring_sizes
from .base import CouplingPredictor

# 小于这个值的 H-H 耦合在常规线宽下不可分辨，不参与裂分
MIN_HH_J = 1.0


def _heavy_set(sig: Signal) -> list[int]:
    """组内去重后的重原子（heavy_indices 是按核逐个列的，会有重复）。"""
    return sorted(set(sig.heavy_indices))


def _is_double_to_carbon(atom: Chem.Atom) -> bool:
    for bond in atom.GetBonds():
        if (
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(atom).GetSymbol() == "C"
        ):
            return True
    return False


def _is_carbonyl(atom: Chem.Atom) -> bool:
    for bond in atom.GetBonds():
        if (
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(atom).GetSymbol() == "O"
        ):
            return True
    return False


class HHCouplingPredictor(CouplingPredictor):
    """按耦合路径长度与几何环境给 J。"""

    def __init__(self) -> None:
        self.t = load_toml("j_hh.toml")

    # ------------------------------------------------------------------
    def predict(self, prep: PreparedMol, signals: list[Signal]) -> None:
        """按「重原子对」枚举，再按 J 值归并等价核个数。

        关键点：等价核个数必须**按耦合到的那个原子**来数，不能直接用
        partner 组的总核数。正辛烷的 -CH3 组含 6 个 H（C1 + C8），
        但其中任一个甲基只跟自己邻位的 2 个 H 耦合 —— 是三重峰，
        不是七重峰。
        """
        mol = prep.heavy
        rings = ring_sizes(mol)
        # H 原子 -> 它在所属重原子上的位置序号（0/1/2…）。
        # 用于在 2D 拓扑下区分双键上的 cis / trans：
        # 两端 H 的序号同奇偶 -> 同侧（cis），否则异侧（trans）。
        self._h_slot = self._h_slot_map(prep)
        # 组 -> {重原子: [该组落在它上面的 H 原子索引]}
        # 等价核个数必须由这张表来数：=CH2 被拆成两个非对映异位组时，
        # 每个组在同一个碳上只贡献 1 个 H，不能拿"碳上总氢数"重复计数。
        h_map: dict[int, dict[int, list[int]]] = {}
        for s in signals:
            per: dict[int, list[int]] = {}
            for h_idx, heavy_idx in zip(s.atom_indices, s.heavy_indices):
                per.setdefault(heavy_idx, []).append(h_idx)
            h_map[s.group_id] = per

        for si in signals:
            si.couplings = []
            if si.is_labile or not si.heavy_indices:
                continue
            # 取组内一个代表重原子：等价核看到的分裂完全相同
            rep = si.heavy_indices[0]
            rep_h = si.atom_indices[0] if si.atom_indices else None
            for sj in signals:
                if sj.group_id == si.group_id or sj.is_labile:
                    continue
                sjb = _heavy_set(sj)
                if not sjb:
                    continue
                per_j = h_map.get(sj.group_id, {})
                # 同一碳上的非对映异位对 -> 2J 偕偶
                if rep in sjb:
                    gem = self._geminal(mol, rep, rings)
                    if abs(gem) > MIN_HH_J:
                        si.couplings.append(
                            Coupling(
                                partner_group=sj.group_id, j_hz=gem,
                                n_equiv=len(per_j.get(rep, [])) or sj.n_nuclei,
                                partner_nucleus="1H",
                                kind=self._last_kind,
                            )
                        )
                    continue
                buckets: dict[float, int] = {}
                kinds: dict[float, str] = {}
                for b in sjb:
                    hb = (per_j.get(b) or [None])[0]
                    j_hz = self._pair_j(mol, rep, b, rings, ha=rep_h, hb=hb)
                    if j_hz is None or abs(j_hz) < MIN_HH_J:
                        continue
                    key = round(j_hz, 1)
                    buckets[key] = buckets.get(key, 0) + (len(per_j.get(b, [])) or 1)
                    kinds.setdefault(key, self._last_kind)
                for key, n in sorted(buckets.items(), key=lambda kv: -abs(kv[0])):
                    si.couplings.append(
                        Coupling(
                            partner_group=sj.group_id, j_hz=float(key),
                            n_equiv=n, partner_nucleus="1H",
                            kind=kinds[key],
                        )
                    )

    # ------------------------------------------------------------------
    @staticmethod
    def _h_slot_map(prep: PreparedMol) -> dict[int, int]:
        """H 原子索引 -> 它在所连重原子上的位置序号。

        索引来自 prep.mol（带显式氢），与 Signal.atom_indices 一致。
        """
        by_heavy: dict[int, list[int]] = {}
        for atom in prep.mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() > 1:
                    by_heavy.setdefault(nb.GetIdx(), []).append(atom.GetIdx())
        out: dict[int, int] = {}
        for hs in by_heavy.values():
            for slot, h in enumerate(sorted(hs)):
                out[h] = slot
        return out

    @staticmethod
    def _h_per_atom(prep: PreparedMol) -> dict[int, int]:
        """重原子 -> 它带几个氢（含显式氢的分子）。"""
        out: dict[int, int] = {}
        for atom in prep.mol.GetAtoms():
            if atom.GetAtomicNum() > 1:
                out[atom.GetIdx()] = sum(
                    1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1
                )
        return out

    # ------------------------------------------------------------------
    def _pair_j(self, mol: Chem.Mol, a: int, b: int, rings,
                ha: int | None = None, hb: int | None = None) -> float | None:
        self._last_kind = ""
        if a < 0 or b < 0:
            return None
        if a == b:
            return self._geminal(mol, a, rings)

        atom_a = mol.GetAtomWithIdx(a)
        atom_b = mol.GetAtomWithIdx(b)

        # --- 芳环内耦合：按环内距离 ---
        if atom_a.GetIsAromatic() and atom_b.GetIsAromatic():
            ring = aromatic_ring_of(mol, a)
            if ring and b in ring:
                d = ring_distance(ring, a, b)
                return self._aromatic_j(mol, ring, a, b, d)

        # --- 双键上的 cis / trans ---
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
            return self._alkene_j(mol, a, b, ha, hb)

        # --- 路径长度（重原子之间）---
        # H-C-C-H 的两个碳是**直接相连**的，所以 3J 对应重原子距离 1：
        #   距离 1 -> 3J 邻位；距离 2 -> 4J 长程；距离 3 -> 5J（通常忽略）
        d = self._bond_distance(mol, a, b)
        if d == 1:
            return self._vicinal(mol, a, b, rings)
        if d == 2:
            return self._long_range(mol, a, b)
        return None

    # ------------------------------------------------------------------
    def _geminal(self, mol: Chem.Mol, c: int, rings) -> float:
        sec = self.t["geminal_2J"]
        base = float(sec["base"])
        mod = sec["modifier"]
        atom = mol.GetAtomWithIdx(c)

        # 烯烃末端 =CH2：2J 是小的正值（乙烯 2.5 / 苯乙烯 1.0）。
        # 必须在这里拦下 —— 否则会掉进下面 sp2 的 -15 分支，
        # 在 400 MHz 上裂出 375 Hz 的假双峰。
        if (atom.GetHybridization() == Chem.HybridizationType.SP2
                and _is_double_to_carbon(atom)):
            self._last_kind = "2J(alkene)"
            return float(sec.get("alkene_terminal", 2.0))

        heavy_nb = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
        value = base
        for n in heavy_nb:
            sym = n.GetSymbol()
            if _is_carbonyl(n):
                value = min(value, float(mod["carbonyl"]))
            elif n.GetIsAromatic() or _is_double_to_carbon(n):
                value = min(value, float(mod["sp2_or_aryl"]))
            elif sym == "O":
                value = min(value, float(mod["oxygen"]))
            elif sym == "N":
                value = min(value, float(mod["nitrogen"]))
            elif sym == "S":
                value = min(value, float(mod["sulfur"]))
            elif sym in ("F", "Cl", "Br", "I"):
                value = min(value, float(mod["halogen"]))
        for size in rings.get(c, []):
            key = str(size)
            if key in sec.get("ring", {}):
                value = float(sec["ring"][key])
                break
        self._last_kind = "2J"
        return value

    # ------------------------------------------------------------------
    def _vicinal(self, mol: Chem.Mol, a: int, b: int, rings) -> float:
        sec = self.t["vicinal_3J"]
        self._last_kind = "3J"
        # 醛基 / 甲酰基上的 3J 明显偏小（CH3-CHO 只有 ~2.9 Hz）
        if _is_carbonyl(mol.GetAtomWithIdx(a)) or _is_carbonyl(mol.GetAtomWithIdx(b)):
            return float(sec.get("aldehyde", 3.0))
        # 环上的 3J：粗略按环大小区分直立/平伏不可知，取中间值
        if rings.get(a) and set(rings.get(a, [])) & set(rings.get(b, [])):
            common = set(rings.get(a, [])) & set(rings.get(b, []))
            if common:
                return (float(sec["cyclic_aa"]) + float(sec["cyclic_ae"])) / 2.0
        return float(sec["free_rotation"])

    # ------------------------------------------------------------------
    def _alkene_j(self, mol: Chem.Mol, a: int, b: int,
                  ha: int | None = None, hb: int | None = None) -> float:
        """乙烯型 3J：优先用双键立体标记，其次用两端 H 的位置序号判顺反。

        =CH2 末端碳上没有立体标记，但它两个 H 必然一个 cis 一个 trans，
        必须分别给 10.5 / 17.0 —— 都取平均值会把 dd 算成两个等值裂分。
        """
        sec = self.t["vicinal_3J"]
        cis, trans = float(sec["alkene_cis"]), float(sec["alkene_trans"])
        self._last_kind = "3J(alkene)"
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond is None:
            return (cis + trans) / 2.0
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOCIS:
            return cis
        if stereo == Chem.BondStereo.STEREOTRANS:
            return trans
        # 无立体标记：用两端 H 在各自碳上的序号奇偶判同侧/异侧
        slots = getattr(self, "_h_slot", {})
        if ha is not None and hb is not None and ha in slots and hb in slots:
            return cis if (slots[ha] - slots[hb]) % 2 == 0 else trans
        # 还判不了（例如环内双键两侧都不带 H）
        return (cis + trans) / 2.0

    # ------------------------------------------------------------------
    def _long_range(self, mol: Chem.Mol, a: int, b: int) -> float | None:
        sec = self.t["long_range_4J"]
        # 烯丙型：路径中间含 C=C
        for path in bond_paths(mol, a, b, max_len=3):
            if len(path) != 4:
                continue
            if any(
                _is_double_to_carbon(mol.GetAtomWithIdx(i)) for i in path[1:3]
            ):
                self._last_kind = "4J(allylic)"
                return float(sec["allylic"])
        self._last_kind = ""
        return None

    # ------------------------------------------------------------------
    def _aromatic_j(self, mol: Chem.Mol, ring: list[int], a: int, b: int,
                    d: int) -> float | None:
        sec = self.t["aromatic_J"]
        if len(ring) == 5:
            name = self._ring5_name(mol, ring)
            tbl = self.t["ring5_J"].get(name, self.t["ring5_J"]["default"])
            if d == 1:
                self._last_kind = "3J(ring5)"
                return float(tbl["adjacent"])
            if d == 2:
                self._last_kind = "4J(ring5)"
                return float(tbl["skip"])
            return None
        key = str(d)
        if key not in sec:
            return None
        value = float(sec[key])
        if d == 1:
            self._last_kind = "3J(ortho)"
        elif d == 2:
            self._last_kind = "4J(meta)"
        else:
            self._last_kind = "5J(para)"
        return value if abs(value) > 0.05 else None

    @staticmethod
    def _ring5_name(mol: Chem.Mol, ring: list[int]) -> str:
        symbols = sorted(mol.GetAtomWithIdx(i).GetSymbol() for i in ring)
        if "O" in symbols:
            return "furan"
        if "S" in symbols:
            return "thiophene"
        if symbols.count("N") >= 2:
            return "imidazole"
        if "N" in symbols:
            return "pyrrole"
        return "default"

    # ------------------------------------------------------------------
    @staticmethod
    def _bond_distance(mol: Chem.Mol, a: int, b: int, max_len: int = 4) -> int:
        paths = bond_paths(mol, a, b, max_len=max_len)
        if not paths:
            return -1
        return min(len(p) - 1 for p in paths)
