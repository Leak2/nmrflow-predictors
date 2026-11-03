"""异核耦合：C-F / C-P / H-F / H-P / F-F / P-F / P-P。

去耦策略（nucleus-level，可配置）：
  13C : 对 1H 宽带去耦 -> 单峰；但保留与 F / P 的耦合
  31P : 对 1H 去耦；保留 P-F / P-P
  19F : 不去耦 -> 保留 H-F 裂分
  1H  : 不去耦 -> 必须看到 H-F / H-P 裂分
"""

from __future__ import annotations

from rdkit import Chem

from ..core.constants import DECOUPLE_DEFAULT
from ..core.types import Coupling, Signal
from ..data.loader import load_toml
from ..mol.envs import bond_paths
from ..mol.prepare import PreparedMol

# 参与耦合的核：观测核 -> [(伙伴核, 耦合表键)]
PAIRS = {
    "13C": [("19F", "C_F"), ("31P", "C_P")],
    "1H": [("19F", "H_F"), ("31P", "H_P")],
    "19F": [("1H", "H_F"), ("19F", "F_F"), ("31P", "P_F")],
    "31P": [("19F", "P_F"), ("31P", "P_P")],
}

# 距离修正：耦合表按「两个核之间有几根键」索引，而我们量的是**重原子**距离。
# 只要 pair 里有 1H，就少算了一根 H-C 键 -> +1。
#   例：P-O-CH3 的重原子距离是 C..P = 2，但 H-C-O-P 是 3 根键 -> 3J = 11 Hz
BOND_OFFSET = 1

# 小于这个值的耦合在常规线宽下不可分辨，报道里也不写
MIN_REPORTABLE_J = 0.4

# 13C 的天然线宽 1~2 Hz，低于 1 Hz 的耦合根本分辨不出来。
# 不按核区分的话，CF3 碳会额外带上 0.6 Hz 的远端 5J(C-F)，
# 多重峰被算成"qq -> 峰太多 -> m"，反而把关键的 1J(C-F)=320 埋掉。
MIN_REPORTABLE_J_BY_NUCLEUS = {"13C": 1.0}

# 芳环 -> sp3 侧链"绕道"耦合的衰减系数（缺省 1.0 = 不衰减；
# 只有 H_F / H_P 这类"绕出环就消失"的表才在 j_hetero.toml 里改小）
DEFAULT_RING_EXIT_DAMPING = 1.0


class HeteroCouplingPredictor:
    def __init__(self) -> None:
        self.t = load_toml("j_hetero.toml")

    def apply(
        self,
        prep: PreparedMol,
        signals_by_nucleus: dict[str, list[Signal]],
        *,
        decouple: dict[str, tuple[str, ...]] | None = None,
    ) -> None:
        dec = dict(DECOUPLE_DEFAULT)
        if decouple:
            dec.update(decouple)

        for nucleus, signals in signals_by_nucleus.items():
            for signal in signals:
                signal.couplings = [
                    c for c in signal.couplings
                    if c.partner_nucleus not in dec.get(nucleus, ())
                ]
                for partner_nuc, table_key in PAIRS.get(nucleus, []):
                    if partner_nuc in dec.get(nucleus, ()):
                        continue
                    partners = signals_by_nucleus.get(partner_nuc, [])
                    if not partners:
                        continue
                    self._add(prep, signal, partners, partner_nuc, table_key)

    # ------------------------------------------------------------------
    def _add(self, prep: PreparedMol, signal: Signal, partners: list[Signal],
             partner_nuc: str, table_key: str) -> None:
        heavy = prep.heavy
        table = self.t.get(table_key, {})
        overrides = table.get("override", [])
        offset = BOND_OFFSET if "1H" in (signal.nucleus, partner_nuc) else 0

        # 观测核的重原子（1H 信号返回的是它连的那个碳）
        obs_atoms = self._heavy_atoms(prep, signal)
        if not obs_atoms:
            return

        for partner in partners:
            if partner.group_id == signal.group_id:
                continue          # 组内磁等价核之间无可观测裂分
            p_atoms = self._heavy_atoms(prep, partner)
            if not p_atoms:
                continue

            # 与 H-H 同理：等价核个数按「真正耦合到的那个原子」来数。
            # 1,4-二氟苯的 F 只跟邻位 2 个 H 耦合（三重峰），
            # 而不是跟芳环上全部 4 个 H 耦合（五重峰）。
            # 观测核取**一个代表原子**：组内是等价核，看到的分裂完全相同。
            # 不能对 obs_atoms 取最小距离 —— 1,4-二氟苯里每个 H 都邻位于
            # 两个 F 中的一个，取最小会把全部 4 个 H 都算成 3J。
            rep_obs = obs_atoms[0]
            buckets: dict[float, int] = {}
            kinds: dict[float, str] = {}
            for b in p_atoms:
                best_bonds, best_path = self._bond_distance(heavy, rep_obs, b)
                if best_bonds is None or best_bonds < 1:
                    continue
                bonds = best_bonds + offset
                ring_path = self._aromatic_ring_path(heavy, best_path)
                j = self._lookup(
                    table, overrides, heavy, [b], bonds, aromatic=ring_path,
                )
                if j is None:
                    continue
                # 芳香环 -> sp3 侧链的"绕道"耦合（Ar-CF3 的 H-F）几乎为 0
                if self._exits_aromatic_ring(heavy, best_path):
                    j *= float(
                        self.t.get("aromatic_exit_damping", {}).get(
                            table_key, DEFAULT_RING_EXIT_DAMPING
                        )
                    )
                if abs(j) < MIN_REPORTABLE_J_BY_NUCLEUS.get(
                    signal.nucleus, MIN_REPORTABLE_J
                ):
                    continue
                key = round(j, 1)
                buckets[key] = buckets.get(key, 0) + self._nuclei_on(
                    prep, partner, b
                )
                kinds.setdefault(key, f"{bonds}J({table_key.replace('_', '-')})")

            for key, n in sorted(buckets.items(), key=lambda kv: -abs(kv[0])):
                signal.couplings.append(
                    Coupling(
                        partner_group=partner.group_id,
                        j_hz=float(key),
                        n_equiv=n,
                        partner_nucleus=partner_nuc,
                        kind=kinds[key],
                    )
                )

    # ------------------------------------------------------------------
    @staticmethod
    def _aromatic_run(mol: Chem.Mol, path: tuple[int, ...]
                      ) -> tuple[int, int]:
        """路径中最长的一段连续芳香原子，返回 (起, 止)。"""
        atoms = [mol.GetAtomWithIdx(i) for i in path]
        best = (0, 0)
        i = 0
        while i < len(atoms):
            if atoms[i].GetIsAromatic():
                j = i
                while j < len(atoms) and atoms[j].GetIsAromatic():
                    j += 1
                if j - i > best[1] - best[0]:
                    best = (i, j)
                i = j
            else:
                i += 1
        return best

    @classmethod
    def _aromatic_ring_path(cls, mol: Chem.Mol, path: tuple[int, ...]) -> bool:
        """耦合路径是否穿过芳环（>=2 个连续芳香原子）。"""
        lo, hi = cls._aromatic_run(mol, path)
        return hi - lo >= 2

    @classmethod
    def _exits_aromatic_ring(cls, mol: Chem.Mol, path: tuple[int, ...]) -> bool:
        """最短耦合路径是否「芳环 -> 非芳香 sp3 碳」。

        典型例子：Ph-CF3 中 H(邻位)-C(ar)-C(ar)-C(sp3)-F。
        这种绕出环的 H-F 耦合实测 < 0.5 Hz（PhCF3 在 -63 ppm 是单峰），
        而同样走 3 根键的 Ar-F 邻位 H-F 却有 ~9 Hz，必须区分开。
        """
        if len(path) < 3:
            return False
        lo, hi = cls._aromatic_run(mol, path)
        if hi - lo < 2:
            return False
        atoms = [mol.GetAtomWithIdx(i) for i in path]
        for k in (lo - 1, hi):
            if 0 <= k < len(atoms):
                a = atoms[k]
                if a.GetSymbol() == "C" and not a.GetIsAromatic():
                    return True
        return False

    # ------------------------------------------------------------------
    @staticmethod
    def _lookup(table: dict, overrides: list[dict], heavy: Chem.Mol,
                p_atoms: list[int], bonds: int,
                aromatic: bool = False) -> float | None:
        # 走芳香环的耦合有自己的一套 J（例如 Ar-F 的邻/间/对位 H-F）
        if aromatic and isinstance(table.get("aromatic"), dict):
            key = str(bonds)
            if key in table["aromatic"]:
                return float(table["aromatic"][key])
        for ov in overrides:
            if int(ov.get("bonds", 0)) != bonds:
                continue
            patt = Chem.MolFromSmarts(ov["smarts"])
            if patt is None:
                continue
            # 注意：RDKit 对对称结构只返回一个代表匹配（CF3 的 3 个 F
            # 只会出现一次 match[0]），所以要看「伙伴原子是否出现在匹配里」
            for match in heavy.GetSubstructMatches(patt):
                if set(match) & set(p_atoms):
                    return float(ov["value"])
        key = str(bonds)
        if key in table:
            return float(table[key])
        return None

    @staticmethod
    def _nuclei_on(prep: PreparedMol, signal: Signal, heavy_idx: int) -> int:
        """某个重原子上挂着几个被观测的核。"""
        if signal.nucleus == "1H":
            return sum(
                1 for n in prep.mol.GetAtomWithIdx(heavy_idx).GetNeighbors()
                if n.GetAtomicNum() == 1
            )
        return 1

    @staticmethod
    def _heavy_atoms(prep: PreparedMol, signal: Signal) -> list[int]:
        out = []
        for idx in signal.atom_indices:
            atom = prep.mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() > 1:
                out.append(idx)
                continue
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() > 1:
                    out.append(nb.GetIdx())
                    break
        return sorted(set(out))

    @staticmethod
    def _bond_distance(mol: Chem.Mol, a: int, b: int, max_len: int = 5
                       ) -> tuple[int | None, tuple[int, ...]]:
        if a == b:
            return 0, (a,)
        paths = bond_paths(mol, a, b, max_len=max_len)
        if not paths:
            return None, ()
        best = min(paths, key=len)
        return len(best) - 1, tuple(best)
