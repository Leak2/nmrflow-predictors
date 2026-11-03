"""13C 化学位移预测：Grant-Paul 烷烃加和 + 取代基修正 + 芳环增量。"""

from __future__ import annotations

from rdkit import Chem

from ..core.solvent import get_solvent
from ..core.types import Confidence, Provenance, Signal
from ..data.loader import shift_table
from ..mol.envs import (
    SubstituentSet,
    aromatic_ring_of,
    relative_position,
    ring_distance,
    ring_substituents,
    scan_alpha_beta_gamma,
)
from ..mol.equivalence import SpinGroup
from ..mol.prepare import PreparedMol, ring_sizes
from .base import ShiftPredictor


def _is_carbonyl(atom: Chem.Atom) -> bool:
    """这个原子是不是羰基碳（C=O）。"""
    return atom.GetSymbol() == "C" and any(
        b.GetBondType() == Chem.BondType.DOUBLE
        and b.GetOtherAtom(atom).GetSymbol() == "O"
        for b in atom.GetBonds()
    )

def _double_bond_partner(mol: Chem.Mol, idx: int) -> int | None:
    """与 idx 成双键的碳原子；没有（或对侧不是碳）返回 None。"""
    atom = mol.GetAtomWithIdx(idx)
    for bond in atom.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 6:
            return other.GetIdx()
    return None


def _chain_carbon_neighbors(mol: Chem.Mol, center: int, exclude: int) -> int:
    """α 碳上"链继续往外延伸"的碳数，用于 ``beta_step``。

    只数**非环**碳邻居。环己烯的烯丙位 CH₂ 另有两只环内碳邻居，但那是环的一部分、
    不是链延伸；实测环己烯烯碳 127.1 ppm，把环内碳算进去会推到 143 ppm。
    """
    n = 0
    for bond in mol.GetAtomWithIdx(center).GetBonds():
        if bond.IsInRing():
            continue
        j = bond.GetOtherAtomIdx(center)
        if j == exclude:
            continue
        if mol.GetAtomWithIdx(j).GetAtomicNum() == 6:
            n += 1
    return n


def _describe(items: list[tuple[str, float]]) -> str:
    """把一侧的增量列表写成 ``(alkyl+10.1,vinyl+7.8)`` 这样的短串，便于追溯。"""
    if not items:
        return "(none)"
    return "(" + ",".join(f"{name}{value:+.1f}" for name, value in items) + ")"


LEVEL_KEY = {1: "a", 2: "b", 3: "g"}

class CarbonShiftPredictor(ShiftPredictor):
    nucleus = "13C"
    tier = "A"

    def __init__(self) -> None:
        self.table = shift_table("13C")
        self.sp3_subs = SubstituentSet(self.table["sp3"]["subst"])
        self.arom_subs = SubstituentSet(self.table["aromatic"]["subst"])
        self.alkene_own = SubstituentSet(self.table["alkene"]["own"])
        self.alkene_far = SubstituentSet(self.table["alkene"]["far"])
        self.alkene_pi = frozenset(self.table["alkene"].get("pi_classes", []))
        # 只保留带 smarts 的条目（functional 里还混着 n_methyl_aromatic 这类子表）
        self.functional = [
            (k, v) for k, v in self.table["functional"].items()
            if isinstance(v, dict) and "smarts" in v and "value" in v
        ]

    # ------------------------------------------------------------------
    def predict(
        self, prep: PreparedMol, groups: list[SpinGroup], solvent_key: str
    ) -> list[Signal]:
        solvent = get_solvent(solvent_key)
        heavy = prep.heavy
        rings = ring_sizes(heavy)
        out: list[Signal] = []
        for group in groups:
            idx = group.atom_indices[0]
            sig = Signal(
                nucleus="13C", group_id=group.group_id,
                atom_indices=list(group.atom_indices),
                heavy_indices=list(group.heavy_indices),
                n_nuclei=group.n_nuclei,
                provenance=Provenance(tier="A", table="shift_13c.toml"),
            )
            self._fill(sig, heavy, idx, rings)
            sig.shift += solvent.shift_offset
            out.append(sig)
        return out

    # ------------------------------------------------------------------
    def _fill(self, sig: Signal, heavy: Chem.Mol, idx: int, rings) -> None:
        atom = heavy.GetAtomWithIdx(idx)

        # --- 官能团碳直接锚定 ---
        # 注意：必须要求 idx 是 SMARTS 的**第一个**原子。
        # 否则丙酮的酮羰基模式会把两个甲基碳也包进 match，导致甲基被判成 208 ppm。
        for name, spec in self.functional:
            patt = Chem.MolFromSmarts(spec["smarts"])
            if patt is None:
                continue
            for match in heavy.GetSubstructMatches(patt, uniquify=False,
                                                  maxMatches=20000):
                if match and match[0] == idx:
                    sig.shift = float(spec["value"])
                    sig.confidence = Confidence.HIGH
                    sig.provenance.key = f"functional.{name}"
                    return

        # --- 芳香碳 ---
        if atom.GetIsAromatic():
            self._aromatic(sig, heavy, idx)
            return

        # --- sp2 烯碳（未被 functional 命中） ---
        if atom.GetHybridization() == Chem.HybridizationType.SP2:
            self._alkene(sig, heavy, idx)
            return

        # --- sp3：Grant-Paul + 取代基 ---
        self._sp3(sig, heavy, idx, rings)

        # 芳环 n-CH3 用同一套分档覆盖（加和法分不开酰亚胺/咪唑的 N-甲基）
        if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 3:
            tbl = self.table["functional"].get("n_methyl_aromatic")
            if tbl:
                for nb in atom.GetNeighbors():
                    if nb.GetSymbol() != "N" or not nb.GetIsAromatic():
                        continue
                    n_co = sum(
                        1 for m in nb.GetNeighbors()
                        if m.GetSymbol() == "C" and _is_carbonyl(m)
                    )
                    sig.shift = float(
                        tbl["two_carbonyl"] if n_co >= 2
                        else tbl["one_carbonyl"] if n_co == 1
                        else tbl["none"]
                    )
                    sig.confidence = Confidence.MEDIUM
                    sig.provenance.key = f"functional.n_methyl_aromatic({n_co} C=O)"
                    break

    # ------------------------------------------------------------------
    def _alkene(self, sig: Signal, heavy: Chem.Mol, idx: int) -> None:
        """烯碳：δ = base + Σ own + Σ far，见 [alkene] 表的说明。

        置信度上限是 MEDIUM，且只要有一个取代基没被识别就降到 LOW。
        原先那条"一律 125.0 + high"的锚定，错得最离谱（38.5 ppm）的原子恰恰
        标得最自信，这是不能接受的。
        """
        tbl = self.table["alkene"]
        base = float(tbl["base"])
        partner = _double_bond_partner(heavy, idx)

        if partner is None:
            sig.shift = base
            sig.confidence = Confidence.LOW
            sig.provenance.key = "alkene.base"
            sig.provenance.note = "非典型双键（对侧不是碳），按乙烯粗估"
            return

        own, unknown_own = self._alkene_side(heavy, idx, partner, self.alkene_own)
        far, unknown_far = self._alkene_side(heavy, partner, idx, self.alkene_far)

        total = base + sum(v for _, v in own) + sum(v for _, v in far)
        parts = [f"own{_describe(own)}", f"far{_describe(far)}"]
        unknown = unknown_own + unknown_far

        sig.shift = total
        sig.provenance.key = "alkene " + " ".join(parts)
        sig.confidence = Confidence.MEDIUM if unknown == 0 else Confidence.LOW
        sig.provenance.note = (
            "烯碳两端加和"
            if unknown == 0
            else f"烯碳两端加和；{unknown} 个取代基未识别，只按已知部分估算"
        )

    def _alkene_side(
        self, heavy: Chem.Mol, center: int, other: int, subs: SubstituentSet
    ) -> tuple[list[tuple[str, float]], int]:
        """一侧的取代基增量列表，以及没能识别的取代基个数。

        同一碳上有多个取代基时，第二个及以后的 π/杂原子项按 ``pi_attenuation``
        打折——两个 π 体系抢同一套轨道，加不满（异戊二烯 vs 异丁烯，见 TOML 注释）。
        """
        # 同一碳上挂多个取代基时，π 项按 ``pi_attenuation`` 打折。两个 σ 给体能线性
        # 相加（异丁烯两个甲基 143.4 = 123.3 + 2×10.1），两个 π 体系抢同一套轨道就不行
        # （异戊二烯 甲基+乙烯基实测 142.5，线性相加会得到 147.6）。
        att = float(self.table["alkene"]["pi_attenuation"])
        raw: list[tuple[str, float, bool]] = []
        unknown = 0
        for nb in heavy.GetAtomWithIdx(center).GetNeighbors():
            j = nb.GetIdx()
            if j == other or nb.GetAtomicNum() == 1:
                continue
            hit = subs.match_at(heavy, j)
            if hit is None:
                unknown += 1
                continue
            name, spec = hit
            value = float(spec["value"])
            step = float(spec.get("beta_step", 0.0))
            if step:
                # 只算第一个"链继续"的碳：1-丁烯与 1-己烯的烯碳实测都是 ~139.5，
                # 说明这个增量是饱和的，按碳数线性累加会把长链算飞
                beta_carbons = min(_chain_carbon_neighbors(heavy, j, center), 1)
                value += step * beta_carbons
            raw.append((name, value, name in self.alkene_pi))

        crowded = len(raw) >= 2
        return [
            (name, value * att if crowded and is_pi else value)
            for name, value, is_pi in raw
        ], unknown

    # ------------------------------------------------------------------
    def _sp3(self, sig: Signal, heavy: Chem.Mol, idx: int, rings) -> None:
        gp = self.table["sp3"]["grant_paul"]
        atom = heavy.GetAtomWithIdx(idx)

        # 统计 α/β/γ/δ 位的碳原子数。
        # 芳环碳照样计入（甲苯甲基的 α 就是芳碳），但**不再向外展开**，
        # 否则整个苯环的碳都会被当成 β/γ，甲基会被算出 35 ppm。
        counts = {1: 0, 2: 0, 3: 0, 4: 0}
        visited = {idx}
        frontier = [idx]
        for level in range(1, 5):
            nxt: list[int] = []
            for cur in frontier:
                for nb in heavy.GetAtomWithIdx(cur).GetNeighbors():
                    j = nb.GetIdx()
                    if j in visited or nb.GetAtomicNum() == 1:
                        continue
                    visited.add(j)
                    if nb.GetSymbol() == "C":
                        counts[level] += 1
                    if not nb.GetIsAromatic():
                        nxt.append(j)
            frontier = nxt
            if not frontier:
                break
        alpha_c, beta_c = counts[1], counts[2]
        gamma_c, delta_c = counts[3], counts[4]

        beta_term = (
            float(gp["beta_first"]) * min(beta_c, 2)
            + float(gp["beta_extra"]) * max(0, beta_c - 2)
        )
        total = (
            float(gp["constant"])
            + float(gp["alpha"]) * alpha_c
            + beta_term
            + float(gp["gamma"]) * gamma_c
            + float(gp["delta"]) * delta_c
        )
        parts = [f"GP(C{alpha_c}/{beta_c}/{gamma_c}/{delta_c})"]

        # 同碳多氟的增量是**饱和**的：CH3F 75、CH2F2 109、CHF3 118、CF4 120，
        # 线性累加三个 63 会算出 200 ppm 这种荒谬值。
        n_f = sum(
            1 for nb in atom.GetNeighbors()
            if nb.GetSymbol() == "F" and nb.GetAtomicNum() == 9
        )
        if n_f:
            sat = self.table["sp3"]["fluorine_saturation"]
            contrib = sum(float(sat[k]) for k in list(sat)[: min(n_f, len(sat))])
            total += contrib
            parts.append(f"F×{n_f}(饱和)")

        for name, level, spec in scan_alpha_beta_gamma(heavy, idx, self.sp3_subs):
            if name == "F":
                continue                       # 已在上面按饱和处理
            inc = float(spec.get(LEVEL_KEY[level], 0.0))
            if inc:
                total += inc
                parts.append(f"{name}@{level}")

        # 羰基对 α 碳的额外去屏蔽：只在本碳的邻位确实是羰基碳时加一次
        neighbor_ids = {nb.GetIdx() for nb in atom.GetNeighbors()}
        applied: set[str] = set()
        for name, spec in self.table["sp3"]["carbonyl"].items():
            patt = Chem.MolFromSmarts(spec["smarts"])
            if patt is None:
                continue
            targets = {m[0] for m in heavy.GetSubstructMatches(patt)}
            if neighbor_ids & targets and name not in applied:
                total += float(spec["a"])
                applied.add(name)
                parts.append(f"C=O({name}).a")

        ring_corr, ring_tag = self._ring_correction(heavy, idx, rings.get(idx, []))
        if ring_corr:
            total += ring_corr
            parts.append(f"{ring_tag}({ring_corr:+.1f})")

        sig.shift = total
        sig.confidence = Confidence.MEDIUM
        sig.provenance.key = " + ".join(parts)

    # ------------------------------------------------------------------
    def _aromatic(self, sig: Signal, heavy: Chem.Mol, idx: int) -> None:
        ring = aromatic_ring_of(heavy, idx)
        subs = ring_substituents(heavy, ring, self.arom_subs) if ring else []
        het = self.table["heteroaromatic"]

        if ring is not None and any(
            heavy.GetAtomWithIdx(i).GetSymbol() != "C" for i in ring
        ):
            base = self._hetero_base(heavy, idx, ring, het)
            scale = float(het.get("subst_scale", 0.70))
            sig.provenance.key = "heteroaromatic"
            sig.confidence = Confidence.LOW
            sig.provenance.note = "杂环碳按杂原子效应粗估"
        else:
            base = float(self.table["aromatic"]["base"])
            scale = 1.0
            sig.provenance.key = "aromatic.base"
            sig.confidence = Confidence.MEDIUM
            fused = self._fused_delta(heavy, idx, ring)
            if fused:
                base += fused
                sig.provenance.key = f"aromatic.fused({fused:+.1f})"

        total = base
        parts = []
        for ring_atom, name, spec in subs:
            dist = ring_distance(ring, idx, ring_atom) if ring else -1
            if dist == 0:
                pos = "i"          # 取代基本身所在的碳（ipso）也要加 i 增量
            else:
                pos = relative_position(dist, len(ring) if ring else 6)
            if not pos:
                continue
            inc = float(spec.get(pos, 0.0)) * scale
            total += inc
            parts.append(f"{name}.{pos}")
        if parts:
            sig.provenance.key += " + " + ", ".join(parts)

        sig.shift = total

    @staticmethod
    def _smallest_ring(ring_list: list[int]) -> int | None:
        """该原子所在的、修正表覆盖的最小环大小；没有则 None。"""
        sizes = [s for s in ring_list if 3 <= s <= 8]
        return min(sizes) if sizes else None

    def _ring_correction(
        self, heavy: Chem.Mol, idx: int, ring_list: list[int]
    ) -> tuple[float, str]:
        """环修正。GP 不含环项，算环烷会系统性偏高 1~9 ppm。

        含 C=C 的环另用一张表：环己烯是半椅式、环戊烯接近平面，饱和环的常数不适用，
        而且两者的偏差**方向相反**（环己烯烯丙位偏高 3.7、环戊烯偏低 5.0）。
        """
        size = self._smallest_ring(ring_list)
        if size is None:
            return 0.0, ""
        role = self._ring_alkene_role(heavy, idx, size)
        corr, tag = 0.0, ""
        if role is not None:
            sub = self.table["sp3"].get("ring_alkene", {}).get(role, {})
            if str(size) in sub:
                corr, tag = float(sub[str(size)]), f"ring{size}-{role}"
        if tag:
            return corr, tag
        table = self.table["sp3"].get("ring", {})
        return float(table.get(str(size), 0.0)), f"ring{size}"

    @staticmethod
    def _ring_alkene_role(heavy: Chem.Mol, idx: int, size: int) -> str | None:
        """这个 sp³ 碳处在"含 C=C 的环"里吗？是的话返回角色。

        - ``"allylic"``：同一环内的邻居里有 sp² 碳（直接连在双键上）
        - ``"other"``：在含双键的环里但不邻接双键
        - ``None``：所在环是饱和的（或用的是通用的 ``[sp3.ring]``）

        只认**环内两个碳之间**的双键。环己酮的 C=O 不算——氧不在环里，
        遍历环内相邻原子对时自然看不到它。
        """
        info = heavy.GetRingInfo()
        for ring in info.AtomRings():
            if len(ring) != size or idx not in ring:
                continue
            has_double = False
            for a, b in zip(ring, list(ring[1:]) + [ring[0]]):
                bond = heavy.GetBondBetweenAtoms(a, b)
                if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                if (
                    heavy.GetAtomWithIdx(a).GetSymbol() == "C"
                    and heavy.GetAtomWithIdx(b).GetSymbol() == "C"
                ):
                    has_double = True
                    break
            if not has_double:
                return None
            allylic = any(
                nb.GetIdx() in ring
                and nb.GetHybridization() == Chem.HybridizationType.SP2
                for nb in heavy.GetAtomWithIdx(idx).GetNeighbors()
            )
            return "allylic" if allylic else "other"
        return None

    # ------------------------------------------------------------------
    def _fused_delta(self, heavy: Chem.Mol, idx: int, ring: list[int] | None) -> float:
        """稠合芳环的桥头 / α / β / meso 修正，非稠合位置返回 0。

        ``aromatic.base`` 只对单环成立：萘的桥头实测 133.5、β 位 125.9，
        用同一个 128.5 会让桥头错 5 ppm。分档见 [aromatic] 的 fused_* 常量。

        判"是否处于稠合体系"要看**整个环**有没有桥头原子——萘的 β 位自己只属于
        一个环（``NumAtomRings == 1``），按自身判断会把它当成单环苯碳。
        """
        tbl = self.table["aromatic"]
        ring_info = heavy.GetRingInfo()
        if not ring or max(ring_info.NumAtomRings(a) for a in ring) < 2:
            return 0.0
        atom = heavy.GetAtomWithIdx(idx)
        if sum(1 for b in atom.GetBonds() if b.IsInRing()) >= 3:
            return float(tbl["fused_bridgehead"])

        # 数一数有几个邻居本身也是桥头碳
        fused_neighbors = sum(
            1
            for nb in atom.GetNeighbors()
            if nb.GetAtomicNum() == 6 and ring_info.NumAtomRings(nb.GetIdx()) >= 2
        )
        if fused_neighbors >= 2:
            return float(tbl["fused_meso"])        # 蒽的 9/10 位
        if fused_neighbors == 1:
            return float(tbl["fused_alpha"])
        return float(tbl["fused_beta"])

    @staticmethod
    def _hetero_base(heavy: Chem.Mol, idx: int, ring: list[int], het: dict) -> float:
        base = float(het.get("base", 128.5))
        total = base
        for nb in heavy.GetAtomWithIdx(idx).GetNeighbors():
            if not nb.GetIsAromatic():
                continue
            sym = nb.GetSymbol()
            if sym == "N":
                key = "alpha_ring_NH" if nb.GetTotalNumHs() > 0 else "alpha_ring_N"
            elif sym == "O":
                key = "alpha_ring_O"
            elif sym == "S":
                key = "alpha_ring_S"
            else:
                continue
            total += float(het.get(key, 0.0))
        for i in ring:
            atom = heavy.GetAtomWithIdx(i)
            if atom.GetSymbol() == "N" and ring_distance(ring, idx, i) == 2:
                key = "beta_ring_NH" if atom.GetTotalNumHs() > 0 else "beta_ring_N"
                total += float(het.get(key, 0.0))
            elif atom.GetSymbol() == "O" and ring_distance(ring, idx, i) == 2:
                total += float(het.get("beta_ring_O", 0.0))
            elif atom.GetSymbol() == "S" and ring_distance(ring, idx, i) == 2:
                total += float(het.get("beta_ring_S", 0.0))
        return total
