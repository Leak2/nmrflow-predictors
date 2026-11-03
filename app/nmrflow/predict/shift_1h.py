"""1H 化学位移预测（Tier A：经验加和 + 整环表 + PMS 烯烃规则）。"""

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
from ..mol.labile import find_labile
from ..mol.prepare import PreparedMol, ring_sizes
from .base import ShiftPredictor

LEVEL_KEY = {1: "a", 2: "b", 3: "g"}


def _is_carbonyl_carbon(atom: Chem.Atom) -> bool:
    for bond in atom.GetBonds():
        if (
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        ):
            return True
    return False


class HeteroRingTable:
    """杂环整环基准值表：按 SMARTS 匹配环，按书写顺序定位点。"""

    def __init__(self, specs: list[dict]):
        # 环越大越具体，优先匹配（吲哚优先于吡咯，否则稠环会被误判成单环）
        ordered = sorted(specs, key=lambda s: -len(s["smarts"]))
        self.specs = []
        for spec in ordered:
            patt = Chem.MolFromSmarts(spec["smarts"])
            if patt is None:
                continue
            self.specs.append((spec["name"], patt, list(spec["shifts"])))

    def lookup(self, mol: Chem.Mol, atom_idx: int) -> tuple[str, float] | None:
        for name, patt, shifts in self.specs:
            for match in mol.GetSubstructMatches(patt):
                if atom_idx in match:
                    pos = match.index(atom_idx)
                    if pos < len(shifts) and shifts[pos] > 0.1:
                        return name, float(shifts[pos])
        return None


class ProtonShiftPredictor(ShiftPredictor):
    nucleus = "1H"
    tier = "A"

    def __init__(self) -> None:
        self.table = shift_table("1H")
        self.aliph_subs = SubstituentSet(self.table["aliphatic"]["subst"])
        self.arom_subs = SubstituentSet(self.table["aromatic"]["subst"])
        self.pms = SubstituentSet(self.table["olefin"]["pms"])
        self.hetero = HeteroRingTable(self.table.get("heteroaromatic", []))
        self.penalty = self.table["aromatic"].get("penalty", {})
        self.fallback = self.table.get("fallback", {})

    # ------------------------------------------------------------------
    def predict(
        self, prep: PreparedMol, groups: list[SpinGroup], solvent_key: str
    ) -> list[Signal]:
        solvent = get_solvent(solvent_key)
        labile = find_labile(prep.mol)
        rings = ring_sizes(prep.heavy)
        signals: list[Signal] = []

        for group in groups:
            h_idx = group.atom_indices[0]
            info = labile.get(h_idx)
            if info is not None:
                sig = self._labile_signal(group, info.cls, solvent_key, prep)
                signals.append(sig)
                continue
            signals.append(
                self._carbon_bound_signal(group, prep, solvent_key, rings)
            )
        return signals

    # ------------------------------------------------------------------
    def _labile_signal(self, group: SpinGroup, cls: str, solvent_key: str,
                       prep: PreparedMol) -> Signal:
        solvent = get_solvent(solvent_key)
        spec = solvent.labile_spec(cls) or solvent.labile_spec("ROH")
        shift = spec.center if spec else 2.5
        width = spec.width_hz if spec else 5.0

        sig = Signal(
            nucleus="1H",
            group_id=group.group_id,
            atom_indices=list(group.atom_indices),
            heavy_indices=list(group.heavy_indices),
            shift=shift,
            n_nuclei=group.n_nuclei,
            fwhm_hz=width,
            is_labile=True,
            is_exchangeable=True,
            confidence=Confidence.LOW,
            provenance=Provenance(
                tier="A", table="solvents.toml", key=f"{solvent.key}.labile.{cls}",
                note="活泼氢：溶剂依赖，实测值波动大",
            ),
        )
        if solvent.drop_exchangeable:
            sig.is_exchangeable = True
            sig.provenance.note += f"；{solvent.display} 中通常被交换掉"
        return sig

    # ------------------------------------------------------------------
    def _carbon_bound_signal(
        self, group: SpinGroup, prep: PreparedMol, solvent_key: str, rings
    ) -> Signal:
        mol = prep.mol
        heavy = prep.heavy
        heavy_idx = group.heavy_indices[0]
        atom = heavy.GetAtomWithIdx(heavy_idx)
        solvent = get_solvent(solvent_key)

        sig = Signal(
            nucleus="1H",
            group_id=group.group_id,
            atom_indices=list(group.atom_indices),
            heavy_indices=list(group.heavy_indices),
            n_nuclei=group.n_nuclei,
            provenance=Provenance(tier="A", table="shift_1h.toml"),
        )

        # --- 醛基氢 ---
        if atom.GetSymbol() == "C" and _is_carbonyl_carbon(atom):
            special = self.table["special"]
            # 只看单键邻居：双键的那个 O 是羰基氧本身，不能算取代基
            neighbors = [
                n for n in atom.GetNeighbors()
                if n.GetAtomicNum() > 1
                and atom.GetOwningMol().GetBondBetweenAtoms(
                    atom.GetIdx(), n.GetIdx()
                ).GetBondType() == Chem.BondType.SINGLE
            ]
            if any(n.GetSymbol() == "N" for n in neighbors):
                # 甲酰胺 / 甲酸酯类的甲酰氢远在高场（~8.0），不能按醛用 9.7
                sig.shift = special["formamide_H"]
                sig.provenance.key = "special.formamide_H"
            elif any(n.GetSymbol() == "O" for n in neighbors):
                sig.shift = special["formate_H"]
                sig.provenance.key = "special.formate_H"
            elif any(n.GetIsAromatic() for n in neighbors):
                sig.shift = special["aldehyde_Ar"]
                sig.provenance.key = "special.aldehyde_Ar"
            else:
                sig.shift = special["aldehyde"]
                sig.provenance.key = "special.aldehyde"
            sig.confidence = Confidence.MEDIUM
            sig.shift += solvent.shift_offset
            return sig

        # --- 芳香氢 ---
        if atom.GetIsAromatic():
            self._aromatic(sig, heavy, heavy_idx, solvent)
            return sig

        # --- 端炔氢 ---
        if atom.GetSymbol() == "C" and atom.GetHybridization() == Chem.HybridizationType.SP:
            sig.shift = self.table["special"]["alkyne_CH"]
            sig.provenance.key = "special.alkyne_CH"
            sig.confidence = Confidence.MEDIUM
            sig.shift += solvent.shift_offset
            return sig

        # --- 烯氢（PMS） ---
        alkene_partner = None
        for bond in atom.GetBonds():
            if (
                bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetOtherAtom(atom).GetSymbol() == "C"
            ):
                alkene_partner = bond.GetOtherAtom(atom).GetIdx()
                break
        if alkene_partner is not None:
            self._olefin(sig, heavy, heavy_idx, alkene_partner, group, solvent)
            return sig

        # --- 芳环上的 N-CH3 ---
        # 加和法把 N 当成普通取代基，会把酰亚胺的 N-CH3 和咪唑的 N-CH3
        # 算成同一个值；按"邻接几个羰基"分档处理。
        if atom.GetSymbol() == "C" and atom.GetTotalNumHs() == 3:
            nmethyl = self.table["special"]["n_methyl_aromatic"]
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() != "N" or not nb.GetIsAromatic():
                    continue
                n_co = sum(
                    1 for m in nb.GetNeighbors()
                    if m.GetSymbol() == "C" and _is_carbonyl_carbon(m)
                )
                sig.shift = float(
                    nmethyl["two_carbonyl"] if n_co >= 2
                    else nmethyl["one_carbonyl"] if n_co == 1
                    else nmethyl["none"]
                )
                sig.provenance.key = f"special.n_methyl_aromatic({n_co} C=O)"
                sig.confidence = Confidence.MEDIUM
                sig.shift += solvent.shift_offset
                return sig

        # --- 脂肪氢 ---
        self._aliphatic(sig, heavy, heavy_idx, group, rings, solvent)
        return sig

    # ------------------------------------------------------------------
    def _aliphatic(self, sig: Signal, heavy, idx: int, group: SpinGroup,
                   rings, solvent) -> None:
        atom = heavy.GetAtomWithIdx(idx)
        n_h = atom.GetTotalNumHs()
        base_tbl = self.table["aliphatic"]["base"]
        if 3 in rings.get(idx, []):
            base = base_tbl["CH_ring3"]
            key = "aliphatic.base.CH_ring3"
        else:
            base = {3: base_tbl["CH3"], 2: base_tbl["CH2"], 1: base_tbl["CH"]}.get(
                n_h, base_tbl["CH"]
            )
            key = f"aliphatic.base.{'CH3' if n_h == 3 else 'CH2' if n_h == 2 else 'CH'}"

        total = base
        used: list[str] = []
        for name, level, spec in scan_alpha_beta_gamma(heavy, idx, self.aliph_subs):
            inc = float(spec.get(LEVEL_KEY[level], 0.0))
            total += inc
            used.append(f"{name}@{level}")

        sig.shift = total + solvent.shift_offset
        sig.provenance.key = f"{key} + {', '.join(used) if used else 'none'}"

        # 非对映异位的 CH2：加和法给不出差异，给一个小分裂并降级
        if group.diastereotopic_pair >= 0:
            sig.shift += 0.10 if group.group_id % 2 == 0 else -0.10
            sig.confidence = Confidence.LOW
            sig.provenance.note = "非对映异位氢：加和法无法区分，已做 ±0.1 ppm 经验分裂"
        elif not used:
            sig.confidence = Confidence.MEDIUM
        else:
            sig.confidence = Confidence.MEDIUM

    # ------------------------------------------------------------------
    def _aromatic(self, sig: Signal, heavy, idx: int, solvent) -> None:
        ring = aromatic_ring_of(heavy, idx)
        hit = self.hetero.lookup(heavy, idx)
        subs = ring_substituents(heavy, ring, self.arom_subs) if ring else []
        fb = self.table.get("heteroaromatic_fallback", {})

        if hit is not None:
            name, base = hit
            scale = fb.get("subst_scale", 0.70)
            sig.provenance.key = f"heteroaromatic.{name}"
            sig.provenance.table = "shift_1h.toml"
            if subs:
                sig.confidence = Confidence.LOW
                sig.provenance.note = f"杂环 {name} 有取代基，增量按苯环 ×{scale} 近似"
            else:
                sig.confidence = Confidence.MEDIUM
        elif ring is not None and self._ring_has_hetero(heavy, ring):
            base = self._hetero_fallback(heavy, idx, ring, fb, sig)
            scale = fb.get("subst_scale", 0.70)
        else:
            base = self.table["aromatic"]["base"]
            scale = 1.0
            sig.provenance.key = "aromatic.base"
            self._apply_aromatic_penalty(sig, subs, ring)

        total = base
        parts = []
        for ring_atom, name, spec in subs:
            dist = ring_distance(ring, idx, ring_atom) if ring else -1
            pos = relative_position(dist, len(ring) if ring else 6)
            if not pos:
                continue
            inc = float(spec.get(pos, 0.0)) * scale
            total += inc
            parts.append(f"{name}.{pos}")
        if parts:
            sig.provenance.key += " + " + ", ".join(parts)

        sig.shift = total + solvent.shift_offset

    @staticmethod
    def _ring_has_hetero(mol: Chem.Mol, ring: list[int]) -> bool:
        return any(mol.GetAtomWithIdx(i).GetSymbol() not in ("C", "H") for i in ring)

    def _hetero_fallback(self, mol: Chem.Mol, idx: int, ring: list[int],
                         fb: dict, sig: Signal) -> float:
        """未命中整环表时，按环内杂原子的邻/间位效应粗估。"""
        base = float(fb.get("base", 7.26))
        n_neighbors, nh_neighbors, o_neighbors, s_neighbors = 0, 0, 0, 0
        meta_n = 0
        for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
            if not nb.GetIsAromatic():
                continue
            sym = nb.GetSymbol()
            if sym == "N":
                if nb.GetTotalNumHs() > 0:
                    nh_neighbors += 1
                else:
                    n_neighbors += 1
            elif sym == "O":
                o_neighbors += 1
            elif sym == "S":
                s_neighbors += 1
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() == "N" and atom.GetTotalNumHs() == 0:
                if ring_distance(ring, idx, i) == 2:
                    meta_n += 1

        total_n = n_neighbors + nh_neighbors
        if total_n >= 2:
            value = float(fb.get("flanked_by_two_N", 7.80))
        else:
            value = base
            value += n_neighbors * float(fb.get("alpha_ring_N", 0.55))
            value += nh_neighbors * float(fb.get("alpha_ring_NH", -0.30))
            value += o_neighbors * float(fb.get("alpha_ring_O", -0.35))
            value += s_neighbors * float(fb.get("alpha_ring_S", -0.10))
            value += meta_n * float(fb.get("meta_ring_N", 0.15))

        sig.provenance.key = "heteroaromatic_fallback"
        sig.provenance.fallback = True
        sig.confidence = Confidence.LOW
        sig.provenance.note = "未命中整环表，按环内杂原子效应粗估"
        return value

    def _apply_aromatic_penalty(self, sig: Signal, subs, ring) -> None:
        if not subs:
            return
        steps = 0
        # 邻位双取代 / 1,2,3-三取代
        if ring and len(subs) >= 2:
            atoms = [s[0] for s in subs]
            if any(
                ring_distance(ring, a, b) == 1
                for i, a in enumerate(atoms)
                for b in atoms[i + 1:]
            ):
                steps += int(self.penalty.get("ortho_disubstituted", 1))
        if len(subs) >= 3:
            steps += int(self.penalty.get("tri_substituted_123", 0))
        push = set(self.penalty.get("push_list", []))
        pull = set(self.penalty.get("pull_list", []))
        names = {s[1] for s in subs}
        if names & push and names & pull:
            steps += int(self.penalty.get("push_pull", 1))
        if steps:
            sig.confidence = sig.confidence.downgrade(steps)
            sig.provenance.note = "多取代/推拉电子体系，加和法误差较大"

    # ------------------------------------------------------------------
    def _olefin(self, sig: Signal, heavy, idx: int, partner: int,
                group: SpinGroup, solvent) -> None:
        base = self.table["olefin"]["base"]
        atom = heavy.GetAtomWithIdx(idx)
        partner_atom = heavy.GetAtomWithIdx(partner)

        # gem：观测碳上除 H 与双键伙伴之外的取代基
        gem_name, gem_spec = None, None
        for nb in atom.GetNeighbors():
            if nb.GetIdx() == partner or nb.GetAtomicNum() == 1:
                continue
            hit = self.pms.match_at(heavy, nb.GetIdx())
            if hit is not None:
                gem_name, gem_spec = hit
                break

        # 伙伴碳上的两个取代基（cis / trans）
        others = [
            nb.GetIdx() for nb in partner_atom.GetNeighbors()
            if nb.GetIdx() != idx and nb.GetAtomicNum() != 1
        ]
        partner_subs = []
        for j in others:
            hit = self.pms.match_at(heavy, j)
            partner_subs.append((j, hit))

        gem = float(gem_spec.get("gem", 0.0)) if gem_spec else 0.0

        def z(name_spec, pos: str) -> float:
            if name_spec is None:
                return 0.0
            return float(name_spec[1].get(pos, 0.0))

        n_h_on_center = atom.GetTotalNumHs()
        if n_h_on_center >= 2:
            # =CH2：两个 H 分别 cis / trans 于伙伴碳上的两个取代基
            s1 = partner_subs[0][1] if len(partner_subs) > 0 else None
            s2 = partner_subs[1][1] if len(partner_subs) > 1 else None
            # 两个 H 一组 cis 于 s1、另一组 cis 于 s2；缺立体信息时无法指认谁是谁
            if group.group_id % 2 == 0:
                total = base + z(s1, "cis") + z(s2, "trans")
            else:
                total = base + z(s2, "cis") + z(s1, "trans")
            sig.confidence = Confidence.LOW
            sig.provenance.note = "末端烯氢：cis/trans 归属未定（缺立体信息时任意指派）"
        else:
            # =CH-：立体未定时取 cis/trans 对称平均，稳健且误差小
            if len(partner_subs) >= 2:
                total = base + gem
                for _, spec in partner_subs:
                    total += (z(spec, "cis") + z(spec, "trans")) / 2.0
                sig.confidence = Confidence.MEDIUM
            elif len(partner_subs) == 1:
                spec = partner_subs[0][1]
                total = base + gem + (z(spec, "cis") + z(spec, "trans")) / 2.0
                sig.confidence = Confidence.MEDIUM
            else:
                total = base + gem
                sig.confidence = Confidence.MEDIUM

        sig.shift = total + solvent.shift_offset
        sig.provenance.key = "olefin.pms"
