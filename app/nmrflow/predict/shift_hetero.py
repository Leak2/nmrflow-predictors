"""19F / 31P 化学位移：SMARTS 环境直接查表（不做加和）。"""

from __future__ import annotations

from rdkit import Chem

from ..core.constants import VAGUE_FWHM
from ..core.types import Confidence, Provenance, Signal
from ..data.loader import shift_table
from ..mol.equivalence import SpinGroup
from ..mol.prepare import PreparedMol
from .base import ShiftPredictor

_NUCLEUS_TABLE = {"19F": "shift_19f.toml", "31P": "shift_31p.toml"}


class HeteroShiftPredictor(ShiftPredictor):
    """按环境 SMARTS 查表；未命中时给宽峰 + low 置信度，绝不假装精确。"""

    def __init__(self, nucleus: str):
        self.nucleus = nucleus
        self.tier = "A"
        self.table = shift_table(nucleus)
        self.envs = []
        for item in self.table.get("env", []):
            patt = Chem.MolFromSmarts(item["smarts"])
            if patt is not None:
                self.envs.append((item, patt))
        self.fb = self.table.get("fallback", {})

    def predict(
        self, prep: PreparedMol, groups: list[SpinGroup], solvent_key: str
    ) -> list[Signal]:
        heavy = prep.heavy
        out: list[Signal] = []
        for group in groups:
            idx = group.atom_indices[0]
            sig = Signal(
                nucleus=self.nucleus, group_id=group.group_id,
                atom_indices=list(group.atom_indices),
                heavy_indices=list(group.heavy_indices),
                n_nuclei=group.n_nuclei,
                provenance=Provenance(tier="A", table=_NUCLEUS_TABLE[self.nucleus]),
            )
            hit = None
            for item, patt in self.envs:
                # match[0] 必须是待查原子本身：表描述的是"这个核处在什么环境"，
                # 而不是"这个核碰巧出现在某个更大的结构里"。后者会让 SO2F 上的
                # F 被 CF2 条目抢走（只要模式里含 F 就算命中）。
                # uniquify=False：否则 CF₃ 这类对称基团只回一个匹配，
                # 组内代表原子恰好不是它时，整组掉进兜底。
                for match in heavy.GetSubstructMatches(
                    patt, uniquify=False, maxMatches=20000
                ):
                    if match[0] == idx:
                        hit = item
                        break
                if hit:
                    break
            if hit is not None:
                sig.shift = float(hit["value"])
                width_ppm = float(hit.get("width", 10.0))
                conf = hit.get("confidence", "medium")
                sig.confidence = Confidence(conf) if conf in list(Confidence) \
                    else Confidence.MEDIUM
                sig.provenance.key = hit["name"]
                sig.width_ppm = width_ppm
                sig.provenance.note = f"环境跨度约 ±{width_ppm:.0f} ppm"
            else:
                sig.shift = float(self.fb.get("value", 0.0))
                sig.confidence = Confidence.LOW
                sig.provenance.fallback = True
                sig.provenance.key = "fallback"
                sig.provenance.note = self.fb.get("note", "")
                sig.width_ppm = float(self.fb.get("width", 40.0))
                sig.fwhm_hz = VAGUE_FWHM
            out.append(sig)
        return out
