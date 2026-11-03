"""预测器抽象接口（Tier A / B / C 统一）。"""

from __future__ import annotations

from abc import ABC, abstractmethod

from ..core.types import Confidence, Signal
from ..mol.equivalence import SpinGroup
from ..mol.prepare import PreparedMol


class ShiftPredictor(ABC):
    """化学位移预测器。"""

    nucleus: str = "1H"
    tier: str = "A"

    @abstractmethod
    def predict(
        self, prep: PreparedMol, groups: list[SpinGroup], solvent_key: str
    ) -> list[Signal] | None:
        """返回与 groups 一一对应的 Signal 列表；返回 None 表示该层未命中。"""


class CouplingPredictor(ABC):
    """J 耦合预测器。"""

    @abstractmethod
    def predict(self, prep: PreparedMol, signals: list[Signal]) -> None:
        """就地给每个 Signal 填 couplings。"""


def merge_tier(base: Signal, override: Signal) -> Signal:
    """Tier B/C 覆盖 Tier A：只覆盖位移，保留归属与耦合。"""
    base.shift = override.shift
    base.confidence = override.confidence
    base.provenance = override.provenance
    if override.fwhm_hz is not None:
        base.fwhm_hz = override.fwhm_hz
    return base


def downgrade_if(signal: Signal, condition: bool, steps: int = 1) -> Signal:
    if condition:
        signal.confidence = signal.confidence.downgrade(steps)
    return signal


def mark_fallback(signal: Signal, table: str, key: str = "", note: str = "") -> Signal:
    signal.provenance.fallback = True
    signal.provenance.table = table
    signal.provenance.key = key or signal.provenance.key
    signal.provenance.note = note
    signal.confidence = Confidence.LOW
    return signal
