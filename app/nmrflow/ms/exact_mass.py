"""精确质量与 m/z 计算。"""

from __future__ import annotations

from dataclasses import dataclass

from ..core.constants import ELECTRON_MASS
from ..core.types import AdductResult
from .adducts import Adduct
from .formula import exact_mass, hill_notation
from .masses import monoisotopic_mass


def adduct_mz(counts: dict[str, int], adduct: Adduct) -> float:
    """计算加合物的质荷比。

    mass = mult * M_neutral + Σ(add) - Σ(sub) - charge * m_e
    m/z  = mass / |charge|
    """
    neutral = exact_mass(counts)
    delta = 0.0
    for sym, n in adduct.add.items():
        delta += monoisotopic_mass(sym) * n
    for sym, n in adduct.sub.items():
        delta -= monoisotopic_mass(sym) * n
    ion_mass = adduct.mult * neutral + delta - adduct.charge * ELECTRON_MASS
    return ion_mass / adduct.abs_charge


def make_adduct_result(counts: dict[str, int], adduct: Adduct) -> AdductResult:
    """生成一条加合物结果；只给 Calculated，不生成 Found。"""
    calc = adduct_mz(counts, adduct)
    return AdductResult(
        adduct=adduct.name,
        formula=hill_notation(adduct.formula_counts(counts)),
        calculated=round(calc, 4),
        charge=adduct.charge,
    )


@dataclass
class MassSummary:
    formula: str
    neutral: float
    average: float
    rdb: float
