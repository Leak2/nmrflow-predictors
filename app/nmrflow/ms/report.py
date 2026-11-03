"""质谱预测编排：分子式 -> 精确质量 + 加合物 m/z + 同位素分布。"""

from __future__ import annotations

from ..core.types import MsResult
from . import adducts as adducts_mod
from . import formula as formula_mod
from .exact_mass import make_adduct_result
from .isotope import isotope_pattern


def predict_ms(
    counts: dict[str, int],
    *,
    method: str = "DART",
    adduct_names: list[str] | None = None,
    polarity: str = "positive",
    with_isotopes: bool = True,
    warnings: list[str] | None = None,
) -> MsResult:
    """由元素计数生成完整质谱预测结果。

    counts       : 中性分子的元素计数
    method       : 电离方式（DART / ESI / EI / ...）
    adduct_names : 指定加合物名；None 时用该电离方式的默认加合物
    """
    warnings = warnings if warnings is not None else []

    unsupported = formula_mod.check_elements(counts)
    if unsupported:
        warnings.append(
            "同位素表中缺少元素 " + ", ".join(unsupported) + "，精确质量不可用"
        )
        return MsResult(
            formula=formula_mod.hill_notation(counts) if counts else "",
            neutral_mass=0.0,
            method=method,
        )

    if adduct_names:
        chosen = [adducts_mod.get_adduct(n) for n in adduct_names]
    else:
        chosen = adducts_mod.adducts_for_method(method, polarity)

    results = [
        make_adduct_result(counts, ad)
        for ad in chosen
    ]

    pattern = []
    if with_isotopes and chosen:
        primary = chosen[0]
        ion_counts = primary.formula_counts(counts)
        pattern = isotope_pattern(ion_counts, charge=primary.charge)

    return MsResult(
        formula=formula_mod.hill_notation(counts),
        neutral_mass=round(formula_mod.exact_mass(counts), 6),
        method=method,
        adducts=results,
        isotope_pattern=pattern,
        rdb=formula_mod.rdb(counts),
    )
