"""质谱：单同位素质量、加合物、Found 抖动、卤素同位素特征比。"""

from __future__ import annotations

import pytest

from nmrflow.ms.adducts import adducts_for_method
from nmrflow.ms.exact_mass import make_adduct_result
from nmrflow.ms.formula import exact_mass, formula_from_mol, parse_formula
from nmrflow.ms.isotope import isotope_pattern
from nmrflow.mol.prepare import prepare


# 用户给的参考分子：C27H26F3N5S，[M+H]+ = 510.1934
REF = {"C": 27, "H": 26, "F": 3, "N": 5, "S": 1}


def test_reference_exact_mass():
    assert exact_mass(REF) == pytest.approx(509.1861, abs=5e-4)


def test_reference_mh_plus():
    ad = adducts_for_method("DART")[0]
    res = make_adduct_result(REF, ad)
    assert res.calculated == pytest.approx(510.1934, abs=1e-4)
    assert res.adduct == "[M+H]+"
    # 报道行里 "for C27H26F3N5S" 用的是**中性**分子式；离子式多一个 H
    assert res.formula == "C27H27F3N5S"


def test_calculated_only_has_no_found_field():
    """0.3.0 起 AdductResult 不再携带 found / found_is_simulated。"""
    ad = adducts_for_method("DART")[0]
    res = make_adduct_result(REF, ad)
    assert not hasattr(res, "found")
    assert not hasattr(res, "found_is_simulated")


def test_calculated_only_is_reproducible():
    """同一输入两次计算结果一致，且仍没有 Found。"""
    ad = adducts_for_method("DART")[0]
    a = make_adduct_result(REF, ad)
    b = make_adduct_result(REF, ad)
    assert a.calculated == b.calculated
    assert not hasattr(a, "found")


def test_formula_from_smiles():
    prep = prepare("CCO")
    assert formula_from_mol(prep.mol) == {"C": 2, "H": 6, "O": 1}


def test_hill_notation_roundtrip():
    assert parse_formula("C27H26F3N5S") == REF


@pytest.mark.parametrize(
    "smiles, mz",
    [
        ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", 195.0877),   # 咖啡因
        ("CCO", 47.0491),
        ("c1ccccc1", 79.0542),
    ],
)
def test_known_mh_plus(smiles, mz):
    prep = prepare(smiles)
    counts = formula_from_mol(prep.mol)
    ad = adducts_for_method("DART")[0]
    res = make_adduct_result(counts, ad)
    assert res.calculated == pytest.approx(mz, abs=1e-3)


def test_chlorine_isotope_ratio():
    """一氯代物的 M : M+2 应接近 3 : 1。"""
    pat = isotope_pattern({"C": 6, "H": 5, "Cl": 1})
    peaks = {round(p.mz, 0): p.abundance for p in pat}
    m = max(p.abundance for p in pat)
    mono = next(p for p in pat if abs(p.mz - 112.0) < 0.6)
    plus2 = next(p for p in pat if abs(p.mz - 114.0) < 0.6)
    assert 0.30 < plus2.abundance / mono.abundance < 0.34


def test_bromine_isotope_ratio():
    """一溴代物的 M : M+2 应接近 1 : 1。"""
    pat = isotope_pattern({"C": 6, "H": 5, "Br": 1})
    mono = next(p for p in pat if abs(p.mz - 156.0) < 0.6)
    plus2 = next(p for p in pat if abs(p.mz - 158.0) < 0.6)
    assert 0.90 < plus2.abundance / mono.abundance < 1.05
