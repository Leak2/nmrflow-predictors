"""¹³C：带电羧酸根、甲酰类化合物的官能团锚定。

为什么单独一个文件
------------------
漏写电荷的 SMARTS 不会抛异常，只会**静默命中别的规则**，图谱照样出，数值
看着像个正常位移。这类错误只能靠测试钉住：

1. 羧酸根 `[O-]` 匹配不上 `[OX2H]`（那是中性羧酸），于是掉到 ``[alkene]`` 的
   两端加和，乙酸钠的羧基碳得到乙烯基值 123.3 —— 实测 182.0，误差 59 ppm。
   而且它比中性羧酸更隐蔽：中性那边好歹有个 178.0 的锚，阴离子这边直接裸奔。

2. ``aldehyde`` 原来的 SMARTS 是 ``[CX3H1]=[OX1]``，只要"碳上带 H 且连羰基"
   就算醛。于是甲酰胺（连 N）被抢成 200.0（实测 165，误差 35 ppm）、甲酸根
   同样被抢。醛的羰基碳另一个邻居必须是碳或氢，收紧后两者各归各位。

3. SMARTS 里 ``X`` 是**总连接数**（含隐式 H）。羧酸根氧只连 1 个原子、0 个 H，
   所以是 ``[OX1-]``；写成 ``[OX2-]`` 匹配为空。这里也一并钉住，
   免得以后有人"顺手修正"回去。

下面的断言全部走 ``CarbonShiftPredictor`` 直连（纯 Tier A），不经 pipeline，
这样 HOSE 索引命不命中都不会影响结论。
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from nmrflow.mol.equivalence import carbon_groups
from nmrflow.mol.prepare import prepare
from nmrflow.predict.shift_13c import CarbonShiftPredictor


def _signals(smiles: str) -> list:
    prep = prepare(smiles)
    groups = carbon_groups(prep, 0)
    return CarbonShiftPredictor().predict(prep, groups, "CDCl3")


def _carbonyl_signal(smiles: str, smarts: str):
    """按 SMARTS 定位羰基碳，返回它的信号。"""
    prep = prepare(smiles)
    patt = Chem.MolFromSmarts(smarts)
    match = prep.heavy.GetSubstructMatch(patt)
    if not match:
        raise AssertionError(f"{smiles} 里没有匹配 {smarts} 的原子")
    target = match[0]
    groups = carbon_groups(prep, 0)
    for group, sig in zip(groups, CarbonShiftPredictor().predict(prep, groups, "CDCl3")):
        if target in group.heavy_indices:
            return sig
    raise AssertionError(f"{smiles} 的原子 {target} 没有对应信号")


# 实测值来源：SDBS / 文献 CDCl3 或 D2O 谱
CARBOXYLATES = [
    # SMILES, 说明, 实测 ppm
    ("CC(=O)[O-].[Na+]", "乙酸钠", 182.0),
    ("CCC(=O)[O-].[Na+]", "丙酸钠", 183.5),
    ("[O-]C(=O)c1ccccc1.[Na+]", "苯甲酸钠", 175.6),
]


@pytest.mark.parametrize(("smiles", "name", "expected"), CARBOXYLATES)
def test_carboxylate_uses_dedicated_rule(smiles, name, expected):
    """羧酸根必须走 functional.carboxylate*，不许掉到 alkene.base。"""
    sig = _carbonyl_signal(smiles, "[CX3](=[OX1])[OX1-]")
    assert sig.provenance.key.startswith("functional.carboxylate"), (
        f"{name} 的羧基碳没走羧酸根规则，实际命中 {sig.provenance.key}"
    )
    assert sig.shift == pytest.approx(expected, abs=3.0)


def test_carboxylate_never_falls_back_to_alkene_base():
    """回归：乙酸钠曾经得到 123.3（乙烯基值），误差 59 ppm。"""
    sig = _carbonyl_signal("CC(=O)[O-].[Na+]", "[CX3](=[OX1])[OX1-]")
    assert sig.shift > 170.0, f"羧基碳掉到了 {sig.shift:.2f}，疑似回落到烯碳基值"
    assert "alkene" not in sig.provenance.key


def test_formate_has_its_own_rule():
    """甲酸根碳上带 H，比同族低 10 ppm，不能与乙酸根共用 182.0。"""
    sig = _carbonyl_signal("[O-]C(=O).[Na+]", "[CX3H1](=[OX1])[OX1-]")
    assert sig.provenance.key == "functional.carboxylate_formate"
    assert sig.shift == pytest.approx(171.8, abs=2.0)
    # 与乙酸根必须分得开，否则说明 formate 规则没被命中
    acetate = _carbonyl_signal("CC(=O)[O-].[Na+]", "[CX3](=[OX1])[OX1-]")
    assert abs(sig.shift - acetate.shift) > 5.0


def test_conjugated_carboxylate_differs_from_aliphatic():
    """共轭羧酸根（苯甲酸钠）低于脂肪羧酸根（乙酸钠），与中性一侧规律一致。"""
    aliphatic = _carbonyl_signal("CC(=O)[O-].[Na+]", "[CX3](=[OX1])[OX1-]")
    conjugated = _carbonyl_signal("[O-]C(=O)c1ccccc1.[Na+]", "[CX3](=[OX1])[OX1-]")
    assert conjugated.provenance.key == "functional.carboxylate_conjugated"
    assert conjugated.shift < aliphatic.shift


# ---------------------------------------------------------------------------
# 中性羧酸：加了阴离子规则后必须纹丝不动
# ---------------------------------------------------------------------------
def test_neutral_carboxylic_acid_unchanged():
    sig = _carbonyl_signal("CC(=O)O", "[CX3](=[OX1])[OX2H]")
    assert sig.provenance.key == "functional.carboxylic_acid"
    assert sig.shift == pytest.approx(178.0, abs=0.1)


def test_carboxylate_and_neutral_acid_do_not_cross_match():
    """中性酸走 178、阴离子走 182，两条规则互不串门。"""
    acid = _carbonyl_signal("CC(=O)O", "[CX3](=[OX1])[OX2H]")
    salt = _carbonyl_signal("CC(=O)[O-].[Na+]", "[CX3](=[OX1])[OX1-]")
    assert acid.provenance.key != salt.provenance.key
    assert salt.shift > acid.shift


# ---------------------------------------------------------------------------
# 醛 SMARTS 过宽：甲酰胺曾被抢成 200.0
# ---------------------------------------------------------------------------
def test_formamide_not_captured_by_aldehyde():
    """回归：甲酰胺实测 165，曾被 aldehyde 抢成 200（误差 35 ppm）。"""
    sig = _carbonyl_signal("O=CN", "[CX3H1](=[OX1])")
    assert sig.provenance.key != "functional.aldehyde", "甲酰胺又被醛规则抢走了"
    assert sig.shift < 175.0, f"甲酰胺得到 {sig.shift:.2f}，偏高"
    assert sig.shift == pytest.approx(165.0, abs=6.0)


@pytest.mark.parametrize(
    ("smiles", "name", "lo", "hi"),
    [
        ("CC=O", "乙醛", 195.0, 205.0),
        ("CCC=O", "丙醛", 197.0, 207.0),
        # Tier A 的 aldehyde 只有单一锚定值 200.0，芳香醛的 ~192 是 Tier B
        # 查表修正出来的（pipeline 里实测 192.04）。这里只验"仍然是醛"，
        # 不验精度——否则测的是 HOSE 索引而非本文件关心的规则匹配。
        ("O=Cc1ccccc1", "苯甲醛", 198.0, 202.0),
    ],
)
def test_aldehyde_still_matches(smiles, name, lo, hi):
    """收紧 SMARTS 后，真正的醛必须照旧命中。"""
    sig = _carbonyl_signal(smiles, "[CX3H1](=[OX1])[#6,#1]")
    assert sig.provenance.key == "functional.aldehyde"
    assert lo <= sig.shift <= hi, f"{name} 得到 {sig.shift:.2f}，超出 {lo}~{hi}"


def test_methyl_formate_not_captured_by_aldehyde():
    """甲酸甲酯同样是甲酰碳，不能被醛抢走（实测约 160）。"""
    sig = _carbonyl_signal("O=COC", "[CX3H1](=[OX1])")
    assert sig.provenance.key != "functional.aldehyde"
    assert sig.shift < 175.0, f"甲酸甲酯得到 {sig.shift:.2f}，偏高"
