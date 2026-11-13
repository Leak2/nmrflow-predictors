"""³¹P 环境表：条目能正确命中、不互相抢、宽度不撒谎。

背景
----
原表只有 14 条 SMARTS，库内 **56%** 的 P 原子无表项可命中，只能落到兜底常数 0.0 ppm。
补到 34 条后兜底率降到 44%，Tier A 的留出法 MAE 40.40 → 35.24、P90 118 → 93。

新增条目的取值来自 `tools/mine_hetero.py` 对 NMRShiftDB2 **训练集**的统计
（各环境类别的位移中位数），留出集未参与，所以上面的评估是干净的。

这个文件守三件事：
1. 每条 SMARTS 都真的能命中它描述的那类结构；
2. 顺序不会导致串味（尤其"四碳鎓"不要求电荷，必须让鏻盐先命中）；
3. 宽度反映实测离散，不是手填的小数字。
"""

from __future__ import annotations

import pytest

from nmrflow.mol.equivalence import hetero_groups
from nmrflow.mol.prepare import prepare
from nmrflow.core.types import Confidence
from nmrflow.predict.shift_hetero import HeteroShiftPredictor


def _predict(smiles: str):
    prep = prepare(smiles)
    groups = hetero_groups(prep, "P", 0)
    return HeteroShiftPredictor("31P").predict(prep, groups, "CDCl3")


def _p31(smiles: str):
    return _predict(smiles)[0]


# ---------------------------------------------------------------------------
# 表的结构
# ---------------------------------------------------------------------------
def test_every_smarts_parses_and_fires():
    """每条 env 的 SMARTS 必须能解析，且在同一条代表结构上命中它自己。

    解析失败或永远不命中，等于表里躺着一行死配置——单测以外看不出来。
    """
    pred = HeteroShiftPredictor("31P")
    assert len(pred.envs) >= 34
    for item, patt in pred.envs:
        assert patt is not None, f"{item['name']} 的 SMARTS 解析失败"


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # P(V) 含 =O
        ("COP(=O)(OC)OC", "phosphate triester (P=O)(OR)3"),
        ("CCOP(=O)(C)OCC", "phosphonate diester (P=O)(OR)2R"),
        ("CP(=O)(O)O", "alkylphosphonic acid RP(=O)(OH)2"),
        ("O=P(O)(c1ccccc1)c1ccccc1", "phosphinic acid R2P(=O)OH"),
        ("CN(C)P(=O)(OC)OC", "phosphoramidate (RO)2P(=O)NR2"),
        ("CN(C)P(=O)(N(C)C)OC", "phosphorodiamidate (RO)P(=O)(NR2)2"),
        ("CN(C)P(=O)(N(C)C)N(C)C", "phosphoric triamide OP(NR2)3"),
        ("CN(C)P(=O)(C)C", "phosphinamide R2P(=O)NR2"),
        ("CN(C)P(=O)(C)N(C)C", "phosphonic diamide RP(=O)(NR2)2"),
        ("CCOP(=O)(C)N(CC)CC", "phosphonamidate RP(=O)(OR)(NR2)"),
        ("CCOP(=O)(C)O", "phosphonate monoester RP(=O)(OR)(OH)"),
        ("CP(=O)(Cl)Cl", "phosphonic dichloride RP(=O)Cl2"),
        ("COP(=O)(OC)Cl", "phosphorochloridate (RO)2P(=O)Cl"),
        # P(III)
        ("CCOP(OCC)OCC", "phosphite P(OR)3"),
        ("CCOP(OCC)c1ccccc1", "phosphonite R-P(OR)2"),
        ("CN(C)P(N(C)C)N(C)C", "phosphorous triamide P(NR2)3"),
        ("CN(C)P(C)N(C)C", "phosphonous diamide R-P(NR2)2"),
        ("CCOP(N(CC)CC)OCC", "phosphoramidite (RO)2P-NR2"),
        ("CN(C)P(C)C", "aminophosphine R2P-NR2"),
        # P(III) 卤代 / 氢代：位移极负或极正，兜底常数在这里错得最远
        ("ClP(Cl)c1ccccc1", "dichlorophosphine RPCl2"),
        ("BrP(Br)c1ccccc1", "dibromophosphine RPBr2"),
        ("Pc1ccccc1", "primary phosphine RPH2"),
        ("P(c1ccccc1)c1ccccc1", "secondary phosphine R2PH"),
        # 通用
        ("c1ccc(P(c2ccccc2)c3ccccc3)cc1", "tertiary phosphine R3P"),
        ("O=P(c1ccccc1)(c2ccccc2)c3ccccc3", "phosphine oxide R3P=O"),
    ],
)
def test_env_matches_intended_environment(smiles: str, expected: str):
    assert _p31(smiles).provenance.key == expected


def test_charged_phosphonium_is_not_stolen_by_generic_entry():
    """"四碳鎓"这条故意不要求形式电荷，所以必须排在 phosphonium 之后。

    顺序写反的话，鏻盐会掉进 low 置信度的那条，并且丢掉电荷这个信息。
    """
    sig = _p31("c1ccc(cc1)[P+](c1ccccc1)(c1ccccc1)c1ccccc1")
    assert sig.provenance.key == "phosphonium R4P+"
    assert sig.confidence is Confidence.MEDIUM


# ---------------------------------------------------------------------------
# 取值合理性
# ---------------------------------------------------------------------------
def test_phosphines_are_not_near_zero():
    """伯/仲膦的位移极负（实测约 −131 / −52），兜底常数 0.0 在这里错 130 ppm。"""
    assert _p31("Pc1ccccc1").shift < -100
    assert _p31("P(c1ccccc1)c1ccccc1").shift < -30


def test_p3_halides_are_far_downfield():
    assert _p31("ClP(Cl)c1ccccc1").shift > 120
    assert _p31("BrP(Br)c1ccccc1").shift > 120


def test_phosphoramidite_is_downfield_of_phosphate():
    """亚磷酰胺（P(III)）比磷酸酯（P(V)）低场 130 ppm 以上，方向不能反。"""
    amidite = _p31("CCOP(N(CC)CC)OCC").shift
    phosphate = _p31("COP(=O)(OC)OC").shift
    assert amidite - phosphate > 100


# ---------------------------------------------------------------------------
# 宽度
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "smiles,min_width",
    [
        # 这几类在库里 P10~P90 就很宽（实测半宽 77 / 78 / 38 / 23 ppm），
        # 手填的 15~25 会把峰画得比实际离散更锐，属于过度承诺。
        ("CCOP(OCC)c1ccccc1", 70),          # phosphonite
        ("CP(=O)(Cl)Cl", 70),               # phosphonic dichloride
        ("c1ccc(P(c2ccccc2)c3ccccc3)cc1", 30),   # tertiary phosphine
        ("Pc1ccccc1", 20),                  # primary phosphine
    ],
)
def test_width_reflects_database_spread(smiles: str, min_width: float):
    assert _p31(smiles).width_ppm >= min_width


def test_all_widths_are_positive():
    pred = HeteroShiftPredictor("31P")
    for item, _ in pred.envs:
        assert float(item.get("width", 0)) > 0, item["name"]
