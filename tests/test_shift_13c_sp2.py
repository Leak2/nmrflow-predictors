"""¹³C 的 sp² 碳：烯碳两端加和、稠环修正、共轭羰基。

为什么单独一个文件
------------------
原先 `shift_13c.toml` 里 `alkene = { value = 125.0 }` 把**所有**非羰基 sp² 碳
一律锚定到 125.0 ppm 并标 `confidence=high`。烯碳实测跨度 86.5~153.4 ppm，
这条规则最坏错 38.5 ppm —— 而且错得越离谱标得越自信。

它能在 98 个测试全绿的情况下活很久，是因为没有任何一条用例碰过烯碳。
所以这里的用例不是"补充覆盖"，而是这个缺陷的直接护栏。
"""

from __future__ import annotations

import pytest

from nmrflow.core.types import Confidence
from nmrflow.mol.equivalence import carbon_groups
from nmrflow.mol.prepare import prepare
from nmrflow.predict.shift_13c import CarbonShiftPredictor


def _predict(smiles: str):
    prep = prepare(smiles)
    groups = carbon_groups(prep, 0)
    return groups, CarbonShiftPredictor().predict(prep, groups, "CDCl3")


def _signals(smiles: str) -> list:
    return _predict(smiles)[1]


def _at(smiles: str, atom_idx: int):
    """SMILES 第 atom_idx 个重原子所在的那条信号。"""
    groups, sigs = _predict(smiles)
    for group, sig in zip(groups, sigs):
        if atom_idx in group.heavy_indices:
            return sig
    raise AssertionError(f"{smiles} 的原子 {atom_idx} 没有对应信号")


# ---------------------------------------------------------------------------
# 烯碳：两端加和
# ---------------------------------------------------------------------------
def test_enol_ether_ends_are_not_a_constant():
    """甲基乙烯基醚 COC=C：=CH 153.4、=CH₂ 86.5。

    这两端差 67 ppm，任何"一个常数"都不可能同时照顾到——原实现给两端都是 125.0，
    =CH₂ 那一端错 38.5 ppm。
    """
    ch = _at("COC=C", 2)          # C2 与 O 相连，是 =CH
    ch2 = _at("COC=C", 3)         # C3 是 =CH₂
    assert ch.shift == pytest.approx(153.4, abs=3.0)
    assert ch2.shift == pytest.approx(86.5, abs=5.0)
    assert ch.shift - ch2.shift > 50


def test_alkene_confidence_is_never_high():
    """走烯碳加和的那些信号不允许标 high。原实现恰恰标了 high。

    只看 ``alkene`` 开头的 provenance——乙酸乙烯酯的酯羰基有明确的 SMARTS 锚定，
    标 high 是应该的，不在本条约束范围内。
    """
    for smiles in ["C=C", "CC=C", "COC=C", "CC(=O)OC=C", "C=CC#N"]:
        alkene_sigs = [s for s in _signals(smiles) if s.provenance.key.startswith("alkene")]
        assert alkene_sigs, smiles
        for sig in alkene_sigs:
            assert sig.confidence is not Confidence.HIGH


def test_both_ends_contribute():
    """丙烷类烯烃两端各自取到自己的增量，不是一个值复制两遍。"""
    ch = _at("CC=C", 1)           # =CH，实测 133.4
    ch2 = _at("CC=C", 2)          # =CH₂，实测 115.7
    assert ch.shift == pytest.approx(133.4, abs=2.0)
    assert ch2.shift == pytest.approx(115.7, abs=2.0)


def test_chain_length_saturates():
    """beta_step 只算第一个链延伸碳：1-丁烯与 1-己烯的烯碳实测都在 139 附近。

    线性按碳数累加会把长链算到 150 以上。
    """
    for smiles, ch_idx, ch2_idx in [
        ("CCC=C", 2, 3),
        ("CCCCC=C", 4, 5),
    ]:
        assert _at(smiles, ch_idx).shift == pytest.approx(139.5, abs=2.0)
        assert _at(smiles, ch2_idx).shift == pytest.approx(114.0, abs=2.0)


def test_pi_substituents_do_not_add_linearly():
    """异戊二烯（甲基+乙烯基同碳）实测 142.5，线性相加会到 147.6。

    两个 σ 给体能线性加（异丁烯两个甲基 143.4 = 123.3 + 2×10.1），
    两个 π 体系抢同一套轨道就加不满。
    """
    assert _at("CC(=C)C=C", 1).shift == pytest.approx(142.5, abs=2.5)
    assert _at("CC(=C)C", 1).shift == pytest.approx(143.4, abs=2.0)


def test_unclassified_substituent_downgrades_confidence():
    """识别不了的取代基不能装作没事：降级并写明。"""
    sig = _at("C=C[Si](C)(C)C", 0)          # 乙烯基三甲基硅烷的 =CH₂
    assert sig.confidence is Confidence.LOW
    assert "未识别" in sig.provenance.note


def test_alkene_key_is_traceable():
    """provenance 要能看出两端各加了什么。"""
    sig = _at("CC=C", 1)
    assert sig.provenance.key.startswith("alkene ")
    assert "own(" in sig.provenance.key
    assert "far(" in sig.provenance.key


# ---------------------------------------------------------------------------
# 稠合芳环
# ---------------------------------------------------------------------------
def test_fused_bridgehead_is_not_plain_benzene():
    """萘的桥头碳实测 133.5，用 128.5 的单环 base 会错 5 ppm。"""
    prep = prepare("c1ccc2ccccc2c1")
    ring_info = prep.heavy.GetRingInfo()
    bridgeheads = [
        a.GetIdx()
        for a in prep.heavy.GetAtoms()
        if ring_info.NumAtomRings(a.GetIdx()) >= 2
        and sum(1 for b in a.GetBonds() if b.IsInRing()) >= 3
    ]
    assert len(bridgeheads) == 2
    for idx in bridgeheads:
        assert _at("c1ccc2ccccc2c1", idx).shift == pytest.approx(133.5, abs=1.5)


def test_single_ring_gets_no_fused_correction():
    """苯不该被稠环修正碰到。"""
    for sig in _signals("c1ccccc1"):
        assert sig.shift == pytest.approx(128.5, abs=0.01)
        assert "fused" not in sig.provenance.key


# ---------------------------------------------------------------------------
# 共轭羰基
# ---------------------------------------------------------------------------
def test_conjugated_carbonyl_is_lower():
    """乙酸甲酯 171.0 但苯甲酸甲酯 166.9、丙烯酸甲酯 167.0——共轭压低约 4 ppm。"""
    benzoate = _at("COC(=O)c1ccccc1", 2)
    acetate = _at("CC(=O)OC", 1)
    acrylate = _at("COC(=O)C=C", 2)
    assert benzoate.shift == pytest.approx(167.0, abs=1.5)
    assert acrylate.shift == pytest.approx(167.0, abs=1.5)
    assert acetate.shift == pytest.approx(171.0, abs=1.0)
    assert "conjugated" in benzoate.provenance.key


def test_conjugated_acid_is_lower():
    """乙酸 178.0，苯甲酸 172.5。"""
    assert _at("OC(=O)c1ccccc1", 1).shift == pytest.approx(172.5, abs=1.5)
    assert _at("CC(=O)O", 1).shift == pytest.approx(178.0, abs=1.0)
