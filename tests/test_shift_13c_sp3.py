"""¹³C 的 sp³ 碳：环烷修正、羧基重复计数。

为什么单独一个文件
------------------
Grant-Paul 是为**开链**烷烃拟合的，公式里没有任何环项。直接拿去算环烷会系统性偏高：
环己烷 32.2（实测 27.1）、环戊烷 34.7（实测 25.9）、环丙烷 10.9（实测 5.7）。

这个偏差会污染所有"环上的 sp³ 碳"，而且容易被误诊。环己烯的烯丙位 CH₂ 原来是 35.2、
实测 25.4，看着像"烯烃增量给错了"，其实 5.1 ppm 是环基线带进来的——剩下那部分才是
烯烃增量的问题，且至今没修好（见 README 已知局限）。

另一处是羧基被算两遍：C-C(=O)-O 里的 O 与 α 位的羰基同属一个官能团，
`OCOR@2` 和 `C=O(ester_acid_amide).a` 各加了一次，丙酸 C2 因此多算 6 ppm。
"""

from __future__ import annotations

import pytest

from nmrflow.data.loader import shift_table
from nmrflow.mol.equivalence import carbon_groups
from nmrflow.mol.prepare import prepare
from nmrflow.predict.shift_13c import CarbonShiftPredictor


def _predict(smiles: str):
    prep = prepare(smiles)
    groups = carbon_groups(prep, 0)
    return groups, CarbonShiftPredictor().predict(prep, groups, "CDCl3")


def _signals(smiles: str) -> list:
    return _predict(smiles)[1]


def _sp3_signals(smiles: str) -> list:
    """只留 sp³ 碳的信号——环烯的分子里还混着 sp² 碳。"""
    prep = prepare(smiles)
    return [
        s for s in _signals(smiles)
        if prep.heavy.GetAtomWithIdx(s.heavy_indices[0]).GetHybridization().name == "SP3"
    ]


def _at(smiles: str, atom_idx: int):
    groups, sigs = _predict(smiles)
    for group, sig in zip(groups, sigs):
        if atom_idx in group.heavy_indices:
            return sig
    raise AssertionError(f"{smiles} 的原子 {atom_idx} 没有对应信号")


# ---------------------------------------------------------------------------
# 环烷修正
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "smiles,expected",
    [
        ("C1CC1", 5.7),          # 环丙烷
        ("C1CCC1", 22.8),        # 环丁烷
        ("C1CCCC1", 25.9),       # 环戊烷
        ("C1CCCCC1", 27.1),      # 环己烷
        ("C1CCCCCC1", 28.7),     # 环庚烷
        ("C1CCCCCCC1", 26.9),    # 环辛烷
    ],
)
def test_cycloalkane_shift(smiles: str, expected: float):
    """六种环烷逐一对上。修正值是"实测 − GP"，所以这条用例等于把标定表钉住。"""
    for sig in _signals(smiles):
        assert sig.shift == pytest.approx(expected, abs=0.2)


def test_ring_correction_is_traceable():
    """provenance 要写明用了哪个环的修正。"""
    sig = _signals("C1CCCCC1")[0]
    assert "ring6(" in sig.provenance.key


def test_acyclic_gets_no_ring_correction():
    """开链不该被环项碰到——否则整条链都会平移。"""
    for sig in _signals("CCCCCC"):
        assert "ring" not in sig.provenance.key
    # 正己烷的三个信号：14.1 / 22.7 / 31.9
    assert sorted(s.shift for s in _signals("CCCCCC")) == pytest.approx(
        [14.1, 22.7, 31.9], abs=0.5
    )


def test_ring_correction_composes_with_substituent():
    """环修正必须和杂原子增量叠加，不是互相覆盖。

    环己胺 C1 实测 50.1：GP(2/2/1/0)=32.2 + NH2@1 22.5 − ring6 5.1 ≈ 51。
    如果环修正写完把取代基增量冲掉，这里会掉到 27 附近。
    """
    sig = _at("NC1CCCCC1", 1)   # N 是 0 号，环上 ipso 碳是 1 号
    assert sig.shift == pytest.approx(50.1, abs=2.5)
    assert "ring6(" in sig.provenance.key
    assert "NH2@1" in sig.provenance.key


def test_smallest_ring_wins():
    """同时属于多个环时取最小环的修正（小环张力主导）。

    双环[2.2.1]庚烷每个桥头碳同时在一个 5 元环和一个 6 元环里，
    取 5 元环的 −8.8。
    """
    sig = _at("C1CC2CCC1C2", 0)
    assert "ring5(" in sig.provenance.key


# ---------------------------------------------------------------------------
# 羧基重复计数
# ---------------------------------------------------------------------------
def test_carboxyl_oxygen_is_not_counted_twice():
    """C-C(=O)-O 的 O 与 α 羰基同属一个羧基，不能各加一次。

    丙酸 C2 实测 27.8。加两遍会得到 35.9（原先的实现），只加 α 羰基得 29.9。
    """
    sig = _at("CCC(=O)O", 1)
    assert sig.shift == pytest.approx(27.8, abs=2.5)
    assert "OCOR@2" not in sig.provenance.key


def test_carboxyl_beta_carbon_is_not_double_counted():
    """丁酸 C3（β 位）同样不该多吃一条 OCOR@2。"""
    sig = _at("CCCC(=O)O", 2)
    assert "OCOR@2" not in sig.provenance.key


def test_acetate_ester_oxygen_still_counted():
    """乙酸乙酯的 OCH₂ 走的是 OCOR@1（a=45.0），这条不能因为上面的修改被误伤。"""
    sig = _at("CC(=O)OCC", 4)   # OCH₂ 是 4 号
    assert "OCOR@1" in sig.provenance.key


# ---------------------------------------------------------------------------
# 仍未修好的：环烯烯丙位
# ---------------------------------------------------------------------------
def test_ring_allylic_residual_is_known():
    """环己烯烯丙位实测 25.4，当前约 29——环修正对饱和环标定，管不到烯丙位。

    这条用例的作用是**记录残差**，不是宣称正确。
    修好之后应当把 abs 收紧到 1.5 以内。
    """
    sig = _at("C1CC=CCC1", 1)
    assert sig.shift == pytest.approx(25.4, abs=5.0)


# ---------------------------------------------------------------------------
# β 位的芳基/烯基不产生去屏蔽
# ---------------------------------------------------------------------------
def test_beta_aryl_increment_is_zero():
    """``aryl`` / ``vinyl`` 的 β 位增量（b）必须是 0。

    原先填的是 6.0 / 3.0，理由是"β 位该有去屏蔽"——与实测相反：
    乙苯的甲基 GP 已经给出 16.2（实测 15.6），再加 6.0 就变成 22.2，错 6.6 ppm；
    异丙苯的甲基同理（25.6 + 6.0 = 31.6，实测 24.1）。归零后两者分别降到 0.6 / 1.5。
    """
    subst = shift_table("13C")["sp3"]["subst"]
    assert float(subst["aryl"]["b"]) == 0.0
    assert float(subst["vinyl"]["b"]) == 0.0


def test_ethylbenzene_methyl():
    """乙苯甲基实测 15.6。这是 β-芳基增量归零的主要依据。"""
    sig = _at("CCc1ccccc1", 0)
    assert sig.shift == pytest.approx(15.6, abs=1.5)
    assert "aryl@2" not in sig.provenance.key


def test_isopropylbenzene_methyl():
    """异丙苯的两个甲基等价，实测 24.1。"""
    sig = _at("CC(C)c1ccccc1", 0)
    assert sig.shift == pytest.approx(24.1, abs=2.0)


def test_ethylbenzene_benzylic_ch2_still_correct():
    """归零 β 项不能把 α 位一起弄坏：苄位 CH₂ 实测 28.9。"""
    sig = _at("CCc1ccccc1", 1)
    assert sig.shift == pytest.approx(28.9, abs=1.0)
    assert "aryl@1" in sig.provenance.key


# ---------------------------------------------------------------------------
# 支链烷烃本来是对的——别再去"修"它
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "smiles,lit",
    [
        # 异丁烷 / 新戊烷：GP 的 α/β 线性计数在支化点并没有失效。
        # 曾误判为"支链模型失效"，实测误差只有 1.3 / 2.4 ppm；
        # 这条用例防止后来者按直觉去调 beta_extra 而弄坏本来就准的地方。
        ("CC(C)C", [24.3, 24.3, 24.3, 25.3]),
        ("CC(C)(C)C", [27.9, 31.7, 31.7, 31.7, 31.7]),
    ],
)
def test_branched_alkanes_are_already_accurate(smiles: str, lit: list[float]):
    got = sorted(s.shift for s in _signals(smiles) for _ in range(s.n_nuclei))
    assert got == pytest.approx(sorted(lit), abs=2.5)


def test_neopentane_quaternary_is_the_known_outlier():
    """新戊烷的季碳实测 27.9，预测 34.1（错 6.2）。

    这条**故意写成记录残差**而不是断言正确：GP 的 ``beta_extra`` 对 3 个 β 碳
    按 9.4+4.0+4.0 计，季碳上给得太满。同时它也是 diag_sp3 那个
    "排序配对会摊平误差"局限的实例——诊断只报 2.4，真实最差是 6.2。
    """
    got = sorted(s.shift for s in _signals("CC(C)(C)C") for _ in range(s.n_nuclei))
    assert got[-1] == pytest.approx(34.1, abs=0.5)      # 预测
    assert abs(got[-1] - 27.9) > 5.0                     # 与实测的差


# ---------------------------------------------------------------------------
# 含 C=C 的环：单独一张修正表
# ---------------------------------------------------------------------------
def test_cyclohexene_allylic_and_other():
    """环己烯：烯丙位 25.4、其余 CH₂ 22.9。

    饱和环的 ring6(−5.1) 在这两个位置都不适用——它是在环己烷上标定的。
    """
    shifts = sorted(s.shift for s in _sp3_signals("C1CC=CCC1") for _ in range(s.n_nuclei))
    assert shifts[:2] == pytest.approx([22.9, 22.9], abs=0.4)
    assert shifts[2:] == pytest.approx([25.4, 25.4], abs=0.4)


def test_cyclopentene_allylic_and_other():
    """环戊烯：烯丙位 32.9（比环己烯**高**）、其余 23.6。

    这一条是分表的原因：环己烯烯丙位需要 −8.8 的修正，环戊烯需要 −3.8，
    方向相同但量级差一倍以上；若只看环己烯会以为环戊烯也该往下修。
    """
    sigs = _sp3_signals("C1CC=CC1")
    allylic = [s for s in sigs if "allylic" in s.provenance.key][0]
    other = [s for s in sigs if "other" in s.provenance.key][0]
    assert allylic.shift == pytest.approx(32.9, abs=0.4)
    assert other.shift == pytest.approx(23.6, abs=0.4)


def test_saturated_ring_keeps_old_correction():
    """饱和环不许被环烯表碰到。"""
    sig = _signals("C1CCCCC1")[0]
    assert sig.shift == pytest.approx(27.1, abs=0.2)
    assert "ring6(" in sig.provenance.key
    assert "allylic" not in sig.provenance.key


def test_ring_ketone_is_not_treated_as_cycloalkene():
    """环己酮的 C=O 不是环内 C=C，不能触发环烯表。

    氧不在环里，所以遍历"环内相邻原子对"时看不到这条双键——
    这条用例守的就是这个判断，写错会把整个酮的环都平移 4 ppm。
    """
    sigs = _signals("O=C1CCCCC1")
    for sig in sigs:
        if "GP" in sig.provenance.key:
            assert "ring6(" in sig.provenance.key
            assert "allylic" not in sig.provenance.key


def test_cyclohexadiene_out_of_sample():
    """1,3-环己二烯**没有**参与环烯表的标定，可作为弱验证。

    改之前 sp³ 碳预测 29.1（实测约 23.0，差 6.1），改后 25.4（差 2.4）。
    它同时有两个烯丙位且共轭，严格说和环己烯不是一个环境，所以只当支持性证据。
    """
    sig = _sp3_signals("C1C=CC=CC1")[0]
    assert sig.shift == pytest.approx(25.4, abs=0.5)
    assert abs(sig.shift - 23.0) < 3.0
