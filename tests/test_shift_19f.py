"""¹⁹F 环境表：CF₃ 接在 sp² 碳上、CF₃X，以及"每个氟都要能匹配上"。

背景
----
原表对 CF₃ 的两条规则（`Alkyl-CF3` / `Rf-CF3`）都要求邻居是 sp³ 碳 `[CX4]`，
于是 `FC(F)(F)C=C(Cl)Cl` 这类 **CF₃ 接乙烯基/亚胺碳** 的结构两条都不命中，
掉进兜底 −120 ppm（实测 −61.8，差 58 ppm）。库内该环境有 147 条训练记录。

另一件事：`GetSubstructMatches` 默认 `uniquify=True` 会按分子自同构折叠等价匹配，
CF₃ 的三个氟只回一个。这个坑第一次出现在挖掘工具里，把"未命中"重复计了 3 次
（缺口从真实 16% 虚报成 59%）。所以这里专门钉一条"组内每个原子都要能匹配"。
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from nmrflow.core.types import Confidence
from nmrflow.data.loader import shift_table
from nmrflow.mol.envs import SubstituentSet
from nmrflow.mol.equivalence import hetero_groups
from nmrflow.mol.prepare import prepare
from nmrflow.pipeline import RunOptions, run
from nmrflow.predict.shift_hetero import HeteroShiftPredictor


def _f19(smiles: str):
    """取核数最多的那个 ¹⁹F 组。

    像 ``FC(F)=N-C(F)(F)F`` 这种分子同时有 CF₂ 和 CF₃ 两组，
    按下标取第 0 个会拿到 CF₂，容易把"测错了组"误当成"条目没命中"。
    """
    prep = prepare(smiles)
    groups = hetero_groups(prep, "F", 0)
    sigs = HeteroShiftPredictor("19F").predict(prep, groups, "CDCl3")
    return max(sigs, key=lambda s: s.n_nuclei)


# ---------------------------------------------------------------------------
# 新增条目
# ---------------------------------------------------------------------------
def test_cf3_on_sp2_carbon():
    """CF₃ 接乙烯基碳：原来两个条目都要求 sp³ 邻居，结果是兜底 −120。"""
    sig = _f19("FC(F)(F)C=C(Cl)Cl")
    assert sig.provenance.key == "CF3-sp2 (CF3 接乙烯基/亚胺碳)"
    assert sig.shift == pytest.approx(-61.5, abs=0.1)
    assert not sig.provenance.fallback


def test_cf3_on_imine_nitrogen_is_covered():
    """CF₃ 接在亚胺氮（NX2）上：原 `N-CF3` 只写了 `[NX3]`，这类会漏。

    库内该环境 15 条记录、中位 −59.3，与 `N-CF3` 的 −55.0 只差 4 ppm，
    所以把 SMARTS 放宽成 `[NX2,NX3]`，不另立条目。
    """
    assert _f19("FC(F)=N-C(F)(F)F").provenance.key == "N-CF3"


@pytest.mark.parametrize("smiles", ["FC(F)(F)Cl", "FC(F)(F)Br"])
def test_cf3_halide(smiles: str):
    sig = _f19(smiles)
    assert sig.provenance.key == "CF3-halogen (CF3Cl / CF3Br)"
    assert sig.confidence is Confidence.LOW      # 记录少、离散大


@pytest.mark.parametrize(
    "smiles,expected",
    [
        # 新增条目不能抢走原有条目
        ("FC(F)(F)c1ccccc1", "Ar-CF3"),
        ("FC(F)(F)CC", "Alkyl-CF3"),
        ("FC(F)(F)S(=O)(=O)c1ccccc1", "CF3-SO2- (三氟甲磺酰 / triflyl)"),
        ("FC(F)(F)Oc1ccccc1", "ArO-CF3 (三氟甲氧基)"),
        ("Fc1ccccc1", "Ar-F"),
    ],
)
def test_new_entries_do_not_hijack_old_ones(smiles: str, expected: str):
    assert _f19(smiles).provenance.key == expected


# ---------------------------------------------------------------------------
# 对称原子的匹配（uniquify 回归护栏）
# ---------------------------------------------------------------------------
def test_every_fluorine_in_cf3_resolves():
    """CF₃ 的三个氟必须**每一个**都能匹配到条目。

    RDKit 默认 `uniquify=True` 只回一个等价匹配，于是 ``match[0] == idx``
    只对其中一个原子成立。小组代表原子恰好不是它时，整组就会掉进兜底——
    这条用例守的是 ``SubstituentSet.match_at`` 里的 ``uniquify=False``。
    """
    table = shift_table("19F")
    subs = SubstituentSet({e["name"]: e for e in table["env"]})
    for smiles in ["FC(F)(F)c1ccccc1", "FC(F)(F)C=C(Cl)Cl", "CCOC(F)(F)F"]:
        mol = Chem.MolFromSmiles(smiles)
        fluorines = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "F"]
        assert len(fluorines) >= 3
        for idx in fluorines:
            assert subs.match_at(mol, idx) is not None, (
                f"{smiles} 的 F{idx} 没有匹配到任何环境"
            )


# ---------------------------------------------------------------------------
# 兜底行为
# ---------------------------------------------------------------------------
def test_uncovered_environment_falls_back_low_and_wide():
    """没覆盖到的环境应该老实兜底：low 置信度 + 宽峰，而不是假装精确。

    XeF₄ / NF₃ 这类无机氟化物不在有机 NMR 表的范围内，
    与其硬编一个条目，不如让它带着"未知"的标记出现。
    """
    sig = _f19("F[Xe](F)(F)F")
    assert sig.provenance.fallback
    assert sig.confidence is Confidence.LOW
    assert sig.width_ppm >= 20


def test_fallback_is_not_silently_zero_width():
    """兜底必须给宽峰；窄峰 + low 是自相矛盾的组合。"""
    result = run("F[Xe](F)(F)F", RunOptions(tier_b=False, with_ms=False))
    sig = result.signals["19F"][0]
    assert sig.provenance.fallback
    assert sig.width_ppm > 10


# ---------------------------------------------------------------------------
# 卤代乙烷（XCF2-CFY）：δ 由两个碳的卤素共同决定
# ---------------------------------------------------------------------------
# 实测值来自 2026-11-12 的 376 MHz 谱（CDCl3），文献值来自
# JFC 1995, 74, 261–267（p.262 正文 / p.266 实验）。两路独立，差值 < 0.5 ppm。
#
# 修复前的表现（本组用例就是为它写的）：
#   I  的 CFClI → 兜底 −120.00（错 48.1 ppm）
#   II 的 CF2I  → HOSE 2 球 −65.50（错 10.8 ppm，2 层球把 CF2Cl/CF2I 混在一起）

CHLOROFLUOROETHANES = [
    # (SMILES, 按位移降序的 (δ, 核数), 说明)
    ("IC(Cl)(F)C(F)(Cl)F", [(-63.40, 2), (-71.90, 1)], "I: CF2Cl / CFClI"),
    ("ClC(Cl)(F)C(F)(I)F", [(-54.80, 2), (-67.30, 1)], "II: CF2I / CFCl2"),
]


@pytest.mark.parametrize("smiles,expected,label", CHLOROFLUOROETHANES)
def test_chlorofluoroethane_shifts(smiles, expected, label):
    """两个异构体都要落在实测值 ±0.6 ppm 内，且不得掉进兜底。"""
    sigs = sorted(run(smiles).signals["19F"], key=lambda s: -s.shift)
    assert len(sigs) == len(expected), f"{label}: 信号数 {len(sigs)}"
    for sig, (want, n) in zip(sigs, expected):
        assert sig.n_nuclei == n, f"{label}: 核数 {sig.n_nuclei} != {n}"
        assert not sig.provenance.fallback, f"{label}: 掉进兜底了"
        assert sig.shift == pytest.approx(want, abs=0.6), f"{label}: {sig.shift:.2f}"


def test_cfcli_no_longer_falls_back_to_minus_120():
    """只带一个 F、连 Cl/I 的 sp³ 碳原先是 −120 兜底（实测 −71.9，差 48 ppm）。

    守的是 `R-CF(X)2` 通用条目：它把"整类缺失"变成 ±12 ppm 的估计，
    而不是一个化学上不可能的 −120。
    """
    sigs = run("IC(Cl)(F)C(F)(Cl)F").signals["19F"]
    cfcli = min(sigs, key=lambda s: s.shift)  # 最靠高场的那个
    assert not cfcli.provenance.fallback
    assert cfcli.shift > -90


def test_pinned_entry_is_not_overridden_by_shallow_hose():
    """文献钉死的环境不能被 2 层球 HOSE 覆盖。

    ¹⁹F 的 HOSE 索引在含卤乙烷上最多只有 2 层球、σ≈8 ppm——那种均值把
    CF2Cl / CF2I / CF2Br 混在一起。表里标了 ``protect = true`` 的条目跳过覆盖，
    这里用 tier == "A" 把行为钉住，否则以后一改 HoseConfig 就会悄悄失效。
    """
    from nmrflow.predict.shift_hose import _pinned_keys

    assert "CF2I-CFCl2 (CF2I 端, 1,2-二氯-1-碘三氟乙烷异构体)" in _pinned_keys("19F")
    sigs = run("ClC(Cl)(F)C(F)(I)F").signals["19F"]
    assert all(s.provenance.tier == "A" for s in sigs)


def test_byproduct_reference_values():
    """X = Cl 的副产物 ClCF2CFCl2：−68.1 / −72.4（JFC 1995 §10.2）。"""
    shifts = sorted(s.shift for s in run("ClC(F)(F)C(F)(Cl)Cl").signals["19F"])
    assert shifts == pytest.approx([-72.4, -68.1], abs=0.3)


# ---------------------------------------------------------------------------
# 非对映异位：CF₂ 紧邻手性中心时不该静默合成一个 2F 峰
# ---------------------------------------------------------------------------
def test_potential_stereocenter_seen_without_explicit_chirality():
    """SMILES 没写 @ 不代表没有手性中心。

    CF₂Cl–CFClI 的 CFClI 碳连 F/Cl/I/CF₂Cl 四个不同基团，样品里就是手性中心；
    ``has_stereochemistry`` 只认已标注的 @，判断两个 F 是否非对映异位必须用
    ``potential_stereocenters``。异构体 II 的两个 Cl 相同，不是中心。
    """
    from nmrflow.mol.equivalence import potential_stereocenters

    mol_i = Chem.MolFromSmiles("IC(Cl)(F)C(F)(Cl)F")
    assert potential_stereocenters(mol_i) == {1}
    assert potential_stereocenters(Chem.MolFromSmiles("ClC(Cl)(F)C(F)(I)F")) == set()


def test_diastereotopic_cf2_raises_warning():
    """有手性中心的那个异构体要提示 AB 分裂；另一个不该提示。

    这条不是"功能"测试，是"别撒谎"测试：模型给不出 AB 分裂，那就必须说清楚，
    否则一条 2F 单峰会被当成预测正确（实测是 −61.82 / −64.76 两条线）。
    """
    warn_i = " ".join(run("IC(Cl)(F)C(F)(Cl)F").report.warnings)
    assert "非对映异位" in warn_i and "AB" in warn_i
    warn_ii = " ".join(run("ClC(Cl)(F)C(F)(I)F").report.warnings)
    assert "非对映异位" not in warn_ii

def test_cf3_so2f_is_not_overridden_by_shallow_hose():
    """CF3-SO2- 表值 -78 不能被 2 层球 HOSE 拉到 -58。"""
    sigs = run("FC(F)(F)S(=O)(=O)F").signals["19F"]
    cf3 = [s for s in sigs if s.n_nuclei == 3]
    assert len(cf3) == 1
    assert abs(cf3[0].shift + 78.0) <= 1.0
    assert cf3[0].provenance.tier == "A"