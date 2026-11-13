"""¹⁹F / ³¹P 走 HOSE：门槛、置信度、以及"浅球层不要乱覆盖"。

背景
----
¹⁹F/³¹P 的 HOSE 支持原本默认关闭，理由是"没做留出法验证"。补完验证后结论是
**该开**，但两个核的门槛不能一样，而且原因不是直觉上的：

    Tier A 兜底子集上的 MAE（兜底常数本身 ³¹P 61.31 / ¹⁹F 48.26）
        ³¹P   L1 22.01  L2 10.87  L3 11.76  L4 2.99    ← 每层都强于兜底
        ¹⁹F   L1 96.06  L2 20.26  L3  6.03  L4 3.20    ← L1 比兜底还差

所以 ³¹P 取 min_level=1、¹⁹F 取 2。这个文件把这两个决定钉住。
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from nmrflow.core.types import Confidence
from nmrflow.pipeline import RunOptions, run
from nmrflow.predict.hose import hose_codes
from nmrflow.predict.shift_hose import (
    DEFAULT_MIN_LEVEL,
    confidence_for,
    default_config,
    measured_mae,
)
from tests.test_shift_hose import _make_index, _row


# ---------------------------------------------------------------------------
# 配置与标定表
# ---------------------------------------------------------------------------
def test_hetero_tier_b_is_on_by_default():
    """验证做完了，就不该再默认关着。"""
    assert RunOptions().tier_b_hetero is True


def test_min_level_differs_between_heteronuclei():
    """¹⁹F 的 1 层球比兜底常数还差，必须从 2 层起；³¹P 的 1 层球有用。"""
    assert DEFAULT_MIN_LEVEL["19F"] == 2
    assert DEFAULT_MIN_LEVEL["31P"] == 1
    assert default_config("19F").min_level == 2
    assert default_config("31P").min_level == 1


def test_hetero_mae_table_is_calibrated():
    """实测 MAE 表是置信度与 σ 警告的唯一依据，钉住关键值防止被随手改掉。"""
    assert measured_mae("19F", 4) == pytest.approx(2.21)
    assert measured_mae("31P", 4) == pytest.approx(1.97)
    assert measured_mae("19F", 1) == pytest.approx(52.31)
    assert measured_mae("31P", 1) == pytest.approx(18.85)


def test_hetero_confidence_never_high():
    """这两个核尾部很重（P90 21 / 27 ppm），封顶 medium。"""
    for nucleus in ("19F", "31P"):
        assert confidence_for(nucleus, 4) is Confidence.MEDIUM
        assert confidence_for(nucleus, 3) is Confidence.MEDIUM
        assert confidence_for(nucleus, 2) is Confidence.LOW
        assert confidence_for(nucleus, 1) is Confidence.LOW


# ---------------------------------------------------------------------------
# 端到端：³¹P 覆盖
# ---------------------------------------------------------------------------
def _code_of(smiles: str, atom_idx: int) -> str:
    codes = hose_codes(Chem.MolFromSmiles(smiles), 4)
    return codes[atom_idx]


def test_shallow_31p_hit_yields_to_curated(tmp_path):
    """³¹P 默认 ``protect_curated=True``：浅层命中不再覆盖高置信度表项。

    磷酸三甲酯的 Tier A 表项 ``phosphate triester`` 给 0.0（实测 ≈2.7），
    而 2 层球 HOSE 会给 14.3。0.2.3 起默认保护 curated——这不是“HOSE 更准”
    的问题，是“用一条远离本环境的记录去改写一个已锚定的表项”不划算。
    整体代价实测很小：³¹P MAE 19.45 → 19.49，P90 不变。
    """
    smi = "COP(=O)(OC)OC"
    code = _code_of(smi, 2)              # C0-O1-P2，P 是 2 号
    path = _make_index(tmp_path, [_row(code, "31P", [2.0, 3.0])])

    result = run(smi, RunOptions(hose_index=str(path), with_ms=False,
                                 with_isotopes=False))
    p_sig = result.signals["31P"][0]
    assert p_sig.provenance.tier == "A"          # 保留 Tier A
    assert p_sig.shift == pytest.approx(0.0, abs=0.01)
    assert p_sig.provenance.key.startswith("phosphate triester")


def test_deep_31p_hit_gets_medium(tmp_path):
    """4 层球（实测 MAE 1.97 ppm）才配 medium。三苯基膦的 P 能到 4 层。"""
    smi = "c1ccc(P(c2ccccc2)c3ccccc3)cc1"
    code = _code_of(smi, 4)              # P 是 4 号原子
    path = _make_index(tmp_path, [_row(code, "31P", [-5.0, -6.0])])

    result = run(smi, RunOptions(hose_index=str(path), with_ms=False,
                                 with_isotopes=False))
    p_sig = result.signals["31P"][0]
    assert p_sig.provenance.tier == "B"
    assert p_sig.shift == pytest.approx(-5.5, abs=0.01)
    assert p_sig.confidence is Confidence.MEDIUM


def test_pipeline_leaves_31p_alone_when_disabled(tmp_path):
    smi = "COP(=O)(OC)OC"
    code = _code_of(smi, 2)
    path = _make_index(tmp_path, [_row(code, "31P", [2.0, 3.0])])

    result = run(smi, RunOptions(hose_index=str(path), tier_b_hetero=False,
                                 with_ms=False, with_isotopes=False))
    assert result.signals["31P"][0].provenance.tier == "A"


# ---------------------------------------------------------------------------
# 端到端：¹⁹F 的浅球层必须被拒
# ---------------------------------------------------------------------------
def test_shallow_19f_hit_is_rejected(tmp_path):
    """1 层球的 ¹⁹F 命中（实测 MAE 52 ppm）不许覆盖 Tier A。

    这是 ¹⁹F 门槛取 2 的直接原因：宁可留着 Tier A 的环境值，也不要一个
    比兜底常数还差的浅层命中。
    """
    smi = "Fc1ccccc1"
    code = _code_of(smi, 0)              # F 是 0 号原子
    path = _make_index(tmp_path, [_row(code, "19F", [-114.0, -115.0], level=1)])

    result = run(smi, RunOptions(hose_index=str(path), with_ms=False,
                                 with_isotopes=False))
    assert result.signals["19F"][0].provenance.tier == "A"


def test_deep_19f_hit_is_accepted(tmp_path):
    """同样的环境、同样的值，只要球层够深就应该覆盖。"""
    smi = "Fc1ccccc1"
    code = _code_of(smi, 0)
    path = _make_index(tmp_path, [_row(code, "19F", [-114.0, -115.0], level=3)])

    result = run(smi, RunOptions(hose_index=str(path), with_ms=False,
                                 with_isotopes=False))
    f_sig = result.signals["19F"][0]
    assert f_sig.provenance.tier == "B"
    assert f_sig.shift == pytest.approx(-114.5, abs=0.01)
    assert f_sig.confidence is Confidence.MEDIUM
