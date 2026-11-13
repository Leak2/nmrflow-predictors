"""Tier B（HOSE 覆盖层）的测试：用合成索引，不依赖真实数据库。"""

from __future__ import annotations

import copy
import sqlite3

import pytest
from rdkit import Chem

from nmrflow.core.types import Confidence, Provenance, Signal
from nmrflow.mol.equivalence import SpinGroup, carbon_groups
from nmrflow.mol.prepare import prepare
from nmrflow.pipeline import RunOptions, run
from nmrflow.predict.hose import hose_codes
from nmrflow.predict.hose_index import HoseIndex
from nmrflow.predict.shift_hose import (
    HoseConfig,
    HoseShiftPredictor,
    confidence_for,
    default_config,
    overlay_tier_b,
)

SCHEMA = """
CREATE TABLE hose(
    nucleus  TEXT    NOT NULL,
    code     TEXT    NOT NULL,
    level    INTEGER NOT NULL,
    n        INTEGER NOT NULL,
    sum_shift REAL   NOT NULL,
    sum_sq    REAL   NOT NULL,
    PRIMARY KEY (nucleus, code)
) WITHOUT ROWID
"""


def _make_index(tmp_path, rows):
    path = tmp_path / "idx.sqlite"
    con = sqlite3.connect(path)
    con.execute(SCHEMA)
    con.executemany("INSERT INTO hose VALUES (?,?,?,?,?,?)", rows)
    con.commit()
    con.close()
    return path


def _row(code, nucleus, shifts, level=None):
    return (
        nucleus,
        code,
        code.count(";") if level is None else level,
        len(shifts),
        sum(shifts),
        sum(s * s for s in shifts),
    )


def _base(group_id, shift, conf=Confidence.MEDIUM, tier="A"):
    return Signal(
        nucleus="13C",
        group_id=group_id,
        atom_indices=[group_id],
        heavy_indices=[group_id],
        shift=shift,
        confidence=conf,
        provenance=Provenance(tier=tier, table="tier_a"),
    )


def _cand(shift, level, nucleus="13C", fallback=False, note="", gid=0):
    sig = Signal(
        nucleus=nucleus,
        group_id=gid,
        shift=shift,
        provenance=Provenance(tier="B", table="nmrshiftdb2.hose"),
    )
    if fallback:
        sig.provenance.fallback = True
        sig.provenance.note = note
    else:
        sig.provenance.key = f"{level}球/7"
        sig.provenance.note = f"HOSE {level} 层球, 7 条同环境记录"
        sig.confidence = confidence_for(nucleus, level)
    return sig


def _group(gid, heavy=0, diastereotopic=-1):
    return SpinGroup(gid, "13C", [heavy], [heavy], 1, diastereotopic_pair=diastereotopic)


# ---------------------------------------------------------------------------
# 置信度分档
# ---------------------------------------------------------------------------
def test_confidence_table():
    assert confidence_for("1H", 4) is Confidence.HIGH
    assert confidence_for("1H", 2) is Confidence.MEDIUM
    assert confidence_for("1H", 1) is Confidence.LOW
    assert confidence_for("13C", 4) is Confidence.MEDIUM      # 1.21 ppm 不宜称 high
    assert confidence_for("13C", 2) is Confidence.LOW
    assert confidence_for("19F", 4) is Confidence.MEDIUM      # 未验证，封顶 medium


# ---------------------------------------------------------------------------
# 覆盖层
# ---------------------------------------------------------------------------
def test_level_below_threshold_is_skipped():
    base = [_base(1, 50.0)]
    cfg = HoseConfig(min_level=3)
    stats = overlay_tier_b(base, [_cand(60.0, 2)], [_group(1)], "13C", cfg)
    assert base[0].shift == 50.0
    assert stats.overridden == 0
    assert stats.skipped_level == 1


def test_override_records_level_and_neighbors():
    base = [_base(1, 50.0)]
    stats = overlay_tier_b(
        base, [_cand(52.5, 4)], [_group(1)], "13C", HoseConfig(min_level=1)
    )
    assert base[0].shift == pytest.approx(52.5)
    assert base[0].provenance.tier == "B"
    assert base[0].provenance.table == "nmrshiftdb2.hose"
    assert base[0].provenance.key == "4球/7"
    assert stats.overridden == 1 and stats.levels == {4: 1}


def test_labile_proton_is_not_overridden():
    base = [_base(1, 2.2, conf=Confidence.LOW)]
    cand = _cand(5.0, 4, nucleus="1H", fallback=True, note="活泼氢：不走 HOSE")
    stats = overlay_tier_b(base, [cand], [_group(1)], "1H", HoseConfig(min_level=1))
    assert base[0].shift == 2.2
    assert stats.skipped_labile == 1


def test_confidence_comes_from_level():
    """置信度由命中球层决定，覆盖后按候选值重算。

    球层是唯一经过留出法标定的量；取两者中较低者会让 Tier A 的 fallback(LOW)
    挡住一条 4 层球的高质量命中。
    """
    base = [_base(1, 200.0, conf=Confidence.HIGH)]
    overlay_tier_b(
        base, [_cand(203.0, 4)], [_group(1)], "13C", HoseConfig(min_level=1)
    )
    assert base[0].shift == pytest.approx(203.0)
    assert base[0].confidence is Confidence.MEDIUM

    # 反过来：Tier A 本来是 low，Tier B 拿到 4 层球，可以升到 medium
    base = [_base(1, 50.0, conf=Confidence.LOW)]
    overlay_tier_b(
        base, [_cand(52.0, 4)], [_group(1)], "13C", HoseConfig(min_level=1)
    )
    assert base[0].confidence is Confidence.MEDIUM


def test_wide_sigma_only_downgrades_confidence():
    """离散度只降置信度，不否决覆盖。

    实测（tools/eval_tier_ab.py）否决反而更差：13C MAE 3.47→4.45、P90 8.70→12.20，
    因为被拒的原子退回 Tier A，而 Tier A 的 13C P90 是 19.20 ppm。
    真实案例：三苯基锍的芳香碳在 1 层球下与普通芳香 CH 同码，σ=7.95。
    """
    base = [_base(1, 133.5, conf=Confidence.MEDIUM)]
    cand = _cand(123.7, 4)
    cand.provenance.sigma = 7.95
    stats = overlay_tier_b(base, [cand], [_group(1)], "13C", HoseConfig(min_level=1))
    assert base[0].shift == pytest.approx(123.7)      # 仍然覆盖
    assert base[0].confidence is Confidence.LOW        # 但打了降级
    assert "离散" in base[0].provenance.note
    assert stats.wide_sigma == 1

    # σ 正常时不动置信度
    base = [_base(1, 133.5, conf=Confidence.MEDIUM)]
    cand = _cand(132.0, 4)
    cand.provenance.sigma = 0.5
    overlay_tier_b(base, [cand], [_group(1)], "13C", HoseConfig(min_level=1))
    assert base[0].confidence is Confidence.MEDIUM


def test_sigma_warning_can_be_disabled():
    """``sigma_warn_ratio=None`` 关掉警告，且**不动位移**——消融表末行靠它证明这一点。

    这里刻意让候选值与 Tier A 只差 0.3 ppm，避开"与经验规则相差过大"那条降级，
    这样置信度的变化只可能来自 σ 警告本身。
    """
    cand = _cand(133.2, 4)
    cand.provenance.sigma = 7.95
    for ratio, expect_downgrade in ((2.5, True), (None, False)):
        base = [_base(1, 133.5, conf=Confidence.MEDIUM)]
        cfg = HoseConfig(min_level=1, sigma_warn_ratio=ratio)
        stats = overlay_tier_b(base, [copy.deepcopy(cand)], [_group(1)], "13C", cfg)
        assert base[0].shift == pytest.approx(133.2)
        assert (base[0].confidence is Confidence.LOW) is expect_downgrade
        assert (stats.wide_sigma == 1) is expect_downgrade
        assert ("离散" in base[0].provenance.note) is expect_downgrade


def test_bare_config_matches_default_for_main_nuclei():
    """``HoseConfig()`` 的 min_level 默认值必须与 ``default_config`` 对 ¹H/¹³C 一致。

    绕过 ``default_config`` 直接构造配置的调用方（工具脚本、外部使用者）如果拿到
    更严的门槛，就会悄悄丢掉一部分修正——而实测每一层 B 都优于同层 A。
    """
    for nucleus in ("1H", "13C"):
        assert HoseConfig().min_level == default_config(nucleus).min_level


def test_large_disagreement_notes_but_keeps_confidence():
    """与 Tier A 分歧只写提示，不降 HOSE 的置信度。

    旧逻辑会把“Tier A 错得离谱、HOSE 正确”的情况（咖啡因羰基/脒碳、
    胆固醇多个环碳）标成 low —— 预测越对，置信度越低。分歧只是“两条证据
    不一致”，不等于“HOSE 一定错”，所以只写 note + warning 提示复核。
    """
    base = [_base(1, 150.0, conf=Confidence.MEDIUM)]
    stats = overlay_tier_b(
        base, [_cand(170.0, 4)], [_group(1)], "13C", HoseConfig(min_level=1)
    )
    assert stats.disagreed == 1
    assert base[0].confidence is Confidence.MEDIUM      # 不再降级
    assert "相差" in base[0].provenance.note


def test_31p_protection_can_be_disabled():
    """关掉 protect_curated 后，浅层命中确实会覆盖 Tier A 的锚定值。

    守的是“机制没坏”：³¹P 默认保护 curated，但保护是配置项而非写死。
    """
    base = [_base(1, 0.0, conf=Confidence.HIGH)]
    stats = overlay_tier_b(
        base, [_cand(14.3, 2)], [_group(1)], "31P",
        HoseConfig(min_level=1, protect_curated=False),
    )
    assert base[0].provenance.tier == "B"
    assert base[0].shift == pytest.approx(14.3, abs=0.01)
    assert stats.overridden == 1


def test_diastereotopic_pair_keeps_split():
    """两个非对映异位氢拿到同一个 HOSE 值，但必须保住 Tier A 的相对位移差。"""
    groups = [_group(1, heavy=0, diastereotopic=1), _group(2, heavy=0, diastereotopic=0)]
    base = [_base(1, 3.50), _base(2, 3.70)]
    cands = [
        _cand(4.00, 4, nucleus="1H", gid=1),
        _cand(4.00, 4, nucleus="1H", gid=2),
    ]
    overlay_tier_b(base, cands, groups, "1H", HoseConfig(min_level=1))
    assert base[0].shift == pytest.approx(3.90)
    assert base[1].shift == pytest.approx(4.10)
    assert "非对映异位" in base[0].provenance.note


def test_stats_note_is_readable():
    base = [_base(1, 50.0), _base(2, 60.0)]
    stats = overlay_tier_b(
        base,
        [_cand(52.0, 4), _cand(58.0, 3)],
        [_group(1), _group(2)],
        "13C",
        HoseConfig(min_level=1),
    )
    note = stats.as_note("13C")
    assert "2/2" in note and "4 层" in note and "3 层" in note


# ---------------------------------------------------------------------------
# 预测器本体
# ---------------------------------------------------------------------------
def test_proton_lookup_uses_heavy_atom(tmp_path):
    """¹H 的原子序号在库里指向带氢的重原子，查表必须用重原子的码。"""
    mol = Chem.MolFromSmiles("CCO")
    codes = hose_codes(mol, 4)
    rows = [_row(c, "1H", [1.2]) for c in codes.values()]
    path = _make_index(tmp_path, rows)

    from nmrflow.mol.equivalence import hydrogen_groups
    from nmrflow.mol.prepare import prepare

    prep = prepare("CCO")
    groups = hydrogen_groups(prep)
    with HoseIndex(path) as idx:
        cands = HoseShiftPredictor("1H", idx, HoseConfig(min_level=1)).predict(
            prep, groups, "CDCl3"
        )
    assert len(cands) == len(groups)
    hits = [c for c in cands if not c.provenance.fallback]
    assert hits and all(c.shift == pytest.approx(1.2) for c in hits)


def test_labile_hydroxyl_falls_back(tmp_path):
    """乙醇的 OH 必须被识别为活泼氢，不参与 HOSE 覆盖。"""
    mol = Chem.MolFromSmiles("CCO")
    codes = hose_codes(mol, 4)
    rows = [_row(c, "1H", [1.2]) for c in codes.values()]
    path = _make_index(tmp_path, rows)

    from nmrflow.mol.equivalence import hydrogen_groups
    from nmrflow.mol.prepare import prepare

    prep = prepare("CCO")
    groups = hydrogen_groups(prep)
    with HoseIndex(path) as idx:
        cands = HoseShiftPredictor("1H", idx, HoseConfig(min_level=1)).predict(
            prep, groups, "CDCl3"
        )
    fallback = [c for c in cands if c.provenance.fallback]
    assert fallback, "OH 应被判为活泼氢"
    assert all("活泼氢" in c.provenance.note for c in fallback)


# ---------------------------------------------------------------------------
# 端到端接线
# ---------------------------------------------------------------------------
def test_pipeline_applies_tier_b(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    codes = hose_codes(mol, 4)
    uniq = sorted(set(codes.values()))
    rows = [_row(uniq[0], "13C", [58.0, 60.0]), _row(uniq[1], "13C", [18.0])]
    path = _make_index(tmp_path, rows)

    result = run(
        "CCO",
        RunOptions(tier_b=True, hose_index=str(path), with_ms=False, with_isotopes=False),
    )
    shifts = {round(s.shift) for s in result.signals["13C"]}
    assert shifts == {59, 18}
    assert all(s.provenance.tier == "B" for s in result.signals["13C"])
    assert any("HOSE" in w for w in result.report.warnings)


def test_pipeline_without_tier_b_keeps_tier_a(tmp_path):
    mol = Chem.MolFromSmiles("CCO")
    codes = hose_codes(mol, 4)
    rows = [_row(c, "13C", [58.0]) for c in codes.values()]
    path = _make_index(tmp_path, rows)

    result = run(
        "CCO",
        RunOptions(tier_b=False, hose_index=str(path), with_ms=False, with_isotopes=False),
    )
    assert all(s.provenance.tier == "A" for s in result.signals["13C"])


def test_pipeline_ignores_missing_index():
    """索引缺失时必须静默跑通，而不是抛异常。"""
    result = run(
        "CCO",
        RunOptions(tier_b=True, hose_index="不存在的路径.sqlite", with_ms=False),
    )
    assert result.signals["1H"]


def test_single_record_hit_reports_no_sigma(tmp_path):
    """库里只有一条记录时，σ 必须是 None 而不是 0.0。

    std 由 ``sum_sq/n − mean²`` 算出，单点上精确等于 0——但那是"无从谈起离散"，
    不是"库内记录高度一致"。这个值会写进 JSON，报 0.0 会让下游误读。
    """
    prep = prepare("CCO")
    groups = carbon_groups(prep, 0)
    code = sorted(hose_codes(Chem.MolFromSmiles("CCO"), 4).values())[0]

    # _make_index 用固定文件名，同一个目录建两次会撞表，所以各用一个子目录
    (tmp_path / "one").mkdir()
    (tmp_path / "two").mkdir()
    one = _make_index(tmp_path / "one", [_row(code, "13C", [58.0])])
    with HoseIndex(one) as idx:
        cand = HoseShiftPredictor("13C", idx, HoseConfig(min_level=1)).predict(
            prep, groups, "CDCl3"
        )
    hit = [c for c in cand if not c.provenance.fallback]
    assert hit, "单条记录也应当命中"
    assert hit[0].provenance.sigma is None

    two = _make_index(tmp_path / "two", [_row(code, "13C", [58.0, 60.0])])
    with HoseIndex(two) as idx:
        cand2 = HoseShiftPredictor("13C", idx, HoseConfig(min_level=1)).predict(
            prep, groups, "CDCl3"
        )
    hit2 = [c for c in cand2 if not c.provenance.fallback]
    assert hit2 and hit2[0].provenance.sigma is not None
