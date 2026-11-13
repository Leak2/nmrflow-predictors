"""HOSE 索引查询层的测试（用合成索引，不依赖真实数据）。"""

from __future__ import annotations

import sqlite3

import pytest
from rdkit import Chem

from nmrflow.predict.hose import hose_codes, truncate
from nmrflow.predict.hose_index import HoseHit, HoseIndex

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
    n = len(shifts)
    return (
        nucleus,
        code,
        code.count(";") if level is None else level,
        n,
        sum(shifts),
        sum(s * s for s in shifts),
    )


CODE3 = "C-4;(-C);(-C);(-C)"
CODE4 = "C-4;(-C);(-C);(-C);(-C)"


def test_exact_level4(tmp_path):
    path = _make_index(tmp_path, [_row(CODE4, "13C", [10.0, 12.0, 14.0])])
    with HoseIndex(path) as idx:
        hit = idx.lookup(CODE4, "13C")
    assert hit is not None
    assert hit.level == 4
    assert hit.n == 3
    assert hit.mean == pytest.approx(12.0)
    assert hit.std == pytest.approx(1.632993, rel=1e-5)


def test_prefix_query_merges_longer_codes(tmp_path):
    """4 层球的码在 3 层查询里必须被算进同一环境。"""
    path = _make_index(
        tmp_path,
        [_row(CODE3, "13C", [10.0]), _row(CODE4, "13C", [20.0])],
    )
    with HoseIndex(path) as idx:
        hit = idx.lookup(CODE4, "13C", max_level=3)
    assert hit is not None
    assert hit.level == 3
    assert hit.n == 2
    assert hit.mean == pytest.approx(15.0)


def test_degrade_when_neighbors_below_threshold(tmp_path):
    """最深层邻居数不够时退一层，并如实报告退化后的球层。"""
    path = _make_index(
        tmp_path,
        [_row(CODE4, "13C", [10.0]), _row(CODE3, "13C", [10.0, 30.0])],
    )
    with HoseIndex(path) as idx:
        strict = idx.lookup(CODE4, "13C", min_neighbors=1)
        loose = idx.lookup(CODE4, "13C", min_neighbors=2)
    assert strict is not None and strict.level == 4 and strict.n == 1
    # 退到 3 层后同时命中 CODE3 与 CODE4 两行
    assert loose is not None and loose.level == 3 and loose.n == 3


def test_sibling_code_not_matched(tmp_path):
    """前缀范围扫描不能把不同分支的码算进来。

    ``sibling`` 在第 3 层球才分叉，所以 3 层及以上的查询必须查不到它；
    退到 2 层时两者环境相同，命中是正确的。
    """
    sibling = "C-4;(-C);(-C);(-X);(-C)"
    path = _make_index(tmp_path, [_row(sibling, "13C", [99.0])])
    with HoseIndex(path) as idx:
        assert idx.lookup(CODE4, "13C", min_level=3) is None
        shallow = idx.lookup(CODE4, "13C", max_level=2)
    assert shallow is not None and shallow.level == 2


def test_nucleus_isolation(tmp_path):
    path = _make_index(tmp_path, [_row(CODE4, "13C", [10.0])])
    with HoseIndex(path) as idx:
        assert idx.lookup(CODE4, "1H") is None


def test_fallback_returns_best_support_then_deepest(tmp_path):
    """邻居数始终不达标时，返回支持数最多的命中；持平时保留更深的球层。"""
    path = _make_index(
        tmp_path,
        [_row(CODE4, "13C", [10.0]), _row(CODE3, "13C", [10.0, 30.0])],
    )
    with HoseIndex(path) as idx:
        # min_level=4 时只允许看第 4 层，支持数 1
        deep = idx.lookup(CODE4, "13C", min_level=4, min_neighbors=5)
        # 允许退到 3 层后，支持数升到 3，应返回 3 层
        shallow = idx.lookup(CODE4, "13C", min_level=3, min_neighbors=5)
    assert deep is not None and deep.level == 4 and deep.n == 1
    assert shallow is not None and shallow.level == 3 and shallow.n == 3


def test_predict_only_returns_matching_atoms(tmp_path):
    """乙醇：¹H 只应落在带氢的碳与氧上，¹³C 只落在碳上。"""
    mol = Chem.MolFromSmiles("CCO")
    codes = hose_codes(mol, 4)
    rows = [_row(c, "1H", [1.2]) for c in codes.values()]
    rows += [_row(c, "13C", [50.0]) for c in codes.values()]
    path = _make_index(tmp_path, rows)

    with HoseIndex(path) as idx:
        h_hits = idx.predict(mol, "1H")
        c_hits = idx.predict(mol, "13C")

    assert set(h_hits) == {0, 1, 2}
    assert set(c_hits) == {0, 1}
    assert all(isinstance(h, HoseHit) for h in h_hits.values())
    assert h_hits[0].mean == pytest.approx(1.2)


def test_missing_index_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        HoseIndex(tmp_path / "nope.sqlite")


def test_truncate_matches_code_at_that_level():
    """降级依赖的前缀性质：截断 k 层必须等于直接按 k 层生成的码。"""
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    full = hose_codes(mol, 5)
    shallow = hose_codes(mol, 3)
    for i, code in full.items():
        assert truncate(code, 3) == shallow[i]
