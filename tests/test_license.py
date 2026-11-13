"""许可与来源声明。

这几条不是形式主义：内置的 HOSE 索引是 NMRShiftDB2 的**衍生数据**，
它的许可（nmrshiftdb2 Database License）要求四件事，代码里丢任何一件都算不合规：

- 4.2/4.3 随数据分发许可证原文，产出物附带来源声明；
- 4.5 依赖它的软件必须用 OSI 认可的开源许可 —— 本项目取 BSD-3-Clause。

以前索引只在开发期本机用，这些都没落地；索引随包发布之后就必须钉住。
"""

from __future__ import annotations

import json
import tomllib
from pathlib import Path

import pytest

from nmrflow.pipeline import RunOptions, dump_json, run

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "app" / "nmrflow" / "data"

#: 许可 4.3 要求逐字附带的声明
EXPECTED = (
    "Contains information from nmrshiftdb2 (www.nmrshiftdb.org),"
)
EXPECTED_TAIL = "which is available here under the nmrshiftdb2 Database License."


# ---------------------------------------------------------------------------
# 项目自身的许可
# ---------------------------------------------------------------------------
def test_license_file_is_bsd_3_clause():
    text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert text.startswith("BSD 3-Clause License")
    # 三条条款都在（少了第 3 条就成了 BSD-2-Clause，与声明不符）
    assert "1. Redistributions of source code" in text
    assert "2. Redistributions in binary form" in text
    assert "3. Neither the name of the copyright holder" in text
    assert "AS IS" in text


def test_pyproject_declares_bsd():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = cfg["project"]
    assert project["license"] == "BSD-3-Clause"
    assert "LICENSE" in project["license-files"]


def test_index_and_license_ship_as_package_data():
    """``pip install .`` 出来的副本也要带索引和许可原文。

    漏了 ``*.sqlite``，装出来的包会静默退化成纯 Tier A；
    漏了 ``*.txt``，则违反 4.2（公开传输须附许可证副本）。
    """
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    patterns = cfg["tool"]["setuptools"]["package-data"]["nmrflow.data"]
    assert "*.sqlite" in patterns
    assert "*.txt" in patterns


def test_upstream_license_text_is_present():
    license_txt = DATA / "nmrshiftdb2_database_license.txt"
    notice = DATA / "HOSE_INDEX_NOTICE.txt"
    assert license_txt.exists(), "索引的许可证原文缺失"
    assert notice.exists(), "来源声明缺失"
    assert "nmrshiftdb2" in license_txt.read_text(encoding="utf-8").lower()


# ---------------------------------------------------------------------------
# 产出物的来源声明
# ---------------------------------------------------------------------------
def test_attribution_only_when_index_is_used():
    on = run("CCO", RunOptions(with_ms=False))
    off = run("CCO", RunOptions(tier_b=False, with_ms=False))
    assert on.report.hose_used is True
    assert on.report.attribution.startswith(EXPECTED)
    assert on.report.attribution.rstrip().endswith(EXPECTED_TAIL)
    # 没用索引就不能挂声明，否则是虚假归属
    assert off.report.hose_used is False
    assert off.report.attribution == ""


def test_attribution_stays_out_of_the_paper_text():
    """``report.text`` 是给用户粘进论文的，不能插英文声明。

    声明放在导出文件末尾和界面小字里；正文必须仍以 HRMS 行结束。
    """
    report = run("CCO", RunOptions()).report
    assert report.hose_used is True
    assert "nmrshiftdb" not in report.text
    assert report.text.rstrip().endswith("[M+H]+.")


def test_json_carries_attribution(tmp_path: Path):
    result = run("CCO", RunOptions(with_ms=False))
    path = dump_json(result, tmp_path / "out.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["attribution"].startswith(EXPECTED)

    result_off = run("CCO", RunOptions(tier_b=False, with_ms=False))
    path_off = dump_json(result_off, tmp_path / "off.json")
    assert json.loads(path_off.read_text(encoding="utf-8"))["attribution"] is None


@pytest.mark.parametrize("smiles", ["CCO", "C1CC=CCC1", "O=P(c1ccccc1)(c2ccccc2)c3ccccc3"])
def test_attribution_is_single_line_per_sentence(smiles: str):
    """声明是两行英文，不能出现三行以上或被折成一段——许可证要求逐字。"""
    attr = run(smiles, RunOptions(with_ms=False)).report.attribution
    assert attr.splitlines() == [EXPECTED, EXPECTED_TAIL]
