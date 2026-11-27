"""版本号一致性：`nmrflow.__version__` 是唯一真源，其余文件必须跟着它。

为什么需要这个文件
------------------
版本号目前在 5 个文件里出现，各自的守护情况并不一致：

| 文件 | 有何守护 |
|---|---|
| ``app/nmrflow/__init__.py`` | 真源本身 |
| ``app/tools/version_info.py`` | 已有（读 ``__version__``，不抄） |
| ``app/tools/make_release.py`` | 已有（``test_packaging`` 断言 ``mr.VERSION``） |
| ``pyproject.toml`` | **此前无守护** |
| ``CITATION.cff`` | **此前无守护** |
| ``.zenodo.json`` | **此前无守护** |

后三个漏改的后果比前三个更重：

- ``pyproject.toml`` 错了，``pip install`` 出来的包版本和 ``nmrflow --version``
  对不上，报 bug 时无法判断用户装的是哪一版；
- ``CITATION.cff`` / ``.zenodo.json`` 错了，Zenodo 归档元数据会永久带着错误版本号，
  而 DOI 记录**不可修改**——只能发新版本覆盖。

所以这里把它们钉死。CITATION.cff 是 YAML，但只取两个标量，用正则足够，
不为此引入 PyYAML 依赖（这个项目连 pandas 都只是可选）。
"""

from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path

import nmrflow

ROOT = Path(__file__).resolve().parent.parent
VERSION = nmrflow.__version__


def _cff_scalar(text: str, key: str) -> str:
    """从 CITATION.cff 里取一个顶层标量的值（去掉可能的引号）。"""
    match = re.search(
        r"^" + re.escape(key) + r":\s*(.+?)\s*$",
        text,
        re.MULTILINE,
    )
    assert match is not None, f"CITATION.cff 里找不到 {key}"
    return match.group(1).strip().strip('"').strip("'")


# ---------------------------------------------------------------------------
# pyproject.toml
# ---------------------------------------------------------------------------
def test_pyproject_version_matches():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert cfg["project"]["version"] == VERSION, (
        f"pyproject.toml 是 {cfg['project']['version']}，"
        f"但 nmrflow.__version__ 是 {VERSION}"
    )


# ---------------------------------------------------------------------------
# CITATION.cff
# ---------------------------------------------------------------------------
def test_citation_version_matches():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert _cff_scalar(text, "version") == VERSION


def test_citation_url_matches_version():
    """url 指向 releases/tag/vX.Y.Z，tag 名必须与版本号一致。"""
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    url = _cff_scalar(text, "url")
    assert url.endswith(f"/tag/v{VERSION}"), f"CITATION.cff 的 url 是 {url}"


# ---------------------------------------------------------------------------
# .zenodo.json
# ---------------------------------------------------------------------------
def test_zenodo_version_matches():
    meta = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    assert meta["version"] == VERSION


def test_zenodo_version_identifier_matches():
    """related_identifiers 里 isVersionOf 那条指向 tag，必须与版本号一致。"""
    meta = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    versioned = [
        item["identifier"]
        for item in meta.get("related_identifiers", [])
        if item.get("relation") == "isVersionOf"
    ]
    assert versioned, ".zenodo.json 缺少 relation=isVersionOf 的条目"
    for url in versioned:
        assert url.endswith(f"/tag/v{VERSION}"), f".zenodo.json 的版本链接是 {url}"


# ---------------------------------------------------------------------------
# 兜底：全仓扫描，防止又冒出一个手抄版本号的地方
# ---------------------------------------------------------------------------
def test_no_stray_version_literals_in_tracked_configs():
    """除已知位置外，配置类文件里不应再出现独立的版本号字面量。

    只扫"配置/元数据"这类容易被忽略的文件；代码与数据表不扫（数据表有自己的
    ``version`` 字段，语义是表结构版本，不是软件版本）。
    """
    watched = ["pyproject.toml", "CITATION.cff", ".zenodo.json"]
    for name in watched:
        text = (ROOT / name).read_text(encoding="utf-8")
        # `cff-version` 是引文文件格式规范的版本（当前 1.2.0），不是软件版本
        text = re.sub(r"^cff-version:.*$", "", text, flags=re.MULTILINE)
        found = set(re.findall(r"\bv?(\d+\.\d+\.\d+)\b", text))
        # 允许出现当前版本；出现别的版本号说明有没跟着改的地方
        unexpected = {v for v in found if v != VERSION}
        assert not unexpected, f"{name} 里出现与当前版本不符的版本号：{unexpected}"
