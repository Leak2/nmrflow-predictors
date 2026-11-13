"""exe 版本资源生成器（``tools/version_info.py``）的护栏。

守三件容易悄悄坏掉的事：
1. 版本号只有**一个**来源（``nmrflow.__version__``），工具不许自己抄一份；
2. 生成的文件必须是**纯 ASCII**——中文直接写进去，PyInstaller 在某些版本/区域设置下
   会按 cp1252/cp936 解码，表现是打包时抛 UnicodeDecodeError 或者属性页乱码；
3. ``FileVersion`` 必须是四个整数，``0.2.7`` 要补成 ``0.2.7.0``。
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _load():
    spec = importlib.util.spec_from_file_location(
        "version_info_tool", ROOT / "app" / "tools" / "version_info.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def vinfo():
    return _load()


def test_version_has_single_source(vinfo):
    import nmrflow

    assert vinfo.read_version() == nmrflow.__version__


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("0.2.7", (0, 2, 7, 0)),
        ("1.10.0", (1, 10, 0, 0)),
        ("0.2.8rc1", (0, 2, 8, 0)),
        ("2.0.0.3", (2, 0, 0, 3)),
    ],
)
def test_version_tuple_pads_to_four(vinfo, raw, expected):
    assert vinfo.version_tuple(raw) == expected


def test_render_is_pure_ascii(vinfo):
    text = vinfo.render("0.2.7")
    assert text.isascii(), "版本资源文件必须纯 ASCII，中文要走 \\uXXXX"
    assert "\u8c31" not in text


def test_render_embeds_version_and_both_languages(vinfo):
    text = vinfo.render("0.2.7")
    assert "filevers=(0, 2, 7, 0)" in text
    assert "prodvers=(0, 2, 7, 0)" in text
    assert "080404b0" in text and "040904b0" in text  # zh-CN / en-US
    assert "StringStruct('FileVersion', '0.2.7.0')" in text
