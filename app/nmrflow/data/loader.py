"""TOML 数据表加载（带缓存与最小 schema 校验）。"""

from __future__ import annotations

import tomllib
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent


@lru_cache(maxsize=32)
def load_toml(name: str) -> dict:
    path = DATA_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"数据表不存在：{path}")
    with path.open("rb") as fh:
        return tomllib.load(fh)


def shift_table(nucleus: str) -> dict:
    fname = {"1H": "shift_1h.toml", "13C": "shift_13c.toml",
             "19F": "shift_19f.toml", "31P": "shift_31p.toml"}[nucleus]
    return load_toml(fname)


def coupling_table(name: str) -> dict:
    return load_toml(f"{name}.toml")


def reload() -> None:
    """清空缓存（改表后热重载）。"""
    load_toml.cache_clear()
