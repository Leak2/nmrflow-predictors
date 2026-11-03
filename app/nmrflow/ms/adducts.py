"""加合物注册表：读取 adducts.toml，按名称/电离方式检索。"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

POSITIVE = "positive"
NEGATIVE = "negative"


@dataclass(frozen=True)
class Adduct:
    name: str            # "[M+H]+"
    mult: int = 1        # 中性分子个数（[2M+H]+ -> 2）
    add: dict[str, int] = field(default_factory=dict)
    sub: dict[str, int] = field(default_factory=dict)
    charge: int = 1
    mode: str = POSITIVE
    ei_only: bool = False

    @property
    def abs_charge(self) -> int:
        return abs(self.charge)

    def formula_counts(self, base: dict[str, int]) -> dict[str, int]:
        """加合物离子的元素组成。"""
        out = {k: v * self.mult for k, v in base.items()}
        for sym, n in self.add.items():
            out[sym] = out.get(sym, 0) + n
        for sym, n in self.sub.items():
            out[sym] = out.get(sym, 0) - n
        return {k: v for k, v in out.items() if v > 0}


@lru_cache(maxsize=1)
def _load() -> tuple[dict[str, Adduct], dict[str, list[str]]]:
    with (DATA_DIR / "adducts.toml").open("rb") as fh:
        raw = tomllib.load(fh)
    table: dict[str, Adduct] = {}
    for mode in (POSITIVE, NEGATIVE):
        for body in raw.get(mode, []):
            ad = Adduct(
                name=body["name"],
                mult=int(body.get("mult", 1)),
                add=dict(body.get("add", {})),
                sub=dict(body.get("sub", {})),
                charge=int(body.get("charge", 1 if mode == POSITIVE else -1)),
                mode=mode,
                ei_only=bool(body.get("ei_only", False)),
            )
            table[ad.name] = ad
    methods = {k.upper(): list(v) for k, v in raw.get("methods", {}).items()}
    return table, methods


def all_adducts() -> dict[str, Adduct]:
    return dict(_load()[0])


def get_adduct(name: str) -> Adduct:
    table, _ = _load()
    if name not in table:
        raise KeyError(f"未注册的加合物 {name!r}，可用：{', '.join(sorted(table))}")
    return table[name]


def adducts_for_method(method: str, polarity: str = "positive") -> list[Adduct]:
    """按电离方式取默认加合物。EI -> M+·；其余按 methods 表。"""
    table, methods = _load()
    names = methods.get(method.upper())
    if names is None:
        names = ["[M+H]+"] if polarity.lower().startswith("pos") else ["[M-H]-"]
    return [table[n] for n in names if n in table]


def available_methods() -> list[str]:
    return sorted(_load()[1])
