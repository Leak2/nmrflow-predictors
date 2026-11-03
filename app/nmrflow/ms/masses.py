"""元素同位素质量与丰度表。"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@dataclass(frozen=True)
class Isotope:
    mass: float
    abundance: float


@dataclass(frozen=True)
class Element:
    symbol: str
    monoisotopic: float
    isotopes: tuple[Isotope, ...]

    @property
    def most_abundant(self) -> Isotope:
        return max(self.isotopes, key=lambda i: i.abundance)


@lru_cache(maxsize=1)
def _load_table() -> dict[str, Element]:
    with (DATA_DIR / "isotopes.toml").open("rb") as fh:
        raw = tomllib.load(fh)
    table: dict[str, Element] = {}
    for sym, body in raw["elements"].items():
        isos = tuple(
            Isotope(float(i["mass"]), float(i["abundance"]))
            for i in body.get("isotopes", [])
        )
        table[sym] = Element(sym, float(body["monoisotopic"]), isos)
    return table


def get_element(symbol: str) -> Element:
    table = _load_table()
    if symbol not in table:
        raise KeyError(f"同位素表中没有元素 {symbol!r}")
    return table[symbol]


def monoisotopic_mass(symbol: str) -> float:
    return get_element(symbol).monoisotopic


def known_elements() -> list[str]:
    return sorted(_load_table())
