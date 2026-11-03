"""溶剂注册表：残余溶剂峰、活泼氢位移、溶剂效应。"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


@dataclass
class ResidualPeak:
    shift: float
    label: str = ""


@dataclass
class LabileSpec:
    center: float
    width_hz: float = 5.0


@dataclass
class Solvent:
    key: str
    display: str
    residual_1h: list[ResidualPeak] = field(default_factory=list)
    residual_13c: list[ResidualPeak] = field(default_factory=list)
    drop_exchangeable: bool = False
    labile_couples: bool = False
    shift_offset: float = 0.0
    labile: dict[str, LabileSpec] = field(default_factory=dict)

    def labile_spec(self, cls: str) -> LabileSpec | None:
        return self.labile.get(cls)


@lru_cache(maxsize=1)
def _load_all() -> tuple[dict[str, Solvent], dict[str, str]]:
    path = DATA_DIR / "solvents.toml"
    with path.open("rb") as fh:
        raw = tomllib.load(fh)

    aliases: dict[str, str] = dict(raw.pop("aliases", {}))
    solvents: dict[str, Solvent] = {}

    for key, body in raw.items():
        labile = {
            k: LabileSpec(center=float(v["center"]), width_hz=float(v.get("width_hz", 5.0)))
            for k, v in (body.get("labile") or {}).items()
        }
        solvents[key] = Solvent(
            key=key,
            display=body.get("display", key),
            residual_1h=[ResidualPeak(float(p["shift"]), p.get("label", ""))
                         for p in body.get("residual_1h", [])],
            residual_13c=[ResidualPeak(float(p["shift"]), p.get("label", ""))
                          for p in body.get("residual_13c", [])],
            drop_exchangeable=bool(body.get("drop_exchangeable", False)),
            labile_couples=bool(body.get("labile_couples", False)),
            shift_offset=float(body.get("shift_offset", 0.0)),
            labile=labile,
        )
    return solvents, aliases


def list_solvents() -> list[str]:
    solvents, _ = _load_all()
    return sorted(solvents)


def get_solvent(name: str) -> Solvent:
    """按名称取溶剂，支持别名。未知溶剂回落 CDCl3 并告警。"""
    solvents, aliases = _load_all()
    if name in solvents:
        return solvents[name]
    key = aliases.get(name)
    if key and key in solvents:
        return solvents[key]
    lowered = {k.lower(): v for k, v in aliases.items()}
    key = lowered.get(name.strip().lower())
    if key and key in solvents:
        return solvents[key]
    raise KeyError(f"未知溶剂 {name!r}，可用：{', '.join(list_solvents())}")
