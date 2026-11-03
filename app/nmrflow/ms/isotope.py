"""同位素分布卷积：计算 M / M+1 / M+2 ... 的相对强度。"""

from __future__ import annotations

from ..core.constants import ELECTRON_MASS
from ..core.types import IsotopePeak
from .masses import get_element

# 分布剪枝：低于最强峰此比例的组分丢弃
_PRUNE_REL = 1e-8
# 最多保留的峰数
_MAX_PEAKS = 60

State = list[tuple[int, float, float]]  # (ordinal, prob, Σ prob × mass_excess)


def _mul(a: State, b: State) -> State:
    """卷积两个分布。excess 累加：sa*pb + sb*pa。"""
    out: dict[int, tuple[float, float]] = {}
    for oa, pa, sa in a:
        if pa <= 0:
            continue
        for ob, pb, sb in b:
            if pb <= 0:
                continue
            o = oa + ob
            p = pa * pb
            s = sa * pb + sb * pa
            if o in out:
                op, os_ = out[o]
                out[o] = (op + p, os_ + s)
            else:
                out[o] = (p, s)
    return [(o, p, s) for o, (p, s) in sorted(out.items())]


def _prune(state: State, max_peaks: int = _MAX_PEAKS, rel: float = _PRUNE_REL) -> State:
    """丢弃极弱组分，防止大分子爆炸。"""
    if not state:
        return state
    peak = max(p for _, p, _ in state)
    if peak <= 0:
        return state
    kept = [t for t in state if t[1] >= peak * rel]
    if len(kept) > max_peaks:
        kept = sorted(kept, key=lambda t: -t[1])[:max_peaks]
        kept.sort()
    return kept


def _element_state(symbol: str, n: int) -> State:
    """n 个某元素的同位素分布（多项式幂，二分加速 + 剪枝）。"""
    elem = get_element(symbol)
    mono = elem.monoisotopic
    base: State = [
        (round(iso.mass) - round(mono), iso.abundance, iso.abundance * (iso.mass - mono))
        for iso in elem.isotopes
        if iso.abundance > 0
    ]
    result: State = [(0, 1.0, 0.0)]
    power = base
    e = int(n)
    while e > 0:
        if e & 1:
            result = _prune(_mul(result, power))
        e >>= 1
        if e:
            power = _prune(_mul(power, power))
    return result


def isotope_pattern(
    counts: dict[str, int],
    *,
    charge: int = 0,
    min_rel_abundance: float = 1e-4,
    max_peaks: int = _MAX_PEAKS,
) -> list[IsotopePeak]:
    """给定离子元素组成，返回同位素峰簇。

    charge 用于扣电子（正离子扣、负离子加）；0 表示不修正。
    """
    state: State = [(0, 1.0, 0.0)]
    for sym, n in sorted(counts.items()):
        if n <= 0:
            continue
        state = _prune(_mul(state, _element_state(sym, n)))

    base = sum(get_element(sym).monoisotopic * n for sym, n in counts.items())
    base -= charge * ELECTRON_MASS

    peaks: list[IsotopePeak] = []
    for ordinal, prob, excess_sum in state:
        if prob <= 0:
            continue
        mean_excess = excess_sum / prob
        peaks.append(IsotopePeak(mz=base + mean_excess, abundance=prob, ordinal=ordinal))

    if not peaks:
        return []

    top = max(p.abundance for p in peaks)
    out = [IsotopePeak(p.mz, p.abundance / top, p.ordinal) for p in peaks]
    out = [p for p in out if p.abundance >= min_rel_abundance]
    out.sort(key=lambda p: p.mz)
    return out[:max_peaks]
