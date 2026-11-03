"""谱图数字化：谱线 -> (δ ppm, 强度) 网格 + 积分曲线。"""

from __future__ import annotations

import numpy as np

from ..core.constants import (
    DEFAULT_FWHM,
    LABILE_FWHM,
    LOW_CONF_FWHM,
    MARGIN_MIN,
    MARGIN_PPM,
    POINTS_PER_PPM,
)
from ..core.types import Confidence, Signal, Spectrum1D, Transition


# 采样分辨率目标：每个点跨 0.1 Hz（远小于 0.5~2 Hz 的线宽）
TARGET_RESOLUTION_HZ = 0.1
MAX_POINTS = 60_000


def _fwhm_for(signal: Signal | None, nucleus: str, resolution_hz: float) -> float:
    """线宽 = 本征线宽 ⊕ 数字分辨率（平方和开方，近似卷积）。

    低置信度确实应该画得“虚”一点，但 12 Hz 的 VAGUE_FWHM 只适合 ¹⁹F/³¹P
    的未知环境；对 400 MHz 的 ¹H 来说 12 Hz = 0.03 ppm，会把 J=7 Hz 的四重峰
    糊成单包。这里按核给上限，详见 ``core.constants.LOW_CONF_FWHM``。
    """
    if signal is not None and signal.fwhm_hz is not None:
        base = signal.fwhm_hz
    else:
        base = DEFAULT_FWHM.get(nucleus, 1.0)
    if signal is not None:
        if signal.is_labile:
            base = max(base, LABILE_FWHM)
        elif signal.confidence is Confidence.LOW:
            base = max(base, LOW_CONF_FWHM.get(nucleus, LABILE_FWHM))
    return float(np.sqrt(base ** 2 + resolution_hz ** 2))


def digitize(
    transitions: list[Transition],
    signals: list[Signal],
    nucleus: str,
    freq_mhz: float,
    solvent: str,
    *,
    points_per_ppm: float = POINTS_PER_PPM,
    margin_ppm: float | None = None,
    shape: str = "lorentzian",
    noise: float = 0.0,
    rng: np.random.Generator | None = None,
) -> Spectrum1D:
    from .lineshape import render

    if not transitions:
        x = np.linspace(-1.0, 12.0, 200)
        return Spectrum1D(nucleus, freq_mhz, solvent, x, np.zeros_like(x))

    freqs = np.array([t.freq_hz for t in transitions])
    signal_span = (freqs.max() - freqs.min()) / freq_mhz
    if margin_ppm is None:
        margin_ppm = min(
            MARGIN_PPM.get(nucleus, 0.6),
            max(MARGIN_MIN.get(nucleus, 0.2), 0.12 * signal_span),
        )
    lo_ppm = freqs.min() / freq_mhz - margin_ppm
    hi_ppm = freqs.max() / freq_mhz + margin_ppm
    span_ppm = hi_ppm - lo_ppm
    # 采样密度必须比线宽细得多，否则谱线会被数字分辨率抹平，
    # 出现「单峰比多重峰高一大截」的假象。
    n_points = int(
        max(
            800,
            span_ppm * points_per_ppm,
            span_ppm * freq_mhz / TARGET_RESOLUTION_HZ,
        )
    )
    n_points = min(n_points, MAX_POINTS)
    x_ppm = np.linspace(lo_ppm, hi_ppm, n_points)
    x_hz = x_ppm * freq_mhz
    resolution_hz = span_ppm / n_points * freq_mhz

    by_group = {s.group_id: s for s in signals}
    y = np.zeros_like(x_hz)
    # 按组渲染：每组的线宽可能不同（活泼氢更宽）
    groups: dict[int, list[Transition]] = {}
    for tr in transitions:
        groups.setdefault(tr.group_id, []).append(tr)
    for gid, trs in groups.items():
        fwhm = _fwhm_for(by_group.get(gid), nucleus, resolution_hz)
        y += render(x_hz, trs, fwhm, shape)

    if noise > 0:
        gen = rng if rng is not None else np.random.default_rng(0)
        y = y + gen.normal(0.0, noise * max(y.max(), 1e-9), size=y.shape)

    return Spectrum1D(
        nucleus=nucleus, freq_mhz=freq_mhz, solvent=solvent,
        x_ppm=x_ppm, y=y, transitions=transitions, signals=signals,
    )


def integral_curve(spec: Spectrum1D) -> np.ndarray:
    """1H 积分曲线：归一化到总质子数。"""
    if spec.x_ppm.size == 0:
        return np.array([])
    dx = np.abs(spec.x_ppm[1] - spec.x_ppm[0])
    # δ 轴是反向的，积分要从高场往低场累加
    cum = np.cumsum(spec.y[::-1])[::-1] * dx
    total = cum[0] if cum.size else 0.0
    if total <= 0:
        return cum
    n_h = sum(s.n_nuclei for s in spec.signals) or 1
    return cum / total * n_h
