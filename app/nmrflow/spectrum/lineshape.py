"""谱线线型：洛伦兹 / 高斯 / Voigt 近似。"""

from __future__ import annotations

import numpy as np

from ..core.types import Transition


def lorentzian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    g = fwhm / 2.0
    return (g * g) / ((x - center) ** 2 + g * g)


def gaussian(x: np.ndarray, center: float, fwhm: float) -> np.ndarray:
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return np.exp(-0.5 * ((x - center) / sigma) ** 2)


def voigt(x: np.ndarray, center: float, fwhm: float, eta: float = 0.5) -> np.ndarray:
    """eta=0 纯高斯，eta=1 纯洛伦兹。"""
    return eta * lorentzian(x, center, fwhm) + (1 - eta) * gaussian(x, center, fwhm)


SHAPES = {"lorentzian": lorentzian, "gaussian": gaussian, "voigt": voigt}


def render(
    x_hz: np.ndarray,
    transitions: list[Transition],
    fwhm_hz: float,
    shape: str = "lorentzian",
    cutoff: float = 8.0,
) -> np.ndarray:
    """把谱线列表渲染成强度数组。

    只在线中心 ±cutoff·FWHM 的窗口内累加：洛伦兹尾巴在 8 倍半高宽外
    已经降到 1/256，肉眼不可见，但省掉了 95% 以上的无效计算
    （对 2^10 个跃迁 × 数万个采样点来说这是决定性的）。
    """
    fn = SHAPES.get(shape, lorentzian)
    y = np.zeros_like(x_hz)
    if x_hz.size == 0 or not transitions:
        return y
    dx = abs(float(x_hz[1] - x_hz[0])) or 1.0
    half = max(cutoff * fwhm_hz, 2.0 * dx)
    for tr in transitions:
        if tr.intensity <= 0:
            continue
        lo = int(np.searchsorted(x_hz, tr.freq_hz - half, side="left"))
        hi = int(np.searchsorted(x_hz, tr.freq_hz + half, side="right"))
        if hi <= lo:
            continue
        seg = x_hz[lo:hi]
        y[lo:hi] += tr.intensity * fn(seg, tr.freq_hz, fwhm_hz)
    return y
