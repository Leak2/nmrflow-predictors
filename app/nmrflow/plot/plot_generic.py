"""13C / 19F / 31P 谱图（无积分曲线）。"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..core.types import Spectrum1D
from .style import apply, color, mark_predicted, save_predicted

TITLES = {
    "13C": "13C NMR 预测谱",
    "19F": "19F NMR 预测谱",
    "31P": "31P NMR 预测谱",
}

# 峰顶标注里的积分后缀（13C 不积分，故不列）
_NUCLEUS_SYMBOL = {"1H": "H", "19F": "F", "31P": "P"}


def _peak_text(peak, nucleus: str) -> str:
    """`-109.0  2F`：位移 + 该峰对应的等价核个数。"""
    text = f"{peak.shift:.1f}"
    symbol = _NUCLEUS_SYMBOL.get(nucleus)
    if symbol and getattr(peak, "integration", 0):
        text += f"  {int(round(peak.integration))}{symbol}"
    return text


def plot_generic(spec: Spectrum1D, out_path: str | Path, *,
                 annotate: bool = True, peaks=None) -> Path:
    apply()
    out_path = Path(out_path)
    fig, ax = plt.subplots(figsize=(11, 3.6), dpi=160)

    ax.plot(spec.x_ppm, spec.y, color=color("curve"), lw=1.1)
    ax.fill_between(spec.x_ppm, spec.y, color=color("curve"), alpha=0.08)

    lo, hi = float(spec.x_ppm.min()), float(spec.x_ppm.max())
    ax.set_xlim(hi, lo)
    ax.set_ylim(-0.02 * max(spec.y.max(), 1e-9), spec.y.max() * 1.22)

    if annotate:
        if peaks:
            for p in peaks:
                ax.annotate(
                    _peak_text(p, spec.nucleus),
                    xy=(p.shift, _y_at(spec, p.shift)),
                    xytext=(0, 5), textcoords="offset points",
                    color=color("annot"), fontsize=7.5, ha="center", rotation=90,
                )
        elif spec.nucleus != "13C":
            for p in _tops(spec):
                ax.annotate(
                    f"{p:.1f}", xy=(p, _y_at(spec, p)), xytext=(0, 6),
                    textcoords="offset points", color=color("annot"),
                    fontsize=7.5, ha="center", rotation=90,
                )

    ax.set_xlabel("化学位移 (ppm)")
    ax.set_ylabel("强度 (相对)")
    ax.set_title(
        f"{TITLES.get(spec.nucleus, spec.nucleus + ' NMR')}  ·  "
        f"{spec.freq_mhz:.0f} MHz  ·  {spec.solvent}"
    )
    ax.set_yticks([])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    mark_predicted(fig)
    save_predicted(fig, out_path)
    plt.close(fig)
    return out_path


def _y_at(spec: Spectrum1D, ppm: float) -> float:
    idx = int(np.argmin(np.abs(spec.x_ppm - ppm)))
    return float(spec.y[idx])


def _tops(spec: Spectrum1D, min_prominence: float = 0.05) -> list[float]:
    y = spec.y
    if y.size < 3:
        return []
    thr = y.max() * min_prominence
    tops = []
    for i in range(1, len(y) - 1):
        if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] > thr:
            tops.append(float(spec.x_ppm[i]))
    out: list[float] = []
    for v in sorted(tops, reverse=True):
        if not out or abs(out[-1] - v) > 0.15:
            out.append(v)
    return out
