"""1H 谱图：积分曲线 + 峰顶 δ 标注 + 残余溶剂峰。"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from ..core.solvent import get_solvent
from ..core.types import Spectrum1D
from .style import apply, color, mark_predicted, save_predicted


def plot_1h(spec: Spectrum1D, out_path: str | Path, *,
            show_integral: bool = True, show_solvent: bool = True,
            annotate: bool = True, peaks=None) -> Path:
    apply()
    out_path = Path(out_path)
    fig, ax = plt.subplots(figsize=(11, 4.6), dpi=160)

    ax.plot(spec.x_ppm, spec.y, color=color("curve"), lw=1.0)
    ax.fill_between(spec.x_ppm, spec.y, color=color("curve"), alpha=0.06)

    lo, hi = float(spec.x_ppm.min()), float(spec.x_ppm.max())
    ax.set_xlim(hi, lo)                      # δ 轴反向
    ymax = spec.y.max() if spec.y.size else 1.0
    # 顶部留 22% 余量，否则峰顶的竖排标注会顶到标题上
    ax.set_ylim(-0.02 * max(ymax, 1e-9), ymax * 1.22)

    # 积分曲线
    if show_integral:
        from ..spectrum.digitize import integral_curve
        integ = integral_curve(spec)
        if integ.size:
            total = integ[0] if integ[0] > 0 else 1.0
            scaled = integ / total * ymax * 0.32
            ax.plot(spec.x_ppm, scaled, color=color("integral"), lw=0.8)

    # 残余溶剂峰
    if show_solvent:
        try:
            sol = get_solvent(spec.solvent)
        except KeyError:
            sol = None
        if sol is not None:
            for rp in sol.residual_1h:
                if lo <= rp.shift <= hi:
                    ax.axvline(rp.shift, color=color("residual"), lw=0.7, ls=":")
                    ax.annotate(
                        rp.label or "solvent", xy=(rp.shift, ymax * 0.97),
                        color=color("residual"), fontsize=7, ha="center", va="top",
                    )

    # 峰标注：优先用报告峰（每个多重峰只标中心 δ + 质子数）
    if annotate:
        if peaks:
            for p in peaks:
                _label(ax, spec, p.shift, f"{p.shift:.2f}",
                       extra=f"{int(round(p.integration))}H",
                       lo=p.shift_lo, hi=p.shift_hi)
        else:
            for p in _peak_tops(spec):
                _label(ax, spec, p, f"{p:.2f}")

    ax.set_xlabel("化学位移 (ppm)")
    ax.set_ylabel("强度 (相对)")
    ax.set_title(
        f"1H NMR 预测谱  ·  {spec.freq_mhz:.0f} MHz  ·  {spec.solvent}"
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


def _label(ax, spec: Spectrum1D, ppm: float, text: str,
           extra: str | None = None, lo: float | None = None,
           hi: float | None = None) -> None:
    """在峰中心立一条竖排标注：δ 值 + 质子数。

    纵坐标取**该多重峰区间内的最高点**，否则多重峰中心的谷底会把
    标注压到谱线里面去。
    """
    if lo is not None and hi is not None:
        a, b = min(lo, hi), max(lo, hi)
        mask = (spec.x_ppm >= a) & (spec.x_ppm <= b)
        y = float(spec.y[mask].max()) if mask.any() else _y_at(spec, ppm)
    else:
        y = _y_at(spec, ppm)
    label = f"{text}  {extra}" if extra else text
    ax.annotate(
        label, xy=(ppm, y), xytext=(0, 5), textcoords="offset points",
        color=color("annot"), fontsize=7.5, ha="center", va="bottom",
        rotation=90,
    )


def _peak_tops(spec: Spectrum1D, min_prominence: float = 0.02) -> list[float]:
    """简单拾峰：局部极大且高于最大强度的 min_prominence。"""
    y = spec.y
    if y.size < 3:
        return []
    thr = y.max() * min_prominence
    tops = []
    for i in range(1, len(y) - 1):
        if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] > thr:
            tops.append(float(spec.x_ppm[i]))
    return _dedupe(tops, 0.06)


def _dedupe(values: list[float], tol: float) -> list[float]:
    out: list[float] = []
    for v in sorted(values, reverse=True):
        if not out or abs(out[-1] - v) > tol:
            out.append(v)
    return out
