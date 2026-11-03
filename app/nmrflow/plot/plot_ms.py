"""质谱同位素分布棒图。"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from ..core.types import MsResult
from .style import apply, color, mark_predicted, save_predicted


def plot_ms(ms: MsResult, out_path: str | Path) -> Path:
    apply()
    out_path = Path(out_path)
    if not ms.isotope_pattern:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "无同位素分布数据", ha="center", va="center")
        mark_predicted(fig)
        save_predicted(fig, out_path)
        plt.close(fig)
        return out_path

    mz = [p.mz for p in ms.isotope_pattern]
    rel = [p.abundance * 100 for p in ms.isotope_pattern]

    fig, ax = plt.subplots(figsize=(8, 3.4), dpi=160)
    # 棒宽固定 0.05 Da：乘 max(mz) 会得到接近 1 Da 的粗块，把整段谱糊掉
    ax.bar(mz, rel, width=0.05, color=color("ms"))
    for x, y, p in zip(mz, rel, ms.isotope_pattern):
        if y >= 1.0:
            ax.annotate(
                f"M+{p.ordinal}" if p.ordinal else "M",
                xy=(x, y), xytext=(0, 3), textcoords="offset points",
                ha="center", fontsize=7.5, color=color("annot"),
            )
    # 只给有峰的位点打刻度，避免一堆空刻度
    ticks = [x for x, y in zip(mz, rel) if y >= 0.5]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{x:.2f}" for x in ticks], fontsize=7.5)
    if len(ticks) > 1:
        ax.set_xlim(min(ticks) - 0.45, max(ticks) + 0.45)
    ax.set_ylim(0, max(rel) * 1.25 + 5)
    ax.set_xlabel("m/z")
    ax.set_ylabel("相对丰度 (%)")
    ax.set_title(f"预测同位素分布  ·  {ms.formula}  ·  {ms.method}")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    mark_predicted(fig)
    save_predicted(fig, out_path)
    plt.close(fig)
    return out_path
