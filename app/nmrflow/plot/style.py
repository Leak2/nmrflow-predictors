"""matplotlib 样式：中文字体、δ 轴反向、配色。"""

from __future__ import annotations

import matplotlib
from matplotlib import font_manager
from pathlib import Path

from .. import __version__

CJK_CANDIDATES = [
    Path("C:/Windows/Fonts/simhei.ttf"),
    Path("C:/Windows/Fonts/msyh.ttc"),
    Path("C:/Windows/Fonts/simsun.ttc"),
    Path("/System/Library/Fonts/PingFang.ttc"),
    Path("/usr/share/fonts/truetype/noto/NotoSansCJK-Regular.ttc"),
]

_COLOR = {
    "bg": "#ffffff",
    "axes": "#1c1c1c",
    "grid": "#e6e6e6",
    "curve": "#123f6b",
    "integral": "#7a7a7a",
    "annot": "#333333",
    "residual": "#c8c8c8",
    "ms": "#8a3b3b",
}


def color(name: str) -> str:
    return _COLOR.get(name, "#000000")


_SETUP_DONE: str | None = None


def setup() -> str:
    """注册中文字体并返回实际使用的字体名（只做一次）。"""
    global _SETUP_DONE
    if _SETUP_DONE is not None:
        return _SETUP_DONE
    matplotlib.use("Agg")
    # 关掉 mathtext 解析：装了中文字体后 mathtext 每次都要重扫字体表，
    # 一条 $^{1}$H 标题就能让单张图从 0.2 s 涨到 3 s。
    matplotlib.rcParams["text.parse_math"] = False
    for path in CJK_CANDIDATES:
        if path.exists():
            try:
                font_manager.fontManager.addfont(str(path))
                name = font_manager.FontProperties(fname=str(path)).get_name()
                matplotlib.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
                matplotlib.rcParams["axes.unicode_minus"] = False
                _SETUP_DONE = name
                return name
            except Exception:
                continue
    matplotlib.rcParams["axes.unicode_minus"] = False
    _SETUP_DONE = "DejaVu Sans"
    return _SETUP_DONE


def apply() -> None:
    setup()
    matplotlib.rcParams.update({
        "figure.facecolor": _COLOR["bg"],
        "axes.facecolor": _COLOR["bg"],
        "axes.edgecolor": _COLOR["axes"],
        "axes.labelcolor": _COLOR["axes"],
        "xtick.color": _COLOR["axes"],
        "ytick.color": _COLOR["axes"],
        "text.color": _COLOR["axes"],
        "axes.grid": False,
        "savefig.dpi": 160,
        "font.size": 11,
    })


def mark_predicted(fig) -> None:
    """右下角固定水印：裁掉标题也还在。"""
    fig.text(
        0.995, 0.006, "PREDICTED - not experimental",
        ha="right", va="bottom", fontsize=7,
        color=color("annot"), alpha=0.55,
    )


def save_predicted(fig, out_path) -> None:
    """保存 PNG，并写入可机读的 Software/Description。"""
    fig.savefig(
        out_path,
        metadata={
            "Software": f"nmrflow {__version__}",
            "Description": "Predicted NMR spectrum - not experimental data",
        },
    )
