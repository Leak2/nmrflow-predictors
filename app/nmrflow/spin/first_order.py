"""一级（弱耦合）裂分：n+1 规则 + 二项强度比，多组耦合做卷积。"""

from __future__ import annotations

from collections.abc import Iterable

from ..core.types import Signal, Transition

# 帕斯卡三角：n 个等价 I=1/2 核的裂分强度比
_PASCAL: dict[int, list[int]] = {}


def pascal(n: int) -> list[int]:
    if n in _PASCAL:
        return _PASCAL[n]
    row = [1]
    for _ in range(n):
        row = [a + b for a, b in zip([0] + row, row + [0])]
    _PASCAL[n] = row
    return row


def split_once(
    lines: list[tuple[float, float]], j_hz: float, n_equiv: int
) -> list[tuple[float, float]]:
    """把每条谱线按 n 个等价核、耦合常数 j 再裂分一次。

    n 个等价 I=1/2 核 -> 2nI+1 条线，强度比取帕斯卡三角第 n 行，
    频率偏移 (k - n/2) * j，k = 0..n。
    """
    weights = pascal(n_equiv)
    out: list[tuple[float, float]] = []
    for freq, inten in lines:
        for k, w in enumerate(weights):
            if w == 0:
                continue
            out.append((freq + (k - n_equiv / 2.0) * j_hz, inten * w))
    return _merge_close(out)


def _merge_close(lines: Iterable[tuple[float, float]], tol: float = 1e-6
                 ) -> list[tuple[float, float]]:
    buckets: dict[float, float] = {}
    for f, i in lines:
        key = round(f / 1e-4) * 1e-4
        buckets[key] = buckets.get(key, 0.0) + i
    return sorted(((k, v) for k, v in buckets.items()), key=lambda t: (t[0], -t[1]))


def first_order_transitions(signal: Signal, freq_mhz: float) -> list[Transition]:
    """一个信号的一级裂分谱线。

    返回 (频率偏移 Hz, 相对强度)，强度总和 = 该组核数（即积分面积）。
    """
    center_hz = signal.shift * freq_mhz
    lines: list[tuple[float, float]] = [(0.0, 1.0)]
    for c in signal.couplings:
        lines = split_once(lines, c.j_hz, c.n_equiv)

    total = sum(i for _, i in lines) or 1.0
    scale = signal.n_nuclei / total
    return [
        Transition(freq_hz=center_hz + f, intensity=i * scale, group_id=signal.group_id)
        for f, i in lines
    ]
