"""报告峰簇构建：把同一信号的谱线合成一个峰（或峰区间），再按需合并。"""

from __future__ import annotations

from ..core.constants import MERGE_SPAN_PPM, MERGE_TOL_PPM
from ..core.types import Confidence, Peak, Signal, Transition
from .multiplicity import multiplicity


def build_peaks(
    signals: list[Signal],
    transitions: list[Transition],
    freq_mhz: float,
    nucleus: str,
    *,
    round_digits: int = 2,
) -> list[Peak]:
    """每个信号生成一个 Peak；谱线来自精确对角化时按频率区间表示。"""
    by_group: dict[int, list[Transition]] = {}
    for tr in transitions:
        by_group.setdefault(tr.group_id, []).append(tr)

    peaks: list[Peak] = []
    for sig in signals:
        trs = by_group.get(sig.group_id, [])
        if trs:
            freqs = sorted(t.freq_hz for t in trs)
            lo, hi = freqs[0] / freq_mhz, freqs[-1] / freq_mhz
            n_lines = len({round(f, 4) for f in freqs})
            # 峰簇中心取强度加权平均
            wsum = sum(t.intensity for t in trs) or 1.0
            center = sum(t.freq_hz * t.intensity for t in trs) / wsum / freq_mhz
        else:
            lo = hi = center = sig.shift
            n_lines = 1

        js = [(c.j_hz, c.n_equiv, c.partner_nucleus) for c in sig.couplings]
        mult, reported_j, reported_nuc = multiplicity(
            js, broad=sig.is_labile, n_lines=n_lines,
            drop_j_for_m=(sig.nucleus == "1H"),
        )
        peaks.append(
            Peak(
                shift=round(center, round_digits),
                shift_lo=round(lo, round_digits),
                shift_hi=round(hi, round_digits),
                intensity=sum(t.intensity for t in trs),
                integration=float(sig.n_nuclei),
                multiplicity=mult,
                j_values=reported_j,
                j_partners=reported_nuc,
                label=sig.label,
                confidence=sig.confidence,
                n_lines=n_lines,
            )
        )
    peaks.sort(key=lambda p: -p.shift)
    return _merge_close(peaks, nucleus, round_digits)


def _merge_close(peaks: list[Peak], nucleus: str,
                 round_digits: int) -> list[Peak]:
    """δ 容差内且归属相同的峰合并成一个区间，积分求和。

    ``MERGE_SPAN_PPM`` 限制合并后的总跨度，避免相邻 CH₂ 峰链式合并成
    一个覆盖整段脂肪区的假宽峰（胆固醇 8H / 0.35~1.68 ppm）。
    """
    tol = MERGE_TOL_PPM.get(nucleus, 0.05)
    span_cap = MERGE_SPAN_PPM.get(nucleus, 1.0)
    out: list[Peak] = []
    for p in peaks:
        if out:
            prev = out[-1]
            same_label = prev.label == p.label
            close = abs(prev.shift - p.shift) <= tol
            lo = min(prev.shift_lo, p.shift_lo)
            hi = max(prev.shift_hi, p.shift_hi)
            if close and same_label and (hi - lo) <= span_cap:
                merged = Peak(
                    shift=round((lo + hi) / 2, round_digits),
                    shift_lo=round(lo, round_digits),
                    shift_hi=round(hi, round_digits),
                    intensity=prev.intensity + p.intensity,
                    integration=prev.integration + p.integration,
                    multiplicity="m",
                    j_values=[],
                    j_partners=[],
                    label=prev.label,
                    confidence=_worst(prev.confidence, p.confidence),
                    n_lines=prev.n_lines + p.n_lines,
                )
                out[-1] = merged
                continue
        out.append(p)
    return out


def _worst(a: Confidence, b: Confidence) -> Confidence:
    return a if a.rank >= b.rank else b
