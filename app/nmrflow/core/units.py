"""单位换算与频率推导。"""

from __future__ import annotations

from .constants import GAMMA_MHZ_T, REFERENCE_NUCLEUS, STRONG_COUPLING_RATIO


def freq_of(nucleus: str, proton_mhz: float) -> float:
    """由 1H 观测频率推导其它核的观测频率（MHz）。

    400 MHz 的 1H 谱仪 -> 19F 376.28 / 13C 100.58 / 31P 161.92 MHz
    """
    return proton_mhz * GAMMA_MHZ_T[nucleus] / GAMMA_MHZ_T[REFERENCE_NUCLEUS]


def round_freq(nucleus: str, proton_mhz: float) -> float:
    """文献写法：频率取整到整数 MHz（400 -> 376 / 100 / 162）。"""
    return round(freq_of(nucleus, proton_mhz))


def ppm_to_hz(delta_ppm: float, freq_mhz: float) -> float:
    """化学位移差（ppm）-> 频率差（Hz）。"""
    return delta_ppm * freq_mhz


def hz_to_ppm(delta_hz: float, freq_mhz: float) -> float:
    """频率差（Hz）-> 化学位移差（ppm）。"""
    return delta_hz / freq_mhz


def coupling_ratio(
    delta_a_ppm: float, delta_b_ppm: float, j_hz: float, freq_mhz: float
) -> float:
    """Δν/|J|，判断是否强耦合。J 为 0 时返回 inf（视为弱耦合）。"""
    if abs(j_hz) < 1e-9:
        return float("inf")
    return abs(ppm_to_hz(delta_a_ppm - delta_b_ppm, freq_mhz) / j_hz)


def is_strongly_coupled(
    delta_a_ppm: float,
    delta_b_ppm: float,
    j_hz: float,
    freq_mhz: float,
    ratio: float = STRONG_COUPLING_RATIO,
) -> bool:
    """Δν/|J| < ratio 即强耦合（必须按 Hz 算，随场强变化）。"""
    return coupling_ratio(delta_a_ppm, delta_b_ppm, j_hz, freq_mhz) < ratio


def ppm_error(measured: float, calculated: float) -> float:
    """质量误差，单位 ppm（质谱用，与 NMR 的 ppm 无关）。"""
    if calculated == 0:
        return 0.0
    return (measured - calculated) / calculated * 1e6
