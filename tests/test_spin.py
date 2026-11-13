"""自旋体系：一级裂分、帕斯卡三角、AB 体系精确解与解析解对拍。"""

from __future__ import annotations

import numpy as np
import pytest

from nmrflow.core.units import freq_of
from nmrflow.spin.first_order import first_order_transitions, pascal
from nmrflow.spin.hamiltonian import (
    _DENSE_MAX_SPINS,
    _solve_blocked,
    _solve_dense,
    exact_transitions,
    solve_spin_system,
)
from nmrflow.core.types import Signal


def test_pascal():
    assert pascal(0) == [1]
    assert pascal(1) == [1, 1]
    assert pascal(2) == [1, 2, 1]
    assert pascal(3) == [1, 3, 3, 1]
    assert pascal(4) == [1, 4, 6, 4, 1]


def test_nucleus_frequencies():
    """旋磁比锚定：400 MHz 质子 -> 19F 376.28 / 13C 100.58 / 31P 161.92"""
    assert freq_of("19F", 400) == pytest.approx(376.28, abs=0.01)
    assert freq_of("13C", 400) == pytest.approx(100.60, abs=0.01)
    assert freq_of("31P", 400) == pytest.approx(161.92, abs=0.01)


def test_quartet_from_three_equivalent_h():
    sig = Signal(nucleus="1H", group_id=1, shift=1.20, n_nuclei=3)
    sig.couplings = []
    from nmrflow.core.types import Coupling
    target = Signal(nucleus="1H", group_id=2, shift=3.50, n_nuclei=2)
    target.couplings = [Coupling(partner_group=1, j_hz=7.0, n_equiv=3)]
    trs = first_order_transitions(target, 400.0)
    assert len(trs) == 4
    js = np.diff(sorted(t.freq_hz for t in trs))
    assert all(abs(j - 7.0) < 1e-6 for j in js)
    # 强度总和 = 该组核数（2H），所以是帕斯卡比 × 2/8
    weights = [t.intensity for t in sorted(trs, key=lambda t: t.freq_hz)]
    assert weights == pytest.approx([0.25, 0.75, 0.75, 0.25])


def test_ab_system_matches_analytic():
    """AB：Δν=6 Hz、J=-12 Hz、60 MHz。数值解必须与解析解逐位吻合。"""
    freq = 60.0
    dnu, j = 6.0, -12.0
    nu_a = 100.0
    nu_b = nu_a + dnu
    packed = solve_spin_system(
        shifts_hz=[nu_a, nu_b],
        j_matrix=[[0.0, j], [j, 0.0]],
    )
    lines = [f for f, _ in packed]
    inten = [i for _, i in packed]
    # 解析解：四条线在 νc ± D/2 ± J/2，D = sqrt(Δν² + J²)
    d = np.hypot(dnu, j)
    centre = (nu_a + nu_b) / 2
    expect = [
        centre - d / 2 - j / 2,
        centre - d / 2 + j / 2,
        centre + d / 2 - j / 2,
        centre + d / 2 + j / 2,
    ]
    order = np.argsort(lines)
    got = np.array(lines)[order]
    assert got == pytest.approx(sorted(expect), abs=1e-9)

    # 强度：外侧 (1 - |J|/D)/2、内侧 (1 + |J|/D)/2，内外各两条
    outer = (1 - abs(j) / d) / 2
    inner = (1 + abs(j) / d) / 2
    w = np.array(inten)[order]
    assert w == pytest.approx([outer, inner, inner, outer], abs=1e-9)
    assert w.sum() == pytest.approx(2.0, abs=1e-9)


def test_first_order_becomes_second_order_at_low_field():
    """同一对核，60 MHz 是二级谱、800 MHz 近似一级谱。"""
    from nmrflow.core.units import is_strongly_coupled

    # Δδ = 0.1 ppm：60 MHz 下 Δν = 6 Hz，Δν/J = 0.86 -> 强耦合（二级谱）
    #                 800 MHz 下 Δν = 80 Hz，Δν/J = 11.4 -> 一级谱
    assert is_strongly_coupled(0.0, 0.1, 7.0, 60)
    assert not is_strongly_coupled(0.0, 0.1, 7.0, 800)


def test_exact_transitions_counts():
    """两个等价核 + 一个异核：ΔM=±1 的跃迁条数有限且强度归一。"""
    trs = exact_transitions(
        shifts_ppm=[1.0, 3.0], j_matrix=[[0.0, 12.0], [12.0, 0.0]],
        freq_mhz=400.0, group_ids=[0, 1], group_sizes=[1, 1],
    )
    assert len(trs) == 4
    total = sum(t.intensity for t in trs)
    assert total == pytest.approx(2.0, abs=1e-6)   # 2 个核，总强度 2


def _sorted_peaks(trs: list[tuple[float, float]]):
    arr = np.array(sorted(trs), dtype=float)
    return arr[:, 0], arr[:, 1]


@pytest.mark.parametrize(
    "shifts,jmat",
    [
        # AB：Δν = 6 Hz、J = -12 Hz
        ([100.0, 106.0], [[0.0, -12.0], [-12.0, 0.0]]),
        # ABX：两两均有耦合
        (
            [100.0, 118.0, 460.0],
            [[0.0, -12.0, 7.0], [-12.0, 0.0, 3.0], [7.0, 3.0, 0.0]],
        ),
        # AA'BB' 近似：4 自旋、位移两档、邻位/对位 J 不同
        (
            [700.0, 702.0, 690.0, 692.0],
            [
                [0.0, 0.0, 8.0, 2.5],
                [0.0, 0.0, 2.5, 8.0],
                [8.0, 2.5, 0.0, 0.0],
                [2.5, 8.0, 0.0, 0.0],
            ],
        ),
        # 5 自旋链：邻位强耦合、隔位弱耦合
        (
            [300.0, 320.0, 345.0, 372.0, 400.0],
            [
                [0.0, 7.5, 1.2, 0.4, 0.1],
                [7.5, 0.0, 7.0, 1.1, 0.3],
                [1.2, 7.0, 0.0, 6.8, 1.0],
                [0.4, 1.1, 6.8, 0.0, 7.2],
                [0.1, 0.3, 1.0, 7.2, 0.0],
            ],
        ),
    ],
)
def test_blocked_matches_dense(shifts, jmat):
    """按 M 分块对角化必须与整块对角化给出同一组谱线。

    直接调内部实现，避免被 solve_spin_system 的小体系阈值短路成同一条路径。
    """
    sh = np.asarray(shifts, dtype=float)
    jm = np.asarray(jmat, dtype=float)
    blocked = _solve_blocked(sh, jm, None)
    dense = _solve_dense(sh, jm, None)

    fb, ib = _sorted_peaks(blocked)
    fd, idn = _sorted_peaks(dense)

    assert len(fb) == len(fd)
    assert fb == pytest.approx(fd, abs=1e-6)
    assert ib == pytest.approx(idn, abs=1e-6)
    assert ib.sum() == pytest.approx(len(shifts), abs=1e-6)


def test_solver_path_switch_consistent():
    """n 越过阈值后默认走分块，结果仍与强制整块一致。"""
    n = _DENSE_MAX_SPINS + 1
    shifts = np.array([100.0 + 60.0 * i for i in range(n)])
    jm = np.zeros((n, n))
    for i in range(n - 1):
        jm[i, i + 1] = jm[i + 1, i] = 7.5 - 0.3 * i

    _fa, ia = _sorted_peaks(solve_spin_system(shifts, jm))
    _fd, idn = _sorted_peaks(solve_spin_system(shifts, jm, force_dense=True))
    assert ia == pytest.approx(idn, abs=1e-6)


def test_blocked_single_spin():
    """单自旋无耦合：一条谱线，强度归一为 1。"""
    out = solve_spin_system([1234.5], [[0.0]])
    assert len(out) == 1
    assert out[0][0] == pytest.approx(1234.5)
    assert out[0][1] == pytest.approx(1.0)


def test_strongly_coupled_ab_keeps_four_lines():
    """AB 体系必须保留 4 条线。

    0.2.4 曾把“归属距离 ≤ Σ|J|/2”当过滤条件，而 AB 的四条线正好落在
    νc ± (D/2 ± J/2)，离各自中心约 |J|/2，二阶修正后略超就被丢掉——
    AX/AB 都只剩 2 条。强度随后被重新归一化，总积分仍然正确，所以
    光看积分发现不了。预算必须是**该自旋自己的** Σ|J|。
    """
    trs = exact_transitions(
        shifts_ppm=[1.000, 1.020], j_matrix=[[0.0, 12.0], [12.0, 0.0]],
        freq_mhz=400.0, group_ids=[0, 1], group_sizes=[1, 1],
    )
    assert len(trs) == 4
    assert sum(t.intensity for t in trs) == pytest.approx(2.0, abs=1e-6)
    # 四条线两两关于各自中心对称：390.79 / 402.79 / 405.21 / 417.21 Hz
    freqs = sorted(t.freq_hz for t in trs)
    assert len({round(f, 1) for f in freqs}) == 4
