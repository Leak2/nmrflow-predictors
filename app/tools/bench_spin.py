"""分块对角化的收益基准：耗时与峰值内存。

用法：
    python tools/bench_spin.py [dense上限] [blocked上限]

dense 路径内存按 8·4^n 字节估算，blocked 按最大块 C(n,n/2)² 估算。
两条路径都直接调内部实现，避免被 solve_spin_system 的小体系阈值掩盖。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import comb

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from nmrflow.spin.hamiltonian import (  # noqa: E402
    _DENSE_MAX_SPINS,
    _solve_blocked,
    _solve_dense,
)


def _case(n: int, seed: int = 7):
    rng = np.random.default_rng(seed + n)
    shifts = np.sort(rng.uniform(100.0, 900.0, n)) + rng.uniform(-3.0, 3.0, n)
    j = np.triu(rng.uniform(0.5, 12.0, (n, n)), 1)
    return shifts, j + j.T


def main() -> None:
    max_dense = int(sys.argv[1]) if len(sys.argv) > 1 else 11
    max_block = int(sys.argv[2]) if len(sys.argv) > 2 else 13

    print(f"小体系阈值 _DENSE_MAX_SPINS = {_DENSE_MAX_SPINS}\n")
    print(f"{'n':>3} {'dense/s':>9} {'blocked/s':>10} {'加速':>7} "
          f"{'dense内存':>10} {'block内存':>10} {'谱线':>7}")
    print("-" * 62)

    for n in range(4, max_block + 1):
        shifts, j = _case(n)

        t0 = time.perf_counter()
        blocked = _solve_blocked(shifts, j, None)
        t_block = time.perf_counter() - t0

        if n <= max_dense:
            t0 = time.perf_counter()
            dense = _solve_dense(shifts, j, None)
            t_dense = time.perf_counter() - t0
            assert len(dense) == len(blocked), "两条路径谱线数不一致"
            cell_d, cell_x = f"{t_dense:>9.3f}", f"{t_dense / max(t_block, 1e-9):>6.1f}x"
        else:
            cell_d, cell_x = f"{'—':>9}", f"{'—':>7}"

        mem_dense = 8.0 * (4.0 ** n) / 2**20
        mem_block = 8.0 * float(comb(n, n // 2)) ** 2 / 2**20
        print(f"{n:>3} {cell_d} {t_block:>10.3f} {cell_x} "
              f"{mem_dense:>8.0f}MB {mem_block:>8.1f}MB {len(blocked):>7}")


if __name__ == "__main__":
    main()
