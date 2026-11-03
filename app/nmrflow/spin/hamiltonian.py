"""自旋哈密顿量精确对角化（处理 AB / AA'BB' 等二级效应）。

H = Σ ν_i Iz_i + Σ_{i<j} J_ij (I_i · I_j)          （单位 Hz）

要点：
1. [H, Iz_total] = 0，所以 M = Σ m_i 是好量子数；
2. eigh 在简并子空间里返回任意组合，不能按特征向量下标配对，
   必须用 ΔM = ±1 的选择定则筛选跃迁；
3. 强度 ∝ |<p|Ix_total|q>|²，Ix_total 在积基下只有"单自旋翻转"矩阵元 1/2。

第 1 点意味着 H 天然按 M 分块对角，每块维度只有 C(n, k)。对整块 2^n 做
eigh 等于在大量零块上白算：12 自旋时 4096² 稠密矩阵 134 MB，而最大的块
只有 C(12,6)=924 维、6.8 MB，eigh 的运算量降一个量级。

主路径 `_solve_blocked` 因此逐块对角化，Ix 矩阵元只在相邻块之间计算
（ΔM=±1 恰好落在块间）。`force_dense=True` 保留原来的整块路径作对拍基准。
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

from ..core.types import Transition

# 跃迁强度阈值：低于此值不记为一条独立谱线
_INTENSITY_CUTOFF = 1e-10

# 小体系仍走整块路径：分块要在 n+1 个小矩阵上重复承担 Python 层开销，
# n <= 6 时这笔开销超过省下来的 eigh 运算量（实测 n=4 反而慢约 40%）。
# 交叉点实测在 n = 7 附近，收益随 n 迅速放大（n=10 约 10x）。
_DENSE_MAX_SPINS = 6


def _bit_matrices(n: int):
    """返回 (m_table, flip_index)。

    m_table[s] = 长度 n 的数组，元素为 ±0.5
    flip_index[k][s] = s ^ (1 << k)
    """
    size = 1 << n
    idx = np.arange(size)
    m_table = np.empty((size, n))
    flips = []
    for k in range(n):
        bit = (idx >> k) & 1
        m_table[:, k] = np.where(bit == 1, 0.5, -0.5)
        flips.append(idx ^ (1 << k))
    return m_table, np.array(flips)


def _m_groups(n: int) -> list[np.ndarray]:
    """按 M 分块：返回长度 n+1 的列表，第 k 项是"自旋向上数为 k"的状态索引。

    M = k - n/2，k 递增即 M 递增，所以相邻下标正好对应 ΔM = ±1。
    """
    size = 1 << n
    idx = np.arange(size)
    popcount = np.zeros(size, dtype=np.int64)
    for k in range(n):
        popcount += (idx >> k) & 1
    return [np.nonzero(popcount == k)[0] for k in range(n + 1)]


def _block_matrix(
    shifts_hz: np.ndarray,
    j_matrix: np.ndarray,
    states: np.ndarray,
    inverse: np.ndarray,
) -> np.ndarray:
    """states 张成的 M 子空间内的哈密顿量（块内稠密，维度 C(n, k)）。"""
    n = len(shifts_hz)
    m = np.empty((len(states), n))
    for k in range(n):
        m[:, k] = np.where((states >> k) & 1, 0.5, -0.5)

    diag = m @ shifts_hz
    pairs: list[tuple[np.ndarray, np.ndarray, float]] = []
    for i in range(n):
        for k in range(i + 1, n):
            j = j_matrix[i, k]
            if j == 0:
                continue
            diag = diag + j * m[:, i] * m[:, k]
            # flip-flop 是 i 与 k 一上一下互换，M 不变，目标仍留在本块。
            # 反平行的一对态互为源与目标，遍历 src 时两个方向自然各写一次。
            src = np.nonzero((m[:, i] * m[:, k]) < 0)[0]
            if src.size:
                dst = inverse[states[src] ^ (1 << i) ^ (1 << k)]
                pairs.append((dst, src, 0.5 * j))

    h = np.zeros((len(states), len(states)))
    np.fill_diagonal(h, diag)
    for dst, src, v in pairs:
        h[dst, src] += v
    return h


def build_hamiltonian(shifts_hz: np.ndarray, j_matrix: np.ndarray) -> np.ndarray:
    """构造整块稠密哈密顿量（2^n × 2^n）。

    主流程不再调用（见 `_solve_blocked`），保留给 force_dense 对拍与外部调试。
    """
    shifts_hz = np.asarray(shifts_hz, dtype=float)
    j_matrix = np.asarray(j_matrix, dtype=float)
    n = len(shifts_hz)
    size = 1 << n
    m_table, flips = _bit_matrices(n)

    h = np.zeros((size, size))
    # 对角项：Σ ν_i m_i + Σ_{i<k} J_ik m_i m_k
    diag = m_table @ shifts_hz
    for i in range(n):
        for k in range(i + 1, n):
            if j_matrix[i, k] == 0:
                continue
            diag = diag + j_matrix[i, k] * m_table[:, i] * m_table[:, k]
    np.fill_diagonal(h, diag)

    # 非对角（flip-flop）项：J/2 连接相差一次自旋翻转的态。
    # 只写 h[dst, src]：遍历所有 src 时反向的一对会自然补上 h[src, dst]，
    # 若此处同时写两个方向会重复计入一倍。
    for i in range(n):
        for k in range(i + 1, n):
            j = j_matrix[i, k]
            if j == 0:
                continue
            mask = (m_table[:, i] * m_table[:, k]) < 0
            src = np.nonzero(mask)[0]
            dst = flips[i][flips[k][src]]
            h[dst, src] += 0.5 * j
    return h


def _normalise(
    transitions: list[tuple[float, float]],
    n: int,
    weights: np.ndarray | None,
) -> list[tuple[float, float]]:
    """归一化到"每组核数之和"。"""
    if not transitions:
        return []
    total = sum(i for _, i in transitions) or 1.0
    nspin = float(np.sum(weights)) if weights is not None else float(n)
    scale = nspin / total
    return [(f, i * scale) for f, i in transitions]


def _solve_dense(
    shifts_hz: np.ndarray,
    j_matrix: np.ndarray,
    weights: np.ndarray | None,
) -> list[tuple[float, float]]:
    """整块对角化（原实现，作为对拍基准保留）。"""
    n = len(shifts_hz)
    size = 1 << n
    h = build_hamiltonian(shifts_hz, j_matrix)
    energies, u = eigh(h)

    m_table, flips = _bit_matrices(n)
    m_product = m_table.sum(axis=1)
    # 每个本征态的 M（M 是好量子数，用概率加权平均即精确值）
    m_eig = (u * u).T @ m_product
    m_eig = np.round(m_eig * 2) / 2.0

    # Ix_total 在本征基下的矩阵 A = 0.5 Σ_k U^T P_k U
    a = np.zeros((size, size))
    for k in range(n):
        a += 0.5 * (u.T @ u[flips[k]])

    transitions: list[tuple[float, float]] = []
    levels = sorted(set(m_eig.tolist()))
    level_set = set(levels)
    for m in levels:
        if m + 1 not in level_set:
            continue
        low = np.nonzero(np.isclose(m_eig, m, atol=1e-6))[0]
        high = np.nonzero(np.isclose(m_eig, m + 1, atol=1e-6))[0]
        if low.size == 0 or high.size == 0:
            continue
        block = a[np.ix_(low, high)]
        de = energies[high][None, :] - energies[low][:, None]
        inten = block * block
        keep = inten > _INTENSITY_CUTOFF
        for (i, jj) in zip(*np.nonzero(keep)):
            freq = float(de[i, jj])
            # 频率可为负：δ < 0（如 19F 的 -122 ppm）对应负的绝对频率，
            # 早期版本在这里丢掉 freq<=0，导致只有正位移的核出峰。
            if abs(freq) < 1e-9:
                continue
            transitions.append((freq, float(inten[i, jj])))
    return _normalise(transitions, n, weights)


def _solve_blocked(
    shifts_hz: np.ndarray,
    j_matrix: np.ndarray,
    weights: np.ndarray | None,
) -> list[tuple[float, float]]:
    """按 M 子空间分块对角化（主路径）。"""
    n = len(shifts_hz)
    groups = _m_groups(n)
    inverse = np.full(1 << n, -1, dtype=np.int64)
    for states in groups:
        inverse[states] = np.arange(len(states))

    energies: list[np.ndarray] = []
    vectors: list[np.ndarray] = []
    for states in groups:
        e, u = eigh(_block_matrix(shifts_hz, j_matrix, states, inverse))
        energies.append(e)
        vectors.append(u)

    transitions: list[tuple[float, float]] = []
    for k in range(n):
        lo, hi = groups[k], groups[k + 1]
        u_lo, u_hi = vectors[k], vectors[k + 1]

        # Ix 只在 ΔM = ±1 的块之间有矩阵元：第 i 位为 0（m = -1/2）的态
        # 翻转后进入 M+1 块，贡献 1/2。
        a = np.zeros((len(lo), len(hi)))
        for i in range(n):
            sel = np.nonzero(((lo >> i) & 1) == 0)[0]
            if sel.size == 0:
                continue
            dst = inverse[lo[sel] ^ (1 << i)]
            a += 0.5 * (u_lo[sel, :].T @ u_hi[dst, :])

        de = energies[k + 1][None, :] - energies[k][:, None]
        inten = a * a
        keep = inten > _INTENSITY_CUTOFF
        for (p, q) in zip(*np.nonzero(keep)):
            freq = float(de[p, q])
            if abs(freq) < 1e-9:
                continue
            transitions.append((freq, float(inten[p, q])))
    return _normalise(transitions, n, weights)


def solve_spin_system(
    shifts_hz: np.ndarray,
    j_matrix: np.ndarray,
    weights: np.ndarray | None = None,
    *,
    force_dense: bool = False,
) -> list[tuple[float, float]]:
    """对角化并返回 [(跃迁频率 Hz, 相对强度)]。

    自旋数超过 `_DENSE_MAX_SPINS` 时按 M 分块对角化；小体系走整块（更快）。
    force_dense=True 强制走整块路径，仅用于验证两条路径等价。
    """
    shifts_hz = np.asarray(shifts_hz, dtype=float)
    j_matrix = np.asarray(j_matrix, dtype=float)
    if shifts_hz.size == 0:
        return []
    if force_dense or shifts_hz.size <= _DENSE_MAX_SPINS:
        return _solve_dense(shifts_hz, j_matrix, weights)
    return _solve_blocked(shifts_hz, j_matrix, weights)


def exact_transitions(
    shifts_ppm: list[float],
    j_matrix: list[list[float]],
    freq_mhz: float,
    group_ids: list[int],
    group_sizes: list[int],
    center_nucleus: int = 0,
) -> list[Transition]:
    """把一个自旋块整体对角化，输出谱线。

    shifts_ppm : 每个自旋的化学位移（磁等价核重复列出）
    j_matrix   : 自旋-自旋耦合矩阵（磁等价核之间为 0）
    group_ids  : 每个自旋所属的信号组（用于归属）
    """
    shifts_hz = np.array(shifts_ppm) * freq_mhz
    jm = np.array(j_matrix, dtype=float)
    raw = solve_spin_system(shifts_hz, jm, np.array(group_sizes, dtype=float))

    # 分块对角化会同时给出 M→M+1 与反向的“发射/镜像”跃迁。对 7 个以上、
    # 位移又接近的 ¹H 块，这些镜像线的能量差可以是负的，画出来就是负 ppm
    # 的假峰（胆固醇的 "-0.29-1.68 (m, 8H)" 就是这么来的）。
    #
    # 全正/全负位移的块只保留物理吸收方向：¹H/¹³C 是全正，¹⁹F 通常是全负。
    # 混合符号的块（少数 ³¹P）不做过滤，避免把真实谱线误删；它们如果出现
    # 镜像线，会由调用方的置信度/告警机制兜底。
    if shifts_ppm and all(s >= 0 for s in shifts_ppm):
        raw = [(f, i) for f, i in raw if f > 0]
    elif shifts_ppm and all(s <= 0 for s in shifts_ppm):
        raw = [(f, i) for f, i in raw if f < 0]

    # 再做一次物理窗口过滤：强耦合跃迁会偏离中心位移，偏离量的上界取“该自旋
    # 参与的全部耦合绝对值之和”Σ_k|J_ik|，而不是全局 Σ|J|/2。
    # ★ 用 Σ|J|/2 会把 AX / AB 体系直接砍掉一半谱线：这两类体系的四条线正好
    #   落在 νc ± (D/2 ± J/2)，离各自中心约 |J|/2，二阶修正后略超 |J|/2 就被
    #   判成“离得太远”。而强度随后又被重新归一化，总积分看起来仍然正确——
    #   这类错误不会自己报出来，只能靠 exact_transitions 的条数回归抓。
    n = len(shifts_hz)
    # 每个组自己的位移预算（Hz）。单组的单量子跃迁偏离中心不会超过它。
    limits_hz = [
        sum(abs(float(jm[i, k])) for k in range(n) if k != i) for i in range(n)
    ]
    if n > 1:
        pad = max(limits_hz)
        lo = float(shifts_hz.min()) - pad
        hi = float(shifts_hz.max()) + pad
        raw = [(f, i) for f, i in raw if lo <= f <= hi]

    # 归属：把每条谱线归到位移最接近的组；同时要求它离该组中心的距离不超过
    # **该组自己**的耦合预算（见上）。这一层把强耦合块里仍落在“全块窗口”内的
    # 组合/伪线剔掉（胆固醇 2.86 ppm 的假线）。
    out: list[Transition] = []
    for freq_hz, inten in raw:
        best, best_d, best_idx = group_ids[0], 1e18, 0
        for idx, (gid, s) in enumerate(zip(group_ids, shifts_ppm)):
            d = abs(freq_hz - s * freq_mhz)
            if d < best_d:
                best, best_d, best_idx = gid, d, idx
        if best_d > limits_hz[best_idx] + 1e-9:
            continue
        out.append(Transition(freq_hz=freq_hz, intensity=inten, group_id=best))

    # 最后再归一化：过滤掉线以后强度总和会小于“自旋数”，按存活下来的线重新
    # 分配到 target，保证 ¹H 积分不缩水。
    # ★ 必须放在**所有**过滤之后。0.2.4 把它放在归属过滤之前，被归属过滤掉的线
    #   会连带把总强度带走（胆固醇 46.0 → 45.99），积分对不上，而峰形看起来正常。
    total = sum(t.intensity for t in out)
    target = float(sum(group_sizes)) if group_sizes else float(len(shifts_ppm))
    if total > 0 and target > 0:
        scale = target / total
        out = [
            Transition(freq_hz=t.freq_hz, intensity=t.intensity * scale,
                       group_id=t.group_id)
            for t in out
        ]
    return out
