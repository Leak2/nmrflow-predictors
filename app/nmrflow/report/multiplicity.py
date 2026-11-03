"""多重峰命名：由 J 的集合推出 s / d / t / q / dd / td / m / br s。"""

from __future__ import annotations

_NAMES = {1: "d", 2: "t", 3: "q", 4: "quint", 5: "sext", 6: "sept", 7: "oct"}

# 单个耦合表项：(J / Hz, 等价核个数, 伙伴核)
Item = tuple[float, int, str]


def _normalize(j_values) -> list[Item]:
    """接受 [(j, n, partner)] / [(j, n)] / [j] 三种写法。"""
    out: list[Item] = []
    for item in j_values or []:
        if isinstance(item, (tuple, list)):
            j = float(item[0])
            n = int(item[1]) if len(item) > 1 else 1
            partner = str(item[2]) if len(item) > 2 else "1H"
        else:
            j, n, partner = float(item), 1, "1H"
        if abs(j) > 0.05:
            out.append((abs(j), n, partner))
    return out


def _cluster(items: list[Item], tol: float = 0.35) -> list[tuple[float, int, str]]:
    """把相近的 J 聚成一类，返回 [(代表值, 等价核总数, 伙伴核)]，按 |J| 降序。

    count 是**等价核的个数**之和，不是耦合条目的条数：
    -CH2- 被一个 -CH3 (3H) 裂分 -> [(7.0, 3, '1H')] -> q（而不是 d）。
    """
    if not items:
        return []
    ordered = sorted(items, key=lambda t: -t[0])
    clusters: list[list[Item]] = [[ordered[0]]]
    for item in ordered[1:]:
        if abs(clusters[-1][0][0] - item[0]) <= tol:
            clusters[-1].append(item)
        else:
            clusters.append([item])
    out = [
        (
            sum(j for j, _, _ in c) / len(c),
            sum(n for _, n, _ in c),
            c[0][2],
        )
        for c in clusters
    ]
    return sorted(out, key=lambda t: -abs(t[0]))


def multiplicity(j_values, *, broad: bool = False, n_lines: int = 1,
                 drop_j_for_m: bool = False) -> tuple[str, list[float], list[str]]:
    """返回 (多重峰名, 报告里要写的 J 列表, 每个 J 的伙伴核)。

    drop_j_for_m：1H 谱里 m 不写 J（谱线太密，报出来没有意义）；
    但 13C/19F/31P 的 J 本身就是关键数据（如 CF3 的 JC-F=320），必须保留。
    """
    if broad:
        return "br s", [], []
    items = _normalize(j_values)
    if not items:
        return "s", [], []

    clusters = _cluster(items)
    distinct = [round(v, 1) for v, _, _ in clusters]
    partners = [p for _, _, p in clusters]

    if len(clusters) == 1:
        _, count, _ = clusters[0]
        return _NAMES.get(count, "m"), distinct, partners

    if len(clusters) == 2:
        (_, c1, _), (_, c2, _) = clusters
        first, second = _NAMES.get(c1, "m"), _NAMES.get(c2, "m")
        if n_lines <= 8:
            return first + second, distinct, partners
        return ("m", [], partners) if drop_j_for_m else ("m", distinct, partners)

    if len(clusters) == 3:
        if n_lines <= 8:
            return "ddd", distinct, partners
        return ("m", [], partners) if drop_j_for_m else ("m", distinct, partners)

    return ("m", [], partners) if drop_j_for_m else ("m", distinct, partners)
