"""自旋体系分块：耦合图连通分量 + 强耦合判定。

注意：Δν/J 必须按 **Hz** 计算，所以 freq_mhz 是必需参数
—— 同一个分子在 60 MHz 是二级谱，在 800 MHz 可能就成了一级谱。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..core.constants import MAX_EXACT_SPINS, STRONG_COUPLING_RATIO
from ..core.types import Signal
from ..core.units import ppm_to_hz


@dataclass
class SpinBlock:
    group_ids: list[int] = field(default_factory=list)
    spin_group: list[int] = field(default_factory=list)   # 每个自旋属于哪个组
    shifts: list[float] = field(default_factory=list)     # 每个自旋的位移 ppm
    j: list[list[float]] = field(default_factory=list)
    n_per_group: list[int] = field(default_factory=list)
    exact: bool = False
    note: str = ""

    @property
    def n_spins(self) -> int:
        return len(self.shifts)


def _components(signals: list[Signal]) -> tuple[list[list[int]], dict[tuple[int, int], float]]:
    """耦合图连通分量 + 边上的 J。"""
    jmap: dict[tuple[int, int], float] = {}
    adj: dict[int, set[int]] = {s.group_id: set() for s in signals}
    for s in signals:
        for c in s.couplings:
            a, b = s.group_id, c.partner_group
            if b not in adj:
                continue
            key = (min(a, b), max(a, b))
            if key not in jmap:
                jmap[key] = c.j_hz
            adj[a].add(b)
            adj[b].add(a)

    seen: set[int] = set()
    comps: list[list[int]] = []
    for gid in adj:
        if gid in seen:
            continue
        stack, comp = [gid], []
        seen.add(gid)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        comps.append(sorted(comp))
    return comps, jmap


def build_blocks(
    signals: list[Signal], freq_mhz: float, max_spins: int = MAX_EXACT_SPINS
) -> tuple[list[SpinBlock], list[str]]:
    """把全部信号切成自旋块，并决定每块走精确对角化还是一级近似。"""
    by_id = {s.group_id: s for s in signals}
    comps, jmap = _components(signals)
    warnings: list[str] = []
    blocks: list[SpinBlock] = []

    for comp in comps:
        if len(comp) == 1:
            gid = comp[0]
            sig = by_id[gid]
            blocks.append(
                SpinBlock(
                    group_ids=[gid], spin_group=[gid] * sig.n_nuclei,
                    shifts=[sig.shift] * sig.n_nuclei,
                    j=[[0.0] * sig.n_nuclei for _ in range(sig.n_nuclei)],
                    n_per_group=[sig.n_nuclei], exact=False,
                )
            )
            continue

        # 是否存在强耦合对
        strong = False
        for i, a in enumerate(comp):
            for b in comp[i + 1:]:
                j = jmap.get((min(a, b), max(a, b)), 0.0)
                if abs(j) < 1e-6:
                    continue
                dnu = abs(ppm_to_hz(by_id[a].shift - by_id[b].shift, freq_mhz))
                if dnu / abs(j) < STRONG_COUPLING_RATIO:
                    strong = True
                    break
            if strong:
                break

        total_spins = sum(by_id[g].n_nuclei for g in comp)
        spin_group: list[int] = []
        shifts: list[float] = []
        for g in comp:
            s = by_id[g]
            spin_group.extend([g] * s.n_nuclei)
            shifts.extend([s.shift] * s.n_nuclei)
        n = len(shifts)
        jmat = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for k in range(i + 1, n):
                ga, gb = spin_group[i], spin_group[k]
                if ga == gb:
                    continue                      # 磁等价核之间无可观测裂分
                j = jmap.get((min(ga, gb), max(ga, gb)), 0.0)
                jmat[i][k] = jmat[k][i] = j

        if strong and total_spins <= max_spins:
            blocks.append(
                SpinBlock(comp, spin_group, shifts, jmat,
                          [by_id[g].n_nuclei for g in comp],
                          exact=True, note="二级（强耦合）精确对角化")
            )
        else:
            note = ""
            if strong:
                note = f"自旋数 {total_spins} 超过上限 {max_spins}，强耦合效应未处理"
                warnings.append(note)
            blocks.append(
                SpinBlock(comp, spin_group, shifts, jmat,
                          [by_id[g].n_nuclei for g in comp],
                          exact=False, note=note)
            )
    return blocks, warnings
