"""把自旋块转成谱线：一级裂分或精确对角化。"""

from __future__ import annotations

from ..core.types import Signal, Transition
from .first_order import first_order_transitions
from .hamiltonian import exact_transitions
from .system import SpinBlock


def block_transitions(block: SpinBlock, signals: list[Signal], freq_mhz: float
                      ) -> list[Transition]:
    by_id = {s.group_id: s for s in signals}

    if block.exact:
        return exact_transitions(
            shifts_ppm=block.shifts,
            j_matrix=block.j,
            freq_mhz=freq_mhz,
            group_ids=block.spin_group,
            group_sizes=[1] * len(block.shifts),
        )

    out: list[Transition] = []
    for gid in block.group_ids:
        out.extend(first_order_transitions(by_id[gid], freq_mhz))
    return out


def all_transitions(
    blocks: list[SpinBlock], signals: list[Signal], freq_mhz: float
) -> list[Transition]:
    out: list[Transition] = []
    for block in blocks:
        out.extend(block_transitions(block, signals, freq_mhz))
    return out
