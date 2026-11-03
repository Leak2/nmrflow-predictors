"""Tier B 化学位移：用 NMRShiftDB2 的 HOSE 索引覆盖 Tier A 的结果。

定位
----
Tier A 是"经验加和 + 环境 SMARTS 查表"，它的问题不是平均精度，而是**结构性盲区**：
多环芳烃、全取代杂环、N-甲基分档这类环境在规则表里没有对应条目，或条目过粗。
HOSE 恰好补这一块——它不依赖任何人为规则，只看"数据库中同环境的原子曾经是多少"。

为什么 Tier B 不是 `ShiftPredictor` 的完整替代
--------------------------------------------
因为它需要 Tier A 的结果才能工作：
1. 球层不足时要退回 Tier A（HOSE 1 层球的 ¹H MAE 0.64 ppm，比 Tier A 还差）；
2. 非对映异位的一对氢在本实现里码相同（没有立体化学扩展），必须借 Tier A 的
   相对位移差把两者拉开；
3. 活泼氢的位移由溶剂决定，库内 78% 记录未标溶剂，混起来会污染。
因此它是一个**覆盖层**：``predict`` 只负责给出候选值，``overlay_tier_b`` 负责合并。

置信度
------
以**球层**为主，n=1 再降一级（单条记录、σ 无从谈起）。实测（``tools/eval_hose_index.py``，留出法）：
邻居数 n 与精度无关甚至负相关（n≥16 时 ¹³C MAE 3.90，n=1 时 2.47），
因为 n 大意味着环境大众化——``CH3-4`` 这种 1 层球码能命中几千条，恰恰最不特异。
把 n 当“越大越可信”会系统性地把最差的预测标成最可信；只对 n=1 做折扣，
是因为另一个极端（单条记录）同样不可信。
"""

from __future__ import annotations

from dataclasses import dataclass

from ..core.types import Confidence, Provenance, Signal
from ..mol.equivalence import SpinGroup
from ..mol.labile import find_labile
from ..mol.prepare import PreparedMol
from .base import ShiftPredictor, merge_tier
from .hose_index import HoseIndex

# 留出法实测的 MAE（ppm），来自 tools/eval_hose_index.py 的分层统计。
# 既用于给置信度定档，也用于判断"Tier A 与 HOSE 差多少才算异常"。
_MAE = {
    "1H": {1: 0.64, 2: 0.37, 3: 0.24, 4: 0.13},
    "13C": {1: 6.96, 2: 3.33, 3: 2.04, 4: 1.21},
    # 19F/31P：tools/eval_tier_ab.py --nuclei 31P,19F，留出 20% 化合物，
    # 只统计 Tier A 与 HOSE 都能给出值的原子。
    "19F": {1: 52.31, 2: 10.73, 3: 4.26, 4: 2.21},
    "31P": {1: 18.85, 2: 8.44, 3: 6.66, 4: 1.97},
}


def measured_mae(nucleus: str, level: int) -> float | None:
    """该核在给定球层上的实测 MAE；未标定的核返回 None。"""
    table = _MAE.get(nucleus)
    if not table:
        return None
    return table.get(level)


def confidence_for(
    nucleus: str, level: int, n: int | None = None
) -> Confidence:
    """球层 + 同码支持数 -> 置信度。

    ¹H 4 层球 0.13 ppm 算 high；¹³C 4 层球 1.21 ppm 只是常态水平，
    与 Tier A 的 medium 同档，因此封顶 medium——把这个值标成 high 是过度承诺。
    ¹⁹F/³¹P 现在已做留出法验证（tools/eval_tier_ab.py --nuclei 31P,19F），
    4 层球 2.21 / 1.97 ppm 相对这两个核的跨度很小，但绝对误差仍是 2 ppm 量级，
    而且尾部很重（P90 21 / 27 ppm，最大 430 / 289），所以同样封顶 medium。

    ``n`` 是同码记录数。球层是主信号，但 **n=1 的 4 层球命中只是“一条记录”**，
    不应当和 n=10 的同样球层同级；实测里 n=1 的尾部明显更重。这里对 n<=1
    降一级，作为“单条记录”的折扣（球层 1~2 本来就是 low，保持不变）。
    """
    table = _CONFIDENCE.get(nucleus)
    if table is not None:
        conf = table.get(level, Confidence.LOW)
    else:
        conf = Confidence.MEDIUM if level >= 3 else Confidence.LOW
    if n is not None and n <= 1:
        conf = conf.downgrade(1)
    return conf


_CONFIDENCE: dict[str, dict[int, Confidence]] = {
    "1H": {
        4: Confidence.HIGH,
        3: Confidence.MEDIUM,
        2: Confidence.MEDIUM,   # 实测 0.37 ppm，仍优于 Tier A 的 0.48
        1: Confidence.LOW,      # 实测 0.64 ppm，不如 Tier A，仅作兜底
    },
    "13C": {
        4: Confidence.MEDIUM,
        3: Confidence.MEDIUM,
        2: Confidence.LOW,
        1: Confidence.LOW,
    },
    # ¹⁹F/³¹P：3 层球 4.26 / 6.66 ppm 仍明显优于 Tier A 有规则子集的 15.96 / 14.10，
    # 记 medium；1~2 层球（52.31 / 10.73、18.85 / 8.44）只当"比兜底常数强"的兜底，记 low。
    "19F": {
        4: Confidence.MEDIUM,
        3: Confidence.MEDIUM,
        2: Confidence.LOW,
        1: Confidence.LOW,
    },
    "31P": {
        4: Confidence.MEDIUM,
        3: Confidence.MEDIUM,
        2: Confidence.LOW,
        1: Confidence.LOW,
    },
}

#: 各核的最低可接受球层。低于它保留 Tier A。
#: 依据：留出法实测"整条合并链路的 MAE/P90 是否变好"——注意判据是**合并后**，
#: 不是"本层 MAE 是否低于 Tier A"。门抬得越高，被拒的原子越多，它们退回 Tier A
#: 反而更差：³¹P 合并 MAE 在 ≥1/2/3/4 球下是 19.71 / 28.05 / 35.84 / 37.95，
#: ¹⁹F 是 8.76 / 10.11 / 13.94 / 19.51，**都是 ≥1 球最好**。
#: ¹⁹F 例外，取 2：它 1 层球的 MAE（全样本 52.31、兜底子集 96.06）比它要替换的东西
#: 还差——Tier A 兜底常数是 48.26、有规则子集是 15.96，两头都输。2 层球起才划算
#: （兜底子集 20.26、有规则子集 10.73）。³¹P 的 1 层球（兜底子集 22.01 vs 兜底 61.31）
#: 明显有用，所以仍取 1。
DEFAULT_MIN_LEVEL = {"1H": 1, "13C": 1, "19F": 2, "31P": 1}


@dataclass
class HoseConfig:
    """覆盖策略。默认值由留出法实测选出，见 tools/eval_tier_ab.py。"""

    max_spheres: int = 4
    #: 低于该球层视为未命中，保留 Tier A。默认值因核而异，见 ``DEFAULT_MIN_LEVEL``——
    #: 四个核实测最优都是 1，所以直接构造 ``HoseConfig()`` 也拿到 1；
    #: 保留 ``default_config()`` 是为了将来某个核需要单独收紧时有地方改。
    min_level: int = 1
    min_neighbors: int = 1
    #: Tier A 已用 SMARTS 锚定到具体官能团（confidence=high）时是否保护它。
    #: ¹³C/¹H 实测 **False 更好**：¹³C 合并 MAE ≥1球 3.41 vs ≥1球+护锚定 3.68
    #: （tools/eval_tier_ab.py，留出 20%）。
    #: ³¹P 例外，默认 True：磷酸三甲酯的 curated 0.0 比 2 层球 HOSE 的 14.35
    #: 更接近实测 2.7；浅层 HOSE 不该覆盖已锚定的高置信度表项。
    protect_curated: bool = False
    #: 同环境记录的离散度超过 ``k × 该核最深球层实测 MAE`` 时，降一级置信度并写明。
    #: **不能用它来否决覆盖**——实测（tools/eval_tier_ab.py，留出法）否决会更差。
    #: 当初的数字是 13C MAE 3.47 → 4.45、P90 8.70 → 12.20，因为被拒的原子只能退回
    #: Tier A，而 Tier A 的 P90 是 19.20 ppm。**这些数字早于后来的 Tier A 修正**
    #: （现在 Tier A 是 6.94 / P90 16.40），结论方向不变但幅度已收窄。
    #: 离散只说明"别太相信"，不说明"别用"。
    sigma_warn_ratio: float | None = 2.5
    #: 与 Tier A 相差超过 k 倍实测 MAE 时**只写提示、不降 HOSE 置信度**。
    #: 理由：分歧是“两条证据不一致”，不是“HOSE 一定错”；咖啡因的羰基/脒碳
    #: 就是 Tier A 错得离谱、HOSE 正确，旧逻辑反而把正确值降成 low。
    #: 分歧仍会出现在 provenance.note 和 report warnings 里，提示人工复核。
    disagree_sigmas: float = 3.0


def default_config(nucleus: str) -> HoseConfig:
    """按核给出默认策略。

    ³¹P 额外打开 ``protect_curated``：磷酸三甲酯这类高置信度环境表项
    （0.0 ppm）比 2 层球 HOSE 的 14.35 ppm 更接近实测，不能被浅层命中覆盖。
    """
    return HoseConfig(
        min_level=DEFAULT_MIN_LEVEL.get(nucleus, 3),
        protect_curated=(nucleus == "31P"),
    )


@dataclass
class OverlayStats:
    """一次覆盖的统计，用于在报告里说明"这条预测到底来自哪一层"。"""

    total: int = 0
    overridden: int = 0
    skipped_labile: int = 0
    skipped_level: int = 0
    wide_sigma: int = 0
    skipped_curated: int = 0
    skipped_no_hit: int = 0
    disagreed: int = 0
    levels: dict[int, int] | None = None

    def __post_init__(self) -> None:
        if self.levels is None:
            self.levels = {}

    def as_note(self, nucleus: str) -> str:
        if not self.overridden:
            return ""
        dist = "、".join(f"{lv} 层 {n} 个" for lv, n in sorted(self.levels.items()))
        text = f"{nucleus} 有 {self.overridden}/{self.total} 个信号改用 HOSE 查表（{dist}）"
        if self.disagreed:
            text += f"，其中 {self.disagreed} 个与经验规则相差较大，建议人工复核"
        return text


class HoseShiftPredictor(ShiftPredictor):
    """逐组给出 Tier B 候选值；未达标的组返回 ``fallback=True`` 的占位信号。"""

    tier = "B"

    def __init__(
        self,
        nucleus: str,
        index: HoseIndex,
        config: HoseConfig | None = None,
    ) -> None:
        self.nucleus = nucleus
        self.index = index
        self.config = config or HoseConfig()

    # ------------------------------------------------------------------
    def predict(
        self, prep: PreparedMol, groups: list[SpinGroup], solvent_key: str = ""
    ) -> list[Signal]:
        """返回与 groups 一一对应的候选信号。

        命中可用时：``shift`` 为库内均值，``provenance.tier='B'``，``fallback=False``。
        不可用时：``fallback=True``，其余字段无意义——``overlay_tier_b`` 会跳过它。

        ``solvent_key`` 不参与计算：库内 78% 记录未标溶剂，按溶剂过滤会把覆盖率
        从 98% 打到个位数。溶剂带来的系统偏差已经包含在实测 MAE 里。
        """
        cfg = self.config
        from .hose import hose_codes  # 局部导入，避免顶层循环

        # 索引构建时用的是**不含显式氢**的 molfile，查询必须同源，
        # 否则 `GetNumImplicitHs()` 与 `GetDegree()` 的口径不一致，码永远对不上。
        codes = hose_codes(prep.heavy, cfg.max_spheres)
        labile = find_labile(prep.mol) if self.nucleus == "1H" else {}

        out: list[Signal] = []
        for group in groups:
            heavy_idx = group.heavy_indices[0] if group.heavy_indices else group.atom_indices[0]
            sig = Signal(
                nucleus=self.nucleus,
                group_id=group.group_id,
                atom_indices=list(group.atom_indices),
                heavy_indices=list(group.heavy_indices),
                n_nuclei=group.n_nuclei,
                provenance=Provenance(tier="B", table="nmrshiftdb2.hose"),
            )

            # 活泼氢：位移由溶剂与交换速率决定，库内混合溶剂会把 0.5~5 ppm 的
            # 各种取值平均成一个没有任何物理意义的中间值。
            if any(h in labile for h in group.atom_indices):
                sig.provenance.fallback = True
                sig.provenance.note = "活泼氢：不走 HOSE"
                out.append(sig)
                continue

            min_level = cfg.min_level

            code = codes.get(heavy_idx)
            if code is None:
                sig.provenance.fallback = True
                sig.provenance.note = "无 HOSE 码"
                out.append(sig)
                continue

            hit = self.index.lookup(
                code,
                self.nucleus,
                max_level=cfg.max_spheres,
                min_level=min_level,
                min_neighbors=cfg.min_neighbors,
            )
            if hit is None:
                sig.provenance.fallback = True
                sig.provenance.note = "库内无同环境记录"
                out.append(sig)
                continue

            sig.shift = hit.mean
            sig.confidence = confidence_for(self.nucleus, hit.level, hit.n)
            sig.provenance.key = f"{hit.level}球/{hit.n}"
            sig.provenance.note = hit.note
            # n < 2 时 std 恒为 0（sum_sq/n − mean² 在单点上是精确的 0），
            # 但那是"只有一条记录、无从谈起离散"，不是"库内高度一致"。
            # 报 0.0 会写进 JSON 让下游误读，所以单条记录一律给 None。
            # 注意：这不改变任何判定——σ 警告本来就不会在 σ=0 时触发。
            sig.provenance.sigma = hit.std if hit.n >= 2 else None
            out.append(sig)
        return out


def _curated(sig: Signal) -> bool:
    """Tier A 是否已把该信号锚定到具体官能团（而不是加和/兜底）。"""
    return (
        sig.provenance.tier == "A"
        and sig.confidence is Confidence.HIGH
        and not sig.provenance.fallback
    )


def _level_of(sig: Signal) -> int:
    """从 ``provenance.key``（形如 ``"4球/12"``）取出命中球层。"""
    key = sig.provenance.key
    head = key.split("球", 1)[0]
    try:
        return int(head)
    except ValueError:
        return 0


def _sigma_limit(cfg: HoseConfig, nucleus: str) -> float | None:
    """超过该离散度就降置信度；未标定的核或显式关闭时返回 None（不检查）。"""
    if cfg.sigma_warn_ratio is None:
        return None
    best = measured_mae(nucleus, cfg.max_spheres)
    return None if best is None else cfg.sigma_warn_ratio * best


def overlay_tier_b(
    base: list[Signal],
    candidates: list[Signal],
    groups: list[SpinGroup],
    nucleus: str,
    config: HoseConfig | None = None,
) -> OverlayStats:
    """把 Tier B 候选值合并进 Tier A 的 ``base``（就地修改），返回统计。

    合并规则
    --------
    - 候选是 fallback（球层不足 / 无码 / 活泼氢）→ 保留 Tier A
    - 球层低于 ``min_level`` → 保留 Tier A
    - Tier A 已是 SMARTS 锚定的 high，而候选不足 4 层球 → 保留 Tier A
    - 否则覆盖位移，置信度取候选值（由球层定档，那是唯一经留出法标定的量），
      并在库内记录离散过大时再降一级
    - 非对映异位的一对氢：候选码相同（本实现无立体化学扩展），
      把 Tier A 已有的相对位移差平移过来，保留分裂但不采信其绝对值
    """
    cfg = config or HoseConfig()
    by_group = {g.group_id: g for g in groups}
    stats = OverlayStats(total=len(base))

    # 先处理非对映异位对：同一重原子上的两个组共享一个候选值
    pair_shift = _diastereotopic_adjust(base, candidates, groups)

    for sig, cand in zip(base, candidates):
        group = by_group.get(sig.group_id)

        if cand.provenance.fallback:
            note = cand.provenance.note
            if "活泼氢" in note:
                stats.skipped_labile += 1
            elif "球层" in note or "环境" in note:
                stats.skipped_level += 1
            else:
                stats.skipped_no_hit += 1
            continue

        # 球层门槛在合并层再判一次：候选可能由另一个配置算出，
        # 也可能被调用方直接喂进来，不能假设它已经按本配置过滤过。
        level = _level_of(cand)
        if level < cfg.min_level:
            stats.skipped_level += 1
            continue

        if cfg.protect_curated and _curated(sig) and level < 4:
            stats.skipped_curated += 1
            continue

        old_shift = sig.shift
        merge_tier(sig, cand)
        # 置信度直接采用候选值：球层是唯一经过留出法标定的量，
        # 取 min 会让 Tier A 的 fallback(LOW) 挡住一条 4 层球的高质量命中。

        # 离散度只降置信度、不否决覆盖：库内记录自己分歧大，说明这条均值不宜轻信，
        # 但退回 Tier A 通常更差（见 HoseConfig.sigma_warn_ratio 的注释）。
        sigma = cand.provenance.sigma
        limit = _sigma_limit(cfg, nucleus)
        if sigma is not None and limit is not None and sigma > limit:
            sig.confidence = sig.confidence.downgrade(1)
            sig.provenance.note += (
                f"；同环境记录离散较大（σ={sigma:.2f} ppm），宜人工复核"
            )
            stats.wide_sigma += 1

        pair = group is not None and group.group_id in pair_shift
        if pair:
            sig.shift = pair_shift[group.group_id]
            sig.confidence = sig.confidence.downgrade(1)
            sig.provenance.note += "；非对映异位：保留经验规则的相对位移差"

        stats.levels[level] = stats.levels.get(level, 0) + 1
        stats.overridden += 1

        if pair:
            continue        # 已按非对映异位处理并降级，不再叠加差异告警

        mae = measured_mae(nucleus, level)
        if mae and abs(sig.shift - old_shift) > cfg.disagree_sigmas * mae:
            # 只提示，不降置信度：Tier A 错得离谱时（咖啡因羰基/脒碳），
            # 降级会把正确的 HOSE 值标成 low。分歧本身不等于 HOSE 错。
            sig.provenance.note += (
                f"；经验规则给 {old_shift:.2f}，相差 {abs(sig.shift - old_shift):.2f} ppm，"
                f"超过 {cfg.disagree_sigmas:.0f}×实测 MAE，建议人工复核"
            )
            stats.disagreed += 1

    return stats


def _diastereotopic_adjust(
    base: list[Signal], candidates: list[Signal], groups: list[SpinGroup]
) -> dict[int, float]:
    """非对映异位对：用 HOSE 值当中心，套上 Tier A 的相对位移差。"""
    by_group = {g.group_id: g for g in groups}
    buckets: dict[int, list[Signal]] = {}
    for sig in base:
        group = by_group.get(sig.group_id)
        if group is None or group.diastereotopic_pair < 0:
            continue
        heavy = group.heavy_indices[0] if group.heavy_indices else -1
        buckets.setdefault(heavy, []).append(sig)

    cand_by_group = {c.group_id: c for c in candidates}
    out: dict[int, float] = {}
    for sigs in buckets.values():
        usable = []
        for s in sigs:
            cand = cand_by_group.get(s.group_id)
            if cand is not None and not cand.provenance.fallback:
                usable.append(s)
        if len(usable) < 2:
            continue
        base_mean = sum(s.shift for s in usable) / len(usable)
        cand_mean = sum(cand_by_group[s.group_id].shift for s in usable) / len(usable)
        for s in usable:
            out[s.group_id] = cand_mean + (s.shift - base_mean)
    return out
