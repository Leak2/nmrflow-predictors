"""贯穿全流程的数据类型定义。"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np


class Confidence(str, Enum):
    """信号置信度。影响报告里的提示与谱图上的线宽。"""

    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"

    @property
    def rank(self) -> int:
        return {"high": 0, "medium": 1, "low": 2}[self.value]

    def downgrade(self, steps: int = 1) -> "Confidence":
        order = [Confidence.HIGH, Confidence.MEDIUM, Confidence.LOW]
        idx = min(len(order) - 1, self.rank + steps)
        return order[idx]


@dataclass
class Provenance:
    """一个预测值的溯源：哪一层、哪张表、哪条规则给的。"""

    tier: str = "A"                 # A=经验规则 B=指纹检索 C=ML
    table: str = ""                 # 数据文件名
    key: str = ""                   # 命中的表项
    fallback: bool = False          # 是否回落到了通用值
    note: str = ""
    sigma: float | None = None      # Tier B：同环境记录的位移标准差（ppm）

    def __str__(self) -> str:
        tail = " (fallback)" if self.fallback else ""
        base = f"{self.table}:{self.key}" if self.key else self.table
        return f"[Tier {self.tier}] {base}{tail}"


@dataclass
class Coupling:
    """一条 J 耦合：本信号被 partner 组以 j_hz 裂分，partner 有 n_equiv 个等价核。"""

    partner_group: int              # 耦合伙伴自旋组的 id
    j_hz: float
    n_equiv: int = 1                # 伙伴组内等价核个数（决定 n+1 裂分与强度比）
    partner_nucleus: str = "1H"     # 伙伴核素（异核耦合用）
    kind: str = ""                  # "2J" "3J" "4J" "ortho" "meta" "para" 等

    def __post_init__(self) -> None:
        self.j_hz = float(self.j_hz)


@dataclass
class Signal:
    """一个自旋组（化学等价的一组核）的预测结果。"""

    nucleus: str                    # "1H" / "13C" / ...
    group_id: int
    atom_indices: list[int] = field(default_factory=list)
    heavy_indices: list[int] = field(default_factory=list)  # 核所连的重原子
    shift: float = 0.0              # ppm，已相对本核参考物质
    n_nuclei: int = 1               # 该组等价核个数（积分）
    couplings: list[Coupling] = field(default_factory=list)
    fwhm_hz: float | None = None    # None -> 用核默认值
    width_ppm: float | None = None  # 19F/31P 查表得到的是 ppm 跨度，绘制时换算成 Hz
    label: str = ""                 # 归属标签，例如 "-CH3" / "Ph-H" / "-NH"
    confidence: Confidence = Confidence.MEDIUM
    provenance: Provenance = field(default_factory=Provenance)
    is_labile: bool = False         # 活泼氢
    is_exchangeable: bool = False   # 可在 D2O 中消失

    def __post_init__(self) -> None:
        self.shift = float(self.shift)

    @property
    def is_broad(self) -> bool:
        return self.is_labile or self.confidence is Confidence.LOW


@dataclass
class Transition:
    """一条谱线：频率（Hz，相对参考）+ 相对强度。"""

    freq_hz: float
    intensity: float
    group_id: int = -1


@dataclass
class Spectrum1D:
    """一张一维谱：数字化的 (δ, I) 数组 + 积分曲线 + 峰表。"""

    nucleus: str
    freq_mhz: float                 # 该核的观测频率
    solvent: str
    x_ppm: np.ndarray
    y: np.ndarray
    integral: np.ndarray | None = None
    transitions: list[Transition] = field(default_factory=list)
    signals: list[Signal] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ppm_range(self) -> tuple[float, float]:
        return float(self.x_ppm.min()), float(self.x_ppm.max())


@dataclass
class Peak:
    """拾峰结果：一个信号在谱图上呈现的峰（或峰簇）。"""

    shift: float                    # 峰簇中心（或区间中点）
    shift_lo: float                 # 区间下界
    shift_hi: float                 # 区间上界
    intensity: float
    integration: float              # 归一化质子数
    multiplicity: str = "s"         # s/d/t/q/dd/m/br s ...
    j_values: list[float] = field(default_factory=list)
    j_partners: list[str] = field(default_factory=list)  # 每个 J 对应的伙伴核
    label: str = ""
    confidence: Confidence = Confidence.MEDIUM
    n_lines: int = 1                # 峰簇内谱线条数


@dataclass
class AdductResult:
    """一个加合物的质荷比结果。"""

    adduct: str                     # "[M+H]+"
    formula: str                    # 加合物总式（可选）
    calculated: float               # 理论 m/z
    charge: int = 1


@dataclass
class IsotopePeak:
    mz: float
    abundance: float                # 相对强度，基峰 = 1.0
    ordinal: int                    # M=0, M+1=1 ...


@dataclass
class MsResult:
    """质谱预测结果。"""

    formula: str                    # Hill 记法中性式
    neutral_mass: float             # 单同位素中性精确质量
    method: str = "DART"
    adducts: list[AdductResult] = field(default_factory=list)
    isotope_pattern: list[IsotopePeak] = field(default_factory=list)
    rdb: float | None = None        # 不饱和度


@dataclass
class NmrBlock:
    """一个核的报道块（文本 + 峰表）。"""

    nucleus: str
    freq_mhz: float
    solvent: str
    peaks: list[Peak] = field(default_factory=list)
    text: str = ""


#: HOSE 索引的数据来源声明。nmrshiftdb2 Database License 4.3 要求：公开使用由该库
#: 产生的成果时须附带此声明。所以它出现在导出的 report.txt 与 JSON 里，但**不进**
#: ``FullReport.text``——那段文本是给用户粘进论文的，插一行英文声明会碍事。
NMR_SHIFTDB_ATTRIBUTION = (
    "Contains information from nmrshiftdb2 (www.nmrshiftdb.org),\n"
    "which is available here under the nmrshiftdb2 Database License."
)

PREDICTED_DATA_NOTE_TEMPLATE = (
    "Predicted by nmrflow {version} - not experimental data."
)


@dataclass
class FullReport:
    """一次完整运行的全部产物。"""

    smiles: str
    solvent: str
    proton_mhz: float
    formula: str = ""
    blocks: list[NmrBlock] = field(default_factory=list)   # 顺序：1H, 19F, 31P, 13C
    ms: MsResult | None = None
    extra_ms: list[MsResult] = field(default_factory=list)  # 盐的其余离子
    warnings: list[str] = field(default_factory=list)
    #: 本次结果里是否有信号被 HOSE 索引改写过。决定要不要附来源声明。
    hose_used: bool = False

    @property
    def attribution(self) -> str:
        """需要随成果一并提供的来源声明；没用 Tier B 时为空字符串。"""
        return NMR_SHIFTDB_ATTRIBUTION if self.hose_used else ""

    @property
    def text(self) -> str:
        """完整报道值文本（NMR 各块 + 末尾 HRMS）。"""
        parts = [f"{self.smiles}:"]
        for blk in self.blocks:
            parts.append(blk.text)
        if self.ms is not None:
            parts.append(self.ms_text)
        for extra in self.extra_ms:
            parts.append(format_hrms_line(extra))
        return "".join(parts)

    @property
    def ms_text(self) -> str:
        return format_hrms_line(self.ms) if self.ms is not None else ""


def format_hrms_line(ms: "MsResult | None") -> str:
    """生成 HRMS 行，例如：

    HRMS (DART) for C27H26F3N5S m/z: Calculated, 510.1934 [M+H]+.
    """
    if ms is None or not ms.adducts:
        return ""
    chunks = []
    for ad in ms.adducts:
        chunk = f"Calculated, {ad.calculated:.4f} {ad.adduct}"
        chunks.append(chunk)
    line = f"HRMS ({ms.method}) for {ms.formula} m/z: " + ", ".join(chunks)
    # EI 的加合物是自由基阳离子 M⁺·，自带结尾符号，再补句点会变成 "M⁺·."
    return line if line.endswith("·") else line + "."
