"""单一编排点：SMILES -> 谱图 + 报道值 + 质谱。

中间态可以 dump 成 JSON，调试新化合物时这是救命的。
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem

from . import __version__
from .core.constants import MAX_EXACT_SPINS, NUCLEUS_ORDER, POINTS_PER_PPM
from .core.solvent import get_solvent
from .core.types import FullReport, MsResult, NmrBlock, Spectrum1D
from .core.units import freq_of, round_freq
from .mol.equivalence import (
    carbon_groups, hetero_groups, hydrogen_groups, potential_stereocenters,
)
from .mol.prepare import prepare
from .mol.labile import find_labile
from .ms.formula import formula_from_mol
from .ms.report import predict_ms
from .predict.coupling_hh import HHCouplingPredictor
from .predict.coupling_hetero import HeteroCouplingPredictor
from .predict.shift_1h import ProtonShiftPredictor
from .predict.shift_13c import CarbonShiftPredictor
from .predict.shift_hetero import HeteroShiftPredictor
from .predict.shift_hose import (
    HoseShiftPredictor,
    default_config,
    overlay_tier_b,
)
from .predict.hose_index import HoseIndex, resolve_index_path
from .report.formatter import format_block
from .report.label import assign_labels
from .report.merge import build_peaks
from .spin.system import build_blocks
from .spin.transitions import all_transitions
from .spectrum.digitize import digitize, integral_curve

ELEMENT_TO_NUCLEUS = {"F": "19F", "P": "31P"}

# 每个核的 group_id 段。**必须全局唯一**：13C 的 group 0 与 19F 的 group 0
# 若撞号，分块时会被当成同一个自旋，把整个碳骨架误判成一个强耦合体系。
GROUP_ID_BASE = {"1H": 0, "13C": 100_000, "19F": 200_000, "31P": 300_000}


@dataclass
class RunOptions:
    solvent: str = "CDCl3"
    proton_mhz: float = 400.0
    nucleus_digits: dict[str, int] = field(default_factory=dict)
    ms_method: str = "DART"
    ms_adducts: list[str] | None = None
    with_ms: bool = True
    with_isotopes: bool = True
    decouple: dict[str, tuple[str, ...]] = field(default_factory=dict)
    points_per_ppm: float = POINTS_PER_PPM
    noise: float = 0.0
    lineshape: str = "lorentzian"
    max_exact_spins: int = MAX_EXACT_SPINS
    all_fragments: bool = True         # 盐类：合并所有含碳碎片，而非只取最大的
    tier_b: bool = True                # 用 HOSE 索引覆盖 1H/13C；无索引则自动跳过
    #: 19F/31P 也走 HOSE。已做留出法验证（tools/eval_tier_ab.py --nuclei 31P,19F）：
    #: 合并后 ³¹P MAE 40.40→19.71、¹⁹F 23.62→8.76。主要收益来自 Tier A 的兜底常数
    #: （³¹P 有 56% 的原子退到固定 0.0 ppm），HOSE 把它们换成同环境的实测均值。
    tier_b_hetero: bool = True
    hose_index: str | None = None      # 显式指定索引路径
    hose_min_level: dict[str, int] = field(default_factory=dict)  # 每核覆盖球层门槛


@dataclass
class RunResult:
    report: FullReport
    spectra: dict[str, Spectrum1D]
    ms: MsResult | None
    prepared: object
    signals: dict[str, list]
    options: RunOptions


def detect_nuclei(formula: dict[str, int]) -> list[str]:
    """按元素决定要出哪些核的谱。

    只有分子里真的有 H / C 时才出 ¹H / ¹³C。否则水、氨、四氯化碳会得到
    ``13C NMR ... δ:`` 这种空块——空块不是“预测”，只是排版噪声。
    """
    out = []
    if formula.get("H", 0) > 0:
        out.append("1H")
    if formula.get("C", 0) > 0:
        out.append("13C")
    for element, nucleus in ELEMENT_TO_NUCLEUS.items():
        if formula.get(element, 0) > 0:
            out.append(nucleus)
    return [n for n in NUCLEUS_ORDER if n in out]


def run(smiles: str, options: RunOptions | None = None) -> RunResult:
    opt = options or RunOptions()
    prep = prepare(smiles, all_fragments=opt.all_fragments)
    solvent = get_solvent(opt.solvent)

    nuclei = detect_nuclei(prep.formula)

    # --- 1H ---
    signals: dict[str, list] = {}
    groups_by_nucleus: dict[str, list] = {}
    p1h = ProtonShiftPredictor()
    h_groups = hydrogen_groups(prep, start=GROUP_ID_BASE["1H"])
    h_signals = p1h.predict(prep, h_groups, opt.solvent)
    if solvent.drop_exchangeable:
        # 组与信号必须同步过滤，否则 Tier B 的对齐会错位
        kept = [
            (g, s) for g, s in zip(h_groups, h_signals) if not s.is_exchangeable
        ]
        h_groups = [g for g, _ in kept]
        h_signals = [s for _, s in kept]
    assign_labels(prep.mol, h_signals, "1H")
    signals["1H"] = h_signals
    groups_by_nucleus["1H"] = h_groups

    # --- 13C ---
    p13 = CarbonShiftPredictor()
    c_groups = carbon_groups(prep, start=GROUP_ID_BASE["13C"])
    c_signals = p13.predict(prep, c_groups, opt.solvent)
    assign_labels(prep.heavy, c_signals, "13C")
    signals["13C"] = c_signals
    groups_by_nucleus["13C"] = c_groups

    # --- 19F / 31P ---
    labile_map = find_labile(prep.mol)
    if "19F" in nuclei:
        f_groups = hetero_groups(prep, "F", start=GROUP_ID_BASE["19F"])
        f_signals = HeteroShiftPredictor("19F").predict(prep, f_groups, opt.solvent)
        for s in f_signals:
            s.label = "-F"
        signals["19F"] = f_signals
        groups_by_nucleus["19F"] = f_groups
    if "31P" in nuclei:
        p_groups = hetero_groups(prep, "P", start=GROUP_ID_BASE["31P"])
        p_signals = HeteroShiftPredictor("31P").predict(prep, p_groups, opt.solvent)
        for s in p_signals:
            s.label = "-P"
        signals["31P"] = p_signals
        groups_by_nucleus["31P"] = p_groups

    # --- Tier B：HOSE 查表覆盖 ---
    tier_b_notes = _apply_tier_b(opt, prep, signals, groups_by_nucleus)

    # --- 耦合 ---
    HHCouplingPredictor().predict(prep, h_signals)
    HeteroCouplingPredictor().apply(prep, signals, decouple=opt.decouple or None)

    # --- 谱图 ---
    spectra: dict[str, Spectrum1D] = {}
    blocks = []
    block_warnings: list[str] = []
    for nucleus in nuclei:
        freq = freq_of(nucleus, opt.proton_mhz)
        sigs = signals[nucleus]
        # 注意：19F/31P 查表给出的 width_ppm 是**化学位移不确定度**
        # （环境跨度），不是谱线宽度，绝不能拿来当 FWHM 用。
        # 线宽统一交给 digitize 按核素默认值 + 置信度处理。
        blk, spin_warnings = build_blocks(
            sigs, freq, max_spins=opt.max_exact_spins
        )
        block_warnings.extend(f"{nucleus}: {w}" for w in spin_warnings)
        if nucleus == "19F":
            block_warnings.extend(_diastereotopic_f_warnings(prep, sigs))
        trans = all_transitions(blk, sigs, freq)
        spec = digitize(
            trans, sigs, nucleus, freq, solvent.display,
            points_per_ppm=opt.points_per_ppm,
            shape=opt.lineshape, noise=opt.noise,
        )
        if nucleus == "1H":
            spec.integral = integral_curve(spec)
        spec.warnings.extend(_nucleus_warnings(sigs))
        spec.warnings.extend(_degenerate_shift_warnings(nucleus, sigs))
        spec.warnings.extend(spin_warnings)
        spectra[nucleus] = spec

        peaks = build_peaks(sigs, trans, freq, nucleus)
        block = NmrBlock(
            nucleus=nucleus, freq_mhz=round_freq(nucleus, opt.proton_mhz),
            solvent=solvent.display, peaks=peaks,
        )
        block.text = format_block(
            block,
            nucleus_digits=opt.nucleus_digits or None,
        )
        blocks.append(block)

    # --- 质谱 ---
    ms = None
    extra_ms: list = []
    if opt.with_ms:
        species = _ms_species(opt, prep, smiles)
        if species:
            # 离子化合物/带电物种：每个含碳的带电碎片单独报 [M]+ / [M]-。
            # 不能把阴离子母体再加一个质子——乙酸根 59.013 会被算成乙酸 60.021，
            # 四甲基铵 74.096 会被算成 75.104。
            for frag_counts, charge, label in species:
                name = "[M]+" if charge > 0 else "[M]-"
                one = predict_ms(
                    frag_counts,
                    method=opt.ms_method,
                    adduct_names=[name],
                    with_isotopes=opt.with_isotopes and not extra_ms,
                )
                one.formula = one.formula + ("+" if charge > 0 else "-")
                if ms is None:
                    ms = one
                else:
                    extra_ms.append(one)
            prep.warnings.append(
                "离子化合物/带电物种：HRMS 按 [M]+ / [M]- 报告，未加质子；"
                "不含碳的对离子已剔除"
            )
        else:
            counts = formula_from_mol(prep.mol)
            ms = predict_ms(
                counts,
                method=opt.ms_method,
                adduct_names=opt.ms_adducts,
                with_isotopes=opt.with_isotopes,
            )

    report = FullReport(
        smiles=smiles, solvent=solvent.display,
        proton_mhz=opt.proton_mhz,
        formula=ms.formula if ms else "",
        blocks=blocks, ms=ms, extra_ms=extra_ms,
        warnings=list(prep.warnings) + block_warnings + tier_b_notes,
        # 只要有任何一个信号被索引改写过就要附来源声明（许可 4.3）
        hose_used=any(
            sig.provenance.tier == "B"
            for sigs in signals.values() for sig in sigs
        ),
    )
    return RunResult(report, spectra, ms, prep, signals, opt)


def _ms_species(opt, prep, smiles: str) -> list[tuple[dict[str, int], int, str]]:
    """返回需要按完整离子式报告 [M]+/[M]- 的物种。

    优先用原始碎片里“带电且含碳”的那些：盐的阳离子/阴离子各报一条，
    Na+/Cl-/BF4-/PF6- 这类不含碳的对离子被排除。单分子没有碎片信息时，
    若母体净电荷非零，就把母体本身当作一个离子（例如乙酸根、四甲基铵）。
    """
    if opt.ms_adducts:
        return []        # 用户显式指定加合物时，尊重用户选择
    if opt.all_fragments and prep.fragments:
        ions = [
            (f.counts, f.charge, f.smiles)
            for f in prep.fragments
            if f.charge != 0 and f.counts.get("C", 0) > 0
        ]
        if ions:
            return ions
    charge = Chem.GetFormalCharge(prep.mol)
    if charge != 0:
        return [(formula_from_mol(prep.mol), charge, smiles)]
    return []


def _apply_tier_b(opt, prep, signals: dict, groups_by_nucleus: dict) -> list[str]:
    """用 HOSE 索引覆盖 1H/13C 的位移，返回给报告的说明。

    索引不存在时静默跳过（便携包里索引是可选的），不报错——
    否则用户拿到一个没有 build/ 的源码树，连最基本的谱都出不来。
    """
    if not opt.tier_b:
        return []
    path = resolve_index_path(opt.hose_index)
    if path is None:
        return []

    nuclei = ["1H", "13C"] + (["19F", "31P"] if opt.tier_b_hetero else [])
    notes: list[str] = []
    with HoseIndex(path) as index:
        for nucleus in nuclei:
            sigs = signals.get(nucleus)
            groups = groups_by_nucleus.get(nucleus)
            if not sigs or not groups:
                continue
            cfg = default_config(nucleus)
            if nucleus in opt.hose_min_level:
                cfg.min_level = opt.hose_min_level[nucleus]
            candidates = HoseShiftPredictor(nucleus, index, cfg).predict(
                prep, groups, opt.solvent
            )
            stats = overlay_tier_b(sigs, candidates, groups, nucleus, cfg)
            note = stats.as_note(nucleus)
            if note:
                notes.append(note)
    return notes


def _degenerate_shift_warnings(nucleus, signals) -> list[str]:
    """异核查表自检：两个化学不等价的环境被查成同一个位移 → 规则没区分开。

    19F/31P 走的是"环境 SMARTS -> 值"的查表，一旦某条规则的 SMARTS 写得过松
    （例如把"链内 CF2"写成只要有一个碳邻居），不同环境会一起命中它，得到完全
    相同的位移，谱图上就叠成一个峰，看起来像"这个分子只有一个氟环境"。
    这种情况必须显式告警，而不是安静地输出一张漂亮的单峰图。
    """
    if nucleus == "1H":
        return []
    buckets: dict[float, list] = {}
    for sig in signals:
        if sig.provenance.fallback:      # 兜底值本来就都是同一个，不重复告警
            continue
        buckets.setdefault(round(sig.shift, 3), []).append(sig)
    out: list[str] = []
    for shift, group in sorted(buckets.items()):
        if len(group) < 2:
            continue
        keys = "、".join(sorted({s.provenance.key for s in group}))
        nuclei = sum(s.n_nuclei for s in group)
        out.append(
            f"{len(group)} 个化学不等价的环境被查成同一化学位移 "
            f"δ {shift:.2f}（共 {nuclei} 个核，规则「{keys}」）："
            f"查表按环境类型给值，同类型环境共用一个位移，谱图上会叠成一个峰；"
            f"真实谱中它们通常是一组靠得很近的多重峰"
        )
    return out


def _diastereotopic_f_warnings(prep, signals) -> list[str]:
    """CF₂ 紧邻手性中心时，两个 F 非对映异位——真实谱是 AB，不是 2F 单峰。

    本项目的 HOSE 实现没有立体化学扩展（见 ``predict/hose.py`` 的"已知偏离"第 2 条），
    等价类按拓扑对称合并，所以这两个 F 会被合成一条 2F 线。位移本身（中心）是对的，
    但多重度会少掉一个 AB 分裂——实测 CF₂Cl–CFClI 是 −61.82 / −64.76（Δδ 2.94 ppm，
    ²J = 165 Hz）。与其编一个 Δδ 硬拆，不如把这件事写在提示里。
    """
    centers = potential_stereocenters(prep.mol)
    if not centers:
        return []
    out: list[str] = []
    for sig in signals:
        if sig.n_nuclei < 2:
            continue
        carbons = set()
        for idx in sig.atom_indices:
            atom = prep.mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 9:
                continue
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 6:
                    carbons.add(nb.GetIdx())
        for c in carbons:
            atom_c = prep.mol.GetAtomWithIdx(c)
            if any(nb.GetIdx() in centers for nb in atom_c.GetNeighbors()):
                out.append(
                    f"δ {sig.shift:.2f} 的 {sig.n_nuclei} 个 F 同处一个碳，该碳紧邻手性中心"
                    f" → 两个 F 非对映异位，真实谱应为 AB 分裂（同类环境实测 Δδ≈2.9 ppm、"
                    f"²J≈165 Hz；本模型无立体化学扩展，按单峰输出）"
                )
                break
    return out


def _nucleus_warnings(signals) -> list[str]:
    out = []
    lows = sum(1 for s in signals if s.confidence.value == "low")
    if lows:
        out.append(f"{lows} 个信号置信度为 low（未命中经验表或为活泼氢）")
    return out


def dump_json(result: RunResult, path: str | Path) -> Path:
    """把中间态 dump 成 JSON，便于调试新化合物。"""
    path = Path(path)
    payload = {
        "schema_version": 1,
        "data_type": "predicted",
        "generator": {"name": "nmrflow", "version": __version__},
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {
            "smiles": result.report.smiles,
            "solvent": result.report.solvent,
            "proton_mhz": result.report.proton_mhz,
        },
        "options": {
            "tier_b": result.options.tier_b,
            "tier_b_hetero": result.options.tier_b_hetero,
            "with_ms": result.options.with_ms,
            "with_isotopes": result.options.with_isotopes,
            "ms_method": result.options.ms_method,
            "points_per_ppm": result.options.points_per_ppm,
            "noise": result.options.noise,
            "lineshape": result.options.lineshape,
        },
        "smiles": result.report.smiles,
        "solvent": result.report.solvent,
        "proton_mhz": result.report.proton_mhz,
        "formula": result.report.formula,
        "warnings": result.report.warnings,
        # 用到了 HOSE 索引就得带上来源声明（nmrshiftdb2 Database License 4.3）
        "attribution": result.report.attribution or None,
        "signals": {
            nuc: [
                {
                    "group_id": s.group_id,
                    "shift": round(s.shift, 3),
                    "n": s.n_nuclei,
                    "label": s.label,
                    "confidence": s.confidence.value,
                    "provenance": str(s.provenance),
                    "tier": s.provenance.tier,
                    "sigma": s.provenance.sigma,
                    "note": s.provenance.note,
                    "couplings": [
                        {"partner": c.partner_group, "J": round(c.j_hz, 2),
                         "n": c.n_equiv, "kind": c.kind}
                        for c in s.couplings
                    ],
                }
                for s in sigs
            ]
            for nuc, sigs in result.signals.items()
        },
        "blocks": [
            {"nucleus": b.nucleus, "freq_mhz": b.freq_mhz, "text": b.text}
            for b in result.report.blocks
        ],
        "ms": {
            "formula": result.ms.formula if result.ms else "",
            "neutral_mass": result.ms.neutral_mass if result.ms else 0.0,
            "adducts": [
                {"adduct": a.adduct, "calculated": a.calculated}
                for a in (result.ms.adducts if result.ms else [])
            ],
            "isotope_pattern": [
                {"mz": round(p.mz, 4), "rel": round(p.abundance, 4),
                 "ordinal": p.ordinal}
                for p in (result.ms.isotope_pattern if result.ms else [])
            ],
        } if result.ms else None,
        "extra_ms": [
            {
                "formula": extra.formula,
                "adducts": [
                    {"adduct": a.adduct, "calculated": a.calculated}
                    for a in extra.adducts
                ],
            }
            for extra in result.report.extra_ms
        ],
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return path


def write_manifest(
    result: RunResult,
    outdir: str | Path,
    outputs: list[str | Path],
) -> Path:
    """在结果目录写 manifest.json：输入、选项、输出版本与哈希。

    不把 manifest 自己列进 outputs（自哈希没有意义）；调用方应在
    所有输出文件写完之后再调用。
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_hashes: dict[str, str] = {}
    for item in outputs:
        path = Path(item)
        if not path.is_file():
            continue
        try:
            rel = path.resolve().relative_to(outdir.resolve()).as_posix()
        except ValueError:
            rel = path.name
        output_hashes[rel] = hashlib.sha256(path.read_bytes()).hexdigest()

    payload = {
        "schema_version": 1,
        "data_type": "predicted",
        "generator": {"name": "nmrflow", "version": __version__},
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {
            "smiles": result.report.smiles,
            "solvent": result.report.solvent,
            "proton_mhz": result.report.proton_mhz,
        },
        "options": {
            "tier_b": result.options.tier_b,
            "tier_b_hetero": result.options.tier_b_hetero,
            "with_ms": result.options.with_ms,
            "with_isotopes": result.options.with_isotopes,
            "ms_method": result.options.ms_method,
            "points_per_ppm": result.options.points_per_ppm,
            "noise": result.options.noise,
            "lineshape": result.options.lineshape,
        },
        "attribution": result.report.attribution or None,
        "outputs": output_hashes,
    }
    path = outdir / "manifest.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                    encoding="utf-8")
    return path
