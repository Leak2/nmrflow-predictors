"""命令行入口：SMILES -> 谱图 PNG + 报道值文本（+ 可选 JSON / 质谱图）。

用法：
    python -m nmrflow.cli "CCO" -o out/ethanol
    python -m nmrflow.cli "FC(F)(F)c1ccccc1" --solvent DMSO-d6 --mhz 600
    python -m nmrflow.cli examples/smiles.txt --outdir results

输出（以 out/ethanol 为前缀）：
    ethanol_1H.png / ethanol_13C.png [ethanol_19F.png] [ethanol_31P.png]
    ethanol_report.txt         报道值（含 HRMS）
    ethanol.json               全部中间态（可选）
    ethanol_MS.png             同位素分布（可选）
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .core.constants import MAX_EXACT_SPINS
from .core.types import PREDICTED_DATA_NOTE_TEMPLATE
from .pipeline import RunOptions, dump_json, run, write_manifest
from .plot.plot_1h import plot_1h
from .plot.plot_generic import plot_generic
from .plot.plot_ms import plot_ms

FILE_SUFFIX = {"1H": "1H", "13C": "13C", "19F": "19F", "31P": "31P"}


def _safe_name(text: str, fallback: str = "mol") -> str:
    """把 SMILES 变成文件名（去掉 / \\ : 等非法字符）。"""
    keep = []
    for ch in text.strip():
        if ch.isalnum() or ch in "-_()[]#=+.,":
            keep.append(ch)
        elif ch in " /\\:":
            keep.append("_")
    name = "".join(keep).strip("_")
    return name[:60] or fallback


def _read_inputs(args: argparse.Namespace) -> list[str]:
    if args.smiles:
        return [args.smiles]
    src = Path(args.file)
    if not src.exists():
        raise SystemExit(f"找不到输入文件：{src}")
    out = []
    for line in src.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def run_one(smiles: str, args: argparse.Namespace) -> list[Path]:
    opt = RunOptions(
        solvent=args.solvent,
        proton_mhz=args.mhz,
        ms_method=args.ms,
        with_ms=not args.no_ms,
        with_isotopes=not args.no_isotopes,
        points_per_ppm=args.points_per_ppm,
        noise=args.noise,
        lineshape=args.lineshape,
        max_exact_spins=args.max_spins,
        tier_b=not args.no_tier_b,
        tier_b_hetero=args.hose_hetero,
        hose_index=args.hose_index,
    )
    result = run(smiles, opt)

    outdir = Path(args.outdir or ".")
    outdir.mkdir(parents=True, exist_ok=True)
    base = args.prefix or _safe_name(smiles)
    written: list[Path] = []

    # --- 谱图 PNG ---
    peaks_of = {b.nucleus: b.peaks for b in result.report.blocks}
    for nucleus, spec in result.spectra.items():
        path = outdir / f"{base}_{FILE_SUFFIX[nucleus]}.png"
        peaks = peaks_of.get(nucleus, [])
        if nucleus == "1H":
            plot_1h(spec, path, show_integral=not args.no_integral, peaks=peaks)
        else:
            plot_generic(spec, path, peaks=peaks)
        written.append(path)

    # --- 质谱图 ---
    if result.ms is not None and not args.no_ms_plot:
        path = outdir / f"{base}_MS.png"
        plot_ms(result.ms, path)
        written.append(path)

    # --- 报道值 ---
    report_path = outdir / f"{base}_report.txt"
    text = result.report.text           # 已含末尾的 HRMS 行
    # 用了 HOSE 索引就附来源声明，并加“这是预测”的页脚。
    tails = [
        PREDICTED_DATA_NOTE_TEMPLATE.format(version=__version__),
    ]
    attribution = result.report.attribution
    if attribution:
        tails.append(attribution)
    report_path.write_text(
        text + ("\n\n" + "\n\n".join(tails) if tails else "") + "\n",
        encoding="utf-8",
    )
    written.append(report_path)

    # --- 中间态 JSON ---
    if args.json:
        written.append(dump_json(result, outdir / f"{base}.json"))

    # --- 机器可读清单（最后写，包含上面所有输出） ---
    written.append(write_manifest(result, outdir, written))

    if not args.quiet:
        print("=" * 72)
        print(smiles)
        print("-" * 72)
        print(text)
        print()
        print(PREDICTED_DATA_NOTE_TEMPLATE.format(version=__version__))
        if result.report.warnings:
            print("\n[提示] " + "；".join(result.report.warnings))
        print("-" * 72)
        for p in written:
            print("  ->", p)
    return written


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="nmrflow",
        description="由 SMILES 预测 1H / 13C（+ 自动检测 19F / 31P）NMR 谱图与报道值",
    )
    p.add_argument("--version", action="version",
                   version=f"nmrflow {__version__}")
    p.add_argument("smiles", nargs="?", help="单个 SMILES")
    p.add_argument("-f", "--file", help="包含若干 SMILES 的文本文件（每行一个）")
    p.add_argument("-o", "--outdir", default="results", help="输出目录")
    p.add_argument("-p", "--prefix", help="输出文件名前缀（默认由 SMILES 生成）")
    p.add_argument("--solvent", default="CDCl3",
                   help="氘代溶剂：CDCl3 / DMSO-d6 / CD3OD / D2O / acetone-d6 / C6D6 / CD3CN")
    p.add_argument("--mhz", type=float, default=400.0, help="质子观测频率（MHz）")
    p.add_argument("--ms", default="DART",
                   help="质谱方法：DART / ESI / APCI / MALDI / CI / EI")
    p.add_argument("--no-ms", action="store_true", help="不做质谱预测")
    p.add_argument("--no-ms-plot", action="store_true", help="不画同位素分布图")
    p.add_argument("--no-isotopes", action="store_true", help="不算同位素分布")
    p.add_argument("--no-integral", action="store_true", help="1H 谱不画积分曲线")
    p.add_argument("--json", action="store_true", help="额外 dump 中间态 JSON")
    p.add_argument("--points-per-ppm", type=float, default=40.0, help="谱图采样密度")
    p.add_argument("--noise", type=float, default=0.0, help="谱图噪声水平（0~1）")
    p.add_argument("--lineshape", choices=["lorentzian", "gaussian", "voigt"],
                   default="lorentzian")
    p.add_argument("--max-spins", type=int, default=MAX_EXACT_SPINS,
                   help="精确对角化的自旋数上限")
    p.add_argument("--no-tier-b", action="store_true",
                   help="不做 HOSE 查表覆盖，只用 Tier A 经验规则")
    p.add_argument("--hose-index", default=None,
                   help="HOSE 索引路径（默认自动查找 nmrflow/data 或 build/）")
    p.add_argument("--no-hose-hetero", dest="hose_hetero", action="store_false",
                   help="19F/31P 不走 HOSE 查表，只用 Tier A 环境表")
    p.add_argument("-q", "--quiet", action="store_true", help="只写文件，不打印")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.smiles and not args.file:
        build_parser().print_help()
        return 2

    ok = 0
    failed = 0
    for smi in _read_inputs(args):
        try:
            run_one(smi, args)
            ok += 1
        except (ValueError, KeyError) as exc:
            failed += 1
            print(f"[跳过] {smi}: {exc}", file=sys.stderr)
    if failed:
        # 部分失败也返回非零：批处理脚本里"3 个成功 1 个失败"不能被当成成功，
        # 否则脚本会继续拿一份不完整的结果往下走。
        print(f"[汇总] 成功 {ok}，失败 {failed}", file=sys.stderr)
        return 1
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
