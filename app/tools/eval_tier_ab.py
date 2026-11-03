"""Tier A（经验加和）vs Tier B（HOSE）vs 合并结果，同一批分子上逐原子比较。

为什么必须做这个实验
--------------------
HOSE 索引单独评估的 MAE（¹H 0.29、¹³C 3.08）只说明"这个索引本身有多准"，
不能说明"接进来是否让输出变好"。Tier A 的平均误差其实很小，
问题在结构盲区；如果盲目用 HOSE 覆盖，很可能在 Tier A 本来很准的地方把它弄坏。

因此这里在同一批**留出化合物**上，对**同一个原子**同时算三个预测值：

    Tier A      —— 现有经验规则
    HOSE        —— 索引查表（按球层截断）
    合并        —— overlay_tier_b 的最终结果，含球层门槛与 σ 警告

只有"A 和 B 都能给出值"的原子才进入对照，避免覆盖率差异冒充精度差异。
另外按"芳香 / 非芳香"拆分——这一项直接决定了球层门槛是否该分类型设置。

训练索引默认保留在 ``build/hose_index_train.sqlite``，重复实验不必重建（约 6 分钟）。

用法
----
    python tools/eval_tier_ab.py --sd /path/to/nmrshiftdb2withsignals.sd
"""

from __future__ import annotations

import argparse
import copy
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "tools"))

from build_hose_index import CALC_PROPS, SUPPORTED, _atom_matches_nucleus, _iter_shift_entries, build  # noqa: E402
from eval_hose_index import _inchikey, _is_holdout  # noqa: E402
from nmrflow.mol.equivalence import (  # noqa: E402
    carbon_groups,
    hetero_groups,
    hydrogen_groups,
)
from nmrflow.mol.prepare import PreparedMol  # noqa: E402
from nmrflow.predict.hose_index import HoseIndex  # noqa: E402
from nmrflow.predict.shift_13c import CarbonShiftPredictor  # noqa: E402
from nmrflow.predict.shift_1h import ProtonShiftPredictor  # noqa: E402
from nmrflow.predict.shift_13c import CarbonShiftPredictor  # noqa: E402,F811
from nmrflow.predict.shift_hetero import HeteroShiftPredictor  # noqa: E402
from nmrflow.predict.shift_hose import HoseConfig, HoseShiftPredictor, overlay_tier_b  # noqa: E402

SOLVENT = "CDCl3"

CONFIGS = [
    ("≥1球", HoseConfig(min_level=1)),
    ("≥2球", HoseConfig(min_level=2)),
    ("≥3球", HoseConfig(min_level=3)),
    ("≥4球", HoseConfig(min_level=4)),
    ("≥1球/护锚定", HoseConfig(min_level=1, protect_curated=True)),
    # σ 只降置信度、不改位移，因此这一行与"≥1球"**应当完全一致**。
    # 保留它是为了让回归者一眼看到"σ 警告没有偷偷动数值"。
    ("≥1球/关σ警告", HoseConfig(min_level=1, sigma_warn_ratio=None)),
]


def _experimental(mol: Chem.Mol, nucleus: str) -> dict[int, float]:
    """``{原子索引: 位移}``；同一原子出现多次时取平均。"""
    acc: dict[int, list[float]] = defaultdict(list)
    n_atoms = mol.GetNumAtoms()
    prefix = f"Spectrum {nucleus} "
    for prop in mol.GetPropNames():
        if not prop.startswith(prefix):
            continue
        for shift, idx in _iter_shift_entries(mol.GetProp(prop)):
            if 0 <= idx < n_atoms and _atom_matches_nucleus(
                mol.GetAtomWithIdx(idx), nucleus
            ):
                acc[idx].append(shift)
    return {k: sum(v) / len(v) for k, v in acc.items()}


def _prep_from_mol(mol: Chem.Mol) -> PreparedMol:
    """把数据库里的 molfile 包成 PreparedMol。

    molfile 本身不含显式氢，``AddHs`` 保留重原子顺序、把氢追加在末尾，
    因此重原子索引与数据库里的归属序号一一对应。
    """
    heavy = Chem.Mol(mol)
    return PreparedMol(
        smiles_in=Chem.MolToSmiles(mol),
        mol=Chem.AddHs(Chem.Mol(mol)),
        heavy=heavy,
        smiles_out=Chem.MolToSmiles(mol),
    )


def _groups(prep: PreparedMol, nucleus: str) -> list:
    """按核取等价分组：¹H 按带氢重原子，¹³C 按碳，¹⁹F/³¹P 每个原子一组。"""
    if nucleus == "1H":
        return hydrogen_groups(prep, 0)
    if nucleus == "13C":
        return carbon_groups(prep, 0)
    return hetero_groups(prep, "F" if nucleus == "19F" else "P", 0)


def _mean(values: list[float]) -> float:
    return statistics.mean(values) if values else float("nan")


def _pct(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * q))]


def main() -> int:
    ap = argparse.ArgumentParser(description="Tier A vs HOSE 逐原子对照")
    ap.add_argument("--sd", required=True, type=Path)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--salt", default="nmrflow-hose-v1")
    ap.add_argument("--spheres", type=int, default=4)
    ap.add_argument("--nuclei", default="1H,13C")
    ap.add_argument("--limit", type=int, default=None, help="只跑前 N 条记录")
    ap.add_argument(
        "--train-index",
        type=Path,
        default=_ROOT / "build" / "hose_index_train.sqlite",
        help="训练集索引的缓存路径；已存在则直接复用",
    )
    ap.add_argument("--rebuild", action="store_true", help="强制重建训练索引")
    args = ap.parse_args()

    nuclei = tuple(n.strip() for n in args.nuclei.split(",") if n.strip())
    bad = [n for n in nuclei if n not in SUPPORTED]
    if bad:
        print(f"不支持的核: {bad}", file=sys.stderr)
        return 2

    t0 = time.perf_counter()
    if args.rebuild or not args.train_index.exists():
        print(f"① 用 {100 * (1 - args.holdout):.0f}% 化合物建训练索引 ...", flush=True)
        build(
            args.sd,
            args.train_index,
            args.spheres,
            nuclei,
            args.limit,
            measure=False,
            keep=lambda mol: (
                (key := _inchikey(mol)) is not None
                and not _is_holdout(key, args.holdout, args.salt)
            ),
            verbose=False,
        )
    else:
        print(f"① 复用训练索引 {args.train_index}", flush=True)
    print(f"   {time.perf_counter() - t0:.0f}s", flush=True)

    print("② 逐原子对照 ...", flush=True)
    err_a: dict[str, list[float]] = defaultdict(list)
    # Tier A 退化到兜底常数的那部分必须单独统计：早期版本 ³¹P 兜底恒为 0.0 ppm、
    # ¹⁹F 恒为 −120.0，实测却有 −457 / +477 这种值，混在一起会把"我们没有这条规则"
    # 算成"我们的化学预测错了很多"，从而误判 HOSE 的价值。
    err_a_fallback: dict[str, list[float]] = defaultdict(list)
    err_a_real: dict[str, list[float]] = defaultdict(list)
    err_b_on_real: dict[str, list[float]] = defaultdict(list)
    n_a_fallback: dict[str, int] = defaultdict(int)
    #: 逐原子留档 (target, a_err, b_err|None, a_fallback, is_carbon_fluorine)
    #: 用来回答"把物理上不可能的实验值剔掉之后，误差还剩多少"——
    #: 这跟模型好坏是两件事，混在一起会既高估问题、也看不清该改哪一层。
    atom_rows: dict[str, list[tuple[float, float, float | None, bool, bool]]] = defaultdict(list)
    #: Tier A 兜底子集上，HOSE 按球层分别有多准——决定"该不该用它替掉兜底常数"
    fb_by_level: dict[str, dict[int, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    err_b: dict[str, list[float]] = defaultdict(list)
    per_cfg: dict[str, dict[str, list[float]]] = {n: defaultdict(list) for n, _ in CONFIGS}
    wins: dict[str, int] = defaultdict(int)
    # (核, 芳香/非芳香) -> 球层 -> (A 误差列表, B 误差列表)
    by_level: dict[tuple[str, str], dict[int, tuple[list[float], list[float]]]] = defaultdict(dict)
    counts: dict[str, int] = defaultdict(int)

    predictors = {
        "1H": ProtonShiftPredictor(),
        "13C": CarbonShiftPredictor(),
        "19F": HeteroShiftPredictor("19F"),
        "31P": HeteroShiftPredictor("31P"),
    }

    with HoseIndex(args.train_index) as index:
        hose = {
            nuc: HoseShiftPredictor(nuc, index, HoseConfig(min_level=1))
            for nuc in nuclei
        }
        supp = Chem.SDMolSupplier(str(args.sd), removeHs=False, sanitize=False)
        seen = 0
        for mol in supp:
            seen += 1
            if args.limit and seen > args.limit:
                break
            if mol is None or any(mol.HasProp(p) for p in CALC_PROPS):
                continue
            key = _inchikey(mol)
            if key is None or not _is_holdout(key, args.holdout, args.salt):
                continue
            try:
                Chem.SanitizeMol(mol)
                prep = _prep_from_mol(mol)
            except Exception:
                continue

            for nucleus in nuclei:
                exp = _experimental(mol, nucleus)
                if not exp:
                    continue
                groups = _groups(prep, nucleus)
                base = predictors[nucleus].predict(prep, groups, SOLVENT)
                cand = hose[nucleus].predict(prep, groups, SOLVENT)
                merged = {name: copy.deepcopy(base) for name, _ in CONFIGS}

                for group, sig, cd in zip(groups, base, cand):
                    keys = [i for i in group.heavy_indices if i in exp]
                    if not keys:
                        continue
                    target = sum(exp[i] for i in keys) / len(keys)
                    a_err = abs(sig.shift - target)
                    err_a[nucleus].append(a_err)
                    counts[nucleus] += 1
                    a_fb = bool(sig.provenance.fallback)
                    if a_fb:
                        err_a_fallback[nucleus].append(a_err)
                        n_a_fallback[nucleus] += 1
                    else:
                        err_a_real[nucleus].append(a_err)

                    is_cf = nucleus == "19F" and all(
                        nb.GetAtomicNum() == 6
                        for i in group.heavy_indices
                        for nb in prep.heavy.GetAtomWithIdx(i).GetNeighbors()
                    )
                    atom_rows[nucleus].append(
                        (target, a_err, None, a_fb, is_cf)
                    )
                    if cd.provenance.fallback:
                        continue
                    b_err = abs(cd.shift - target)
                    atom_rows[nucleus][-1] = (target, a_err, b_err, a_fb, is_cf)
                    err_b[nucleus].append(b_err)
                    if not a_fb:
                        err_b_on_real[nucleus].append(b_err)
                    elif not cd.provenance.fallback:
                        fb_by_level[nucleus][
                            int(cd.provenance.key.split("球")[0])
                        ].append(b_err)
                    if b_err < a_err:
                        wins[nucleus] += 1
                    level = int(cd.provenance.key.split("球")[0])
                    heavy_idx = group.heavy_indices[0] if group.heavy_indices else -1
                    kind = (
                        "芳香"
                        if heavy_idx >= 0
                        and prep.heavy.GetAtomWithIdx(heavy_idx).GetIsAromatic()
                        else "非芳香"
                    )
                    pair = by_level[(nucleus, kind)].setdefault(level, ([], []))
                    pair[0].append(a_err)
                    pair[1].append(b_err)

                for name, cfg in CONFIGS:
                    overlay_tier_b(merged[name], cand, groups, nucleus, cfg)
                for name, _ in CONFIGS:
                    for group, sig in zip(groups, merged[name]):
                        keys = [i for i in group.heavy_indices if i in exp]
                        if not keys:
                            continue
                        target = sum(exp[i] for i in keys) / len(keys)
                        per_cfg[name][nucleus].append(abs(sig.shift - target))

            if seen % 5000 == 0:
                print(f"   {seen:,} 条记录 ... {time.perf_counter() - t0:.0f}s", flush=True)

    print("\n" + "=" * 84)
    print(f"留出 {args.holdout:.0%} 化合物 | 只统计 A 与 B 都能给出值的原子")
    print("=" * 84)
    print(f"{'核':<5}{'A 可预测':>10}{'B 覆盖':>9}{'A MAE':>9}{'B MAE':>9}{'B 更优':>9}")
    print("-" * 60)
    for nucleus in nuclei:
        total = counts[nucleus]
        if not total:
            continue
        print(
            f"{nucleus:<5}{total:>10,}{len(err_b[nucleus]) / total:>8.0%}"
            f"{_mean(err_a[nucleus]):>9.2f}{_mean(err_b[nucleus]):>9.2f}"
            f"{wins[nucleus] / max(1, len(err_b[nucleus])):>8.0%}"
        )

    print("\nTier A 的覆盖质量（同环境有无 SMARTS 规则，还是退到兜底常数）")
    head2 = (
        f"{'核':<6}{'原子':>8}{'A 兜底':>10}{'A MAE':>9}{'兜底子集':>11}"
        f"{'有规则子集':>12}{'B 在规则子集':>14}{'B P90':>10}{'B 最大':>10}"
    )
    print(head2)
    print("-" * len(head2))
    for nucleus in nuclei:
        if not counts[nucleus]:
            continue
        n_fb = n_a_fallback[nucleus]
        print(
            f"{nucleus:<6}{counts[nucleus]:>8,}{n_fb:>9,}({100 * n_fb / counts[nucleus]:.0f}%)"
            f"{_mean(err_a[nucleus]):>9.2f}{_mean(err_a_fallback[nucleus]):>11.2f}"
            f"{_mean(err_a_real[nucleus]):>12.2f}{_mean(err_b_on_real[nucleus]):>14.2f}"
            f"{_pct(err_b[nucleus], 0.9):>10.1f}{max(err_b[nucleus]):>10.1f}"
        )

    print("\nTier A 走到兜底常数时，HOSE 逐球层能否替掉它")
    for nucleus in nuclei:
        if n_a_fallback[nucleus] == 0:
            continue
        print(f"  {nucleus}（兜底子集 {n_a_fallback[nucleus]:,} 个原子，"
              f"Tier A 兜底值 MAE {_mean(err_a_fallback[nucleus]):.2f}）")
        print(f"    {'球层':<6}{'原子':>8}{'HOSE MAE':>11}{'比兜底强':>10}")
        for level in sorted(fb_by_level[nucleus]):
            vals = fb_by_level[nucleus][level]
            print(f"    {level:<6}{len(vals):>8,}{_mean(vals):>11.2f}"
                  f"{'是' if _mean(vals) < _mean(err_a_fallback[nucleus]) else '否':>10}")

    print("\n实验值的物理合理性过滤：误差里有多少来自『数据本身不像真的』")
    print("  论文口径 = 剔除 δ>300 或 δ<−400（Li et al. 2024 的规则）")
    print("  有机C–F  = F 的邻居必须是碳，且 −250 ≤ δ ≤ 60")
    head3 = f"{'核':<5}{'口径':<12}{'原子':>8}{'A MAE':>9}{'B MAE':>9}{'B 最大':>9}"
    print(head3)
    print("-" * len(head3))
    for nucleus, rows in atom_rows.items():
        if not rows:
            continue
        for label, keep in (
            ("全部", lambda r: True),
            ("论文口径", lambda r: -400.0 <= r[0] <= 300.0),
            ("有机C–F", lambda r: r[4] and -250.0 <= r[0] <= 60.0),
        ):
            sel = [r for r in rows if keep(r)]
            if not sel:
                continue
            a = _mean([r[1] for r in sel])
            bb = [r[2] for r in sel if r[2] is not None]
            b = _mean(bb) if bb else float("nan")
            mx = max(bb) if bb else float("nan")
            print(f"{nucleus:<5}{label:<12}{len(sel):>8,}{a:>9.2f}{b:>9.2f}{mx:>9.1f}")

    print("\n按球层 × 原子类型：A 与 B 在同一批原子上的 MAE（ppm）")
    for nucleus in nuclei:
        for kind in ("芳香", "非芳香"):
            rows = by_level.get((nucleus, kind))
            if not rows:
                continue
            print(f"  {nucleus} {kind}")
            print(f"    {'球层':<5}{'条目':>8}{'A MAE':>9}{'B MAE':>9}{'B 更优':>8}")
            for level in sorted(rows):
                a_list, b_list = rows[level]
                w = sum(1 for a, b in zip(a_list, b_list) if b < a)
                print(
                    f"    {level:<5}{len(a_list):>8,}{_mean(a_list):>9.2f}"
                    f"{_mean(b_list):>9.2f}{w / max(1, len(a_list)):>7.0%}"
                )

    # 只看均值会掩盖"个别离谱覆盖"，因此同时报 90 分位与最大值
    print("\n合并策略（覆盖全部可预测原子；未覆盖的保留 Tier A）")
    head = f"{'策略':<18}" + "".join(
        f"{n + ' MAE':>9}{n + ' P90':>9}{n + ' 最大':>9}" for n in nuclei
    )
    print(head)
    print("-" * len(head))
    print(
        f"{'仅 Tier A':<18}"
        + "".join(
            f"{_mean(err_a[n]):>9.2f}{_pct(err_a[n], 0.9):>9.2f}{max(err_a[n]):>9.1f}"
            if err_a[n] else f"{'—':>9}" * 3
            for n in nuclei
        )
    )
    for name, _ in CONFIGS:
        cells = []
        for n in nuclei:
            vals = per_cfg[name][n]
            cells.append(
                f"{_mean(vals):>9.2f}{_pct(vals, 0.9):>9.2f}{max(vals):>9.1f}"
                if vals else f"{'—':>9}" * 3
            )
        print(f"{name:<18}" + "".join(cells))

    print(f"\n总耗时 {time.perf_counter() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
