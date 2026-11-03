"""用留出法评估 HOSE 索引的真实预测精度。

为什么按 InChIKey 划分
----------------------
同一化合物在 NMRShiftDB2 里常有多条记录（不同溶剂、不同场强）。若按记录随机划分，
同一分子的另一条谱会留在训练集里，测出来的是"背答案"的精度。
因此以 InChIKey 为单位划分，保证训练集里没有测试集的任何化合物。

流程
----
1. 用 80% 化合物建一份临时索引（``--holdout`` 可调）；
2. 遍历全部记录，只取剩下 20% 的化合物，逐条比对
   「预测值 vs 数据库实测值」，按核统计 MAE / 中位误差 / 90 分位；
3. 同时统计球层分布与覆盖率——MAE 必须和"能预测多少"一起看，
   只报 MAE 不报覆盖率会掩盖大量降级。

用法
----
    python tools/eval_hose_index.py --sd /path/to/nmrshiftdb2withsignals.sd
"""

from __future__ import annotations

import argparse
import hashlib
import statistics
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "tools"))

from build_hose_index import CALC_PROPS, SUPPORTED, _atom_matches_nucleus, _iter_shift_entries, build  # noqa: E402
from nmrflow.predict.hose import hose_codes  # noqa: E402
from nmrflow.predict.hose_index import HoseIndex  # noqa: E402


def _inchikey(mol: Chem.Mol) -> str | None:
    for prop in ("INChI key", "INChIKey", "INCHIKEY"):
        if mol.HasProp(prop):
            value = mol.GetProp(prop).strip()
            if value:
                return value
    return None


def _is_holdout(key: str, holdout: float, salt: str) -> bool:
    digest = hashlib.blake2b(f"{salt}:{key}".encode(), digest_size=8).digest()
    return (int.from_bytes(digest, "big") % 1_000_000) < holdout * 1_000_000


def _pct(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * q))]


def _n_bucket(n: int) -> str:
    if n <= 1:
        return "1"
    if n == 2:
        return "2"
    if n <= 5:
        return "3-5"
    if n <= 15:
        return "6-15"
    return "16+"


N_BUCKETS = ("1", "2", "3-5", "6-15", "16+")


def main() -> int:
    ap = argparse.ArgumentParser(description="留出法评估 HOSE 索引精度")
    ap.add_argument("--sd", required=True, type=Path)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--salt", default="nmrflow-hose-v1")
    ap.add_argument("--spheres", type=int, default=4)
    ap.add_argument("--min-neighbors", type=int, default=1)
    ap.add_argument("--nuclei", default="1H,13C")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    nuclei = tuple(n.strip() for n in args.nuclei.split(",") if n.strip())
    bad = [n for n in nuclei if n not in SUPPORTED]
    if bad:
        print(f"不支持的核: {bad}", file=sys.stderr)
        return 2

    t0 = time.perf_counter()

    print(f"① 用 {100 * (1 - args.holdout):.0f}% 化合物建临时索引 ...", flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        idx_path = Path(tmp) / "train.sqlite"
        r = build(
            args.sd,
            idx_path,
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
        print(f"   训练集分子 {r['stats']['mol'] - r['stats']['skip_filter']:,}，"
              f"索引行 {sum(x[1] for x in r['rows']):,}", flush=True)

        print("② 在留出的化合物上比对预测值与实测值 ...", flush=True)
        errors: dict[str, list[float]] = defaultdict(list)
        level_of: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        by_level: dict[str, dict[int, list[float]]] = defaultdict(lambda: defaultdict(list))
        by_n: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
        entries: dict[str, int] = defaultdict(int)
        missed: dict[str, int] = defaultdict(int)

        with HoseIndex(idx_path) as index:
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
                except Exception:
                    continue

                codes = hose_codes(mol, args.spheres)
                n_atoms = mol.GetNumAtoms()
                for nucleus in nuclei:
                    prefix = f"Spectrum {nucleus} "
                    for prop in mol.GetPropNames():
                        if not prop.startswith(prefix):
                            continue
                        for shift, idx in _iter_shift_entries(mol.GetProp(prop)):
                            if not (0 <= idx < n_atoms):
                                continue
                            if not _atom_matches_nucleus(mol.GetAtomWithIdx(idx), nucleus):
                                continue
                            entries[nucleus] += 1
                            code = codes.get(idx)
                            hit = (
                                index.lookup(
                                    code,
                                    nucleus,
                                    max_level=args.spheres,
                                    min_neighbors=args.min_neighbors,
                                )
                                if code
                                else None
                            )
                            if hit is None:
                                missed[nucleus] += 1
                                continue
                            err = abs(hit.mean - shift)
                            errors[nucleus].append(err)
                            level_of[nucleus][hit.level] += 1
                            by_level[nucleus][hit.level].append(err)
                            by_n[nucleus][_n_bucket(hit.n)].append(err)
                if seen % 5000 == 0:
                    print(f"   {seen:,} 条记录 ... {time.perf_counter() - t0:.0f}s", flush=True)

    print("\n" + "=" * 66)
    print(f"留出 {args.holdout:.0%} 化合物 | {args.spheres} 层球 | min_neighbors={args.min_neighbors}")
    print("=" * 66)
    header = f"{'核':<5}{'实测条目':>9}{'已预测':>9}{'覆盖率':>8}{'MAE':>8}{'中位':>8}{'90分位':>9}"
    print(header)
    print("-" * len(header))
    for nucleus in nuclei:
        errs = errors[nucleus]
        total = entries[nucleus]
        if not total:
            print(f"{nucleus:<5}{0:>9}")
            continue
        cover = len(errs) / total
        print(
            f"{nucleus:<5}{total:>9,}{len(errs):>9,}{cover:>7.1%}"
            f"{statistics.mean(errs):>8.2f}{statistics.median(errs):>8.2f}{_pct(errs, 0.9):>9.2f}"
        )

    print("\n命中球层分布（占已预测条目）：")
    for nucleus in nuclei:
        dist = level_of[nucleus]
        total = sum(dist.values())
        if not total:
            continue
        parts = [f"球{k}: {v / total:.0%}" for k, v in sorted(dist.items())]
        print(f"  {nucleus:<4} " + "  ".join(parts))

    print("\n按球层拆分（绝对误差，ppm）：")
    for nucleus in nuclei:
        rows = by_level[nucleus]
        if not rows:
            continue
        tot = sum(len(v) for v in rows.values())
        print(f"  {nucleus}")
        print(f"    {'球层':<6}{'条目':>9}{'占比':>7}{'MAE':>8}{'中位':>8}")
        for level in sorted(rows):
            errs = rows[level]
            print(
                f"    {level:<6}{len(errs):>9,}{len(errs) / tot:>7.0%}"
                f"{statistics.mean(errs):>8.2f}{statistics.median(errs):>8.2f}"
            )

    print("\n按同码邻居数拆分（绝对误差，ppm）：")
    for nucleus in nuclei:
        rows = by_n[nucleus]
        if not rows:
            continue
        tot = sum(len(v) for v in rows.values())
        print(f"  {nucleus}")
        print(f"    {'邻居数':<6}{'条目':>9}{'占比':>7}{'MAE':>8}{'中位':>8}")
        for bucket in N_BUCKETS:
            errs = rows.get(bucket)
            if not errs:
                continue
            print(
                f"    {bucket:<6}{len(errs):>9,}{len(errs) / tot:>7.0%}"
                f"{statistics.mean(errs):>8.2f}{statistics.median(errs):>8.2f}"
            )

    print(f"\n总耗时 {time.perf_counter() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
