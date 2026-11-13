"""诊断：curated 表里"钉死"的位移，有多少被浅层 HOSE（<4 球）盖掉了。

为什么值得单独查
----------------
Tier A 的 ¹⁹F / ³¹P 有一张按 SMARTS 写死位移的环境表（``nmrflow/data/shift_19f.toml``
等）。这些条目是"这个环境实测就在这个位置"的结论，比 2 层球的均值可信得多。

但 ¹⁹F 的 ``HoseConfig.protect_curated`` 默认是 **False**（只有 ³¹P 是 True），
于是 2 层球命中就可以把 curated 值覆盖掉。已知的受害者：

    CF₃SO₂F 的 CF₃：表里条目 −78，被 2 球 HOSE 拉到 −58.2（差 20 ppm）

2 层球只编码到中心碳的取代元素，分不出"CF₃ 接 SO₂F 还是接苯环"；
而 ¹⁹F 的 3/4 球命中率本来就不高，所以这类覆盖经常是拿一个粗均值换掉一个精确值。

本工具回答的问题是：**全库有多少个这样的原子，覆盖之后离实测更远还是更近。**
不是"覆盖了就一定坏"——有些 curated 条目也许本身就偏，用实测说话。

用法
----
    python tools/diag_curated_override.py --sd <SD> [--nuclei 19F,31P] [--top 15]

每个核输出：

    钉死值总数 / 被浅层覆盖数 / 覆盖前 MAE / 覆盖后 MAE / 变差数（>1、>5 ppm）
    以及最差的若干条：SMILES、原子、球层、Tier A、覆盖后、实测、Δ

``--nuclei`` 也可以给 ``1H,13C``，但那两个核的 curated 覆盖面窄、跑一轮要十几分钟。
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

from build_hose_index import (  # noqa: E402
    CALC_PROPS,
    SUPPORTED,
    _atom_matches_nucleus,
    _iter_shift_entries,
)
from nmrflow.mol.equivalence import (  # noqa: E402
    carbon_groups,
    hetero_groups,
    hydrogen_groups,
)
from nmrflow.mol.prepare import PreparedMol  # noqa: E402
from nmrflow.predict.hose_index import HoseIndex  # noqa: E402
from nmrflow.predict.shift_13c import CarbonShiftPredictor  # noqa: E402
from nmrflow.predict.shift_1h import ProtonShiftPredictor  # noqa: E402
from nmrflow.predict.shift_hetero import HeteroShiftPredictor  # noqa: E402
from nmrflow.predict.shift_hose import (  # noqa: E402
    HoseConfig,
    HoseShiftPredictor,
    _curated,
    _level_of,
    overlay_tier_b,
)

SOLVENT = "CDCl3"


def _prep_from_mol(mol: Chem.Mol) -> PreparedMol:
    heavy = Chem.Mol(mol)
    return PreparedMol(
        smiles_in=Chem.MolToSmiles(mol),
        mol=Chem.AddHs(Chem.Mol(mol)),
        heavy=heavy,
        smiles_out=Chem.MolToSmiles(mol),
    )


def _groups(prep: PreparedMol, nucleus: str) -> list:
    if nucleus == "1H":
        return hydrogen_groups(prep, 0)
    if nucleus == "13C":
        return carbon_groups(prep, 0)
    return hetero_groups(prep, "F" if nucleus == "19F" else "P", 0)


def _experimental(mol: Chem.Mol, nucleus: str) -> dict[int, float]:
    """``{原子索引: 位移}``；只取第 0 张谱（见 diag_multispectrum.py）。"""
    prefix = f"Spectrum {nucleus} "
    props = [p for p in mol.GetPropNames() if p.startswith(prefix)]
    if not props:
        return {}

    def index(prop: str) -> int:
        tail = prop.rsplit(" ", 1)[-1]
        return int(tail) if tail.isdigit() else 0

    acc: dict[int, list[float]] = defaultdict(list)
    n_atoms = mol.GetNumAtoms()
    for shift, idx in _iter_shift_entries(mol.GetProp(min(props, key=index))):
        if 0 <= idx < n_atoms and _atom_matches_nucleus(
            mol.GetAtomWithIdx(idx), nucleus
        ):
            acc[idx].append(shift)
    return {k: sum(v) / len(v) for k, v in acc.items()}


def _mean(values: list[float]) -> float:
    return statistics.mean(values) if values else float("nan")


def run(sd_path: Path, nuclei: list[str], index_path: Path, top: int, limit: int | None):
    predictors = {
        "1H": ProtonShiftPredictor,
        "13C": CarbonShiftPredictor,
        "19F": lambda: HeteroShiftPredictor("19F"),
        "31P": lambda: HeteroShiftPredictor("31P"),
    }
    tiers = {n: predictors[n]() for n in nuclei}
    cfg = HoseConfig(min_level=1, protect_curated=False)  # 刻意关掉保护：要看"不保护会怎样"

    stats: dict[str, dict] = {
        n: {"curated": 0, "overridden": 0, "err_before": [], "err_after": [], "cases": []}
        for n in nuclei
    }

    with HoseIndex(index_path) as index:
        hose = {n: HoseShiftPredictor(n, index, cfg) for n in nuclei}
        seen = 0
        for mol in Chem.ForwardSDMolSupplier(str(sd_path)):
            seen += 1
            if limit and seen > limit:
                break
            if mol is None or any(mol.HasProp(p) for p in CALC_PROPS):
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
                base = tiers[nucleus].predict(prep, groups, SOLVENT)
                cand = hose[nucleus].predict(prep, groups, SOLVENT)
                merged = copy.deepcopy(base)
                overlay_tier_b(merged, cand, groups, nucleus, cfg)

                s = stats[nucleus]
                for group, a_sig, m_sig in zip(groups, base, merged):
                    if not _curated(a_sig):
                        continue
                    keys = [i for i in group.heavy_indices if i in exp]
                    if not keys:
                        continue
                    target = _mean([exp[i] for i in keys])
                    s["curated"] += 1
                    if abs(a_sig.shift - m_sig.shift) < 1e-9:
                        continue          # 没被覆盖
                    level = _level_of(m_sig)
                    if level >= 4:
                        continue          # 只关心浅层（<4 球）
                    s["overridden"] += 1
                    err_a = abs(a_sig.shift - target)
                    err_m = abs(m_sig.shift - target)
                    s["err_before"].append(err_a)
                    s["err_after"].append(err_m)
                    s["cases"].append(
                        {
                            "smiles": Chem.MolToSmiles(mol),
                            "atom": keys[0],
                            "level": level,
                            "a": round(a_sig.shift, 2),
                            "m": round(m_sig.shift, 2),
                            "exp": round(target, 2),
                            "delta": round(err_m - err_a, 2),
                            "key": a_sig.provenance.key,
                        }
                    )

    report: dict[str, dict] = {}
    for nucleus, s in stats.items():
        before, after = s["err_before"], s["err_after"]
        worse = [c for c in s["cases"] if c["delta"] > 0]
        better = [c for c in s["cases"] if c["delta"] < 0]
        report[nucleus] = {
            "curated_atoms": s["curated"],
            "overridden_shallow": s["overridden"],
            "mae_before": round(_mean(before), 2) if before else None,
            "mae_after": round(_mean(after), 2) if after else None,
            "worse_than_1ppm": sum(1 for c in s["cases"] if c["delta"] > 1),
            "worse_than_5ppm": sum(1 for c in s["cases"] if c["delta"] > 5),
            "improved": len(better),
            # 只有真变差的进 worst、真变好的进 best，否则 top=N 会把两边混在一起
            "worst": sorted(worse, key=lambda c: -c["delta"])[:top],
            "best": sorted(better, key=lambda c: c["delta"])[:top],
        }
    return report


def main() -> int:
    ap = argparse.ArgumentParser(description="查 curated 值被浅层 HOSE 覆盖的情况")
    ap.add_argument("--sd", required=True, type=Path)
    ap.add_argument("--nuclei", default="19F,31P")
    ap.add_argument(
        "--train-index",
        type=Path,
        default=_ROOT / "build" / "hose_index_hetero.sqlite",
        help="HOSE 索引；查 1H/13C 时传 hose_index_train.sqlite",
    )
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()

    nuclei = [n.strip() for n in args.nuclei.split(",") if n.strip()]
    bad = [n for n in nuclei if n not in SUPPORTED]
    if bad:
        print(f"不支持的核: {bad}", file=sys.stderr)
        return 2
    if not args.sd.exists():
        print(f"找不到 {args.sd}", file=sys.stderr)
        return 2

    t0 = time.perf_counter()
    report = run(args.sd, nuclei, args.train_index, args.top, args.limit)

    if args.json:
        import json

        print(json.dumps(report, ensure_ascii=False, indent=1))
        return 0

    for nucleus, r in report.items():
        print(f"=== {nucleus}")
        print(f"  钉死值原子（Tier A curated）  {r['curated_atoms']:,}")
        print(f"  被 <4 球 HOSE 覆盖            {r['overridden_shallow']:,}")
        if not r["overridden_shallow"]:
            print("  （没有浅层覆盖，这个核是安全的）\n")
            continue
        print(f"  MAE 覆盖前 / 覆盖后           {r['mae_before']} / {r['mae_after']}")
        print(
            f"  变差 >1 ppm / >5 ppm          {r['worse_than_1ppm']} / {r['worse_than_5ppm']}"
            f"   （变好 {r['improved']}）"
        )
        print("  --- 覆盖后变差的（Δ>0 表示离实测更远）")
        if r["worst"]:
            for c in r["worst"]:
                print(
                    f"    {c['smiles'][:48]:48s} lv{c['level']} "
                    f"A={c['a']:>7.2f} → {c['m']:>7.2f}  实测={c['exp']:>7.2f}  "
                    f"Δ={c['delta']:+.2f}"
                )
        else:
            print("    （没有变差的案例）")
        if r["best"]:
            print("  --- 覆盖后变好的（Δ<0 表示更接近实测）")
            for c in r["best"][: max(1, args.top // 3)]:
                print(
                    f"    {c['smiles'][:48]:48s} lv{c['level']} "
                    f"A={c['a']:>7.2f} → {c['m']:>7.2f}  实测={c['exp']:>7.2f}  "
                    f"Δ={c['delta']:+.2f}"
                )
        print()
    print(f"总耗时 {time.perf_counter() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
