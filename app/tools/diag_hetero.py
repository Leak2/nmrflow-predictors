"""查 ¹⁹F / ³¹P 的离群值：是实验数据本身脏，还是 Tier A 算错了。

背景
----
`eval_tier_ab.py --nuclei 31P,19F` 给出的数字很可疑：

    31P  A MAE 40.40  最大误差 457 ppm
    19F  A MAE 23.62  最大误差 597 ppm

³¹P 的实际范围大致 −250 ~ +250 ppm（相对 85% H₃PO₄），¹⁹F 大致 −250 ~ +100
（相对 CFCl₃）。457 / 597 ppm 的"误差"不可能来自正常的化学位移差异，
要么库里有脏记录，要么我们的取值方式有问题。

所以这个脚本直接逐条打印最差的样本，把 (实验值, Tier A, HOSE) 摆在一起看。
不看清楚就下结论，等于用脏数据去决定要不要给一个功能默认开启。
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "tools"))

from build_hose_index import CALC_PROPS  # noqa: E402
from eval_hose_index import _inchikey, _is_holdout  # noqa: E402
from eval_tier_ab import SOLVENT, _experimental, _groups, _prep_from_mol  # noqa: E402
from nmrflow.predict.hose_index import HoseIndex  # noqa: E402
from nmrflow.predict.shift_hetero import HeteroShiftPredictor  # noqa: E402
from nmrflow.predict.shift_hose import HoseConfig, HoseShiftPredictor  # noqa: E402


def _pct(values: list[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * q))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sd", required=True, type=Path)
    ap.add_argument("--nuclei", default="31P,19F")
    ap.add_argument("--index", type=Path,
                    default=_ROOT / "build" / "hose_index_hetero.sqlite")
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--salt", default="nmrflow-hose-v1")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()

    nuclei = tuple(n.strip() for n in args.nuclei.split(",") if n.strip())
    tier_a = {n: HeteroShiftPredictor(n) for n in nuclei}
    worst: dict[str, list[tuple[float, str, float, float, float]]] = {n: [] for n in nuclei}
    raw: dict[str, list[float]] = {n: [] for n in nuclei}
    fallback_count: dict[str, Counter] = {n: Counter() for n in nuclei}

    with HoseIndex(args.index) as index:
        hose = {n: HoseShiftPredictor(n, index, HoseConfig(min_level=1)) for n in nuclei}
        supp = Chem.SDMolSupplier(str(args.sd), removeHs=False, sanitize=False)
        seen = 0
        for mol in supp:
            seen += 1
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
                raw[nucleus].extend(exp.values())
                groups = _groups(prep, nucleus)
                a_sigs = tier_a[nucleus].predict(prep, groups, SOLVENT)
                b_sigs = hose[nucleus].predict(prep, groups, SOLVENT)
                for group, a, b in zip(groups, a_sigs, b_sigs):
                    keys = [i for i in group.heavy_indices if i in exp]
                    if not keys:
                        continue
                    target = sum(exp[i] for i in keys) / len(keys)
                    err = abs(a.shift - target)
                    fallback_count[nucleus][a.provenance.key.split("(")[0]] += 1
                    item = (err, _smiles(mol), a.shift, b.shift, target)
                    bucket = worst[nucleus]
                    bucket.append(item)
                    bucket.sort(key=lambda t: -t[0])
                    del bucket[args.top:]

            if seen % 20000 == 0:
                print(f"  {seen:,} 条 ...", flush=True)

    for nucleus in nuclei:
        vals = raw[nucleus]
        print("=" * 96)
        print(f"{nucleus}：库内实验值分布（留出集，{len(vals):,} 个原子）")
        print("-" * 96)
        print(
            f"  最小 {min(vals):8.1f}   1% {_pct(vals, 0.01):8.1f}   "
            f"中位 {_pct(vals, 0.5):8.1f}   99% {_pct(vals, 0.99):8.1f}   "
            f"最大 {max(vals):8.1f}"
        )
        outside = sum(1 for v in vals if abs(v) > 300)
        print(f"  |δ| > 300 ppm 的记录：{outside} 条 ({100 * outside / max(1, len(vals)):.1f}%)")

        print(f"\n  Tier A 的命中来源（前 6）：")
        for name, cnt in fallback_count[nucleus].most_common(6):
            print(f"    {name:<28}{cnt:>6,}")

        print(f"\n  最差的 {len(worst[nucleus])} 个原子：")
        print(f"    {'误差':>8}{'TierA':>9}{'HOSE':>9}{'实测':>9}   SMILES")
        for err, smi, a, b, t in worst[nucleus]:
            print(f"    {err:>8.1f}{a:>9.1f}{b:>9.1f}{t:>9.1f}   {smi[:52]}")
        print()

    return 0


def _smiles(mol: Chem.Mol) -> str:
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return "<unparsable>"


if __name__ == "__main__":
    raise SystemExit(main())
