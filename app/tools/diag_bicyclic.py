"""桥环/稠环诊断：按「环体系签名」统计 Tier A 与 Tier B 的逐原子误差。

动机（评审 P1-1）：`shift_13c.py` 对小双环体系只有三行硬编码
（`_SMALL_BIKE_BRIDGEHEAD` / `_SMALL_BIKE_BRIDGE_CH2`），标定只用了十氢萘与
降冰片烷两个分子。相邻骨架（双环[3.3.0]辛烷、双环[2.2.2]辛烷、甾体……）
要么落在「双环 → 不采信 HOSE」的判定里却没有对应常数，要么干脆走 HOSE。

本脚本从 NMRShiftDB2 原始 SD 里筛出脂肪族多环分子，按原子所属环体系分类，
分别算 Tier A（纯规则）与默认（规则 + HOSE）的误差，回答两件事：
  1. 哪些签名现在错得最厉害；
  2. 对同一签名，Tier A 和 HOSE 应该信谁。

★ 原子对齐方式与 `tools/eval_tier_ab.py` 一致：直接拿 SD 里的 molfile 包成
  ``PreparedMol``（``_prep_from_mol``），**不能**走 ``MolToSmiles`` 往返——
  规范化会重排原子，序号与数据库里的归属对不上，误差会变成 30~60 ppm 的噪声。

用法：
    python tools/diag_bicyclic.py --sd <withsignals.sd> [--set holdout] [--min-count 15]
"""

from __future__ import annotations

import argparse
import copy
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from rdkit import Chem, RDLogger  # noqa: E402

from build_hose_index import CALC_PROPS  # noqa: E402
from eval_tier_ab import (  # noqa: E402
    _experimental,
    _groups,
    _inchikey,
    _is_holdout,
    _prep_from_mol,
)
from nmrflow.predict.hose_index import HoseIndex  # noqa: E402
from nmrflow.predict.shift_13c import CarbonShiftPredictor  # noqa: E402
from nmrflow.predict.shift_hose import (  # noqa: E402
    HoseConfig,
    HoseShiftPredictor,
    overlay_tier_b,
)

RDLogger.DisableLog("rdApp.*")

SOLVENT = "CDCl3"


def _is_aliphatic_polycycle(mol: Chem.Mol, max_rings: int = 4,
                            max_ring_atoms: int = 22) -> bool:
    """全碳环（可带杂原子取代基）、无芳香、2~4 个环。"""
    rings = mol.GetRingInfo().AtomRings()
    if not 2 <= len(rings) <= max_rings:
        return False
    ring_atoms = {i for ring in rings for i in ring}
    if len(ring_atoms) > max_ring_atoms:
        return False
    if any(mol.GetAtomWithIdx(i).GetIsAromatic()
           or mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in ring_atoms):
        return False
    return all(a.GetSymbol() in ("C", "H", "O", "F", "Cl", "Br", "I")
               for a in mol.GetAtoms())


def _atom_signature(mol: Chem.Mol, idx: int) -> str:
    """原子视角的签名：所在环的尺寸组合 + 角色（桥头/桥碳/稠合/环上）。"""
    rings = [r for r in mol.GetRingInfo().AtomRings() if idx in r]
    sizes = tuple(sorted(len(r) for r in rings))
    atom = mol.GetAtomWithIdx(idx)
    ring_bonds = sum(1 for b in atom.GetBonds() if b.IsInRing())
    if ring_bonds >= 3:
        role = "桥头"
    elif len(rings) >= 2:
        role = "稠合"
    elif any(sum(1 for b in nb.GetBonds() if b.IsInRing()) >= 3
             for nb in atom.GetNeighbors()):
        role = "桥碳"
    else:
        role = "环上"
    return f"{sizes}{role}"


def _p90(values: list[float]) -> float:
    return statistics.quantiles(values, n=10)[8] if len(values) >= 10 else max(values)


def main() -> int:
    ap = argparse.ArgumentParser(description="桥环/稠环逐原子诊断")
    ap.add_argument("--sd", required=True)
    ap.add_argument("--set", choices=["all", "holdout", "train"], default="holdout")
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--salt", default="nmrflow-hose-v1")
    ap.add_argument("--index", default="build/hose_index_train.sqlite",
                    help="训练集索引；留出集评估必须用它，否则是自评自测")
    ap.add_argument("--min-count", type=int, default=15)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    a_pred = CarbonShiftPredictor()
    stats: dict[str, dict[str, list[float]]] = defaultdict(lambda: {"a": [], "b": []})
    mols = atoms = 0

    with HoseIndex(Path(args.index)) as index:
        hose = HoseShiftPredictor("13C", index, HoseConfig(min_level=1))
        supp = Chem.SDMolSupplier(str(args.sd), removeHs=False, sanitize=False)
        seen = 0
        for mol in supp:
            seen += 1
            if args.limit and seen > args.limit:
                break
            if mol is None or any(mol.HasProp(p) for p in CALC_PROPS):
                continue                      # 跳过计算谱，只用实验谱
            key = _inchikey(mol)
            if key is None:
                continue
            hold = _is_holdout(key, args.holdout, args.salt)
            if args.set == "holdout" and not hold:
                continue
            if args.set == "train" and hold:
                continue
            try:
                Chem.SanitizeMol(mol)
                if not _is_aliphatic_polycycle(mol):
                    continue
                prep = _prep_from_mol(mol)
            except Exception:
                continue

            exp = _experimental(mol, "13C")
            if not exp:
                continue
            groups = _groups(prep, "13C")
            base = a_pred.predict(prep, groups, SOLVENT)
            # 必须 deepcopy：overlay_tier_b 是原地改写 Signal，浅拷贝会让 A/B 变成同一对象
            merged = copy.deepcopy(base)
            cand = hose.predict(prep, groups, SOLVENT)
            if cand:
                overlay_tier_b(merged, cand, groups, "13C", HoseConfig(min_level=1))

            used = 0
            for group, sig_a, sig_b in zip(groups, base, merged):
                keys = [i for i in group.heavy_indices if i in exp]
                if not keys:
                    continue
                target = sum(exp[i] for i in keys) / len(keys)
                signature = _atom_signature(prep.heavy, keys[0])
                stats[signature]["a"].append(abs(sig_a.shift - target))
                stats[signature]["b"].append(abs(sig_b.shift - target))
                used += 1
            if used:
                mols += 1
                atoms += used

    print(f"筛选到 {mols} 个脂肪族多环分子 / {atoms} 个 ¹³C 信号（{args.set}）")
    header = (f"\n{'签名':<16}{'信号':>6}{'A MAE':>9}{'A P90':>8}"
              f"{'B MAE':>9}{'B P90':>8}{'B-A':>8}")
    print(header)
    print("-" * len(header))
    rows = [(s, v) for s, v in stats.items() if len(v["a"]) >= args.min_count]
    for sig, v in sorted(rows, key=lambda kv: -len(kv[1]["a"])):
        a, b = v["a"], v["b"]
        print(f"{sig:<16}{len(a):>6}{statistics.mean(a):>9.2f}{_p90(a):>8.2f}"
              f"{statistics.mean(b):>9.2f}{_p90(b):>8.2f}"
              f"{statistics.mean(b) - statistics.mean(a):>+8.2f}")

    all_a = [x for v in stats.values() for x in v["a"]]
    all_b = [x for v in stats.values() for x in v["b"]]
    if all_a:
        print("-" * len(header))
        print(f"{'合计':<16}{len(all_a):>6}{statistics.mean(all_a):>9.2f}"
              f"{_p90(all_a):>8.2f}{statistics.mean(all_b):>9.2f}{_p90(all_b):>8.2f}"
              f"{statistics.mean(all_b) - statistics.mean(all_a):>+8.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
