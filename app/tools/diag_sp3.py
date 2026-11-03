"""sp³ 碳（环烷 / 烯丙位 / 苄位）的 Tier A 诊断。

为什么单独写这个
----------------
Grant-Paul 是给开链烷烃拟合的，公式里没有任何环项。直接拿去算环烷，
环己烷会算出 32.2（实测 27.1）、环戊烷 34.7（实测 25.9）。

这个偏差会污染所有"环上的 sp³ 碳"，包括烯丙位：环己烯的烯丙位 CH₂ 预测 35.2、
实测 25.4，看着像"烯烃增量给错了"，其实 5.1 ppm 是环基线带进来的。

所以先单独把环烷基线量出来，再谈烯丙位/苄位增量。

数据约定
--------
和 diag_sp2 一致：给出"逐 sp³ 原子展开"的文献值列表，与预测值各自排序后按位相减。
不写 `{原子序号: 位移}`——SMILES 序号与文献编号的对应太容易写错（diag_sp2 里我
就在联苯上栽过一次）；也不做"最近邻配对"——挑的过程会把错藏起来。

★ **排序配对会摊平误差，这是它的固有局限**：它给出的是"在不知道原子归属的前提下
  能达到的最小总误差"，对**汇总**指标是对的，但会把单个原子的离谱误差分散掉。
  典型例子是新戊烷：季碳预测 34.1、实测 27.9（错 6.2），但预测的多重集是 4×29.6 + 34.1、
  文献是 27.9 + 4×31.7，排序配对后最大只报 2.4。
  所以看汇总没问题，**要定位单个错得离谱的原子必须回到逐原子对照**。

用法
----
    python tools/diag_sp3.py               # 全部
    python tools/diag_sp3.py --only ring   # 只看环烷
    python tools/diag_sp3.py --set holdout # 只看未参与标定的分子
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from nmrflow.mol.equivalence import carbon_groups  # noqa: E402
from nmrflow.mol.prepare import prepare  # noqa: E402
from nmrflow.predict.shift_13c import CarbonShiftPredictor  # noqa: E402

#: (名称, SMILES, 类别, 文献 ¹³C 逐 sp³ 原子展开，CDCl₃）
#: 对称分子必须按原子数展开：苯环上的 2 个等价 CH₃ 是两条，不是一条。
#: 展开错了会直接报"原子数不符"，不会悄悄放过——这是这一版比手写原子序号可靠的地方。
CASES: list[tuple[str, str, str, list[float]]] = [
    # ---- 开链烷烃：GP 本来就该准的基线，用来确认改动没把好东西弄坏 ----
    ("丙烷", "CCC", "chain", [15.8, 15.8, 16.1]),
    ("正丁烷", "CCCC", "chain", [13.4, 13.4, 25.2, 25.2]),
    ("正戊烷", "CCCCC", "chain", [13.9, 13.9, 22.6, 22.6, 34.5]),
    ("正己烷", "CCCCCC", "chain", [14.1, 14.1, 22.7, 22.7, 31.9, 31.9]),
    ("乙烷", "CC", "chain", [5.7, 5.7]),
    # ---- 支链烷烃：GP 的 α/β 线性计数在支化点失效的暴露点 ----
    # 异丁烷 CH3 实测 24.3，GP 只给 16.2（少 8）
    ("异丁烷", "CC(C)C", "branch", [24.3, 24.3, 24.3, 25.3]),
    # 新戊烷：季碳 + 3 个等价甲基
    ("新戊烷", "CC(C)(C)C", "branch", [27.9, 31.7, 31.7, 31.7, 31.7]),
    # ---- 支链醇：杂原子与支化叠加 ----
    ("2-丙醇", "CC(C)O", "branch", [25.2, 25.2, 64.5]),
    ("叔丁醇", "CC(C)(C)O", "branch", [31.3, 31.3, 31.3, 69.0]),
    # ---- 环烷：**环修正表的唯一标定依据**。GP 公式里没有任何环项，
    #      环己烷会算成 32.2（实测 27.1）、环戊烷 34.7（实测 25.9）。
    ("环丙烷", "C1CC1", "ring", [5.7, 5.7, 5.7]),
    ("环丁烷", "C1CCC1", "ring", [22.8] * 4),
    ("环戊烷", "C1CCCC1", "ring", [25.9] * 5),
    ("环己烷", "C1CCCCC1", "ring", [27.1] * 6),
    ("环庚烷", "C1CCCCCC1", "ring", [28.7] * 7),
    ("环辛烷", "C1CCCCCCC1", "ring", [26.9] * 8),
    # ---- 烯丙位 / 苄位：`vinyl` / `aryl` 增量的标定依据 ----
    ("丙烯", "CC=C", "allylic", [19.6]),
    ("1-丁烯", "CCC=C", "allylic", [26.5, 13.5]),
    ("1-戊烯", "CCCC=C", "allylic", [33.9, 22.3, 14.1]),
    ("1-己烯", "CCCCC=C", "allylic", [33.8, 31.4, 22.3, 14.1]),
    ("3-甲基-1-丁烯", "CC(C)C=C", "allylic", [31.0, 22.4, 22.4]),
    ("甲苯", "Cc1ccccc1", "benzylic", [21.5]),
    ("乙苯", "CCc1ccccc1", "benzylic", [28.9, 15.6]),
    ("异丙苯", "CC(C)c1ccccc1", "benzylic", [34.3, 24.1, 24.1]),
    # 环烯的 sp³ 碳：环修正与烯烃增量必须同时成立才能对上
    ("环己烯", "C1CC=CCC1", "allylic", [25.4, 25.4, 22.9, 22.9]),
    ("环戊烯", "C1CC=CC1", "allylic", [32.9, 32.9, 23.6]),
]

#: 留出集：**环修正与取代基增量都没有用它们标定**。
#: 唯一检验的是"环修正 + 杂原子/羧基增量能不能正确叠加"——这是最容易被忽略的接缝。
HELDOUT: list[tuple[str, str, str, list[float]]] = [
    ("环己醇", "OC1CCCCC1", "ring_subst", [70.0, 35.4, 35.4, 24.2, 24.2, 25.5]),
    ("环己胺", "NC1CCCCC1", "ring_subst", [50.1, 36.9, 36.9, 25.2, 25.2, 24.8]),
    ("环己甲酸", "OC(=O)C1CCCCC1", "ring_subst",
     [43.3, 29.1, 29.1, 25.9, 25.9, 25.6]),
    ("对二甲苯", "Cc1ccc(C)cc1", "benzylic", [20.9, 20.9]),
    # 羰基/羧基的 α 位：检验"同一个官能团会不会被算两遍"
    ("乙酸甲酯", "CC(=O)OC", "carbonyl", [20.8, 51.7]),
    ("丙酸", "CCC(=O)O", "carbonyl", [27.8, 9.0]),
    ("丁酸", "CCCC(=O)O", "carbonyl", [36.0, 18.4, 13.5]),
    ("2-丁酮", "CCC(C)=O", "carbonyl", [36.8, 29.2, 8.0]),
    ("环己酮", "O=C1CCCCC1", "carbonyl", [42.0, 42.0, 27.0, 27.0, 25.0]),
]


def _run_case(name: str, smiles: str, kind: str, lit: list[float]) -> dict:
    prep = prepare(smiles)
    groups = carbon_groups(prep, 0)
    sigs = CarbonShiftPredictor().predict(prep, groups, "CDCl3")

    pred: list[float] = []
    keys: list[str] = []
    for group, sig in zip(groups, sigs):
        atom = prep.heavy.GetAtomWithIdx(group.heavy_indices[0])
        if atom.GetIsAromatic():
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        for _ in group.heavy_indices:
            pred.append(sig.shift)
            keys.append(sig.provenance.key)

    if len(pred) != len(lit):
        return {"name": name, "kind": kind, "skip": True,
                "why": f"sp³ 原子数不符：预测 {len(pred)} / 文献 {len(lit)}"}

    order = sorted(range(len(pred)), key=lambda i: pred[i])
    errs = [abs(pred[i] - l) for i, l in zip(order, sorted(lit))]
    worst = max(range(len(errs)), key=lambda i: errs[i])
    return {
        "name": name, "kind": kind, "skip": False,
        "pred": sorted(pred), "lit": sorted(lit), "errs": errs,
        "max_err": max(errs), "mean_err": sum(errs) / len(errs),
        "worst_key": keys[order[worst]],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="sp³ 碳诊断")
    ap.add_argument("--set", choices=["calib", "holdout", "all"], default="all")
    ap.add_argument("--only", default=None,
                    help="chain / ring / allylic / benzylic / ring_subst")
    ap.add_argument("--tol", type=float, default=3.0)
    args = ap.parse_args()

    groups = []
    if args.set in ("calib", "all"):
        groups.append(("标定集（3 ppm 以上才算问题）", CASES, args.tol))
    if args.set in ("holdout", "all"):
        groups.append(("留出集（未参与标定）", HELDOUT, args.tol))

    all_errs: list[float] = []
    for title, cases, tol in groups:
        cases = [c for c in cases if args.only is None or c[2] == args.only]
        if not cases:
            continue
        results = [_run_case(*c) for c in cases]
        print("=" * 92)
        print(title)
        print("-" * 92)
        print(f"{'分子':<16}{'类别':<12}{'原子':>5}{'最大':>8}{'平均':>8}   最差原子的来源")
        print("-" * 92)
        for r in results:
            if r.get("skip"):
                print(f"{r['name']:<16}{r['kind']:<12}  {r['why']}")
                continue
            flag = "  <<<" if r["max_err"] > tol else ""
            print(
                f"{r['name']:<16}{r['kind']:<12}{len(r['errs']):>5}"
                f"{r['max_err']:>8.1f}{r['mean_err']:>8.1f}   "
                f"{r['worst_key'][:34]}{flag}"
            )
            all_errs.extend(r["errs"])

        bad = [r for r in results if not r.get("skip") and r["max_err"] > tol]
        if bad:
            print(f"\n误差 > {tol:.0f} ppm（排序后按位对照）")
            for r in bad:
                print(f"\n  {r['name']}  ({r['kind']})")
                print("    预测  " + " ".join(f"{v:>6.1f}" for v in r["pred"]))
                print("    文献  " + " ".join(f"{v:>6.1f}" for v in r["lit"]))
                print("    误差  " + " ".join(f"{v:>6.1f}" for v in r["errs"]))
        print()

    print("-" * 92)
    if all_errs:
        print(
            f"合计 {len(all_errs)} 个 sp³ 原子：平均 {sum(all_errs) / len(all_errs):.2f} ppm，"
            f"最大 {max(all_errs):.1f} ppm，超 3 ppm 占 "
            f"{sum(1 for e in all_errs if e > 3) / len(all_errs):.0%}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
