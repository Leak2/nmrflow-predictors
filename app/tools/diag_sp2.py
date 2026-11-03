"""sp² 碳（烯碳 / 稠合芳碳）的 Tier A 诊断。

为什么单独写这个
----------------
`shift_13c.toml` 里有一条
``alkene = { value = 125.0, smarts = "[CX3;$(C=[CX3]);!$(C=[OX1])]" }``，
它把**所有**非羰基 sp² 碳一律锚定到 125.0 ppm 并标 ``confidence=high``。

这条规则在乙烯附近碰巧是对的（123.3），但烯碳的真实跨度是 80~160 ppm：
甲基乙烯基醚的 =CH 是 153.4，乙烯基乙酸酯的 =CH₂ 只有 97.7。
一个常数不可能覆盖，而且它拿到的是 **high** —— 错得越离谱标得越自信。

稠合芳环有同类问题：``aromatic.base = 128.5`` 只对单环成立，
萘的桥头碳 133.5、β 位 125.9，全都被算成 128.5。

比较方式
--------
每个分子给出**该分子全部 sp² 碳**的文献位移，展开成"一原子一项"的列表，
然后与预测值各自排序、按位相减。

不用 ``{原子序号: 位移}`` 是因为 SMILES 原子序号与文献编号的对应关系极易写错
（联苯的 ipso 就不是 0 号，RDKit 的环起点也不在取代位），而写错的对照表会
把正确预测报成错误。也不用"谱峰集合配对"——那要给每个预测值挑最近的文献值，
挑的过程会掩盖错误。

排序对齐没有选择余地：多一个错值就会整体错位，误差立刻显形。

用法
----
    python tools/diag_sp2.py                 # 全部
    python tools/diag_sp2.py --only alkene   # 只看烯烃
    python tools/diag_sp2.py --tol 3         # 收紧阈值
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

#: (名称, SMILES, 类别, 文献 ¹³C 逐 sp² 原子展开，CDCl₃）
CASES: list[tuple[str, str, str, list[float]]] = [
    # ---------------- 链状烯烃 ----------------
    ("乙烯", "C=C", "alkene", [123.3, 123.3]),
    # =CH 133.4 / =CH2 115.7
    ("丙烯", "CC=C", "alkene", [133.4, 115.7]),
    ("1-丁烯", "CCC=C", "alkene", [139.8, 113.8]),
    ("异丁烯", "CC(=C)C", "alkene", [143.4, 110.4]),
    ("反-2-丁烯", "C/C=C/C", "alkene", [124.5, 124.5]),
    # 末端 117.0 / 内侧 137.5
    ("1,3-丁二烯", "C=CC=C", "alkene", [117.0, 137.5, 137.5, 117.0]),
    ("异戊二烯", "CC(=C)C=C", "alkene", [142.5, 114.6, 137.6, 117.3]),
    # ---------------- 推电子取代：偏移最大 ----------------
    ("甲基乙烯基醚", "COC=C", "alkene", [153.4, 86.5]),
    ("乙基乙烯基醚", "CCOC=C", "alkene", [151.3, 86.4]),
    ("乙烯基乙酸酯", "CC(=O)OC=C", "alkene", [170.4, 141.3, 97.7]),
    # ---------------- 吸电子取代 ----------------
    ("丙烯酸甲酯", "COC(=O)C=C", "alkene", [167.0, 130.4, 128.2]),
    ("丙烯腈", "C=CC#N", "alkene", [117.9, 132.6]),
    # ---------------- 环烯 ----------------
    ("环己烯", "C1CC=CCC1", "alkene", [127.1, 127.1]),
    ("环戊烯", "C1CC=CC1", "alkene", [130.9, 130.9]),
    ("1,3-环己二烯", "C1C=CC=CC1", "alkene", [124.5, 133.0, 124.5, 133.0]),
    # ---------------- 单环芳烃基线 ----------------
    ("苯", "c1ccccc1", "arene", [128.5] * 6),
    # ipso 137.9 / o 129.1×2 / m 128.2×2 / p 125.3
    ("甲苯", "Cc1ccccc1", "arene",
     [137.9, 129.1, 129.1, 128.2, 128.2, 125.3]),
    # ipso 159.9 / o 114.1×2 / m 129.7×2 / p 120.7
    ("苯甲醚", "COc1ccccc1", "arene",
     [159.9, 114.1, 114.1, 129.7, 129.7, 120.7]),
    ("苯酚", "Oc1ccccc1", "arene",
     [154.9, 115.4, 115.4, 129.8, 129.8, 121.1]),
    # 乙烯基 113.8 / 136.9；环 ipso 137.6 / o 126.1×2 / m 128.4×2 / p 127.5
    ("苯乙烯", "C=Cc1ccccc1", "arene",
     [113.8, 136.9, 137.6, 126.1, 126.1, 128.4, 128.4, 127.5]),
    # 联苯：ipso 141.2×2 / m 128.8×4 / o 127.2×4 / p 127.2×2
    ("联苯", "c1ccc(-c2ccccc2)cc1", "arene",
     [141.2, 141.2, 128.8, 128.8, 128.8, 128.8,
      127.2, 127.2, 127.2, 127.2, 127.2, 127.2]),
    # ---------------- 稠合芳烃 ----------------
    # 萘：桥头 133.5×2 / α 128.0×4 / β 125.9×4
    ("萘", "c1ccc2ccccc2c1", "fused",
     [133.5, 133.5, 128.0, 128.0, 128.0, 128.0,
      125.9, 125.9, 125.9, 125.9]),
    # 蒽：桥头 131.9×4 / meso 125.3×2 / α 128.4×4 / β 126.5×4
    ("蒽", "c1ccc2cc3ccccc3cc2c1", "fused",
     [131.9, 131.9, 131.9, 131.9, 125.3, 125.3,
      128.4, 128.4, 128.4, 128.4,
      126.5, 126.5, 126.5, 126.5]),
]


#: 留出集：**没有参与 [alkene] / [aromatic] 增量的标定**，用来检查上面的数字不是拟合出来的。
#: 只是品种太少（手头能确认文献值的就这些），所以它只能否定、不能证明。
HELDOUT: list[tuple[str, str, str, list[float]]] = [
    # 1-戊烯 / 1-己烯：长链末端烯烃，用来验证 beta_step 的饱和假设
    ("1-戊烯", "CCCC=C", "alkene", [139.4, 114.2]),
    ("1-己烯", "CCCCC=C", "alkene", [139.2, 114.1]),
    # 甲基丙烯酸甲酯：α 位有甲基又有酯基，两个取代基同碳
    ("甲基丙烯酸甲酯", "COC(=O)C(=C)C", "alkene", [167.4, 136.4, 125.0]),
    # 苯胺 / 氯苯 / 硝基苯：单环芳烃的给电子/弱吸电子/强吸电子三档
    ("苯胺", "Nc1ccccc1", "arene",
     [146.3, 115.1, 115.1, 129.3, 129.3, 118.5]),
    ("氯苯", "Clc1ccccc1", "arene",
     [134.3, 128.8, 128.8, 129.5, 129.5, 126.9]),
    ("硝基苯", "[O-][N+](=O)c1ccccc1", "arene",
     [148.3, 123.5, 123.5, 129.4, 129.4, 134.7]),
    # 苯甲酸甲酯：验证 COOR 的 i/o/m/p 增量
    ("苯甲酸甲酯", "COC(=O)c1ccccc1", "arene",
     [166.9, 130.6, 129.7, 129.7, 128.4, 128.4, 132.8]),
]


def _run_case(name: str, smiles: str, kind: str, lit: list[float]) -> dict:
    prep = prepare(smiles)
    groups = carbon_groups(prep, 0)
    sigs = CarbonShiftPredictor().predict(prep, groups, "CDCl3")

    # 每个等价组有多少原子，就展开多少项；**只取 sp² 碳**
    # （sp³ 的烯丙位误差是另一类问题，混进来会让结论失焦）
    pred: list[float] = []
    keys: list[str] = []
    for g, s in zip(groups, sigs):
        atom = prep.heavy.GetAtomWithIdx(g.heavy_indices[0])
        if atom.GetHybridization() != Chem.HybridizationType.SP2:
            continue
        for _ in g.heavy_indices:
            pred.append(s.shift)
            keys.append(s.provenance.key)

    if len(pred) != len(lit):
        return {
            "name": name, "kind": kind, "skip": True,
            "why": f"原子数不符：预测 {len(pred)} / 文献 {len(lit)}",
        }

    order = sorted(range(len(pred)), key=lambda i: pred[i])
    errs = [abs(pred[i] - l) for i, l in zip(order, sorted(lit))]
    worst = max(range(len(errs)), key=lambda i: errs[i])
    return {
        "name": name, "kind": kind, "skip": False,
        "pred": sorted(pred), "lit": sorted(lit),
        "errs": errs, "max_err": max(errs),
        "mean_err": sum(errs) / len(errs),
        "worst_key": keys[order[worst]],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="sp² 碳诊断")
    ap.add_argument(
        "--set", choices=["calib", "holdout", "all"], default="all",
        help="calib=参与标定的分子；holdout=未参与标定；all=两者都跑",
    )
    ap.add_argument("--only", choices=["alkene", "arene", "fused"], default=None)
    ap.add_argument("--tol", type=float, default=5.0, help="超过该误差逐条列出")
    args = ap.parse_args()

    groups = []
    if args.set in ("calib", "all"):
        groups.append(("标定集（增量由这批分子拟合，误差偏乐观）", CASES))
    if args.set in ("holdout", "all"):
        groups.append(("留出集（未参与标定，才是可信的数字）", HELDOUT))

    all_errs: list[float] = []
    for title, cases in groups:
        cases = [c for c in cases if args.only is None or c[2] == args.only]
        if not cases:
            continue
        results = [_run_case(*c) for c in cases]
        print("=" * 86)
        print(title)
        print("-" * 86)
        print(f"{'分子':<20}{'类别':<8}{'原子数':>7}{'最大误差':>10}{'平均误差':>10}   最差原子的来源")
        print("-" * 86)
        for r in results:
            if r.get("skip"):
                print(f"{r['name']:<20}{r['kind']:<8}  {r['why']}")
                continue
            flag = "  <<<" if r["max_err"] > args.tol else ""
            print(
                f"{r['name']:<20}{r['kind']:<8}{len(r['errs']):>7}"
                f"{r['max_err']:>10.1f}{r['mean_err']:>10.1f}   "
                f"{r['worst_key'][:28]}{flag}"
            )
            all_errs.extend(r["errs"])

        bad = [r for r in results if not r.get("skip") and r["max_err"] > args.tol]
        if bad:
            print(f"\n误差 > {args.tol:.0f} ppm（排序后按位对照）")
            for r in bad:
                print(f"\n  {r['name']}  ({r['kind']})")
                print("    预测  " + " ".join(f"{v:>6.1f}" for v in r["pred"]))
                print("    文献  " + " ".join(f"{v:>6.1f}" for v in r["lit"]))
                print("    误差  " + " ".join(f"{v:>6.1f}" for v in r["errs"]))
        print()

    print("-" * 86)
    if all_errs:
        print(
            f"合计 {len(all_errs)} 个 sp² 原子：平均 {sum(all_errs) / len(all_errs):.2f} ppm，"
            f"最大 {max(all_errs):.1f} ppm，超 {args.tol:.0f} ppm 占 "
            f"{sum(1 for e in all_errs if e > args.tol) / len(all_errs):.0%}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
