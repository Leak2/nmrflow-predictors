"""从 NMRShiftDB2 里统计 ³¹P / ¹⁹F 的环境-位移分布，用来补 Tier A 的环境表。

为什么要从数据库里挖，而不是照着文献手写
--------------------------------------
上一轮把 Tier A 的兜底率报成 56% / 24% 时，我只看了 `provenance.key == "fallback"`，
没有核对"到底哪些环境没命中"。手写 SMARTS 的问题是：写出来的模式可能看起来合理，
但**没覆盖真实分子里存在的环境**，而这种缺失在单测里看不出来——只有拿整批数据统计
才能发现。

所以这个脚本做两件事：

1. 把每个 P / F 原子按"邻居签名"归类（`P(=O)(OH)(OR)(C)` 这种），
   统计每类的记录数与实测位移分布。**没有对应 SMARTS 的签名就是待补充的条目**。
2. 对已经命中表项的原子，比较表值与库内中位数，看哪些条目的 value 偏了。

两者都只在**训练集**（非 holdout）上做，避免把评估集的信息写进表里。

用法
----
    python tools/mine_hetero.py --sd ... --nucleus 31P
    python tools/mine_hetero.py --sd ... --nucleus 19F --top 30
"""

from __future__ import annotations

import argparse
import statistics
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(_ROOT / "tools"))

from build_hose_index import CALC_PROPS, _atom_matches_nucleus, _iter_shift_entries  # noqa: E402
from eval_hose_index import _inchikey, _is_holdout  # noqa: E402
from nmrflow.data.loader import shift_table  # noqa: E402


def _substituent_tag(mol: Chem.Mol, center: int, nb: Chem.Atom) -> str:
    """把一个邻居写成能直接翻译成 SMARTS 的短标签。

    ★ 氢必须用 ``GetTotalNumHs()`` 判断，不能只看显式邻居：
      ``O=P(O)(c1ccccc1)c1ccccc1`` 里的 -OH 没有显式 H，
      只看邻居会把它和一个真正的 O⁻ 混成同一个签名，把同一类拆成两堆
      （实测这会让次膦酸的 38.6 与 30.0 两个中位数看起来像两个不同的环境）。
    """
    sym = nb.GetSymbol()
    bond = mol.GetBondBetweenAtoms(center, nb.GetIdx())
    order = bond.GetBondTypeAsDouble() if bond is not None else 1.0
    n_h = nb.GetTotalNumHs()
    if sym == "H":
        return "-H"
    if order >= 1.9 and sym in ("O", "S"):
        return f"={sym}"
    if sym == "O":
        if n_h > 0:
            return "-OH"
        rest = [n for n in nb.GetNeighbors() if n.GetIdx() != center]
        if rest and {n.GetSymbol() for n in rest} <= {"C", "c"}:
            return "-OR"
        return "-O"
    if sym == "N":
        # 氨基按取代数分档：-NH2 / -NHR / -NR2
        rest = sum(1 for n in nb.GetNeighbors() if n.GetIdx() != center)
        return "-NH2" if n_h >= 2 else "-NHR" if n_h == 1 and rest else "-NR2"
    if sym == "C":
        # ★ 芳香碳必须用 GetIsAromatic() 判断：GetSymbol() 永远返回 "C"，
        #   不会返回 SMILES 里的小写 "c"。早期写成 sym == "c" 导致 -Ar 永不触发，
        #   芳基和烷基被混成同一类——这正是"二甲基膦酸实测 55、类中位数却是 30"
        #   的原因（30 是被二芳基的 ~30 拉过去的）。
        return "-Ar" if nb.GetIsAromatic() else "-C"
    return f"-{sym}"


def _signature(mol: Chem.Mol, idx: int) -> str:
    """规范签名（标签已排序，原子编号无关）。

    ★ ¹⁹F 这类**端基核**不能只看自己的邻居：氟只连一个原子，
      邻居签名会退化成 ``F(-C)``，把所有含氟环境压成一类。
      所以端基卤素改描述**母原子**的环境：``F-C(-F)(-F)(-F)`` = CF3、
      ``F-C(-F)(-F)(-C)`` = CF2、``F-C(-Ar)(-Ar)`` = 芳环上的氟。
    """
    atom = mol.GetAtomWithIdx(idx)
    charge = "+" if atom.GetFormalCharge() > 0 else ""

    heavy = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1]
    if atom.GetSymbol() in ("F", "Cl", "Br", "I") and len(heavy) == 1:
        parent = heavy[0]
        tags = sorted(
            _substituent_tag(mol, parent.GetIdx(), nb)
            for nb in parent.GetNeighbors()
            if nb.GetIdx() != idx
        )
        return f"{atom.GetSymbol()}{charge}-{parent.GetSymbol()}({''.join(tags)})"

    tags = sorted(_substituent_tag(mol, idx, nb) for nb in atom.GetNeighbors())
    return f"{atom.GetSymbol()}{charge}(" + "".join(tags) + ")"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sd", required=True, type=Path)
    ap.add_argument("--nucleus", default="31P")
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--salt", default="nmrflow-hose-v1")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--min-count", type=int, default=5)
    ap.add_argument("--show", type=int, default=0, help="打印前 N 个签名的样例 SMILES")
    args = ap.parse_args()

    table = shift_table(args.nucleus)
    envs = []
    for item in table.get("env", []):
        patt = Chem.MolFromSmarts(item["smarts"])
        if patt is not None:
            envs.append((item, patt))

    # 签名 -> 位移列表（只统计没被任何 SMARTS 命中的原子）
    missing: dict[str, list[float]] = defaultdict(list)
    examples: dict[str, list[str]] = defaultdict(list)
    # 表项名 -> 位移列表（被命中的原子）
    covered: dict[str, list[float]] = defaultdict(list)
    total = 0

    supp = Chem.SDMolSupplier(str(args.sd), removeHs=False, sanitize=False)
    seen = 0
    for mol in supp:
        seen += 1
        if mol is None or any(mol.HasProp(p) for p in CALC_PROPS):
            continue
        key = _inchikey(mol)
        if key is None or _is_holdout(key, args.holdout, args.salt):
            continue                      # 只用训练集，别把评估集写进表
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue

        n_atoms = mol.GetNumAtoms()
        for prop in mol.GetPropNames():
            if not prop.startswith(f"Spectrum {args.nucleus} "):
                continue
            for shift, idx in _iter_shift_entries(mol.GetProp(prop)):
                if not (0 <= idx < n_atoms):
                    continue
                if not _atom_matches_nucleus(mol.GetAtomWithIdx(idx), args.nucleus):
                    continue
                total += 1
                hit = None
                for item, patt in envs:
                    # 与 shift_hetero 的匹配规则保持一致：match[0] 必须是被观测的那个核。
                    # ★ uniquify=False 同样不可省：默认值会把 CF₃ 三个等价氟折叠成一个
                    #   匹配，导致"未命中"被重复计 3 次，把缺口规模估大。
                    for match in mol.GetSubstructMatches(
                        patt, uniquify=False, maxMatches=20000
                    ):
                        if match and match[0] == idx:
                            hit = item
                            break
                    if hit:
                        break
                if hit is not None:
                    covered[hit["name"]].append(shift)
                else:
                    sig = _signature(mol, idx)
                    missing[sig].append(shift)
                    if len(examples[sig]) < 3:
                        try:
                            examples[sig].append(Chem.MolToSmiles(mol))
                        except Exception:
                            pass

        if seen % 20000 == 0:
            print(f"  {seen:,} 条 ...", flush=True)

    print("=" * 96)
    print(f"{args.nucleus}：训练集共 {total:,} 个原子，"
          f"未命中任何表项 {sum(len(v) for v in missing.values()):,} "
          f"({100 * sum(len(v) for v in missing.values()) / max(1, total):.1f}%)")
    print("=" * 96)

    fallback_value = float(table.get("fallback", {}).get("value", 0.0))

    print(f"\n【待补】未命中表项的邻居签名（记录数 ≥ {args.min_count}）")
    print("  中位MAE = 用该类中位数预测的 MAE；兜底MAE = 用 fallback 值预测的 MAE；")
    print("  **只有中位MAE 明显更小才值得新增条目**——否则等于把兜底常数换个写法。")
    print(f"{'签名':<40}{'记录':>6}{'中位':>8}{'中位MAE':>9}{'兜底MAE':>9}{'改善':>7}")
    print("-" * 96)
    rows = [(sig, v) for sig, v in missing.items() if len(v) >= args.min_count]
    ranked = []
    for sig, vals in rows:
        med = statistics.median(vals)
        mae_med = sum(abs(v - med) for v in vals) / len(vals)
        mae_fb = sum(abs(v - fallback_value) for v in vals) / len(vals)
        ranked.append((mae_fb - mae_med, sig, len(vals), med, mae_med, mae_fb))
    for gain, sig, n, med, mae_med, mae_fb in sorted(ranked, reverse=True)[: args.top]:
        print(f"{sig:<40}{n:>6}{med:>8.1f}{mae_med:>9.1f}{mae_fb:>9.1f}{gain:>7.1f}")

    if args.show:
        print("\n  样例（供人工确认化学类别）：")
        for gain, sig, *_ in sorted(ranked, reverse=True)[: args.show]:
            print(f"    {sig}")
            for smi in examples[sig]:
                print(f"        {smi[:78]}")

    print(f"\n【校准】已命中表项的实际分布 vs 表里的 value")
    print(f"{'表项':<40}{'记录':>7}{'表值':>8}{'中位':>8}{'P10':>8}{'P90':>8}{'半宽':>7}")
    print("-" * 100)
    for name, vals in sorted(covered.items(), key=lambda kv: -len(kv[1])):
        if len(vals) < args.min_count:
            continue
        table_value = next(
            (float(i["value"]) for i in table["env"] if i["name"] == name), float("nan")
        )
        ordered = sorted(vals)
        med = statistics.median(ordered)
        p10 = ordered[int(0.1 * (len(ordered) - 1))]
        p90 = ordered[int(0.9 * (len(ordered) - 1))]
        # 半宽 = (P90-P10)/2，直接对应表里的 width（"环境跨度约 ±width ppm"）。
        # 用 max-min 会被个别离群值撑爆，用平均绝对偏差又会偏窄。
        half = (p90 - p10) / 2
        print(f"{name:<40}{len(vals):>7}{table_value:>8.1f}{med:>8.1f}"
              f"{p10:>8.1f}{p90:>8.1f}{half:>7.0f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
