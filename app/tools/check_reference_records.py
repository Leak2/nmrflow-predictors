"""核对 SD 里的简单双环参考记录，防止"用错记录"悄悄改掉模型。

为什么需要
----------
0.2.3 用十氢萘/降冰片烷标定了三个常数，0.2.6 又因为 SD 里 44.2 的十氢萘记录把这
三个常数删掉，0.2.8 又因为"44.2 不合理"把它们恢复。三次改动都基于**对同一条记录
的不同解读**，而中间没有任何一步去核对文献值。

这个脚本把那一步补上：把 SD 里目标分子的每条记录，按"对称等价碳数"聚成
`(位移, 碳数)`，与文献值逐项比对，明确指出某条记录是 **trans / cis / 无法解释**。

文献值来源（写在这里，便于复核）
--------------------------------
- **降冰片烷**（bicyclo[2.2.1]heptane）：Lippmaa & Paasivirta,
  *Carbon-13 Chemical Shifts of Bicyclic Compounds*, 1970, Table 2 第 1 条——
  C1/C4 37.3、C2/C3/C5/C6 30.6、**C7 39.2**。
  注意 C7 是**最去屏蔽**的碳；写成 24.5（比其它碳都高场）在化学上不成立。
- **十氢萘**：桥头碳是判断顺反的标志——**trans 44.0、cis 36.0**，相差约 8 ppm。
  trans 的完整三段为 44.0(2C, CH) / 34.7(4C) / 27.1(4C)；cis 只锚桥头 36.0，
  其余两段因文献归属不统一，只做区间检查（24~31）。

**一个必须记住的坑**：`InChIKey` 的前 14 位是**骨架哈希**，cis/trans 十氢萘完全相同
（`NNBZCPXTIHJBJL`）。只比这 14 位会看到"同一分子两条记录互相矛盾"，进而把两条
都当成脏数据。本脚本按完整 27 位打印，并单独标注骨架冲突。

用法
----
    python tools/check_reference_records.py --sd <nmrshiftdb2withsignals.sd>
    python tools/check_reference_records.py --sd <SD> --json
    python tools/check_reference_records.py --sd <SD> --warn-only

退出码：0 = 每个目标分子都有记录能对上文献；1 = 有分子对不上或出现无法解释的记录。
`--warn-only` 总是退出 0（只打印，不阻断）。
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

NUCLEUS = "13C"
#: 聚类容差：同一条记录的等价碳之间通常 <0.1 ppm，0.3 足够分开真实的两个环境。
CLUSTER_TOL = 0.3


# --------------------------------------------------------------- 参考表 ---
#: ``pattern`` 用于整体匹配；``anchor`` 只用于判定这个记录属于哪个异构体。
REFS: dict[str, dict] = {
    "decalin": {
        "label": "十氢萘 / decalin",
        "smiles": "C1CCC2CCCCC2C1",
        "n_carbons": 10,
        "tol": 2.0,
        "isomers": {
            "trans": {
                "pattern": [(44.0, 2), (34.7, 4), (27.1, 4)],
                "note": "桥头 44.0（文献值，顺反判据相差约 8 ppm）",
                "required": True,
            },
            "cis": {
                "anchor": (36.0, 2),
                "ch2_range": (24.0, 31.0),
                "note": "只锚桥头 36.0；其余两段文献归属不统一",
                "required": False,
            },
        },
    },
    "norbornane": {
        "label": "降冰片烷 / norbornane",
        "smiles": "C1CC2CCC1C2",
        "n_carbons": 7,
        "tol": 2.0,
        "isomers": {
            "parent": {
                "pattern": [(37.3, 2), (30.6, 4), (39.2, 1)],
                "note": "Lippmaa 1970：C1/C4 37.3、C2/C3/C5/C6 30.6、C7 39.2",
                "required": True,
            }
        },
    },
}


# --------------------------------------------------------------- 解析 ---
def _iter_shift_entries(body: str):
    """解析 ``位移;强度;原子序号|...``。"""
    for chunk in body.split("|"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(";")
        if len(parts) < 3:
            continue
        try:
            yield float(parts[0]), int(float(parts[2]))
        except ValueError:
            continue


def _experimental_carbon(mol: Chem.Mol) -> dict[int, float]:
    """``{原子索引: δ}``；只取碳原子、同一原子多次出现取平均。"""
    acc: dict[int, list[float]] = defaultdict(list)
    n_atoms = mol.GetNumAtoms()
    for prop in mol.GetPropNames():
        if not prop.startswith(f"Spectrum {NUCLEUS} "):
            continue
        for shift, idx in _iter_shift_entries(mol.GetProp(prop)):
            if 0 <= idx < n_atoms and mol.GetAtomWithIdx(idx).GetSymbol() == "C":
                acc[idx].append(shift)
    return {k: sum(v) / len(v) for k, v in acc.items()}


def _clusters(shifts: dict[int, float]) -> list[tuple[float, int]]:
    """按位移把等价碳聚成 ``(均值, 个数)``——个数就是对称性证据。"""
    out: list[list[float]] = []
    for value in sorted(shifts.values()):
        if out and value - out[-1][2] <= CLUSTER_TOL:
            total = out[-1][0] * out[-1][1] + value
            out[-1][1] += 1
            out[-1][2] = value
            out[-1][0] = total / out[-1][1]
        else:
            out.append([float(value), 1.0, float(value)])
    return [(round(v, 2), int(n)) for v, n, _ in out]


def _canonical(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=False) if mol else None


# --------------------------------------------------------------- 判定 ---
def _match_pattern(
    observed: list[tuple[float, int]], expected: list[tuple[float, int]], tol: float
) -> tuple[bool, float, list[str]]:
    """按"碳数相同"配对；返回 (是否全中, 最大偏差, 未配上的描述)。"""
    pool = list(observed)
    worst = 0.0
    misses: list[str] = []
    for want, count in expected:
        candidates = [o for o in pool if o[1] == count]
        if not candidates:
            misses.append(f"缺 {want:g}x{count}")
            continue
        value = min(candidates, key=lambda o: abs(o[0] - want))
        pool.remove(value)
        worst = max(worst, abs(value[0] - want))
    ok = not misses and worst <= tol
    return ok, worst, misses


def _classify(clusters: list[tuple[float, int]], spec: dict) -> str:
    """按参考表判断这条记录像哪个异构体；认不出就返回 ``unknown``。"""
    for name, ref in spec["isomers"].items():
        if "pattern" in ref:
            ok, _, _ = _match_pattern(clusters, ref["pattern"], spec["tol"])
            if ok:
                return name
        else:
            anchor, count = ref["anchor"]
            heads = [c for c in clusters if c[1] == count and abs(c[0] - anchor) <= spec["tol"]]
            rest = [c for c in clusters if c not in heads]
            lo, hi = ref["ch2_range"]
            if heads and rest and all(lo <= c[0] <= hi for c in rest):
                return name
    return "unknown"


def check(sd_path: Path) -> dict:
    targets = {name: _canonical(spec["smiles"]) for name, spec in REFS.items()}
    hits: dict[str, list[dict]] = defaultdict(list)

    supplier = Chem.SDMolSupplier(str(sd_path), removeHs=False, sanitize=False)
    for mol in supplier:
        if mol is None:
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue
        flat = Chem.MolToSmiles(mol, isomericSmiles=False)
        for name, key in targets.items():
            if flat != key:
                continue
            shifts = _experimental_carbon(mol)
            if not shifts:
                continue
            clusters = _clusters(shifts)
            spec = REFS[name]
            row = {
                "clusters": clusters,
                "isomer": _classify(clusters, spec),
                "n_carbons": int(sum(n for _, n in clusters)),
                "ok_shape": sum(n for _, n in clusters) == spec["n_carbons"],
            }
            hits[name].append(row)

    report = {"sd": str(sd_path), "targets": {}}
    for name, spec in REFS.items():
        rows = hits.get(name, [])
        found = [r["isomer"] for r in rows]
        need = {k for k, v in spec["isomers"].items() if v.get("required")}
        good = sorted(need & set(found))
        unexplained = [r for r in rows if r["isomer"] == "unknown"]
        bad_shape = [r for r in rows if not r["ok_shape"]]
        ok = bool(rows) and good == sorted(need) and not unexplained and not bad_shape
        report["targets"][name] = {
            "label": spec["label"],
            "records": len(rows),
            "isomers_found": found,
            "matched": good,
            "missing": sorted(need - set(found)),
            "unexplained": len(unexplained),
            "bad_shape": len(bad_shape),
            "ok": ok,
            "detail": rows,
        }
    report["ok"] = all(t["ok"] for t in report["targets"].values())
    return report


def _print_text(report: dict) -> None:
    print(f"SD: {report['sd']}\n")
    for name, t in report["targets"].items():
        flag = "PASS" if t["ok"] else "FAIL"
        print(f"[{flag}] {t['label']}  记录 {t['records']} 条")
        for i, row in enumerate(t["detail"], 1):
            pattern = "  ".join(f"{v:g}x{n}" for v, n in row["clusters"])
            shape = "" if row["ok_shape"] else "  <-- 碳数不对，记录可能损坏"
            print(f"    #{i} {row['isomer']:<8} {pattern}{shape}")
        if t["missing"]:
            print(f"    缺: {', '.join(t['missing'])}")
        if t["unexplained"]:
            print("    有记录既不像 trans 也不像 cis —— 需要人工看一眼")
        if not t["records"]:
            print("    这个分子在 SD 里一条记录都没有")
        print()


def main() -> int:
    ap = argparse.ArgumentParser(description="核对 SD 里双环参考记录")
    ap.add_argument("--sd", required=True, type=Path)
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--warn-only", action="store_true", help="总是退出 0")
    args = ap.parse_args()

    if not args.sd.exists():
        print(f"找不到 {args.sd}", file=sys.stderr)
        return 2
    report = check(args.sd)
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        _print_text(report)
    if args.warn_only:
        return 0
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
