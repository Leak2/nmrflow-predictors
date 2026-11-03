"""19F 上"学习式聚合"能不能赢过 HOSE 查表的对照实验。

背景：FluoBase 用 GCN 在 3376 条清洗过的实验位移上做到 MAE 3.636 ppm
（Li et al., Artif. Intell. Chem. 2024）。我们自己的 HOSE 在有机 C–F 上是 7.32。
差距来自两处，必须分开看，否则会得出"我们该上深度学习"的错误结论：

1. **数据清洗**——他们的规则是剔除 δ>300 / δ<−400。
2. **模型形态**——HOSE 是"最近邻复制"，GCN 是"平滑插值"。

本脚本只测第 2 条：**沿用同一个 HOSE 码的信息量**（1~4 层球原子类型计数），
换成一个可学习的线性聚合，看能不能压过查表。

- 特征与 HOSE 码同源，信息量对齐 → 差异只来自聚合方式；
- 不引入 torch / dgl / sklearn，只用 numpy，避免为了做实验把项目依赖搞重；
- 划分沿用 eval_tier_ab 的 InChIKey 20% 留出，与 7.32 那组数字可比。

用法::

    python tools/ml_baseline_19f.py --sd /path/to/withsignals.sd
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from eval_tier_ab import (  # noqa: E402
    CALC_PROPS,
    _experimental,
    _inchikey,
    _is_holdout,
    _prep_from_mol,
)

RDLogger.DisableLog("rdApp.*")

MAX_SHELL = 4


def _atom_category(atom: Chem.Atom) -> str:
    """原子的粗类别。刻意做粗——HOSE 码本身也是这种"元素+杂化+芳香性"粒度。"""
    if atom.GetIsAromatic():
        return "ar"
    sym = atom.GetSymbol()
    hyb = atom.GetHybridization().name
    if sym == "C":
        return f"C.{hyb}"
    if sym in ("N", "O", "S", "P", "F", "Cl", "Br", "I"):
        return sym
    return "X"


def _shell_features(mol: Chem.Mol, idx: int) -> dict[str, float]:
    """以 idx 为中心，按拓扑距离 1~MAX_SHELL 统计邻居类别。

    这就是 HOSE 码在做的事（逐层记录球形环境），只是这里落成向量而非字符串。
    """
    seen = {idx: 0}
    frontier = [idx]
    for depth in range(1, MAX_SHELL + 1):
        nxt = []
        for a in frontier:
            for nb in mol.GetAtomWithIdx(a).GetNeighbors():
                j = nb.GetIdx()
                if j in seen:
                    continue
                seen[j] = depth
                nxt.append(j)
        frontier = nxt

    feats: dict[str, float] = {}
    for j, d in seen.items():
        if d == 0:
            continue
        cat = _atom_category(mol.GetAtomWithIdx(j))
        feats[f"d{d}:{cat}"] = feats.get(f"d{d}:{cat}", 0.0) + 1.0

    atom = mol.GetAtomWithIdx(idx)
    feats["n_F_total"] = float(
        sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "F")
    )
    feats["n_atoms"] = float(mol.GetNumAtoms())
    feats["n_rings"] = float(mol.GetRingInfo().NumRings())
    feats["f_charge"] = float(atom.GetFormalCharge())
    feats["neighbor_deg"] = float(
        sum(nb.GetDegree() for nb in atom.GetNeighbors())
    )
    return feats


def _collect(sd_path: Path, holdout: float, salt: str):
    """返回 (train_X, train_y, test_X, test_y, test_is_cf)。"""
    train: list[tuple[dict[str, float], float]] = []
    test: list[tuple[dict[str, float], float, bool]] = []
    skipped = 0

    supp = Chem.SDMolSupplier(str(sd_path), removeHs=False, sanitize=False)
    for mol in supp:
        if mol is None or any(mol.HasProp(p) for p in CALC_PROPS):
            skipped += 1
            continue
        key = _inchikey(mol)
        if key is None:
            skipped += 1
            continue
        try:
            Chem.SanitizeMol(mol)
            _prep_from_mol(mol)                   # 与评测脚本同样的可用性检查
        except Exception:
            skipped += 1
            continue
        exp = _experimental(mol, "19F")
        if not exp:
            continue
        # _is_holdout 为 True 表示属于留出集（与 eval_tier_ab 的用法一致）
        is_test = _is_holdout(key, holdout, salt)
        for idx, shift in exp.items():
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "F":
                continue
            feats = _shell_features(mol, idx)
            if is_test:
                is_cf = all(nb.GetAtomicNum() == 6 for nb in atom.GetNeighbors())
                test.append((feats, shift, is_cf))
            else:
                train.append((feats, shift))
    return train, test, skipped


def _matrix(rows: list[dict[str, float]], vocab: dict[str, int]) -> np.ndarray:
    m = np.zeros((len(rows), len(vocab)), dtype=np.float64)
    for i, feats in enumerate(rows):
        for k, v in feats.items():
            j = vocab.get(k)
            if j is not None:
                m[i, j] = v
    return m


def _ridge(X: np.ndarray, y: np.ndarray, lam: float) -> np.ndarray:
    """带截距的岭回归（闭式解，无需 sklearn）。"""
    X1 = np.hstack([X, np.ones((X.shape[0], 1))])
    reg = np.eye(X1.shape[1]) * lam
    reg[-1, -1] = 0.0                              # 截距不惩罚
    return np.linalg.solve(X1.T @ X1 + reg, X1.T @ y)


def main() -> int:
    ap = argparse.ArgumentParser(description="19F 学习式聚合 vs HOSE 查表")
    ap.add_argument("--sd", type=Path, required=True)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--salt", default="nmrflow-holdout")
    ap.add_argument("--lambda", dest="lam", type=float, default=1.0)
    args = ap.parse_args()

    train, test, skipped = _collect(args.sd, args.holdout, args.salt)
    print(f"训练 {len(train):,} 条 / 留出 {len(test):,} 条（跳过 {skipped:,} 个分子）")
    if not train or not test:
        print("样本不足")
        return 1

    vocab: dict[str, int] = {}
    for feats, _ in train:
        for k in feats:
            vocab.setdefault(k, len(vocab))
    print(f"特征维度 {len(vocab)}")

    Xtr = _matrix([f for f, _ in train], vocab)
    ytr = np.array([y for _, y in train], dtype=np.float64)
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd[sd == 0] = 1.0
    Xtr = (Xtr - mu) / sd

    w = _ridge(Xtr, ytr, args.lam)

    Xte = _matrix([f for f, _, _ in test], vocab)
    Xte = (Xte - mu) / sd
    pred = np.hstack([Xte, np.ones((Xte.shape[0], 1))]) @ w
    err = np.abs(pred - np.array([y for _, y, _ in test]))
    cf = np.array([c for _, _, c in test])

    def report(label: str, mask) -> None:
        e = err[mask]
        print(
            f"  {label:<12}{e.size:>6} 条  MAE {e.mean():>6.2f}  "
            f"P90 {np.percentile(e, 90):>6.1f}  最大 {e.max():>6.1f}"
        )

    print("\n本次实验（线性聚合，特征信息量与 HOSE 码对齐）")
    report("全部", np.ones_like(cf, dtype=bool))
    report("有机 C–F", cf)

    # 训练集残差，用于判断是欠拟合还是数据本身噪声大
    fit = Xtr @ w[:-1] + w[-1]
    fit_err = np.abs(fit - ytr)
    print(
        f"\n训练集 MAE {fit_err.mean():.2f}（留出 {err.mean():.2f}）"
        " —— 两者接近说明瓶颈在信息量/噪声，不在模型容量"
    )
    print("\n对照：HOSE 查表在同一留出集上是 8.76（有机 C–F 7.32）、Tier A 22.37")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
