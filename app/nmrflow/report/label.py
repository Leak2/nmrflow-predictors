"""归属标签：按 SMARTS 优先级给每个信号命名（-CH3 / Ph-H / -NH ...）。"""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem

from ..data.loader import load_toml


@lru_cache(maxsize=2)
def _rules(kind: str) -> list[tuple[str, Chem.Mol]]:
    table = load_toml("labels.toml")
    out = []
    for item in table.get(kind, []):
        patt = Chem.MolFromSmarts(item["smarts"])
        if patt is not None:
            out.append((item["label"], patt))
    return out


def label_for_1h(mol: Chem.Mol, heavy_idx: int) -> str:
    for label, patt in _rules("h"):
        if mol.GetSubstructMatches(patt, useQueryQueryMatches=False) and any(
            heavy_idx in m for m in mol.GetSubstructMatches(patt)
        ):
            return label
    return "-H"


def label_for_13c(mol: Chem.Mol, atom_idx: int) -> str:
    for label, patt in _rules("c"):
        for match in mol.GetSubstructMatches(patt):
            if atom_idx in match:
                return label
    return "-C"


def assign_labels(mol: Chem.Mol, signals, nucleus: str) -> None:
    """就地给信号填 label。"""
    for sig in signals:
        idx = sig.atom_indices[0] if sig.atom_indices else -1
        if idx < 0:
            continue
        if nucleus == "1H":
            heavy = next(
                (n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors()
                 if n.GetAtomicNum() > 1),
                -1,
            )
            sig.label = label_for_1h(mol, heavy) if heavy >= 0 else "-H"
        else:
            sig.label = label_for_13c(mol, idx)
