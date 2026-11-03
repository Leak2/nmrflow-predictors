"""活泼氢（OH / NH / SH / COOH …）分类。

分类结果用于在 solvents.toml 的 labile 表里查"位移中心 + 线宽"，
因为活泼氢的位移强烈依赖溶剂，且不能靠加和法预测。
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

# 分类名 -> 说明
LABILE_CLASSES = {
    "ROH": "醇羟基",
    "ArOH": "酚羟基",
    "COOH": "羧酸",
    "CONH": "酰胺 N-H",
    "NH2": "胺 N-H",
    "SH": "硫醇",
    "ArNH": "杂环 N-H（吡咯/吲哚等）",
    "SO2NH": "磺酰胺 N-H",
    "oxime": "肟 N-OH",
    "enol": "烯醇 O-H",
}


@dataclass
class LabileInfo:
    cls: str
    heavy_idx: int
    h_idx: int
    description: str = ""


def _is_carbonyl_carbon(atom: Chem.Atom) -> bool:
    for bond in atom.GetBonds():
        if (
            bond.GetBondType() == Chem.BondType.DOUBLE
            and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        ):
            return True
    return False


def classify_labile(mol: Chem.Mol, h_idx: int) -> LabileInfo | None:
    """判定一个氢是否为活泼氢并给出分类；非活泼氢返回 None。"""
    h = mol.GetAtomWithIdx(h_idx)
    heavy_neighbors = [n for n in h.GetNeighbors() if n.GetAtomicNum() > 1]
    if len(heavy_neighbors) != 1:
        return None
    x = heavy_neighbors[0]
    xs = x.GetSymbol()

    if xs == "O":
        others = [n for n in x.GetNeighbors() if n.GetIdx() != h_idx]
        if not others:
            return None
        host = others[0]
        if _is_carbonyl_carbon(host):
            cls = "COOH"
        elif host.GetSymbol() == "N":
            cls = "oxime"
        elif host.GetIsAromatic():
            cls = "ArOH"
        elif host.GetHybridization() == Chem.HybridizationType.SP2 and any(
            b.GetBondType() == Chem.BondType.DOUBLE
            and b.GetOtherAtom(host).GetAtomicNum() == 6
            for b in host.GetBonds()
        ):
            cls = "enol"
        else:
            cls = "ROH"
        return LabileInfo(cls, x.GetIdx(), h_idx, LABILE_CLASSES[cls])

    if xs == "N":
        if x.GetIsAromatic():
            cls = "ArNH"
        else:
            neighbors = [n for n in x.GetNeighbors() if n.GetIdx() != h_idx]
            if any(_is_carbonyl_carbon(n) for n in neighbors):
                cls = "CONH"
            elif any(n.GetSymbol() == "S" for n in neighbors):
                cls = "SO2NH"
            else:
                cls = "NH2"
        return LabileInfo(cls, x.GetIdx(), h_idx, LABILE_CLASSES[cls])

    if xs == "S":
        return LabileInfo("SH", x.GetIdx(), h_idx, LABILE_CLASSES["SH"])

    return None


def find_labile(mol: Chem.Mol) -> dict[int, LabileInfo]:
    """返回 {氢原子索引: LabileInfo}。"""
    out: dict[int, LabileInfo] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        info = classify_labile(mol, atom.GetIdx())
        if info is not None:
            out[atom.GetIdx()] = info
    return out
