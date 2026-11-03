"""分子式处理：RDKit -> 元素计数、Hill 记法、解析、精确质量、不饱和度。"""

from __future__ import annotations

import re
from collections import Counter

from .masses import get_element, monoisotopic_mass

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")

HALOGENS = frozenset({"F", "Cl", "Br", "I", "At"})


def formula_from_mol(mol) -> dict[str, int]:
    """RDKit Mol -> 元素计数（不含显式氢以外的隐式氢）。

    调用前需确保已 AddHs，否则隐式 H 不会计入。
    """
    counts: Counter[str] = Counter()
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym == "H" or sym == "*":
            counts[sym] += atom.GetTotalNumHs() if sym != "H" else 1
        else:
            counts[sym] += 1
            counts["H"] += atom.GetTotalNumHs()
    if counts.get("H", 0) == 0:
        counts.pop("H", None)
    return {k: v for k, v in counts.items() if v > 0}


def hill_notation(counts: dict[str, int]) -> str:
    """Hill 记法：有 C 时 C、H 打头，其余按字母序；无 C 时全部按字母序。"""
    items = [(k, v) for k, v in counts.items() if v > 0]
    if "C" in counts:
        head = []
        for sym in ("C", "H"):
            if sym in counts and counts[sym] > 0:
                head.append((sym, counts[sym]))
        rest = sorted(k for k, _ in items if k not in ("C", "H"))
        ordered = head + [(k, counts[k]) for k in rest]
    else:
        ordered = sorted(items, key=lambda kv: kv[0])
    return "".join(sym if n == 1 else f"{sym}{n}" for sym, n in ordered)


def parse_formula(text: str) -> dict[str, int]:
    """解析分子式字符串，例如 'C27H26F3N5S' -> {'C':27,'H':26,'F':3,'N':5,'S':1}。"""
    text = text.strip()
    if not text:
        raise ValueError("分子式为空")
    counts: Counter[str] = Counter()
    pos = 0
    for m in _FORMULA_TOKEN.finditer(text):
        if m.start() != pos:
            raise ValueError(f"分子式 {text!r} 在位置 {pos} 处无法解析")
        pos = m.end()
        sym, num = m.group(1), m.group(2)
        if sym not in _KNOWN_SYMBOLS:
            raise ValueError(f"未知元素符号 {sym!r}")
        counts[sym] += int(num) if num else 1
    if pos != len(text):
        raise ValueError(f"分子式 {text!r} 在位置 {pos} 处无法解析")
    return dict(counts)


def _known_symbols() -> set[str]:
    from .masses import known_elements

    return set(known_elements())


_KNOWN_SYMBOLS: set[str] = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
    "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
}


def exact_mass(counts: dict[str, int]) -> float:
    """中性分子的单同位素精确质量。"""
    return sum(monoisotopic_mass(sym) * n for sym, n in counts.items())


def average_mass(counts: dict[str, int]) -> float:
    """平均分子量（按天然丰度加权）。"""
    total = 0.0
    for sym, n in counts.items():
        elem = get_element(sym)
        total += sum(i.mass * i.abundance for i in elem.isotopes) * n
    return total


def rdb(counts: dict[str, int]) -> float:
    """不饱和度（环加双键数）：C + 1 - (H + X)/2 + (N + P)/2。"""
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    x = sum(counts.get(sym, 0) for sym in HALOGENS & set(counts))
    n = counts.get("N", 0) + counts.get("P", 0)
    return c + 1 - (h + x) / 2 + n / 2


def check_elements(counts: dict[str, int]) -> list[str]:
    """返回分子式中不受支持（同位素表里没有）的元素列表。"""
    from .masses import known_elements

    known = set(known_elements())
    return sorted(sym for sym in counts if sym not in known)
