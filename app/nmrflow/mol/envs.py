"""原子环境识别：取代基 SMARTS 匹配、α/β/γ 扫描、芳环取代位置判定。"""

from __future__ import annotations

from functools import lru_cache

from rdkit import Chem


@lru_cache(maxsize=512)
def _pattern(smarts: str) -> Chem.Mol | None:
    patt = Chem.MolFromSmarts(smarts)
    return patt


class SubstituentSet:
    """一组取代基定义（SMARTS 的第一个原子就是连接点）。

    匹配优先级：SMARTS 越长越具体，越优先（OCOR 优先于 OR）。
    """

    def __init__(self, table: dict[str, dict]):
        items = []
        for name, spec in table.items():
            smarts = spec.get("smarts")
            if not smarts:
                continue
            patt = _pattern(smarts)
            if patt is None:
                continue
            items.append((name, spec, smarts, patt))
        # 优先级：显式 prio > SMARTS 越长越具体 > 名称
        items.sort(key=lambda t: (-float(t[1].get("prio", 0)), -len(t[2]), t[0]))
        self.items = items
        self.by_name = {name: spec for name, spec, _, _ in items}

    def __len__(self) -> int:
        return len(self.items)

    def match_at(self, mol: Chem.Mol, atom_idx: int) -> tuple[str, dict] | None:
        """找出以 atom_idx 为连接点的取代基（按优先级取第一个）。

        ★ ``uniquify=False`` 不可省：RDKit 默认会按分子自同构折叠等价匹配，
          CF₃ 的三个氟只回一个 match，于是 ``match[0] == atom_idx`` 只在其中一个
          原子上成立，另外两个被判成"无环境"。实测这会让 ¹⁹F 的兜底率虚高到 24%。
        """
        for name, spec, _, patt in self.items:
            for match in mol.GetSubstructMatches(patt, uniquify=False, maxMatches=20000):
                if match and match[0] == atom_idx:
                    return name, spec
        return None


# 这些取代基是"链延续"，命中后仍需继续向外扫描
CHAIN_NAMES = frozenset({"alkyl"})


def scan_alpha_beta_gamma(
    mol: Chem.Mol,
    center: int,
    subs: SubstituentSet,
    max_level: int = 3,
    skip: set[int] | None = None,
) -> list[tuple[str, int, dict]]:
    """从观测碳向外扫描 α/β/γ 取代基。

    命中取代基的原子不再向外展开（避免把取代基内部的原子当成 β/γ 重复计入，
    例如 -O-CO-R 的羰基碳不应再作为 β 取代基加一次）。
    "alkyl" 例外：它是碳链延续，命中后继续展开。
    """
    skip = skip or set()
    found: list[tuple[str, int, dict]] = []
    visited = {center} | set(skip)
    frontier = [center]
    for level in range(1, max_level + 1):
        nxt: list[int] = []
        for cur in frontier:
            for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
                j = nb.GetIdx()
                if j in visited or nb.GetAtomicNum() == 1:
                    continue
                visited.add(j)
                hit = subs.match_at(mol, j)
                if hit is None:
                    nxt.append(j)
                    continue
                found.append((hit[0], level, hit[1]))
                if hit[0] in CHAIN_NAMES:
                    nxt.append(j)
        frontier = nxt
        if not frontier:
            break
    return found


def aromatic_ring_of(mol: Chem.Mol, atom_idx: int) -> list[int] | None:
    """返回包含 atom_idx 的芳香环（按连接顺序排列的原子索引）。"""
    for ring in mol.GetRingInfo().AtomRings():
        if atom_idx in ring and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            return list(ring)
    return None


def ring_distance(ring: list[int], a: int, b: int) -> int:
    """环上两个原子之间的最短环内距离。"""
    if a not in ring or b not in ring:
        return -1
    n = len(ring)
    ia, ib = ring.index(a), ring.index(b)
    d = abs(ia - ib) % n
    return min(d, n - d)


def ring_substituents(
    mol: Chem.Mol, ring: list[int], subs: SubstituentSet
) -> list[tuple[int, str, dict]]:
    """找出环上各取代基：(环原子索引, 取代基名, 条目)。

    稠环的"另一侧环原子"不是取代基，必须排除（否则喹啉的苯环部分会被
    当成吡啶环的取代基，凭空多出增量）。
    """
    ring_set = set(ring)
    all_rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    out: list[tuple[int, str, dict]] = []
    for idx in ring:
        for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
            j = nb.GetIdx()
            if j in ring_set or nb.GetAtomicNum() == 1:
                continue
            # 与 idx 同处一个环 -> 稠合原子，跳过
            if any(idx in r and j in r for r in all_rings):
                continue
            hit = subs.match_at(mol, j)
            if hit is not None:
                out.append((idx, hit[0], hit[1]))
    return out


def relative_position(dist: int, ring_size: int = 6) -> str:
    """环内距离 -> o/m/p。"""
    if dist <= 0:
        return ""
    if ring_size == 6:
        if dist == 1:
            return "o"
        if dist == 2 or dist == 4:
            return "m"
        if dist == 3:
            return "p"
    if dist == 1:
        return "o"
    if dist == 2:
        return "m"
    return "p"


def shortest_path_bonds(mol: Chem.Mol, a: int, b: int, max_len: int = 6) -> int:
    """两原子之间的最短键数（用于判断耦合路径长度）。"""
    if a == b:
        return 0
    from collections import deque

    seen = {a}
    q = deque([(a, 0)])
    while q:
        cur, d = q.popleft()
        if d >= max_len:
            break
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nb.GetIdx()
            if j in seen:
                continue
            seen.add(j)
            if j == b:
                return d + 1
            q.append((j, d + 1))
    return -1


def bond_paths(mol: Chem.Mol, a: int, b: int, max_len: int = 4) -> list[list[int]]:
    """枚举 a->b 长度 <= max_len 的所有路径（返回原子索引路径）。"""
    results: list[list[int]] = []

    def dfs(cur: int, path: list[int], visited: set[int]) -> None:
        if len(path) - 1 > max_len:
            return
        if cur == b and len(path) > 1:
            results.append(list(path))
            return
        for nb in mol.GetAtomWithIdx(cur).GetNeighbors():
            j = nb.GetIdx()
            if j in visited:
                continue
            visited.add(j)
            path.append(j)
            dfs(j, path, visited)
            path.pop()
            visited.remove(j)

    dfs(a, [a], {a})
    return results
