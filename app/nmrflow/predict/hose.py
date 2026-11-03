"""HOSE 球形环境码生成（自实现，无新增依赖）。

HOSE（Hierarchically Ordered Spherical Environment）码由 Bremser 于 1978 年提出，
是化学位移查表预测的基础数据结构：把某个原子周围的环境按化学键距离逐层展开，
拼成字符串；数据库中球形环境相同的原子，其位移可以互相平均。

为什么自己写而不用现成实现
--------------------------
`Ratsemaat/HOSE_code_generator` 是 AGPL-3.0，直接 import 会把整个便携包拖进 AGPL；
而 HOSE 算法本身公开、无专利，球形展开是一个 BFS。因此这里自行实现，
索引与查询同源，格式只需内部一致。

码格式（本实现自定义）
----------------------
    <中心描述>-<配位数>[R<环大小>];<球1>;<球2>;...

- 中心描述：元素符号（芳香原子小写）+ 氢数 + 电荷/同位素后缀，如 ``ch``、``CH3``、``N+``
- ``-<配位数>``：该原子连出的全部邻居数（重原子 + 氢），对应 CDK 的 ``C-4`` 记法
- ``R<环大小>``：该原子所在最小环的大小；不在环中则省略
- 球 i 段（i ≥ 1）：按上一层父原子分组，每组 ``(子项,子项,...)``；
  子项为 ``键符 + 该原子描述``；键符 ``-`` 单键、``=`` 双键、``#`` 三键、``:`` 芳香键
- 父原子没有子原子时写作 ``()``

**前缀性质**：球之间用 ``;`` 分隔，且每一层的展开不依赖总层数，
因此 ``hose_code(..., spheres=3)`` 是 ``spheres=5`` 的严格前缀。
索引只需存最深一层，查询时截断到第 k 个 ``;`` 即退化为 k 层球——
这是后面做"命中不足就退一层"的降级预测的前提。

已知与 CDK / 原始论文的偏离
--------------------------
1. 符号集自定义（CDK 用 ``/`` 与 ``,`` 分隔，本实现用 ``,`` 与括号分组），
   不保证与 NMRShiftDB2 现有 HOSE 码逐字兼容。本项目的索引从原始数据自建，不受影响。
2. 不做立体化学扩展（Kuhn & Johnson, ACS Omega 2019）。后果是**对映/非对映异位的原子
   会得到相同码**，¹H 预测需要额外手段区分，见 ``probe_hose.py`` 的说明。
3. 环闭合边不额外编码：展开是生成树，环信息靠多层球自然闭合体现，
   低层球（1~2 层）区分不出"在环上"与"不在环上"，靠中心码的 ``R`` 标记补一部分。
"""

from __future__ import annotations

import hashlib

from rdkit import Chem

_BOND_SYM = {
    Chem.BondType.SINGLE: "-",
    Chem.BondType.DOUBLE: "=",
    Chem.BondType.TRIPLE: "#",
    Chem.BondType.AROMATIC: ":",
    Chem.BondType.DATIVE: "-",
}

# 芳香原子用小写表示（与 SMILES 习惯一致）；这些元素在 RDKit 中才会被标为芳香
_AROMATIC_LOWER = {"C", "N", "O", "S", "P", "B", "Si", "Se", "As"}

DEFAULT_SPHERES = 6


def _hash(text: str) -> int:
    """把字符串压成 64 位整数，只用于排序与去重。"""
    return int.from_bytes(hashlib.blake2b(text.encode(), digest_size=8).digest(), "big")


def _atom_desc(atom: Chem.Atom, h_as_nodes: bool) -> str:
    sym = atom.GetSymbol()
    if atom.GetIsAromatic() and sym in _AROMATIC_LOWER:
        sym = sym.lower()
    out = sym
    if not h_as_nodes:
        n_h = atom.GetTotalNumHs()
        if n_h == 1:
            out += "H"
        elif n_h > 1:
            out += f"H{n_h}"
    chg = atom.GetFormalCharge()
    if chg:
        out += "+" * chg if chg > 0 else "-" * (-chg)
    iso = atom.GetIsotope()
    if iso:
        out += f"[{iso}]"
    return out


class _Graph:
    """一个分子的展开所需静态信息，构造一次供所有原子复用。"""

    __slots__ = ("n", "adj", "desc", "deg", "rank", "rings", "h_as_nodes")

    def __init__(self, mol: Chem.Mol, h_as_nodes: bool = False):
        self.h_as_nodes = h_as_nodes
        self.n = mol.GetNumAtoms()
        self.adj: list[list[tuple[int, str]]] = [[] for _ in range(self.n)]
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            sym = _BOND_SYM.get(bond.GetBondType(), "-")
            self.adj[i].append((j, sym))
            self.adj[j].append((i, sym))
        atoms = [mol.GetAtomWithIdx(i) for i in range(self.n)]
        self.desc = [_atom_desc(a, h_as_nodes) for a in atoms]
        # 配位数 = 图上的邻居数 + 未显式成点的隐式氢
        self.deg = [
            a.GetDegree() + (0 if h_as_nodes else a.GetNumImplicitHs()) for a in atoms
        ]
        self.rings = _smallest_rings(mol)
        self.rank = self._ranks()

    def _ranks(self) -> list[int]:
        """Weisfeiler-Lehman 迭代细化，给每个原子一个与编号无关的序。

        只用于决定兄弟节点的拼接顺序，不进入码字符串。迭代到等价类数目不再增长为止，
        此时对称原子的 rank 必然相同，从而同一子结构在任何编号下都拼出同一字符串。
        """
        rank = [_hash(d) for d in self.desc]
        prev_classes = -1
        for _ in range(max(4, self.n)):
            new = []
            for i in range(self.n):
                parts = sorted(f"{sym}{rank[j]}" for j, sym in self.adj[i])
                new.append(_hash(f"{rank[i]}|{','.join(parts)}"))
            classes = len(set(new))
            rank = new
            if classes == prev_classes:
                break
            prev_classes = classes
        return rank


def _smallest_rings(mol: Chem.Mol) -> list[int]:
    info = mol.GetRingInfo()
    out = [0] * mol.GetNumAtoms()
    if info.NumRings() == 0:
        return out
    for ring in info.AtomRings():
        size = len(ring)
        for idx in ring:
            if out[idx] == 0 or size < out[idx]:
                out[idx] = size
    return out


def _levels(g: _Graph, root: int, spheres: int) -> list[str]:
    """从 root 出发逐层展开，返回球 1..k 的文本段（不含中心）。"""
    seen = {root}
    frontier = [root]
    levels: list[str] = []
    for _ in range(spheres):
        groups: list[list[tuple[str, int, str, int]]] = []
        nxt: list[int] = []
        for parent in frontier:
            kids = []
            for j, sym in g.adj[parent]:
                if j in seen:
                    continue
                seen.add(j)
                kids.append((sym, g.rank[j], g.desc[j], j))
            kids.sort(key=lambda t: (t[0], t[1], t[2]))
            groups.append(kids)
            nxt.extend(t[3] for t in kids)
        if not nxt:
            break
        levels.append("".join("(" + ",".join(s + d for s, _, d, _ in grp) + ")" for grp in groups))
        frontier = nxt
    return levels


def _center(g: _Graph, idx: int, with_ring: bool) -> str:
    text = f"{g.desc[idx]}-{g.deg[idx]}"
    if with_ring and g.rings[idx]:
        text += f"R{g.rings[idx]}"
    return text


def hose_code(
    mol: Chem.Mol,
    atom_idx: int,
    spheres: int = DEFAULT_SPHERES,
    *,
    ring_size: bool = True,
    _graph: _Graph | None = None,
) -> str:
    """单个原子的 HOSE 码。``spheres`` 是展开的球层数。"""
    g = _graph or _Graph(mol)
    parts = [_center(g, atom_idx, ring_size), *_levels(g, atom_idx, spheres)]
    return ";".join(parts)


def hose_codes(
    mol: Chem.Mol,
    spheres: int = DEFAULT_SPHERES,
    *,
    ring_size: bool = True,
    atoms: list[int] | None = None,
) -> dict[int, str]:
    """整分子的 HOSE 码，返回 ``{原子索引: 码}``。默认只算重原子（氢折进描述）。"""
    g = _Graph(mol)
    targets = range(g.n) if atoms is None else atoms
    return {
        i: ";".join([_center(g, i, ring_size), *_levels(g, i, spheres)]) for i in targets
    }
def hose_codes_with_hydrogens(
    mol: Chem.Mol, spheres: int = DEFAULT_SPHERES, *, ring_size: bool = True
) -> dict[int, str]:
    """¹H 用的变体：把氢显式加进图里，返回所有原子（含氢）的码。

    **生产路径不调用它**：¹H 的 HOSE 查表走 ``hose_codes`` + 拓扑等价分组。
    这个变体是给 ``tests/test_hose.py`` 用的——用来钉住"氢作为节点参与展开"
    这条性质，以及记录"非对映异位氢仍会同码"这个已知局限。删它就得连那两个
    测试一起删，等于把一条已声明的不确定性抹掉，所以留着。

    ``Chem.AddHs`` 保留原有原子顺序并把氢追加在末尾，所以重原子的索引与输入一致，
    氢的索引从 ``mol.GetNumAtoms()`` 开始。氢作为节点参与展开，
    因此两个化学不等价的氢（例如手性中心上的两个 H）会得到不同码——
    但**对映/非对映异位**仍需立体化学扩展才能完全区分，本实现没有做。
    """
    with_h = Chem.AddHs(mol)
    g = _Graph(with_h, h_as_nodes=True)
    return {
        i: ";".join([_center(g, i, ring_size), *_levels(g, i, spheres)])
        for i in range(g.n)
    }

def truncate(code: str, spheres: int) -> str:
    """把码截断到 k 层球。第 0 段是中心，所以 ``spheres=k`` 保留 k+1 段。"""
    parts = code.split(";")
    if spheres < 0:
        return ""
    return ";".join(parts[: spheres + 1])
