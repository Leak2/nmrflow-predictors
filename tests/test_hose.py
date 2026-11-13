"""HOSE 码生成器的性质测试。

这些测试守的不是"某个字符串长什么样"，而是四条索引能用的前提：
前缀性质、对称原子同码、编号无关、环境区分能力。
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from nmrflow.predict.hose import (
    hose_code,
    hose_codes,
    hose_codes_with_hydrogens,
    truncate,
)


def mol(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    assert m is not None
    return m


@pytest.mark.parametrize("smiles", ["Cc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "CCO"])
def test_truncation_matches_shorter_generation(smiles):
    """截断到 k 层 == 直接按 k 层生成，索引才能只存最深一层。"""
    m = mol(smiles)
    for idx in range(m.GetNumAtoms()):
        full = hose_code(m, idx, 6)
        for k in range(0, 7):
            assert truncate(full, k) == hose_code(m, idx, k)


def test_symmetric_atoms_share_code():
    """甲苯的邻位一对、间位一对必须同码，邻位与间位必须不同码。"""
    m = mol("Cc1ccccc1")
    codes = hose_codes(m, 4)
    assert codes[2] == codes[6]      # 邻位
    assert codes[3] == codes[5]      # 间位
    assert codes[2] != codes[3]      # 邻 ≠ 间
    assert codes[1] != codes[2]      # 连甲基的碳 ≠ 邻位


def test_renumbering_does_not_change_code():
    """原子编号顺序不影响码，否则索引无法跨分子比对。"""
    m = mol("CC(=O)Oc1ccccc1C(=O)O")
    base = hose_codes(m, 4)
    order = list(range(m.GetNumAtoms()))[::-1]
    m2 = Chem.RenumberAtoms(m, order)
    moved = hose_codes(m2, 4)
    for old_idx, old_code in base.items():
        assert moved[order.index(old_idx)] == old_code


def test_hierarchy_separates_similar_methyls():
    """乙醇与正丙醇的端基甲基第 1 层球相同，第 2 层球必须分开。"""
    eth, pro = mol("CCO"), mol("CCCO")
    eth_code, pro_code = hose_code(eth, 0, 3), hose_code(pro, 0, 3)
    assert truncate(eth_code, 1) == truncate(pro_code, 1)
    assert truncate(eth_code, 2) != truncate(pro_code, 2)


def test_ring_size_marker():
    """环上原子的中心码带 R 标记，链状原子不带。"""
    assert "R6" in hose_code(mol("c1ccccc1"), 0, 1)
    assert "R6" in hose_code(mol("C1CCCCC1"), 0, 1)
    assert "R" not in hose_code(mol("CCCCCC"), 0, 1)


def test_degree_counts_hydrogens():
    """中心码的配位数把氢算进去（CDK 的 C-4 记法）。"""
    m = mol("CCO")
    assert hose_code(m, 0, 1).startswith("CH3-4")
    assert hose_code(m, 2, 1).startswith("OH-2")


def test_hydrogen_variant_splits_groups():
    """显式氢变体：乙醇 6 个氢分成 3 组（CH3 / CH2 / OH）。"""
    eth = mol("CCO")
    codes = hose_codes_with_hydrogens(eth, 4)
    n_heavy = eth.GetNumAtoms()
    h_codes = [codes[i] for i in range(n_heavy, n_heavy + 6)]
    assert len(h_codes) == 6
    assert len(set(h_codes)) == 3


def test_diastereotopic_hydrogens_are_not_split():
    """已知局限：没有立体化学扩展，非对映异位氢会得到同一个码。

    2-丁醇 CH2 上的两个氢在化学上不等价，理想情况下应给出不同位移。
    本实现目前区分不了，所以 ¹H 走 HOSE 时必须依赖现有的拓扑等价判定兜底。
    等立体扩展（Kuhn & Johnson 2019）接进来后，这个测试应当改成 assert not equal。
    """
    m = mol("CCC(C)O")
    with_h = Chem.AddHs(m)
    codes = hose_codes_with_hydrogens(m, 4)
    target = None
    for atom in with_h.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() == 4:
            hs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetSymbol() == "H"]
            if len(hs) == 2:
                target = hs
                break
    assert target is not None
    assert codes[target[0]] == codes[target[1]]


def test_truncate_edge_cases():
    assert truncate("CH3-4;-OH;-CH2", 0) == "CH3-4"
    assert truncate("CH3-4;-OH;-CH2", -1) == ""
    assert truncate("CH3-4;-OH;-CH2", 99) == "CH3-4;-OH;-CH2"
