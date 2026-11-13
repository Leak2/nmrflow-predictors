"""``tools/check_reference_records.py`` 的护栏。

这个工具的作用是"别再用错记录改模型"——0.2.3/0.2.6/0.2.8 三次改动都是对同一条
十氢萘记录的不同解读。所以这里既要钉住"对的记录必须 PASS"，也要钉住
"错的记录必须 FAIL"，否则它自己就成了下一个错误依据。

合成 SD 直接由 RDKit 写，不依赖真实谱库。
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parent.parent

DECALIN = "C1CCC2CCCCC2C1"
NORBORNANE = "C1CC2CCC1C2"

#: 文献值（见 tools/check_reference_records.py 的 docstring）
DECALIN_TRANS = [(44.0, 2), (34.7, 4), (27.1, 4)]
DECALIN_CIS = [(36.0, 2), (26.5, 4), (25.0, 4)]
NORBORNANE_LIT = [(37.3, 2), (30.6, 4), (39.2, 1)]
#: 0.2.8 用的"文献值"（C7 24.5）与 0.2.3 的十氢萘"34.5/26.9/26.3"
NORBORNANE_0_2_8 = [(36.8, 2), (28.5, 4), (24.5, 1)]
DECALIN_0_2_3 = [(34.5, 2), (26.9, 4), (26.3, 4)]


def _load():
    spec = importlib.util.spec_from_file_location(
        "refcheck_tool", ROOT / "app" / "tools" / "check_reference_records.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def refcheck():
    return _load()


def _write_sd(path: Path, entries: list[tuple[str, list[tuple[float, int]]]]) -> Path:
    """按 (SMILES, [(δ, 碳数)...]) 写一个带 ``Spectrum 13C`` 属性的 SD。"""
    writer = Chem.SDWriter(str(path))
    for smiles, pattern in entries:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        flat: list[float] = []
        for value, count in pattern:
            flat.extend([value] * count)
        carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
        assert len(flat) == len(carbons), (smiles, flat, len(carbons))
        body = "|".join(f"{v:.2f};1;{i}" for v, i in zip(flat, carbons))
        mol.SetProp("Spectrum 13C 1", body)
        writer.write(mol)
    writer.close()
    return path


def test_clusters_group_equivalent_carbons(refcheck):
    # 升序返回；同一簇内取均值
    assert refcheck._clusters({0: 44.0, 1: 44.05, 2: 34.7}) == [(34.7, 1), (44.02, 2)]


def test_decalin_trans_record_passes(refcheck, tmp_path):
    sd = _write_sd(tmp_path / "trans.sdf", [(DECALIN, DECALIN_TRANS)])
    t = refcheck.check(sd)["targets"]["decalin"]
    assert t["isomers_found"] == ["trans"], t
    assert t["ok"] is True


def test_decalin_cis_record_is_recognised_as_cis(refcheck, tmp_path):
    """cis/trans 共用骨架 InChIKey 前 14 位；工具必须靠 36.0 的桥头把它们分开。"""
    sd = _write_sd(tmp_path / "cis.sdf", [(DECALIN, DECALIN_CIS)])
    report = refcheck.check(sd)
    assert report["targets"]["decalin"]["isomers_found"] == ["cis"]
    # 只有 cis、没有 trans：trans 是必需项，所以整体不通过
    assert report["ok"] is False


def test_decalin_old_patch_values_are_not_trans(refcheck, tmp_path):
    """0.2.3 那组 34.5/26.9/26.3 落到 cis 区间，不能被当成 trans。"""
    sd = _write_sd(tmp_path / "old.sdf", [(DECALIN, DECALIN_0_2_3)])
    report = refcheck.check(sd)
    assert report["targets"]["decalin"]["isomers_found"] == ["cis"]


def test_norbornane_literature_passes(refcheck, tmp_path):
    sd = _write_sd(
        tmp_path / "both.sdf",
        [(DECALIN, DECALIN_TRANS), (NORBORNANE, NORBORNANE_LIT)],
    )
    report = refcheck.check(sd)
    assert report["targets"]["norbornane"]["isomers_found"] == ["parent"]
    assert report["ok"] is True


def test_norbornane_wrong_c7_is_rejected(refcheck, tmp_path):
    """把 C7 写成 24.5（0.2.8 的依据）必须被判为对不上。"""
    sd = _write_sd(tmp_path / "norb_bad.sdf", [(NORBORNANE, NORBORNANE_0_2_8)])
    report = refcheck.check(sd)
    t = report["targets"]["norbornane"]
    assert t["isomers_found"] == ["unknown"], t
    assert t["ok"] is False


def test_missing_molecule_fails(refcheck, tmp_path):
    sd = _write_sd(tmp_path / "ethanol.sdf", [("CCO", [(58.0, 1), (18.0, 1)])])
    report = refcheck.check(sd)
    assert report["ok"] is False
    for target in report["targets"].values():
        assert target["records"] == 0
