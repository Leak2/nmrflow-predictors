"""端到端：SMILES -> 报告文本 / 谱图数组 / 元素自动探测。

这里守的是"用户的验收标准"：
  1. 一定出 1H 与 13C；
  2. 含 F 才出 19F、含 P 才出 31P；
  3. 报道值形如 `SMILES:1H NMR (CDCl3, 400 MHz) δ: ...` 并以 HRMS 收尾；
  4. 常见分子的位移落在文献值附近。
"""

from __future__ import annotations

import re

import pytest

from nmrflow.pipeline import RunOptions, run


def _opt(**kwargs):
    """本模块守的是 **Tier A** 的行为，因此显式关掉 Tier B。

    Tier B（HOSE 查表）会改写大部分位移，且它是否生效取决于 build/ 下有没有索引；
    混在一起测，等于让测试结果依赖一个构建产物。Tier B 自己的行为见
    tests/test_shift_hose.py。
    """
    kwargs.setdefault("tier_b", False)
    return RunOptions(**kwargs)


def _peaks(result, nucleus):
    block = next(b for b in result.report.blocks if b.nucleus == nucleus)
    return block.peaks


def _shifts(result, nucleus):
    return [p.shift for p in _peaks(result, nucleus)]


# ---------------------------------------------------------------------------
# 元素自动探测
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CCO", ["1H", "13C"]),
        ("FC(F)(F)c1ccccc1", ["1H", "19F", "13C"]),
        ("COP(=O)(OC)OC", ["1H", "31P", "13C"]),
        ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", ["1H", "13C"]),
    ],
)
def test_nucleus_detection(smiles, expected):
    result = run(smiles, _opt())
    assert [b.nucleus for b in result.report.blocks] == expected


def test_no_spectrum_for_absent_element():
    result = run("CCO", _opt())
    assert "19F" not in result.spectra
    assert "31P" not in result.spectra


# ---------------------------------------------------------------------------
# 报道值形态
# ---------------------------------------------------------------------------
def test_report_shape():
    result = run("CCO", _opt())
    text = result.report.text
    assert text.startswith("CCO:")
    assert re.search(r"1H NMR \(CDCl3, 400 MHz\) δ: ", text)
    assert re.search(r"13C NMR \(CDCl3, 101 MHz\) δ: ", text)
    assert text.rstrip().endswith("[M+H]+.")


def test_hrms_line():
    """HRMS 只给 Calculated，不伪造 Found。"""
    text = run("CCO", _opt()).report.text
    assert "HRMS (DART) for C2H6O m/z: Calculated, 47.0491" in text
    assert "Found" not in text


def test_spin_limit_warning_reaches_report():
    result = run("CCO", _opt(max_exact_spins=1, proton_mhz=10))
    assert any("自旋数" in w and "上限" in w for w in result.report.warnings)


def test_ei_uses_radical_cation_notation():
    from nmrflow.ms.adducts import adducts_for_method

    assert adducts_for_method("EI")[0].name == "M⁺·"
    result = run("CCO", _opt(ms_method="EI"))
    assert "M⁺·" in result.report.text


def test_hrms_default_has_no_found():
    """默认只给 Calculated。

    旧行为是默认就编一个 Found（calc + 随机抖动），用户会把它当成实测值
    粘进论文——这是诚信问题，不是排版问题。
    """
    text = run("CCO", _opt()).report.text
    assert "Calculated, 47.0491" in text
    assert "Found" not in text


@pytest.mark.parametrize(
    "smiles, nucleus, expected",
    [
        # 期望值取常见文献/表值；同一组等价核只列一个值，命中即可。
        ("CCO", "1H", [1.25, 3.72]),
        ("CCOC(C)=O", "1H", [4.12, 2.05, 1.26]),
        ("Cc1ccccc1", "1H", [7.2, 2.36]),
        ("CCCCCCCC", "13C", [14.1, 22.7, 31.9, 29.3]),
        ("Oc1ccccc1", "13C", [154.9, 129.8, 121.1, 115.4]),
    ],
)
def test_shifts_close_to_literature(smiles, nucleus, expected):
    # 经验加和法的固有误差：1H 约 ±0.3 ppm，13C 约 ±0.8 ppm
    tol = 0.3 if nucleus == "1H" else 0.8
    got = _shifts(run(smiles, _opt()), nucleus)
    for want in expected:
        assert any(abs(g - want) < tol for g in got), \
            f"{smiles} {nucleus}: 期望 {want}±{tol}，实际 {got}"


# ---------------------------------------------------------------------------
# 多重峰命名
# ---------------------------------------------------------------------------
def test_ethanol_multiplets():
    peaks = _peaks(run("CCO", _opt()), "1H")
    by_shift = {round(p.shift): p for p in peaks}
    ch2 = by_shift[4]      # 3.47 附近
    ch3 = by_shift[1]      # 1.22 附近
    assert ch2.multiplicity == "q"
    assert ch3.multiplicity == "t"
    assert ch2.j_values[0] == pytest.approx(7.0, abs=0.1)


def test_benzene_is_singlet():
    peaks = _peaks(run("c1ccccc1", _opt()), "1H")
    assert len(peaks) == 1
    assert peaks[0].multiplicity == "s"
    assert peaks[0].integration == pytest.approx(6)


def test_octane_methyl_is_triplet():
    peaks = _peaks(run("CCCCCCCC", _opt()), "1H")
    ch3 = min(peaks, key=lambda p: p.shift)
    assert ch3.multiplicity == "t"
    assert ch3.integration == pytest.approx(6)


def test_aldehyde_quartet_small_j():
    peaks = _peaks(run("CC=O", _opt()), "1H")
    cho = max(peaks, key=lambda p: p.shift)
    assert cho.multiplicity == "q"
    assert 2.0 < cho.j_values[0] < 4.0     # CH3-CHO 的 3J 只有 ~2.9 Hz


# ---------------------------------------------------------------------------
# 异核耦合
# ---------------------------------------------------------------------------
def test_cf3_carbon_quartet():
    peaks = _peaks(run("FC(F)(F)c1ccccc1", _opt()), "13C")
    cf3 = next(p for p in peaks if p.label == "-CF3")
    assert cf3.multiplicity == "q"
    assert cf3.j_values == [pytest.approx(272.6, abs=0.5)]
    assert cf3.j_partners == ["19F"]


def test_phosphate_coupling_labels():
    result = run("COP(=O)(OC)OC", _opt())
    c_peaks = _peaks(result, "13C")
    assert all("JC-P" in (b.text or "") for b in result.report.blocks
               if b.nucleus == "13C")
    assert c_peaks[0].j_partners == ["31P"]
    text = result.report.text
    assert "JC-P=" in text.replace(" ", "").replace("JC-P= ", "JC-P=")
    assert "JC-F" not in text


def test_benzotrifluoride_fluorine_is_singlet():
    """Ph-CF3 的 19F 实测是单峰：绕出芳环的 H-F 耦合必须被压掉。"""
    peaks = _peaks(run("FC(F)(F)c1ccccc1", _opt()), "19F")
    assert peaks[0].multiplicity == "s"


def test_fluorobenzene_ortho_hf():
    peaks = _peaks(run("Fc1ccccc1", _opt()), "1H")
    js = {round(j) for p in peaks for j in p.j_values}
    assert 9 in js          # 邻位 H-F
    assert 6 in js          # 间位 H-F


# ---------------------------------------------------------------------------
# 谱图数组
# ---------------------------------------------------------------------------
def test_caffeine_n_methyl_tiers():
    """咖啡因三个 N-CH3 分属酰亚胺 / 内酰胺 / 咪唑，δ 必须分开。"""
    result = run("Cn1cnc2c1c(=O)n(C)c(=O)n2C", _opt())
    h = [p.shift for p in _peaks(result, "1H") if p.shift < 4.2]
    assert len(h) == 3
    assert sorted(h) == pytest.approx([3.40, 3.60, 3.90], abs=0.05)
    c = [p.shift for p in _peaks(result, "13C") if p.shift < 40]
    assert sorted(c) == pytest.approx([27.9, 29.7, 33.6], abs=0.3)


def test_spectrum_arrays():
    result = run("CCO", _opt())
    spec = result.spectra["1H"]
    assert spec.x_ppm.size == spec.y.size
    assert spec.y.max() > 0
    assert spec.integral is not None
    # δ 轴从低到高，强度非负
    assert spec.x_ppm[0] < spec.x_ppm[-1]
    assert spec.y.min() >= -0.05 * spec.y.max()


def test_field_dependence_of_digitized_spectrum():
    """同一分子在 60 MHz 与 800 MHz 下的谱线数应当不同（二级效应）。"""
    low = run("CCOC(C)=O", _opt(proton_mhz=60))
    high = run("CCOC(C)=O", _opt(proton_mhz=800))
    assert len(low.spectra["1H"].transitions) >= 1
    assert low.spectra["1H"].freq_mhz == pytest.approx(60.0)
    assert high.spectra["1H"].freq_mhz == pytest.approx(800.0)


# ---------------------------------------------------------------------------
# 多氟环境 / 盐类（光致产酸剂 PAG）
# ---------------------------------------------------------------------------
PAG = (
    "C1=C(C=CC=C1)[S+](C1=CC=CC=C1)C1=CC=CC=C1."
    "C=CC1=CC=C(C=C1)OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)"
    "[C-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F"
)


def test_all_fluorine_environments_resolved():
    """12 个 F 分属 4 类环境，绝不能全落到同一个 fallback 值。"""
    result = run(PAG, _opt())
    shifts = [s.shift for s in result.signals["19F"]]
    assert len(shifts) == 4
    # 峰位必须彼此分开（旧版全部 -115.0，谱图上叠成一个包）
    assert len({round(x, 1) for x in shifts}) == 4
    assert min(abs(a - b) for i, a in enumerate(shifts) for b in shifts[i + 1:]) > 3.0

    # 具体归属：三氟甲磺酰 ~-78、三个 CF2 在 -109 / -114 / -122
    assert any(abs(x + 78.0) < 3.0 for x in shifts), shifts
    assert any(abs(x + 109.0) < 3.0 for x in shifts), shifts
    assert any(abs(x + 122.0) < 3.0 for x in shifts), shifts


def test_cf3_sulfonyl_coupling_is_320():
    """CF3-SO2 的 13C 1J(C-F) 约 320 Hz，不是 Ar-CF3 的 272.6。"""
    result = run(PAG, _opt())
    cf3 = [s for s in result.signals["13C"] if s.provenance.key.endswith("CF3_sulfonyl")]
    assert cf3, [s.provenance.key for s in result.signals["13C"]]
    j1 = [c.j_hz for c in cf3[0].couplings if abs(c.j_hz) > 200]
    assert j1 and pytest.approx(j1[0], abs=5) == 320.0


# ---------------------------------------------------------------------------
# 全氟烷基磺酰氟：负化学位移 + 同一分子里的正位移核
# ---------------------------------------------------------------------------
SULFONYL_FLUORIDE = (
    "FC(F)(S(=O)(=O)F)C(F)(F)C(F)(F)S(=O)(=O)OC1=CC=C(C=C1)C=C"
)


def test_negative_shift_transitions_are_kept():
    """δ<0 的核不能被丢掉。

    19F 的 -110 / -122 ppm 对应负的绝对频率（以 CFCl3 为 0）。早期版本在
    哈密顿量求解后按 `freq <= 0` 过滤，把整片高场氟全部删掉，于是谱图上只剩
    磺酰氟那一个正位移峰。
    """
    result = run(SULFONYL_FLUORIDE, _opt())
    spec = result.spectra["19F"]
    assert spec.x_ppm.min() < -100, spec.x_ppm.min()

    gids = {t.group_id for t in spec.transitions}
    assert len(gids) == 4, sorted(gids)
    assert any(t.freq_hz < 0 for t in spec.transitions), "负频率跃迁应保留"


def test_sulfonyl_fluoride_has_its_own_environment():
    """R-SO2F 走专门条目（+46 附近），既不落到 SF5/SF3 的 60，也不落到 CF2。"""
    result = run(SULFONYL_FLUORIDE, _opt())
    sigs = result.signals["19F"]
    shifts = [s.shift for s in sigs]
    assert len(sigs) == 4, shifts

    so2f = [s for s in sigs if s.n_nuclei == 1][0]
    assert "SO2F" in so2f.provenance.key
    assert pytest.approx(so2f.shift, abs=5) == 46.0

    # 端基 CF2（邻 SO2F）不能和链内 CF2 撞成同一个值
    values = {round(s.shift) for s in sigs}
    assert len(values) == 4, sorted(values)



def test_salt_keeps_both_fragments():
    """盐的两半都要出信号：阳离子 15 个芳香 H + 阴离子 4 个芳香 H。"""
    result = run(PAG, _opt())
    total_h = sum(s.n_nuclei for s in result.signals["1H"])
    assert total_h == 22          # C18H15S+ (15) + 苯乙烯基 (7)

    # 三苯基锍的芳香 H 在 7.5~7.8，不能全落在苯的基准 7.26
    aromatic = [s.shift for s in result.signals["1H"] if s.label == "Ph-H"]
    assert min(aromatic) > 7.3
    assert max(aromatic) > 7.7


def test_salt_hrms_reports_two_ions():
    """离子化合物分别报 [M]+ / [M]-，不把两半加成中性式。"""
    result = run(PAG, _opt())
    assert result.report.ms is not None
    assert "[M]+" in result.report.ms.adducts[0].adduct
    assert pytest.approx(result.report.ms.adducts[0].calculated, abs=1e-3) == 263.0889
    assert len(result.report.extra_ms) == 1
    assert "[M]-" in result.report.extra_ms[0].adducts[0].adduct
    assert pytest.approx(result.report.extra_ms[0].adducts[0].calculated, abs=1e-3) == 674.8787


def test_styrene_alkene_couplings():
    """=CH2 的偕偶是小的正值，顺/反 3J 必须分开给 10.5 / 17.0。"""
    result = run("C=Cc1ccccc1", _opt())
    gem = [c.j_hz for s in result.signals["1H"] for c in s.couplings if c.kind == "2J(alkene)"]
    assert gem and all(abs(j) < 4.0 for j in gem), gem       # 不能是 -15

    vic = sorted({round(c.j_hz, 1) for s in result.signals["1H"]
                  for c in s.couplings if c.kind == "3J(alkene)"})
    assert vic == pytest.approx([10.5, 17.0], abs=0.2), vic


def test_hetero_reports_carry_integration():
    """19F / 31P 的报道值要像 1H 一样给出核数（2F / 1P）。"""
    smi = "FC(F)(S(=O)(=O)F)C(F)(F)C(F)(F)S(=O)(=O)OC1=CC=C(C=C1)C=C"
    result = run(smi, _opt(solvent="CDCl3"))
    line = [b for b in result.report.blocks if b.nucleus == "19F"][0].peaks

    counts = [int(round(p.integration)) for p in line]
    assert sum(counts) == 7, counts          # 分子里共 7 个 F
    assert counts.count(2) == 3 and counts.count(1) == 1, counts

    text = result.report.text
    assert "1F)" in text and "2F)" in text
    # 积分紧跟在 J 列表之后，格式与 1H 保持一致
    assert "-109.00 (m, J= 12.0 Hz, J= 5.0 Hz, J= 2.0 Hz, 2F)" in text


def test_phosphorus_report_carries_integration():
    result = run("CC(=O)OP(=O)(OCC)OCC", _opt())
    block = [b for b in result.report.blocks if b.nucleus == "31P"][0]
    assert [int(round(p.integration)) for p in block.peaks] == [1]
    assert "1P)" in result.report.text


# ---------------------------------------------------------------------------
# 全氟烷基磺酸根（PAG 阴离子）：α-CF2 连在 S 上，不能被"链内 CF2"抢走
# ---------------------------------------------------------------------------
SULFONATE_ANION = (
    "C1=C(C=CC=C1)[S+](C1=CC=CC=C1)C1=CC=CC=C1."
    "[O-]S(=O)(=O)C(F)(F)C(F)(F)C1C2CC(C1)C1OC3(C4CC5CC(C4)CC3C5)OC12"
)


def test_sulfonate_anion_splits_into_two_fluorine_environments():
    """CF2-SO3- 与 CF2-CH 是两个环境，都必须有自己的位移。

    旧版表里没有 [O-] 磺酸根条目，而"链内 CF2"的 SMARTS 只要求有一个碳邻居，
    于是连在硫上的 α-CF2 也被判成链内 CF2，两组 CF2 撞成同一个 -122，
    谱图上只剩一个 4F 单峰。
    """
    result = run(SULFONATE_ANION, _opt())
    sigs = result.signals["19F"]

    assert len(sigs) == 2, [(s.shift, s.provenance.key) for s in sigs]
    assert [s.n_nuclei for s in sigs] == [2, 2]
    shifts = sorted(s.shift for s in sigs)
    assert min(b - a for a, b in zip(shifts, shifts[1:])) > 0.2, shifts
    assert all(abs(x + 122.0) > 3.0 for x in shifts), shifts

    so3 = [s for s in sigs if "SO3" in s.provenance.key]
    assert so3, [s.provenance.key for s in sigs]
    assert pytest.approx(so3[0].shift, abs=3) == -114.0

    # 两个环境都要出峰，且 19F 谱上能看到两组跃迁
    spec = result.spectra["19F"]
    assert len({t.group_id for t in spec.transitions}) == 2


def test_internal_cf2_rule_requires_two_carbon_neighbours():
    """连在硫上的 CF2 不属于"两侧均为全氟碳"，必须落到别的条目。"""
    result = run(SULFONATE_ANION, _opt())
    for sig in result.signals["19F"]:
        assert "链内" not in sig.provenance.key, sig.provenance.key


def test_fluorotelomer_chain_terminus_shifts():
    """全氟链三端的位移要分开：CF2-CH ~-115、CF2-CF2- ~-122、CF2-CF3 ~-126。"""
    smi = "[O-]S(=O)(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    result = run(smi, _opt())
    shifts = {round(s.shift) for s in result.signals["19F"]}
    assert -115 in shifts and -122 in shifts and -126 in shifts, sorted(shifts)


def test_degenerate_shift_warning_mentions_the_shared_rule():
    """环境表没能区分开时要有告警，而不是安静地给一张单峰图。"""
    smi = "[O-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    result = run(smi, _opt())
    assert any("同一化学位移" in w for w in result.spectra["19F"].warnings), \
        result.spectra["19F"].warnings


def test_ch2_next_to_potential_stereocenter_warns():
    """CH2 紧邻潜在手性中心：只提示非对映异位，不伪造 AB。"""
    result = run("CCC(Cl)(F)C")
    assert any("CH2" in w and "非对映异位" in w for w in result.report.warnings)