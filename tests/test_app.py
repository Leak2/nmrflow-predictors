"""Streamlit 界面的冒烟测试：脚本能跑通、点「生成谱图」不抛异常。"""

from __future__ import annotations

from pathlib import Path

import pytest

st_testing = pytest.importorskip("streamlit.testing.v1")

APP = Path(__file__).resolve().parents[1] / "app" / "nmrflow" / "app.py"


def test_app_renders_and_runs():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()
    assert not at.exception

    # 默认不预填 SMILES（下拉默认「（自定义）」），所以显式给一个分子。
    # 早先这条用例靠“下拉默认选中某个示例、输入框里就有值”跑通，
    # 那个隐含依赖本身就是这次要改掉的东西。
    at.text_input[0].set_value("CCO")
    at.button[0].click().run()
    assert not at.exception
    assert at.code, "报道值代码块没渲染出来"
    assert "1H NMR" in at.code[0].value
    assert "HRMS" in at.code[0].value


def test_example_dropdown_defaults_to_custom():
    """下拉默认「（自定义）」，且输入框不预填示例。

    预填的话，用户会以为框里那个分子是自己输入的；而默认选中某个示例时，
    下拉显示的和框里写的又对不上（早先默认是「乙酸乙酯」）。
    """
    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()
    assert at.selectbox[0].value == "（自定义）"
    assert at.text_input[0].value == ""


def test_picking_example_fills_input_but_custom_keeps_edits():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()

    at.selectbox[0].set_value("咖啡因").run()
    assert at.text_input[0].value == "Cn1cnc2c1c(=O)n(C)c(=O)n2C"

    at.text_input[0].set_value("C1CCCCC1").run()
    # 改别的控件不能把输入打回示例值
    at.selectbox[[s.label for s in at.selectbox].index("氘代溶剂")].set_value("DMSO-d6").run()
    assert at.text_input[0].value == "C1CCCCC1"

    # 切回「（自定义）」也不能清掉手打的内容——这是 value= 写法最容易犯的错
    at.selectbox[0].set_value("（自定义）").run()
    assert at.text_input[0].value == "C1CCCCC1"

    # 再选回同一个示例，应当重新填充
    at.selectbox[0].set_value("咖啡因").run()
    assert at.text_input[0].value == "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


def test_empty_smiles_gives_hint_not_traceback():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()
    at.text_input[0].set_value("")
    _click(at, "生成谱图")
    assert not at.exception
    assert any("SMILES" in e.value for e in at.error)



def test_app_handles_bad_smiles():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()
    # 把 SMILES 输入框改成非法串
    at.text_input[0].set_value("这不是SMILES").run()
    at.button[0].click().run()
    assert any("SMILES" in e.value for e in at.error) or not at.exception


def _click(at, label: str):
    """按标签点按钮。侧栏与主区的元素顺序不稳定，不能靠下标。"""
    labels = [b.label for b in at.button]
    at.button[labels.index(label)].click().run()


def _generate(at, smiles: str):
    at.text_input[0].set_value(smiles)
    _click(at, "生成谱图")


def test_results_survive_rerun_so_each_download_works():
    """点下载后页面不能退回空状态。

    ``st.download_button`` 被点击同样会触发一次 rerun，而 ``st.button`` 在 rerun 时
    返回 False。若结果只活在「点了生成」那一支里，用户就会看到「下载一张图之后
    整页回到首页，第二张图下不了」。结果必须存在 ``session_state`` 里。
    """
    at = st_testing.AppTest.from_file(str(APP), default_timeout=240)
    at.run()
    _generate(at, "CCO")

    downloads = len(at.get("download_button"))
    assert downloads >= 3, "报道值 + 各核谱图应该各有一个下载按钮"
    assert at.code, "报道值没渲染"

    # 连点三次下载 = 三次不按钮的 rerun
    for _ in range(3):
        at.run()
        assert not at.exception
        assert at.code, "rerun 后报道值丢了（结果没有持久化）"
        assert len(at.get("download_button")) == downloads, "rerun 后下载按钮少了"

    # 下载内容本身也要稳定：同一结果两次 rerun 应是同一组按钮
    labels_before = [d.label for d in at.get("download_button")]
    at.run()
    assert [d.label for d in at.get("download_button")] == labels_before


def test_results_are_marked_with_the_molecule_they_belong_to():
    """改 SMILES 但不重新生成时，结果区要说明它属于哪个分子。

    持久化的副作用是「输入框改了、下面还是旧图」，不标注就成了误导。
    """
    at = st_testing.AppTest.from_file(str(APP), default_timeout=240)
    at.run()
    _generate(at, "CCO")
    at.text_input[0].set_value("c1ccccc1").run()   # 只改输入，不点生成

    captions = [c.value for c in at.caption if c.value.startswith("当前显示")]
    assert captions and "CCO" in captions[0]
    assert "c1ccccc1" not in captions[0]


def test_new_generation_replaces_previous_result():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=240)
    at.run()
    _generate(at, "CCO")
    _generate(at, "c1ccccc1")

    captions = [c.value for c in at.caption if c.value.startswith("当前显示")]
    assert "c1ccccc1" in captions[0]
    # 替换而不是叠加：下载按钮数量不变
    assert len(at.get("download_button")) == 5


def test_clear_button_returns_to_empty_state():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=240)
    at.run()
    _generate(at, "CCO")
    _click(at, "清空结果")

    assert not at.get("download_button")
    assert at.info, "清空后应回到「输入 SMILES」的提示"
    assert "清空结果" not in [b.label for b in at.button]


def test_batch_results_survive_rerun():
    at = st_testing.AppTest.from_file(str(APP), default_timeout=300)
    at.run()
    labels = [c.label for c in at.checkbox]
    at.checkbox[next(i for i, l in enumerate(labels) if "批量" in l)].set_value(True).run()
    at.text_area[0].set_value("CCO\nC1CCCCC1")
    _click(at, "批量生成")

    assert at.get("dataframe")
    at.run()
    assert at.get("dataframe"), "批量结果 rerun 后丢了"
    assert at.get("download_button"), "批量 ZIP 下载按钮 rerun 后丢了"


def test_zip_payload_contains_reports_and_images(tmp_path):
    from nmrflow.app import _zip_dir

    (tmp_path / "report.txt").write_text("calculated only", encoding="utf-8")
    (tmp_path / "spectrum.png").write_bytes(b"png")
    data = _zip_dir(tmp_path)
    import io
    import zipfile

    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        assert sorted(zf.namelist()) == ["report.txt", "spectrum.png"]


def test_simulated_found_is_not_exposed_in_ui():
    """「模拟 Found」只留命令行入口，界面上不能有这个开关。

    它是内测期的排版演示功能（Found = Calculated + 随机抖动）。正式版把 UI 撤掉，
    但 ``RunOptions(ms_found_mode="simulate")`` 与 ``--found`` 仍然可用——
    所以这里既要断言界面没有开关，也要断言默认输出里不出现 Found。
    """
    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()

    labels = [c.label for c in at.checkbox]
    assert not any("Found" in label for label in labels), f"界面又出现了 Found 开关：{labels}"

    _generate(at, "CCO")
    assert at.code
    assert "Found" not in at.code[0].value


def test_version_badge_shows_current_version():
    """右上角角标必须与 ``nmrflow.__version__`` 一致。

    这个角标存在的意义就是"用户报问题时能说出自己用的哪一版"，
    所以它一旦和真实版本脱钩，比没有更糟——所以用 __version__ 断言，而不是硬编码。
    """
    import nmrflow

    at = st_testing.AppTest.from_file(str(APP), default_timeout=180)
    at.run()

    badges = [m.value for m in at.markdown if "nmrflow-version" in m.value]
    assert badges, "界面没有渲染版本角标"
    assert f"v{nmrflow.__version__}" in badges[0]
    assert not at.exception
