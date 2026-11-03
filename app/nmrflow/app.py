"""Streamlit 单页界面：输入 SMILES -> 即时出谱图与报道值。

启动：
    streamlit run nmrflow/app.py

页面结构：
    左栏  输入与参数
    右栏  上：报道值（可复制）+ HRMS
          下：各核谱图（1H 带积分曲线）

结果同时写到 ``<包根>/results/<slug>/``：PNG、report.txt。
支持单分子、批量（每行一个 SMILES）、结构式预览、整包 ZIP 下载。
"""

from __future__ import annotations

import io
import os
import sys
import tempfile
import time
import uuid
import zipfile
from pathlib import Path

# `streamlit run nmrflow/app.py` 会把脚本所在目录（nmrflow/）放进 sys.path，
# 用相对 import 会炸；这里手动把项目根加进去，两种启动方式都能用。
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# numpy 2.5 把 f2py 做成懒加载属性。scipy 内部的 `from numpy import *`
# 会在 Streamlit 的脚本线程里首次触发它，偶发 re 解析失败
# （"missing ), unterminated subpattern"）。先导入一次就不会再走这条路径。
import numpy.f2py  # noqa: E402,F401

import streamlit as st  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Draw  # noqa: E402

from nmrflow.core.types import PREDICTED_DATA_NOTE_TEMPLATE  # noqa: E402
from nmrflow import __version__, autostop  # noqa: E402
from nmrflow.pipeline import RunOptions, run, write_manifest  # noqa: E402
from nmrflow.plot.plot_1h import plot_1h  # noqa: E402
from nmrflow.plot.plot_generic import plot_generic  # noqa: E402
from nmrflow.plot.plot_ms import plot_ms  # noqa: E402
from nmrflow.predict.hose_index import resolve_index_path  # noqa: E402

# 便携包裁剪掉索引时，界面上把 Tier B 开关置灰，而不是让用户点了没反应
_HAS_INDEX = resolve_index_path() is not None

EXAMPLES = {
    "乙醇": "CCO",
    "乙酸乙酯": "CCOC(C)=O",
    "咖啡因": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "三氟甲苯（19F）": "FC(F)(F)c1ccccc1",
    "氟苯（邻/间/对 H-F）": "Fc1ccccc1",
    "磷酸三甲酯（31P）": "COP(=O)(OC)OC",
    "2,2-二氟丙酸": "CC(F)(F)C(=O)O",
    "喹啉": "c1ccc2ncccc2c1",
}

BATCH_DEFAULT = "CCO\nCCOC(C)=O\nCn1cnc2c1c(=O)n(C)c(=O)n2C\nFC(F)(F)c1ccccc1"

SOLVENTS = ["CDCl3", "DMSO-d6", "CD3OD", "D2O", "acetone-d6", "C6D6", "CD3CN"]
BATCH_LIMIT = 30

#: 结果存活在会话里而不是「这一次运行」里。
#: st.download_button 被点击同样会触发一次 rerun（所有控件都会），
#: 如果结果只存在于「点了生成」那一支，rerun 后 `st.button` 返回 False，
#: 整页就退回空状态——用户看到的是「一点下载就跳回首页，下载不了第二张图」。
_JOB_KEY = "nmrflow_job"
_AUTOSTOP_KEY = "autostop_enabled"


@st.fragment(run_every=30)
def _heartbeat_fragment(session_id: str) -> None:
    """页面存活心跳；所有页面关闭后由 autostop 守护线程退出。

    为什么不是 ``beforeunload``：浏览器不能执行本地脚本，而且刷新、跳转、
    崩溃都会误触发。心跳只回答"还有没有活着的页面"，宽限期由守护线程负责。
    """
    autostop.touch(session_id)
    st.empty()


def _on_autostop_change() -> None:
    autostop.set_enabled(bool(st.session_state.get(_AUTOSTOP_KEY)))

#: 右上角版本号：固定定位，不占正文布局，字号比正文小一档。
#: 用户报问题时第一句总是"你用的哪一版"，界面上能直接看到比翻目录可靠。
#: 版本号来自 ``nmrflow.__version__``，不在这里再抄一遍。
_VERSION_BADGE = f"""
<style>
.nmrflow-version {{
    position: fixed;
    top: 4.1rem;
    right: 1.2rem;
    z-index: 9999;
    font-size: 0.72rem;
    line-height: 1.2;
    color: rgba(128, 128, 128, 0.9);
    background: rgba(128, 128, 128, 0.08);
    border-radius: 4px;
    padding: 2px 6px;
    pointer-events: none;
    user-select: none;
}}
</style>
<div class="nmrflow-version" title="nmrflow version">v{__version__}</div>
"""

#: 「示例分子」下拉里的自定义项；默认就选它，不让用户以为输入框属于某个示例
CUSTOM = "（自定义）"
DEFAULT_SMILES = "CCO"

#: 用户输入框的会话键。用显式 key 而不是每次传 ``value=``：
#: ``value`` 一变 Streamlit 就重建控件，选「自定义」时会把用户刚打的 SMILES 清掉。
_SMILES_KEY = "smiles_input"
#: 记录「当前输入框里是哪个示例填进去的」，用来避免重复写入
_EXAMPLE_KEY = "example_applied"


def _default_results_root() -> Path:
    """优先写包根目录下的 results/；开发期写仓库根；装成 wheel 时退回 cwd。"""
    here = Path(__file__).resolve()
    for root in (here.parents[2], here.parents[1], Path.cwd()):
        if (root / "runtime").is_dir() or (root / "pyproject.toml").is_file():
            return root / "results"
    return Path.cwd() / "results"


def _prepare_outdir(slug: str) -> tuple[Path, str]:
    """建 results/<slug>/；目录不可写时退回临时目录，保证仍能出图。"""
    outdir = _default_results_root() / slug
    try:
        outdir.mkdir(parents=True, exist_ok=True)
        probe = outdir / ".nmrflow_write_test"
        probe.write_text("", encoding="utf-8")
        probe.unlink()
        return outdir, ""
    except OSError as exc:
        fallback = Path(tempfile.mkdtemp(prefix="nmrflow-"))
        return fallback, f"结果目录不可写（{exc}），本次已改存到临时目录：{fallback}"


def _report_file_text(result) -> str:
    """report.txt 的内容：正文 + 预测声明 + 来源声明。"""
    tails = [
        PREDICTED_DATA_NOTE_TEMPLATE.format(version=__version__),
    ]
    attribution = result.report.attribution
    if attribution:
        tails.append(attribution)
    return result.report.text + ("\n\n" + "\n\n".join(tails) if tails else "") + "\n"


def _slug(smiles: str) -> str:
    keep = "".join(c if c.isalnum() else "_" for c in smiles)
    return keep.strip("_")[:40] or "mol"


def _structure_png(smiles: str) -> bytes | None:
    """输入结构的 2D 图；解析失败返回 None。"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        img = Draw.MolToImage(mol, size=(360, 260))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:
        return None


def _zip_dir(path: Path) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file in sorted(path.rglob("*")):
            if file.is_file():
                zf.write(file, file.relative_to(path).as_posix())
    return buf.getvalue()


def _write_result_files(result, slug: str, outdir: Path) -> dict:
    """把 report.txt / 各核 PNG / MS PNG 写到 outdir，返回路径。"""
    outdir.mkdir(parents=True, exist_ok=True)
    report_path = outdir / f"{slug}_report.txt"
    report_path.write_text(_report_file_text(result), encoding="utf-8")

    peaks_of = {b.nucleus: b.peaks for b in result.report.blocks}
    spectra: dict[str, Path] = {}
    for nucleus, spec in result.spectra.items():
        path = outdir / f"{slug}_{nucleus}.png"
        if nucleus == "1H":
            plot_1h(spec, path, peaks=peaks_of.get(nucleus, []))
        else:
            plot_generic(spec, path, peaks=peaks_of.get(nucleus, []))
        spectra[nucleus] = path

    ms_path = None
    if result.ms is not None and result.ms.isotope_pattern:
        ms_path = outdir / f"{slug}_MS.png"
        plot_ms(result.ms, ms_path)
    output_paths = [report_path, *spectra.values()]
    if ms_path is not None:
        output_paths.append(ms_path)
    manifest_path = write_manifest(result, outdir, output_paths)
    return {
        "report": report_path,
        "spectra": spectra,
        "ms": ms_path,
        "manifest": manifest_path,
    }


def _run_options(
    solvent, mhz, ms_method, show_ms, tier_b, hose_hetero
) -> RunOptions:
    return RunOptions(
        solvent=solvent,
        proton_mhz=float(mhz),
        ms_method=ms_method,
        with_ms=show_ms,
        tier_b=tier_b,
        tier_b_hetero=hose_hetero,
    )


def _render_peak_tables(result) -> None:
    with st.expander("峰表（逐峰）"):
        for block in result.report.blocks:
            st.markdown(f"**{block.nucleus} NMR**（{block.freq_mhz:.0f} MHz）")
            rows = [
                {
                    "δ": f"{p.shift:.2f}",
                    "多重度": p.multiplicity,
                    "J / Hz": ", ".join(f"{j:.1f}" for j in p.j_values),
                    "积分": f"{p.integration:g}",
                    "归属": p.label,
                    "置信度": p.confidence.value,
                }
                for p in block.peaks
            ]
            st.dataframe(rows, width="stretch", hide_index=True)


def _render_job_header(job: dict) -> None:
    """结果区顶部：说明这一屏是哪次生成的 + 清空入口。

    放在主区而不是侧栏，是因为侧栏先于「生成」渲染，刚生成完那一次侧栏还看不到
    job，按钮会等到下一次交互才冒出来。清空后立刻 rerun，避免同一遍里继续渲染旧内容。
    """
    left, right = st.columns([5, 1], vertical_alignment="center")
    with left:
        st.caption(f"当前显示：{job['caption']}（生成于 {job['stamp']}）")
    with right:
        if st.button("清空结果", key="clear_job", width="stretch"):
            st.session_state.pop(_JOB_KEY, None)
            st.rerun()


def _render_single(job: dict) -> None:
    result = job["result"]
    outdir: Path = job["outdir"]
    slug: str = job["slug"]
    files = job["files"]

    _render_job_header(job)

    left, right = st.columns([1.0, 1.4])

    with left:
        st.subheader("报道值")
        st.code(result.report.text, language="text")
        st.caption(PREDICTED_DATA_NOTE_TEMPLATE.format(version=__version__))
        st.download_button(
            "下载 report.txt", data=job["report_bytes"],
            file_name=files["report"].name, mime="text/plain", key="dl_report",
        )
        st.download_button(
            "下载全部结果（ZIP）", data=job["zip_bytes"],
            file_name=f"{slug}_results.zip", mime="application/zip", key="dl_zip",
        )
        if result.report.attribution:
            st.caption(result.report.attribution.replace("\n", " "))

        structure = job["structure"]
        if structure is not None:
            st.subheader("结构")
            st.image(structure, caption="输入 SMILES 的 2D 结构", width=340)

        meta = {
            "分子式": result.report.formula,
            "溶剂": result.report.solvent,
            "1H 频率": f"{result.report.proton_mhz:.0f} MHz",
            "出图核种": " / ".join(b.nucleus for b in result.report.blocks),
            "结果目录": str(outdir),
        }
        st.table({"项": list(meta), "值": [str(v) for v in meta.values()]})
        if hasattr(os, "startfile"):
            if st.button("打开结果目录", width="stretch"):
                try:
                    os.startfile(outdir)  # type: ignore[attr-defined]
                except OSError as exc:
                    st.info(f"无法打开目录：{exc}")

        if result.report.warnings:
            st.warning("\n".join(result.report.warnings))

    with right:
        st.subheader("谱图")
        for nucleus, path in files["spectra"].items():
            data = path.read_bytes()
            st.image(data, caption=f"{nucleus} NMR", width="stretch")
            st.download_button(
                f"下载 {nucleus} 谱图（PNG）", data=data, file_name=path.name,
                mime="image/png", key=f"dl_{nucleus}",
            )
        if files["ms"] is not None:
            data = files["ms"].read_bytes()
            st.image(data, caption="同位素分布", width="stretch")
            st.download_button(
                "下载质谱图（PNG）", data=data, file_name=files["ms"].name,
                mime="image/png", key="dl_ms",
            )

    _render_peak_tables(result)


def _parse_batch(text: str) -> list[str]:
    out = []
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def _build_single(smiles: str, options: RunOptions) -> dict:
    """跑一次预测并落盘，把渲染所需的一切打包成 job。

    单独成函数是为了能放进 ``session_state``：下载按钮被点击时会触发 rerun，
    若结果只活在「点了生成」那一支里，rerun 后整页就退回去了。
    """
    result = run(smiles, options)
    slug = _slug(smiles)
    outdir, warning = _prepare_outdir(slug)
    files = _write_result_files(result, slug, outdir)
    return {
        "mode": "single",
        "result": result,
        "slug": slug,
        "outdir": outdir,
        "files": files,
        "report_bytes": _report_file_text(result).encode("utf-8"),
        # ZIP 在这里算一次就够；每次 rerun 重算是 200ms 级的浪费
        "zip_bytes": _zip_dir(outdir),
        "structure": _structure_png(smiles),
        "warning": warning,
        "stamp": time.strftime("%H:%M:%S"),
        "caption": f"{smiles} · {options.solvent} · {options.proton_mhz:.0f} MHz",
    }


def _render_batch(job: dict) -> None:
    rows = job["rows"]
    successes = job["successes"]
    batch_dir: Path = job["batch_dir"]
    total = job["total"]

    _render_job_header(job)
    if job["note"]:
        st.warning(job["note"])
    ok = len(successes)
    st.success(f"完成 {ok}/{total}；结果目录：{batch_dir}")
    st.dataframe(rows, width="stretch", hide_index=True)

    if ok:
        st.download_button(
            "下载全部结果（ZIP）", data=job["zip_bytes"],
            file_name=job["zip_name"], mime="application/zip",
            key="dl_batch_zip",
        )

    for idx, smi, result, outdir in successes:
        with st.expander(f"{idx:02d} · {smi}"):
            st.code(result.report.text, language="text")
            st.caption(PREDICTED_DATA_NOTE_TEMPLATE.format(version=__version__))
            st.caption(f"文件：{outdir}")
            if result.report.warnings:
                st.warning("\n".join(result.report.warnings))


def _build_batch(text: str, options: RunOptions) -> tuple[dict | None, str]:
    """批量跑完再返回 job；输入为空时返回 (None, 错误说明)。"""
    smiles_list = _parse_batch(text)
    if not smiles_list:
        return None, "批量输入为空：每行一个 SMILES。"
    note = ""
    if len(smiles_list) > BATCH_LIMIT:
        note = f"一次最多 {BATCH_LIMIT} 个，已只取前 {BATCH_LIMIT} 个。"
        smiles_list = smiles_list[:BATCH_LIMIT]

    stamp = time.strftime("%Y%m%d-%H%M%S")
    batch_dir, warning = _prepare_outdir(f"batch_{stamp}")

    rows = []
    successes = []
    progress = st.progress(0.0, text="批量生成中…")
    for idx, smi in enumerate(smiles_list, 1):
        try:
            result = run(smi, options)
            slug = f"{idx:02d}_{_slug(smi)}"
            outdir = batch_dir / slug
            _write_result_files(result, slug, outdir)
            successes.append((idx, smi, result, outdir))
            rows.append({
                "#": idx,
                "SMILES": smi,
                "分子式": result.report.formula,
                "状态": "OK",
                "结果目录": str(outdir),
            })
        except Exception as exc:  # 单个失败不能中断整批
            rows.append({
                "#": idx,
                "SMILES": smi,
                "分子式": "—",
                "状态": f"失败：{exc}",
                "结果目录": "—",
            })
        progress.progress(idx / len(smiles_list), text=f"{idx}/{len(smiles_list)}")

    progress.empty()
    return {
        "mode": "batch",
        "rows": rows,
        "successes": successes,
        "batch_dir": batch_dir,
        "zip_bytes": _zip_dir(batch_dir) if successes else b"",
        "zip_name": f"nmrflow_batch_{stamp}.zip",
        "total": len(smiles_list),
        "warning": warning,
        "note": note,
        "stamp": time.strftime("%H:%M:%S"),
        "caption": f"批量 {len(smiles_list)} 个",
    }, ""


def main() -> None:
    st.set_page_config(page_title="nmrflow · NMR 谱图预测", layout="wide")
    st.markdown(_VERSION_BADGE, unsafe_allow_html=True)
    st.title("nmrflow — 由 SMILES 预测 NMR 谱图与报道值")
    st.caption(
        "输入 SMILES，自动输出 ¹H / ¹³C 谱；检测到 F 追加 ¹⁹F、检测到 P 追加 ³¹P，"
        "并给出文献风格报道值与 HRMS 预测。支持批量与整包下载。"
    )

    # 每个页面每 30 秒 touch 一次；多标签页各自有 session_id，共享同一个最新心跳。
    if "_autostop_sid" not in st.session_state:
        st.session_state["_autostop_sid"] = uuid.uuid4().hex
    _heartbeat_fragment(st.session_state["_autostop_sid"])

    with st.sidebar:
        st.header("输入")
        batch_mode = st.checkbox("批量模式（每行一个 SMILES）", value=False)
        if batch_mode:
            smiles_text = st.text_area(
                "SMILES 列表", value=BATCH_DEFAULT, height=180,
                help="每行一个；# 后为注释；一次最多 30 个。",
            )
            smiles = ""
        else:
            picked = st.selectbox(
                "示例分子", [CUSTOM] + list(EXAMPLES), index=0,
                help="选一个示例会把它的 SMILES 填进下面的输入框；"
                     "选「（自定义）」则保留你自己输入的内容。",
            )
            if picked == CUSTOM:
                # 不覆盖用户已经输入的内容——否则「自定义」反而会把手打的 SMILES 清掉
                st.session_state.setdefault(_SMILES_KEY, "")
                st.session_state.pop(_EXAMPLE_KEY, None)
            elif st.session_state.get(_EXAMPLE_KEY) != picked:
                # 只在「换了一个示例」时写入；同一个示例重复渲染不重写，
                # 免得用户改完示例分子再点别的控件就被打回原值
                st.session_state[_SMILES_KEY] = EXAMPLES[picked]
                st.session_state[_EXAMPLE_KEY] = picked
            smiles = st.text_input(
                "SMILES", key=_SMILES_KEY,
                placeholder=f"例如 {DEFAULT_SMILES}、c1ccccc1",
            )
            smiles_text = ""

        st.header("参数")
        solvent = st.selectbox("氘代溶剂", SOLVENTS, index=0)
        mhz = st.select_slider("质子观测频率 (MHz)",
                               options=[60, 90, 200, 300, 400, 500, 600, 800],
                               value=400)
        ms_method = st.selectbox("质谱方法", ["DART", "ESI", "APCI", "MALDI", "CI", "EI"])
        show_ms = st.checkbox("预测质谱（HRMS）", value=True)
        use_tier_b = st.checkbox(
            "HOSE 数据库校正（Tier B）", value=_HAS_INDEX,
            disabled=not _HAS_INDEX,
            help=(
                "用 NMRShiftDB2 的同环境实测位移覆盖经验值。¹H 4 层球 MAE 0.13 ppm、"
                "¹³C 1.21 ppm；命中球层不足时自动保留经验值。"
                if _HAS_INDEX else "未找到 HOSE 索引，先运行 tools/build_hose_index.py"
            ),
        )
        hose_hetero = st.checkbox(
            "¹⁹F / ³¹P 也走 HOSE 校正", value=True,
            disabled=not (use_tier_b and _HAS_INDEX),
            help="个别 ³¹P 分子（如磷酸三甲酯）会被 HOSE 带偏，可在这里关掉。",
        )
        auto_stop = st.checkbox(
            "关闭所有页面后自动停止服务",
            value=autostop.is_enabled(),
            key=_AUTOSTOP_KEY,
            on_change=_on_autostop_change,
            help=(
                "所有浏览器页面关闭 3 分钟后自动退出（刷新、切标签、短暂断网不会停）。"
                "默认关闭；想立刻退出可以双击 停止 NMR.bat。"
            ),
        )
        if auto_stop:
            st.caption("已启用：所有页面关闭 3 分钟后自动停止；多标签只要还有一个开着就不会停。")
        go = st.button("批量生成" if batch_mode else "生成谱图",
                       type="primary", width="stretch")

    options = _run_options(
        solvent, mhz, ms_method, show_ms, use_tier_b, hose_hetero
    )
    wanted = "batch" if batch_mode else "single"
    job = st.session_state.get(_JOB_KEY)

    if go:
        # 生成新结果：失败时清掉旧结果，避免页面上留着上一次的分子
        st.session_state.pop(_JOB_KEY, None)
        job = None
        if not batch_mode and not smiles.strip():
            st.error("请先输入 SMILES；也可以从左侧「示例分子」里选一个。")
            return
        if batch_mode:
            job, err = _build_batch(smiles_text, options)
            if job is None:
                st.error(err)
                return
        else:
            try:
                job = _build_single(smiles, options)
            except (ValueError, KeyError) as exc:
                st.error(f"无法处理这个 SMILES：{exc}")
                return
        st.session_state[_JOB_KEY] = job

    if job is None:
        st.info("在左侧输入 SMILES 后点「生成谱图」；批量模式下每行一个 SMILES。")
        return

    # 模式切换后旧结果不再适用；提示而不是强行渲染到错误的版式里。
    # 不删 job——切回原模式就能继续用，删了反而要重跑一遍。
    if job["mode"] != wanted:
        st.info("已切换模式，当前结果来自上一次生成。点「生成谱图」刷新。")
        return

    if job["warning"]:
        st.warning(job["warning"])

    if job["mode"] == "batch":
        _render_batch(job)
    else:
        _render_single(job)


if __name__ == "__main__":
    main()
