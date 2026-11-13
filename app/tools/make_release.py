"""Assemble the portable release package and zip it.

Layout of the produced folder::

    nmrflow-<version>-portable/
    +- NMR-Predictor.exe      frozen launcher (double-click me)
    +- 启动 NMR.bat            console version, useful for debugging
    +- 安装.bat                (re)create the desktop shortcut
    +- 卸载.bat                remove the desktop shortcut
    +- 安装使用说明.txt         quick start for end users
    +-     app/                    application source (importable as `nmrflow`)
    +- app/tools/            evaluation kit (re-runs the README numbers)
    +- app/build/            the two training indexes the eval kit uses
    +- runtime/python/         relocatable CPython + every dependency
    +- assets/                 icons
    +- examples/               sample SMILES
    +- results/                generated spectra land here
    +- logs/                   streamlit log when something goes wrong
    +- README.md

Usage::

    python tools/make_release.py                 # full build
    python tools/make_release.py --no-zip        # folder only
    python tools/make_release.py --skip-runtime  # reuse build/runtime
    python tools/make_release.py --skip-exe      # reuse a previous exe
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
APP = ROOT / "app"
BUILD = APP / "build"
DIST = ROOT / "dist"
RUNTIME = BUILD / "runtime" / "python"

sys.path.insert(0, str(APP / "tools"))
import version_info as vinfo  # noqa: E402  (tools/ 内的兄弟模块)

# 与 pyproject.toml 保持一致。0.2.0 起便携包含 HOSE 索引（Tier B 默认开启）；
# 0.2.1 只更新了包内 README（FluoBase 对照、开发环境说明），代码与 0.2.0 相同。
VERSION = "0.3.3"
PKG_NAME = f"nmrflow-{VERSION}-portable"
EXE_NAME = "NMR-Predictor.exe"
PORT = 8901

#: 随包的评估脚本。判定标准只有一条：README 里某个数字靠它产生，或者
#: 某个结论靠它推翻。打包工具（make_release / build_runtime / launcher /
#: make_icon / verify_package / win_ansi）不随包——它们只对维护者有意义。
EVAL_TOOLS = [
    "check_eval_kit.py",
    "check_reference_records.py",
    "build_hose_index.py",
    "eval_hose_index.py",
    "eval_tier_ab.py",
    "diag_sp2.py",
    "diag_sp3.py",
    "diag_bicyclic.py",
    "diag_hetero.py",
    "mine_hetero.py",
    "ml_baseline_19f.py",
    "bench_spin.py",
    "README.md",
]

#: 随包的两份**训练索引**（留出 20% 之后的 80%），是复算 3.41 的必需品：
#: 用运行时那份全量索引去评估等于拿答案考自己。约 27 MB。
EVAL_INDICES = ["hose_index_train.sqlite", "hose_index_hetero.sqlite"]

EXCLUDE_TOP = {
    "build",
    "dist",
    "results",
    ".pytest_cache",
    ".workbuddy",
    ".venv",
    "__pycache__",
    ".git",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


# ------------------------------------------------------------- artefacts ---
INSTALL_BAT = """@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist "{exe}" (
    echo [ERROR] {exe} not found next to this script.
    echo.
    pause
    exit /b 1
)
"{exe}" --install
""".replace("{exe}", EXE_NAME)

#: 停止 / 卸载 / 自检 三个脚本太复杂（含参数解析、校验、临时清理脚本），
#: 内联成 Python 字符串既难读又容易在 %% 与引号上出错，因此放在 packaging/bat/
#: 下面按文件维护；打包时只做 {{PORT}} 替换。
BAT_DIR = ROOT / "packaging" / "bat"


def _read_bat(name: str) -> str:
    path = BAT_DIR / name
    if not path.exists():
        raise SystemExit(f"missing batch template: {path}")
    return path.read_text(encoding="utf-8").replace("{{PORT}}", str(PORT))


UNINSTALL_BAT = _read_bat("卸载.bat")
STOP_BAT = _read_bat("停止 NMR.bat")
SELFCHECK_BAT = _read_bat("自检.bat")

LAUNCH_BAT = r"""@echo off
chcp 65001 >nul
setlocal
cd /d "%~dp0"

set "PY=%~dp0runtime\python\python.exe"
set "APP=%~dp0app\nmrflow\app.py"

if not exist "%PY%" (
    echo [ERROR] runtime\python\python.exe is missing.
    echo         The package was not unpacked correctly.
    echo.
    pause
    exit /b 1
)

set "HOME=%~dp0runtime\home"
if not exist "%HOME%\.streamlit" mkdir "%HOME%\.streamlit"
if not exist "%HOME%\.streamlit\credentials.toml" (
    > "%HOME%\.streamlit\credentials.toml" echo [general]
    >> "%HOME%\.streamlit\credentials.toml" echo email = ""
)

echo Starting NMR Spectrum Predictor...
echo Keep this window open; closing it stops the server.
echo.

set "PYTHONIOENCODING=utf-8"
set "PYTHONUTF8=1"
set "MPLBACKEND=Agg"
set "USERPROFILE=%HOME%"
set "PYTHONPATH=%~dp0app"

"%PY%" "%~dp0app\_win_ansi.py"

"%PY%" -m streamlit run "%APP%" ^
  --server.port 8901 ^
  --server.address 127.0.0.1 ^
  --server.headless true ^
  --browser.gatherUsageStats false ^
  --server.fileWatcherType none ^
  --client.toolbarMode minimal

echo.
echo Server stopped.
pause
"""

README_TXT = r"""NMR 谱图预测 —— 安装使用说明
================================================

一、安装（只需三步）
------------------------------------------------
1. 把整个压缩包解压到任意位置。
   建议放在空间充足的磁盘，例如 D:\NMR\ 。
   解压后请不要只拷贝单个文件，必须保留完整的目录结构，
   尤其是 runtime 文件夹（里面自带 Python 和全部依赖）。

2. 双击目录里的  NMR-Predictor.exe

3. 首次启动会询问“是否在桌面创建快捷方式”，选“是”即可。
   之后直接双击桌面上的「NMR 谱图预测」图标就能打开。

   注意：解压目录不要放在中文或空格过多的深层路径下，
   也不要直接放在需要管理员权限的目录（如 C:\Windows\）。

二、使用
------------------------------------------------
1. 浏览器会自动打开 http://127.0.0.1:{port}
   若没有自动打开，手动在浏览器输入该地址即可。

2. 左侧填写：
   - SMILES 结构式（必填）
   - 溶剂：CDCl3 / DMSO-d6 / CD3OD / D2O 等
   - 核磁频率：400 MHz（也可填 300 / 500 / 600）

   左侧还有一个「HOSE 数据库校正（Tier B）」勾选项，默认已勾上：
   程序用内置的离线索引（约 32 MB，来自 NMRShiftDB2 实验谱）去修正
   经验规则的偏差，1H 与 13C 通常能明显改善。取消勾选可对比纯规则的结果。
   其下的「19F / 31P 也走 HOSE 校正」默认勾上；个别 31P 分子（如磷酸三甲酯）
   被带偏时可取消。

   质谱默认只给 Calculated 精确质量，不伪造实测值 Found。
   界面上不提供「模拟 Found」开关（内测期的排版演示功能，正式版已隐藏入口；
   需要在命令行用 --found 打开时，导出的 report.txt 会带 [SIMULATED] 标注）。

   全部计算都在本机完成，不联网、不上传结构式。

3. 点击运行后，右侧给出：
   - 1H NMR 谱图（含积分曲线与峰顶标注）
   - 13C NMR 谱图
   - 若结构含 F，自动追加 19F NMR
   - 若结构含 P，自动追加 31P NMR
   - 文献风格的报道值文本，以及质谱（HRMS）预测

4. 生成后，谱图 PNG 和报道值 TXT 保存在解压目录下的 results\<SMILES>\ 文件夹；
   页面上每个文件也有下载按钮（谱图不用右键另存），还可以「下载全部结果（ZIP）」。

5. 结构预览与批量：
   - 左侧勾选「批量模式（每行一个 SMILES）」可一次跑最多 30 个分子，
     结果写到 results\batch_时间戳\，并提供一个整包 ZIP。
   - 单分子模式下，报道值左侧会显示输入的 2D 结构图。

三、关闭
------------------------------------------------
直接关闭浏览器标签页不会停止后台服务（这是正常的，重新打开页面就能接着用）。
要立即退出：双击目录里的  停止 NMR.bat  （它结束监听 {port} 的进程，并复查端口是否释放）。

也可以让它自己停：在左侧勾选「关闭所有页面后自动停止服务」，之后所有页面关掉
3 分钟就自动退出。默认不勾——因为刷新、切标签、电脑休眠都可能让人以为"页面还开着"，
宽限期就是留给这些情况的；正在跑批量任务时也建议先别勾。

四、常见问题
------------------------------------------------
Q：双击后没有反应？
A：先确认 runtime 文件夹完整；若目录里出现 logs\streamlit.log，把该文件发给维护者。
   最常见的原因是解压不完整（runtime 文件夹缺失）。

Q：提示端口被占用？
A：说明服务已经在运行，直接在浏览器打开 http://127.0.0.1:{port} 即可。

Q：谱图里的中文显示为方块？
A：程序使用系统自带的黑体（simhei.ttf）。若目标系统缺少该字体，
   安装任意中文字体后重启程序即可。

Q：提示缺少 VC++ 运行库？
A：安装 Microsoft Visual C++ Redistributable (x64) 后重试。

Q：双击 exe 提示“未知发布者”，被 SmartScreen 拦截？
A：exe 没有数字签名（本项目没有代码签名证书）。可以核对 SHA256SUMS.txt 里
   NMR-Predictor.exe 的哈希，或改用  启动 NMR.bat  启动；两条路径跑的是同一套 app。

Q：怎么知道自己用的是哪个版本？
A：右键 NMR-Predictor.exe → 属性 → 详细信息，看“文件版本”（例如 {version}）；
   也可以在命令行执行  NMR-Predictor.exe --version  （窗体程序没有控制台，会弹小窗）。
   报问题时请带上版本号与构建日期，这两项在弹窗里都有。

Q：为什么个别预测和文献差了几个 ppm？
A：1H/13C 的规则来自经验加和，本身有系统误差；勾选 HOSE 校正后多数情况会更好，
   但个别分子（尤其 19F/31P）也可能变差。想要两边对比就取消勾选再跑一次。

Q：想换一台电脑？
A：把整个目录（或原压缩包）复制过去即可，不需要重新安装。

Q：想确认这个包没坏/没被改过？
A：双击  自检.bat  。它会核对随包的评估脚本、三份索引（含逐字节哈希）以及
   训练索引是否真的排除了留出集，最后打印 PASS/FAIL。核不上就是包不完整。

Q：如何卸载？
A：双击  卸载.bat  打开卸载向导：
     [1] 只删桌面快捷方式，程序和 results 都保留；
     [2] 彻底删除（会先把 results 移到桌面，再停服务、删目录）；
     [3] 取消。
   选 2 需要手动输入 DELETE 确认；服务停不掉时会中止删除，不会留下删一半的目录。
   本程序不写注册表，不在系统目录留文件。

五、许可
------------------------------------------------
本程序：BSD 3-Clause License，全文见随包的 LICENSE 文件。

程序内置的 HOSE 位移索引来自 NMRShiftDB2 实验谱数据，按
nmrshiftdb2 Database License 发布，许可证原文见
app\nmrflow\data\nmrshiftdb2_database_license.txt。
按照该许可的要求，使用本程序得到的结果在公开时请附带以下声明：

Contains information from nmrshiftdb2 (www.nmrshiftdb.org),
which is available here under the nmrshiftdb2 Database License.

（导出的 report.txt 与 JSON 末尾已自动带这行，不需要手工补。）

六、目录说明
------------------------------------------------
NMR-Predictor.exe   主程序，双击启动
启动 NMR.bat        命令行版启动器，排错时可看到完整日志
停止 NMR.bat        结束监听 {port} 的服务（可带端口参数：停止 NMR.bat 8902）
安装.bat            重新创建桌面快捷方式
卸载.bat            卸载向导（删快捷方式 / 删整个目录）
自检.bat            核对包完整性与评估套件（打印 PASS/FAIL）
app\                程序源码
app\tools\          评估套件：复算 README 里的留出法误差（见 app\tools\README.md）
app\build\          评估用的训练索引，约 26 MB；程序运行时不需要，删掉不影响使用
runtime\            内置 Python 与全部依赖（体积最大，请勿删除）
assets\             图标
examples\           示例 SMILES 列表
results\            生成的谱图与报道值（首次生成时自动创建）
LICENSE             本程序的 BSD-3-Clause 许可全文

七、想自己核对精度（可跳过）
------------------------------------------------
README 里的误差数字（1H 0.39 / 13C 3.41 / 19F 8.76 / 31P 19.71）不是让你"相信"的，
是可以复算的。先跑自检：

    runtime\python\python.exe app\tools\check_eval_kit.py

九项全 PASS 之后，按 app\tools\README.md 里的命令跑评估。你还需要自己下载
NMRShiftDB2 的 nmrshiftdb2withsignals.sd（151 MB，SourceForge 免费，
不需要注册）——它不随包，原因写在那个文件里。评估用的是随包的"训练索引"
（80% 化合物），不是程序运行时那份全量索引，否则等于拿答案考自己。
""".replace("{port}", str(PORT)).replace("{version}", VERSION)


# ----------------------------------------------------------------- steps ---
def copy_source(pkg: Path) -> None:
    app = pkg / "app"
    app.mkdir(parents=True, exist_ok=True)
    log("copying application source")
    # Overlay rather than delete-then-copy: removing a large tree trips the
    # workspace's bulk-delete guard, and the source layout rarely loses files.
    shutil.copytree(
        APP / "nmrflow",
        app / "nmrflow",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache"),
        dirs_exist_ok=True,
    )
    # LICENSE 必须随包：项目自身是 BSD-3-Clause，而内置 HOSE 索引的许可 4.5
    # 要求依赖它的软件必须是 OSI 认可的开源许可，用户得能看到原文。
    # requirements-dev.txt 也要随包：README 让开发者用 `pip install -r requirements-dev.txt`
    # 跑测试，包里没有这个文件的话那句话就是死的。
    for extra in ("README.md", "requirements.txt", "requirements-dev.txt",
                  "pyproject.toml", "LICENSE", "CHANGELOG.md", "REPRODUCE.md"):
        src = ROOT / extra
        if src.exists():
            shutil.copy2(src, pkg / extra)
    if (ROOT / "examples").exists():
        shutil.copytree(
            ROOT / "examples",
            pkg / "examples",
            ignore=shutil.ignore_patterns("__pycache__"),
            dirs_exist_ok=True,
        )
    if (ROOT / "assets").exists():
        shutil.copytree(ROOT / "assets", pkg / "assets", dirs_exist_ok=True)
    # ship an empty results/ and logs/ -- never leftovers from local test runs
    for folder in ("results", "logs"):
        target = pkg / folder
        target.mkdir(parents=True, exist_ok=True)
        for item in target.iterdir():
            try:
                item.unlink()
            except OSError:
                pass

    # the console launcher needs the ANSI helper too
    shutil.copy2(APP / "tools" / "win_ansi.py", app / "_win_ansi.py")

    # 包内回归只带两个文件：stdlib 入口 + 那 17 条用例。
    # 仓库里其余测试是 pytest 风格（依赖 fixture），便携 runtime 没有 pytest。
    tests_dir = pkg / "tests"
    tests_dir.mkdir(exist_ok=True)
    for name in ("run_tests.py", "test_regressions.py"):
        shutil.copy2(ROOT / "tests" / name, tests_dir / name)


def copy_eval_kit(pkg: Path) -> None:
    """把评估脚本与两份训练索引放进 ``app/`` 下。

    为什么不放包根目录：脚本第一段就是
    ``sys.path.insert(0, Path(__file__).parent.parent)``——只有让 ``tools/`` 与
    ``nmrflow/`` 同级，``import nmrflow`` 才成立。放在 ``app/tools/`` 之后，
    源码仓库与便携包两边**不需要给随包版本打任何补丁**，同一份文件两处可跑。

    训练索引（80% 化合物）必须随包：运行时用的是全量索引，拿它做评估等于
    拿答案考自己；重建一次要 6 分钟，也没有必要。
    """
    tools_dir = pkg / "app" / "tools"
    tools_dir.mkdir(parents=True, exist_ok=True)
    log(f"copying evaluation kit -> {tools_dir}")
    for name in EVAL_TOOLS:
        src = APP / "tools" / name
        if not src.exists():
            raise SystemExit(f"eval tool missing: {src}")
        shutil.copy2(src, tools_dir / name)

    build_dir = pkg / "app" / "build"
    build_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        "# nmrflow " + VERSION + " eval kit — SHA256",
        "# 路径相对 app/；用 tools/check_eval_kit.py 自动核对",
        "",
    ]
    for name in EVAL_INDICES:
        src = APP / "build" / name
        if not src.exists():
            raise SystemExit(f"training index missing: {src} -- rerun eval_tier_ab")
        shutil.copy2(src, build_dir / name)
        lines.append(f"{_sha256(src)}  build/{name}")
    lines.append(f"{_sha256(pkg / 'app' / 'nmrflow' / 'data' / 'hose_index.sqlite')}"
                 "  nmrflow/data/hose_index.sqlite")
    (build_dir / "EVAL_KIT_SHA256.txt").write_text("\n".join(lines) + "\n",
                                                   encoding="utf-8")
    total = sum((build_dir / n).stat().st_size for n in EVAL_INDICES)
    log(f"eval kit ready: {len(EVAL_TOOLS)} scripts, indices {total / 1e6:.1f} MB")


def ensure_runtime(pkg: Path, skip: bool) -> None:
    """``--skip-runtime``：包内已有 runtime 就沿用，没有就从 build/runtime 拷。

    **不能把「跳过」理解成「什么都不做」**——换版本号后 ``dist/<新版本>`` 是空目录，
    原写法会直接产出一个没有解释器的包（17 MB），只有 check_layout 才拦得住。
    和 ``reuse_exe`` 是同一类假设，一并修掉。
    """
    dst = pkg / "runtime" / "python"
    if skip and (dst / "python.exe").exists():
        log(f"reusing runtime already in package: {dst}")
        return
    if skip:
        log("no runtime in package yet -- copying from build/runtime anyway")
    copy_runtime(pkg)


def copy_runtime(pkg: Path) -> None:
    if not (RUNTIME / "python.exe").exists():
        raise SystemExit(
            "build/runtime/python missing -- run tools/build_runtime.py first"
        )
    dst = pkg / "runtime" / "python"
    log(f"copying runtime -> {dst} (this takes a while)")
    # Overlay in place: wiping 8000+ files trips the bulk-delete guard.
    # __pycache__ must travel with the package -- it holds the pre-compiled
    # bytecode that keeps first-run startup fast and race-free.
    shutil.copytree(
        RUNTIME,
        dst,
        ignore=shutil.ignore_patterns("tests"),
        dirs_exist_ok=True,
    )


def bat_bytes(text: str) -> bytes:
    """bat 的最终字节：CRLF + 中文控制台能直接解析的编码。

    两条硬约束都是 0.2.8 踩出来的：

    1. **必须 CRLF。** LF-only 时 cmd 会把多行 ``if``/``for``、标签、``%~1`` 子串
       拼错——实测 ``卸载.bat`` 被拆成 ``'cho'``、``'HOICE"=="3" goto cancel'``
       这类碎片去执行，卸载向导完全不可用。
    2. **含中文时必须是 GBK，且不能带 ``chcp 65001``。** cmd 按**当前代码页**
       读批处理，而 ``chcp`` 本身是第 2 行才生效，之前读进的缓冲区仍按 CP936
       解码：三字节的 UTF-8 汉字被当成 1.5 个 GBK 字符，行就从汉字中间断掉。
       实测 CRLF + UTF-8 + chcp 仍有 ``'会' is not recognized`` 两处；换 GBK
       并把 chcp 去掉之后是 0 处。

    纯 ASCII 的脚本原样按 ASCII 写——它们带 ``chcp 65001`` 是给下面 Python 的
    输出用的，没有解码风险。
    """
    data = text.replace("\r\n", "\n").replace("\n", "\r\n")
    if data.isascii():
        return data.encode("ascii")
    if "chcp 65001" in data.lower():
        raise SystemExit("含非 ASCII 的 bat 不能带 chcp 65001：cmd 会按 CP936 读成乱码")
    try:
        return data.encode("gbk")
    except UnicodeEncodeError as exc:
        raise SystemExit(f"bat 里有 GBK 无法表示的字符：{exc}") from None


def write_scripts(pkg: Path) -> None:
    log("writing batch helpers and docs")
    for name, text in (
        ("安装.bat", INSTALL_BAT),
        ("卸载.bat", UNINSTALL_BAT),
        ("启动 NMR.bat", LAUNCH_BAT),
        ("停止 NMR.bat", STOP_BAT),
        ("自检.bat", SELFCHECK_BAT),
    ):
        (pkg / name).write_bytes(bat_bytes(text))
    # 说明书是给记事本看的，保持 UTF-8（Win10 起自动识别），与 bat 的规则无关
    (pkg / "安装使用说明.txt").write_bytes(README_TXT.encode("utf-8"))


def check_bat_files(pkg: Path) -> list[str]:
    """打包后按字节核对每个 bat——这类问题在 Python 里看不出来，只有字节能证明。"""
    problems: list[str] = []
    for path in sorted(pkg.glob("*.bat")):
        data = path.read_bytes()
        lone_lf = data.count(b"\n") - data.count(b"\r\n")
        if lone_lf:
            problems.append(f"{path.name}: {lone_lf} 处纯 LF")
        if data.isascii():
            continue
        try:
            data.decode("utf-8")
        except UnicodeDecodeError:
            pass
        else:
            problems.append(f"{path.name}: 非 ASCII 但编码是 UTF-8（应为 GBK）")
        if b"chcp 65001" in data.lower():
            problems.append(f"{path.name}: 含非 ASCII 又带 chcp 65001")
        try:
            data.decode("gbk")
        except UnicodeDecodeError:
            problems.append(f"{path.name}: 非 ASCII 也不是合法 GBK")
    return problems


def build_exe(pkg: Path) -> None:
    log("freezing the launcher with PyInstaller")
    spec_dir = BUILD / "pyinstaller"
    spec_dir.mkdir(parents=True, exist_ok=True)
    # 版本资源：右键属性里能看到的那些字段。来源只有一个（nmrflow.__version__），
    # 这里再做一次交叉核对，避免出现"包名 0.2.8、exe 里写 0.2.7"这种静默错位。
    version = vinfo.read_version(ROOT)
    if version != VERSION:
        raise SystemExit(
            f"版本号不一致：make_release.VERSION={VERSION}，"
            f"nmrflow.__version__={version}"
        )
    vfile = vinfo.write(version, spec_dir / "version_info.txt")
    log(f"version resource: {vfile} (version {version})")
    cmd = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--onefile",
        "--windowed",
        "--name",
        "NMR-Predictor",
        "--icon",
        str(ROOT / "assets" / "nmr.ico"),
        "--version-file",
        str(vfile),
        "--hidden-import",
        "pylnk3",
        "--collect-submodules",
        "pylnk3",
        "--distpath",
        str(pkg),
        "--workpath",
        str(spec_dir / "work"),
        "--specpath",
        str(spec_dir),
        str(APP / "tools" / "launcher.py"),
    ]
    proc = subprocess.run(cmd, cwd=str(ROOT))
    if proc.returncode != 0:
        raise SystemExit("PyInstaller failed")
    exe = pkg / EXE_NAME
    if not exe.exists():
        raise SystemExit(f"{exe} was not produced")
    verify_exe_version(exe, version)
    log(f"exe ready: {exe} ({exe.stat().st_size / 1e6:.1f} MB)")


def verify_exe_version(exe: Path, version: str) -> dict[str, str]:
    """核对 exe 里**真的**有版本资源，且版本号等于当前版本。

    ``--version-file`` 拼错时 PyInstaller 可能只警告不报错，产出的 exe 照样能跑，
    只是属性页空白——这种静默失效只能在成品上验证。
    """
    want = ".".join(str(n) for n in vinfo.version_tuple(version))
    info = vinfo.describe(exe)
    if not info:
        raise SystemExit(f"{exe.name} 没有版本资源（--version-file 没生效）")
    got = info.get("FileVersion")
    if got != want:
        raise SystemExit(f"{exe.name} 版本资源为 {got}，期望 {want}")
    log(f"exe version resource ok: {info.get('ProductName')} {got} / "
        f"{info.get('FileDescription')}")
    return info


def reuse_exe(pkg: Path) -> None:
    """``--skip-exe``：回收一个已有的 exe。

    **不能假设包里已经有了**——换版本号后 ``dist/<新版本>`` 是空目录，
    原来那种"跳过就等于沿用"的写法会产出一个没有主程序的包（且不报错）。
    这里按时间取最新的一份，并把来源与时间打出来，方便核对它对应哪次改动。

    回收之后还必须核对**版本资源**：0.2.7 起 exe 里写着版本号，如果沿用上一版的
    exe，属性页会明晃晃地显示旧版本，而包名、CHANGELOG、SHA256SUMS 都是新的——
    这种不一致比"没有版本号"更糟。核不上就要求重新冻结。
    """
    target = pkg / EXE_NAME
    version = vinfo.read_version(ROOT)
    if target.exists():
        verify_exe_version(target, version)
        log(f"reusing exe already in package: {target}")
        return
    candidates = [BUILD / "pyinstaller" / "dist" / EXE_NAME]
    candidates += list(DIST.glob(f"*/{EXE_NAME}"))
    candidates = [c for c in candidates if c.exists()]
    if not candidates:
        raise SystemExit(
            "no exe to reuse -- run without --skip-exe (needs PyInstaller installed)"
        )
    newest = max(candidates, key=lambda c: c.stat().st_mtime)
    try:
        verify_exe_version(newest, version)
    except SystemExit as exc:
        raise SystemExit(
            f"不能沿用 {newest}\n  {exc}\n"
            "  版本资源对不上，请去掉 --skip-exe 重新冻结。"
        ) from None
    shutil.copy2(newest, target)
    log(
        f"reused exe: {newest} "
        f"({time.strftime('%Y-%m-%d %H:%M', time.localtime(newest.stat().st_mtime))}, "
        f"{newest.stat().st_size / 1e6:.1f} MB)"
    )


def write_checksums(pkg: Path) -> None:
    """给 exe / 索引 / 启停脚本生成 SHA256SUMS.txt。

    exe 没有代码签名，用户从网盘拿到后无法判断有没有被替换；
    给一份可核对的清单，比“请相信我们”有用。
    """
    import hashlib

    targets = [
        EXE_NAME,
        "app/nmrflow/data/hose_index.sqlite",
        "app/build/hose_index_train.sqlite",
        "app/build/hose_index_hetero.sqlite",
        "启动 NMR.bat",
        "停止 NMR.bat",
        "卸载.bat",
        "自检.bat",
    ]
    lines = [
        "# nmrflow " + VERSION + " portable — SHA256",
        "# 生成时间 " + time.strftime("%Y-%m-%d"),
        "# 用法：Get-FileHash <文件> -Algorithm SHA256",
        "",
    ]
    for rel in targets:
        f = pkg / rel
        if not f.exists():
            raise SystemExit(f"checksum target missing: {rel}")
        digest = hashlib.sha256(f.read_bytes()).hexdigest().upper()
        lines.append(f"{digest}  {rel}")
    (pkg / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log("wrote SHA256SUMS.txt")


def make_zip(pkg: Path) -> Path:
    archive = DIST / f"{PKG_NAME}.zip"
    log(f"zipping {pkg} -> {archive}")
    if archive.exists():
        archive.unlink()
    files = sorted(p for p in pkg.rglob("*") if p.is_file())
    total = sum(p.stat().st_size for p in files)
    done = 0
    t0 = time.time()
    with zipfile.ZipFile(
        archive, "w", zipfile.ZIP_DEFLATED, compresslevel=6, strict_timestamps=False
    ) as zf:
        for path in files:
            zf.write(path, str(path.relative_to(pkg.parent)).replace("\\", "/"))
            done += path.stat().st_size
            if done % (200 * 1024 * 1024) < path.stat().st_size:
                log(
                    f"  {done / 1e6:.0f}/{total / 1e6:.0f} MB "
                    f"({time.time() - t0:.0f}s)"
                )
    log(f"archive ready: {archive} ({archive.stat().st_size / 1e6:.0f} MB)")
    return archive


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-runtime", action="store_true",
                    help="包内已有 runtime 时跳过拷贝（没有则仍从 build/runtime 拷）")
    ap.add_argument("--skip-exe", action="store_true")
    ap.add_argument("--skip-source", action="store_true")
    ap.add_argument("--no-zip", action="store_true")
    args = ap.parse_args()

    DIST.mkdir(parents=True, exist_ok=True)
    pkg = DIST / PKG_NAME
    pkg.mkdir(parents=True, exist_ok=True)

    if not args.skip_source:
        copy_source(pkg)
    copy_eval_kit(pkg)          # 评估套件跟脚本清单走，不跟 --skip-source 走
    ensure_runtime(pkg, args.skip_runtime)
    write_scripts(pkg)
    bat_problems = check_bat_files(pkg)
    if bat_problems:
        raise SystemExit("bat 编码/换行不合格：\n  " + "\n  ".join(bat_problems))
    if args.skip_exe:
        reuse_exe(pkg)
    else:
        build_exe(pkg)

    # 必须在 exe 与脚本都就位之后才算校验和
    write_checksums(pkg)

    if not args.no_zip:
        make_zip(pkg)

    missing = check_layout(pkg)
    if missing:
        raise SystemExit("package incomplete, missing: " + ", ".join(missing))

    size = sum(f.stat().st_size for f in pkg.rglob("*") if f.is_file())
    log(f"package ready: {pkg} ({size / 1e6:.0f} MB)")
    return 0


def check_layout(pkg: Path) -> list[str]:
    """打包后自检。少一个关键文件，用户拿到手就是双击没反应。"""
    required = [
        EXE_NAME,
        "启动 NMR.bat",
        "安装.bat",
        "卸载.bat",
        "安装使用说明.txt",
        "app/nmrflow/app.py",
        "app/nmrflow/__init__.py",
        "tests/run_tests.py",
        "tests/test_regressions.py",
        "app/tools/check_eval_kit.py",
        "app/tools/eval_tier_ab.py",
        "app/tools/README.md",
        "app/build/hose_index_train.sqlite",
        "app/build/hose_index_hetero.sqlite",
        "app/build/EVAL_KIT_SHA256.txt",
        "CHANGELOG.md",
        "REPRODUCE.md",
        "SHA256SUMS.txt",
        "停止 NMR.bat",
        "卸载.bat",
        "自检.bat",
        "app/tools/check_reference_records.py",
        "app/_win_ansi.py",
        "app/nmrflow/data/hose_index.sqlite",
        "LICENSE",
        "app/nmrflow/data/nmrshiftdb2_database_license.txt",
        "app/nmrflow/data/HOSE_INDEX_NOTICE.txt",
        "runtime/python/python.exe",
        "runtime/python/Lib/site-packages/rdkit/__init__.py",
        "runtime/python/Lib/site-packages/streamlit/__init__.py",
    ]
    return [rel for rel in required if not (pkg / rel).exists()]


if __name__ == "__main__":
    raise SystemExit(main())
