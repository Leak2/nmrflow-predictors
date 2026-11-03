"""0.2.3 regression suite — stdlib unittest (portable runtime has no pytest).

Run from the package root:
    .\\runtime\\python\\python.exe tests\\run_tests.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


from PIL import Image  # noqa: E402

from nmrflow import __version__  # noqa: E402
from nmrflow.cli import main as cli_main  # noqa: E402
from nmrflow.core.types import Confidence  # noqa: E402
from nmrflow.pipeline import dump_json, run  # noqa: E402
from nmrflow.plot.plot_generic import plot_generic  # noqa: E402
from nmrflow.predict.hose_index import resolve_index_path  # noqa: E402

HAS_INDEX = resolve_index_path() is not None
CHOLESTEROL = "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"
DECALIN = "C1CCC2CCCCC2C1"
NORBORNANE = "C1CC2CCC1C2"
TMP = "COP(=O)(OC)OC"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


class OutputTests(unittest.TestCase):
    def test_default_report_file_has_predicted_note(self):
        work = Path(tempfile.mkdtemp())
        self.assertEqual(cli_main(["CCO", "-o", str(work), "-q"]), 0)
        text = (work / "CCO_report.txt").read_text(encoding="utf-8")
        self.assertIn("Predicted by nmrflow", text)
        self.assertIn("not experimental data", text)
        self.assertNotIn("[SIMULATED]", text)

    def test_default_hrms_only_calculated(self):
        text = run("CCO").report.text
        self.assertIn("Calculated", text)
        self.assertNotIn("Found", text)

    def test_png_has_predicted_metadata(self):
        work = Path(tempfile.mkdtemp())
        path = work / "c.png"
        plot_generic(run("CCO").spectra["13C"], path)
        with Image.open(path) as image:
            self.assertIn("nmrflow", image.info.get("Software", ""))
            self.assertIn("Predicted", image.info.get("Description", ""))

    def test_json_has_top_level_provenance(self):
        work = Path(tempfile.mkdtemp())
        path = dump_json(run("CCO"), work / "x.json")
        data = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], 1)
        self.assertEqual(data["data_type"], "predicted")
        self.assertEqual(data["generator"]["version"], __version__)
        self.assertIn("generated_at", data)
        self.assertNotIn("found_mode", data["ms"])

    def test_found_flag_removed(self):
        with self.assertRaises(SystemExit) as ctx:
            cli_main(["CCO", "--found"])
        self.assertEqual(ctx.exception.code, 2)

    def test_manifest_present_and_fields(self):
        work = Path(tempfile.mkdtemp())
        self.assertEqual(cli_main(["CCO", "-o", str(work), "-q"]), 0)
        path = work / "manifest.json"
        self.assertTrue(path.exists(), "每次运行都应写 manifest.json")
        data = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(data["schema_version"], 1)
        self.assertEqual(data["data_type"], "predicted")
        self.assertEqual(data["generator"]["version"], __version__)
        self.assertEqual(data["input"]["smiles"], "CCO")
        self.assertIn("generated_at", data)
        self.assertIn("CCO_report.txt", data["outputs"])
        self.assertRegex(data["outputs"]["CCO_report.txt"], r"^[0-9a-f]{64}$")

    def test_ethics_doc_exists(self):
        root = Path(__file__).resolve().parents[1]
        self.assertTrue((root / "ETHICS.md").is_file())

    def test_citation_cff_and_zenodo(self):
        root = Path(__file__).resolve().parents[1]
        cff = (root / "CITATION.cff").read_text(encoding="utf-8")
        self.assertIn('version: "0.3.0"', cff)
        self.assertIn("BSD-3-Clause", cff)
        zen = json.loads((root / ".zenodo.json").read_text(encoding="utf-8"))
        self.assertEqual(zen["version"], "0.3.0")
        self.assertEqual(zen["license"], "BSD-3-Clause")

    def test_example_predicted_image(self):
        root = Path(__file__).resolve().parents[1]
        path = root / "assets" / "example-predicted-13c.png"
        self.assertTrue(path.is_file())
        with Image.open(path) as image:
            self.assertIn("nmrflow 0.3.0", image.info.get("Software", ""))
            self.assertIn("Predicted", image.info.get("Description", ""))

    def test_readme_license_faq(self):
        root = Path(__file__).resolve().parents[1]
        text = (root / "README.md").read_text(encoding="utf-8")
        self.assertIn("代码与数据许可 FAQ", text)
        self.assertIn("nmrshiftdb2 Database License", text)
        self.assertIn("CITATION.cff", text)


class SpectrumTests(unittest.TestCase):
    def test_ethanol_ch2_quartet(self):
        spec = run("CCO").spectra["1H"]
        mask = (spec.x_ppm > 3.55) & (spec.x_ppm < 3.80)
        x, y = spec.x_ppm[mask], spec.y[mask]
        maxima = [
            i for i in range(1, len(y) - 1)
            if y[i] >= y[i - 1] and y[i] > y[i + 1]
        ]
        self.assertGreaterEqual(len(maxima), 4, [round(float(x[i]), 4) for i in maxima])

    def test_cholesterol_no_negative_or_giant_peak(self):
        result = run(CHOLESTEROL)
        peaks = [p for b in result.report.blocks if b.nucleus == "1H" for p in b.peaks]
        self.assertGreaterEqual(min(p.shift_lo for p in peaks), -0.01)
        self.assertLessEqual(max(p.shift_hi - p.shift_lo for p in peaks), 0.30)
        self.assertAlmostEqual(
            sum(t.intensity for t in result.spectra["1H"].transitions), 46.0, places=2
        )
        self.assertTrue(any("超过上限" in w for w in result.report.warnings))

    def test_19f_negative_shift_preserved(self):
        spec = run("FC(F)(F)c1ccccc1").spectra["19F"]
        values = [t.freq_hz / spec.freq_mhz for t in spec.transitions]
        self.assertTrue(values)
        self.assertLess(max(values), 0)


class HeteroMSTests(unittest.TestCase):
    def test_acetate_is_m_minus(self):
        text = run("CC(=O)[O-]").report.text
        self.assertIn("C2H3O2-", text)
        self.assertIn("[M]-", text)
        self.assertNotIn("[M+H]", text)

    def test_tetramethylammonium_is_m_plus(self):
        text = run("C[N+](C)(C)C").report.text
        self.assertIn("C4H12N+", text)
        self.assertIn("74.0964", text)
        self.assertIn("[M]+", text)

    def test_salt_counterions_not_reported(self):
        text = run("CC(=O)[O-].[Na+]").report.text
        self.assertEqual(text.count("HRMS"), 1)
        self.assertIn("for C2H3O2-", text)
        self.assertNotIn("for Na+", text)

        text = run("C[N+](C)(C)C.[Cl-]").report.text
        self.assertEqual(text.count("HRMS"), 1)
        self.assertNotIn("for Cl-", text)

    def test_two_organic_ions(self):
        text = run("CC(=O)[O-].C[NH3+]").report.text
        self.assertEqual(text.count("HRMS"), 2)
        self.assertIn("C2H3O2-", text)
        self.assertIn("CH6N+", text)


@unittest.skipUnless(HAS_INDEX, "HOSE index missing")
class ConfidenceAndRingTests(unittest.TestCase):
    def test_tmp_31p_keeps_curated_value(self):
        sig = run(TMP).signals["31P"][0]
        self.assertLess(abs(sig.shift), 5.0)
        self.assertEqual(sig.provenance.tier, "A")
        self.assertIs(sig.confidence, Confidence.HIGH)

    def test_caffeine_13c_not_marked_low(self):
        result = run(CAFFEINE)
        confs = [s.confidence for s in result.signals["13C"]]
        self.assertTrue(all(c is not Confidence.LOW for c in confs), confs)
        self.assertTrue(any("建议人工复核" in w for w in result.report.warnings))

    def test_decalin_uses_tier_b_not_a_patch(self):
        """十氢萘走 HOSE，实测 27.1 / 34.7 / 44.2（NMRShiftDB2 两条独立谱一致）。

        0.2.6 之前这里有一条"小双环脂肪体系 → 放弃 HOSE + 三个硬编码常数"的特例，
        常数是照"34.5 / 26.9 / 26.3"标定的——**那组参考值是错的**，照它改出来的
        34.7 / 26.4 / 25.2 反而比 HOSE 的原始预测（42.5 / 32.6 / 25.9）更远离实测。
        留出集上 262 个多环分子 / 2452 个信号：特例 7.19 ppm，放开 HOSE 4.07 ppm。
        """
        result = run(DECALIN)
        signals = result.signals["13C"]
        values = sorted(s.shift for s in signals)
        expected = [27.1, 34.7, 44.2]
        self.assertLessEqual(max(abs(a - b) for a, b in zip(values, expected)), 3.0, values)
        self.assertTrue(all(s.provenance.tier == "B" for s in signals),
                        [(s.provenance.tier, s.shift) for s in signals])

    def test_norbornane_uses_tier_b(self):
        """降冰片烷 实测 29.6(×4) / 36.3(×2) / 38.3。"""
        result = run(NORBORNANE)
        signals = result.signals["13C"]
        values = sorted(s.shift for s in signals)
        expected = [29.6, 36.3, 38.3]
        self.assertLessEqual(max(abs(a - b) for a, b in zip(values, expected)), 3.0, values)


class MiscTests(unittest.TestCase):
    def test_water_has_no_empty_13c_block(self):
        text = run("O").report.text
        self.assertIn("1H NMR", text)
        self.assertNotIn("13C NMR", text)

    def test_unknown_solvent_returns_one(self):
        work = Path(tempfile.mkdtemp())
        code = cli_main(["CCO", "--solvent", "foo", "-o", str(work)])
        self.assertEqual(code, 1)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
#: 源码仓库布局 tools/ 在根目录；便携包布局在 app/tools/（跟着 nmrflow 一起可导入）
TOOLS = ROOT / "tools" if (ROOT / "tools" / "check_eval_kit.py").exists() else (
    ROOT / "app" / "tools"
)
#: nmrflow 包所在目录（源码仓库是根目录，便携包是 app/）——子进程要用它当 PYTHONPATH
APP_ROOT = TOOLS.parent


class EvalKitTests(unittest.TestCase):
    """随包的评估套件必须"能用"，而不只是"文件在"。

    三条守的是三件不同的事：脚本能 import（说明 app/ 的路径假设没写错）、
    自检能过（说明索引没被截断或换掉）、自检**真的会报错**（说明它不是在走过场）。
    """

    def test_eval_scripts_import_under_package_layout(self):
        """``--help`` 走到 argparse 就返回，所以这一步只验证 import 链。"""
        for name in ("eval_tier_ab.py", "check_eval_kit.py"):
            proc = subprocess.run(
                [sys.executable, str(TOOLS / name), "--help"],
                capture_output=True, text=True, timeout=180,
            )
            self.assertEqual(proc.returncode, 0, f"{name}\n{proc.stderr[-1500:]}")

    def test_check_eval_kit_passes_quick(self):
        proc = subprocess.run(
            [sys.executable, str(TOOLS / "check_eval_kit.py"), "--quick", "--json"],
            capture_output=True, text=True, timeout=300,
        )
        self.assertEqual(proc.returncode, 0, proc.stdout[-2000:] + proc.stderr[-800:])
        payload = json.loads(proc.stdout)
        self.assertTrue(payload["ok"])
        failed = [c for c in payload["checks"] if c["status"] != "PASS"]
        self.assertEqual(failed, [])

    def test_check_eval_kit_survives_gbk_console(self):
        """默认中文控制台是 GBK，``¹``(U+00B9) 编不出来会让自检直接崩。

        0.2.8 上就是这个现象：用户双击 自检.bat 看到 `UnicodeEncodeError: 'gbk'`。
        所以核名改 ASCII，输出流再加 backslashreplace 兜底。
        """
        env = dict(os.environ, PYTHONIOENCODING="gbk", PYTHONUTF8="0")
        # 子进程按 GBK 输出，父进程也必须按 GBK 解，否则先死在解码上
        proc = subprocess.run(
            [sys.executable, str(TOOLS / "check_eval_kit.py"), "--quick"],
            capture_output=True, text=True, timeout=300, env=env,
            encoding="gbk", errors="backslashreplace",
        )
        self.assertEqual(proc.returncode, 0, (proc.stdout or "") + (proc.stderr or ""))
        self.assertNotIn("UnicodeEncodeError", proc.stderr)

    def test_check_eval_kit_catches_wrong_train_index(self):
        """把训练索引指向全量索引，自检必须失败——否则它挡不住"背答案"。"""
        code = (
            "import sys, json;"
            f"sys.path.insert(0, {str(TOOLS)!r});"
            "import check_eval_kit as m;"
            "m.TRAIN_INDEX = m.FULL_INDEX;"
            "ok, res = m.check(quick=True);"
            "print(json.dumps({'ok': ok, 'fail': [r['check'] for r in res"
            " if r['status'] != 'PASS']}, ensure_ascii=False))"
        )
        proc = subprocess.run([sys.executable, "-c", code],
                              capture_output=True, text=True, timeout=300)
        self.assertEqual(proc.returncode, 0, proc.stderr[-1500:])
        payload = json.loads(proc.stdout.strip().splitlines()[-1])
        self.assertFalse(payload["ok"], payload)
        self.assertTrue(payload["fail"], payload)


class AutostopTests(unittest.TestCase):
    """「关掉页面 3 分钟后自动退出」——默认关，打开后靠心跳超时退出。

    这条功能的风险是"误杀"：刷新、切标签、休眠都可能让心跳短暂中断。
    所以既要验证它**真的会停**，也要确认**默认不停**、关掉之后不会又自己停。
    """

    def test_default_is_off(self):
        from nmrflow import autostop

        self.assertFalse(autostop.is_enabled())

    def test_touch_advances_heartbeat(self):
        from nmrflow import autostop

        before = autostop.latest_heartbeat()
        autostop.touch("unit-test-session")
        self.assertGreater(autostop.latest_heartbeat(), before)

    def test_disable_keeps_process_alive(self):
        """关掉之后，即使心跳很旧也不许退出。"""
        from nmrflow import autostop

        autostop.set_enabled(True)
        autostop.set_enabled(False)
        self.assertFalse(autostop.is_enabled())

    def test_exits_after_grace_without_heartbeat(self):
        code = (
            "import time;"
            "from nmrflow import autostop;"
            "autostop.set_enabled(True);"
            "time.sleep(30);"
            "print('still alive')"
        )
        env = dict(
            os.environ,
            NMRFLOW_AUTOSTOP_GRACE="2",
            NMRFLOW_AUTOSTOP_CHECK="0.3",
            PYTHONPATH=str(APP_ROOT),
            PYTHONIOENCODING="utf-8",
        )
        started = time.time()
        proc = subprocess.run([sys.executable, "-c", code], capture_output=True,
                              text=True, timeout=60, env=env)
        elapsed = time.time() - started
        self.assertEqual(proc.returncode, 0, proc.stderr[-800:])
        self.assertNotIn("still alive", proc.stdout, "心跳超时后应当直接退出")
        self.assertLess(elapsed, 20, elapsed)


class ReferenceRecordTests(unittest.TestCase):
    """``check_reference_records`` 必须能分出"对"和"错"的记录。

    这工具存在的意义就是别再用错记录改模型（0.2.3 / 0.2.6 / 0.2.8 三次都是同一个坑），
    所以这里既钉"对的记录 PASS"，也钉"错的记录 FAIL"。
    """

    DECALIN = "C1CCC2CCCCC2C1"
    NORBORNANE = "C1CC2CCC1C2"
    #: 工具要求十氢萘与降冰片烷都能对上，所以每份 SD 都带一条正确的降冰片烷，
    #: 这样"整体退出码"只反映十氢萘那条记录的对错。
    NORBORNANE_LIT = [(37.3, 2), (30.6, 4), (39.2, 1)]

    def _write_sd(self, decalin_pattern, norbornane_pattern=None):
        """按 ``[(δ, 碳数)...]`` 写一个带 13C 归属的 SD。"""
        from rdkit import Chem

        work = Path(tempfile.mkdtemp())
        path = work / "ref.sdf"
        writer = Chem.SDWriter(str(path))
        entries = [
            (self.DECALIN, decalin_pattern),
            (self.NORBORNANE, norbornane_pattern or self.NORBORNANE_LIT),
        ]
        for smiles, pattern in entries:
            mol = Chem.MolFromSmiles(smiles)
            flat = []
            for value, count in pattern:
                flat.extend([value] * count)
            carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
            self.assertEqual(len(flat), len(carbons), smiles)
            mol.SetProp("Spectrum 13C 1",
                        "|".join(f"{v:.2f};1;{i}" for v, i in zip(flat, carbons)))
            writer.write(mol)
        writer.close()
        return path

    def _run(self, sd, *extra):
        return subprocess.run(
            [sys.executable, str(TOOLS / "check_reference_records.py"), "--sd", str(sd), *extra],
            capture_output=True, text=True, timeout=180,
        )

    def test_trans_decalin_passes(self):
        proc = self._run(self._write_sd([(44.0, 2), (34.7, 4), (27.1, 4)]))
        self.assertEqual(proc.returncode, 0, proc.stdout[-800:] + proc.stderr[-400:])
        self.assertIn("trans", proc.stdout)

    def test_wrong_decalin_fails(self):
        """0.2.3 那组 34.5/26.9/26.3 不是 trans，工具必须判 FAIL。"""
        proc = self._run(self._write_sd([(34.5, 2), (26.9, 4), (26.3, 4)]))
        self.assertEqual(proc.returncode, 1, proc.stdout[-800:])
        self.assertIn("FAIL", proc.stdout)

    def test_warn_only_never_fails(self):
        proc = self._run(self._write_sd([(34.5, 2), (26.9, 4), (26.3, 4)]), "--warn-only")
        self.assertEqual(proc.returncode, 0, proc.stdout[-500:])


class BatFileTests(unittest.TestCase):
    """bat 的字节级约束——这类问题在 Python 里完全看不出来。

    0.2.8 发出去的 5 个 bat 是 LF-only，双击卸载向导直接把行拆成碎片去执行
    （`'cho'`、`'HOICE"=="3" goto cancel'`）；且含中文的 bat 若存成 UTF-8 并带
    `chcp 65001`，cmd 仍按 CP936 读前几行，汉字被从中间切断。
    """

    def test_bats_are_crlf(self):
        bats = sorted(ROOT.glob("*.bat"))
        self.assertTrue(bats, "包里应当有 bat")
        for path in bats:
            data = path.read_bytes()
            lone_lf = data.count(b"\n") - data.count(b"\r\n")
            self.assertEqual(lone_lf, 0, f"{path.name} 有 {lone_lf} 处纯 LF")

    def test_chinese_bat_is_gbk_without_chcp(self):
        for path in sorted(ROOT.glob("*.bat")):
            data = path.read_bytes()
            if data.isascii():
                continue
            with self.assertRaises(UnicodeDecodeError, msg=f"{path.name} 是 UTF-8，应存 GBK"):
                data.decode("utf-8")
            self.assertNotIn(b"chcp 65001", data.lower(), f"{path.name} 含中文又带 chcp 65001")
            data.decode("gbk")  # 必须能按 GBK 解出


if __name__ == "__main__":
    unittest.main(verbosity=2)


