"""Packaging toolchain checks (P1): bat generator, safety gates, packaging invariants."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
TOOLS = ROOT / "app" / "tools"
sys.path.insert(0, str(TOOLS))

import make_release as mr  # noqa: E402


def test_version_matches_package():
    import nmrflow

    assert mr.VERSION == nmrflow.__version__


def test_bat_bytes_forces_crlf_and_gbk():
    for name, text in (
        ("停止 NMR.bat", mr.STOP_BAT),
        ("卸载.bat", mr.UNINSTALL_BAT),
    ):
        data = mr.bat_bytes(text)
        assert data.count(b"\n") == data.count(b"\r\n"), name
        if not data.isascii():
            data.decode("gbk")
            with pytest.raises(UnicodeDecodeError):
                data.decode("utf-8")
            assert b"chcp 65001" not in data.lower()


def test_bat_bytes_rejects_chcp_with_chinese():
    with pytest.raises(SystemExit):
        mr.bat_bytes("@echo off\nchcp 65001 >nul\necho 中文\n")


def test_check_bat_files_flags_bad_bytes(tmp_path):
    (tmp_path / "good.bat").write_bytes(mr.bat_bytes("@echo off\necho ok\n"))
    (tmp_path / "bad.bat").write_bytes(b"@echo off\necho ok\n")
    problems = mr.check_bat_files(tmp_path)
    assert any("bad.bat" in p and "LF" in p for p in problems)
    assert not any("good.bat" in p for p in problems)


def _run_bat(path: Path, args: list[str], tmp_path: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["cmd", "/c", str(path), *args],
        input=b"\n",
        cwd=tmp_path,
        capture_output=True,
        timeout=30,
    )


def test_uninstall_refuses_non_install_dir(tmp_path):
    bat = tmp_path / "卸载.bat"
    bat.write_bytes(mr.bat_bytes(mr.UNINSTALL_BAT))
    proc = _run_bat(bat, ["2", "DELETE", "N"], tmp_path)
    assert proc.returncode != 0
    assert tmp_path.exists()
    assert "ERROR" in proc.stdout.decode("gbk", errors="replace")


def test_uninstall_refuses_drive_root():
    root = Path("C:/")
    # Do not execute against a real drive; assert the template's guard is present.
    assert "refuse_root" in mr.UNINSTALL_BAT
    assert 'if "%ROOT:~1,1%"==":"' in mr.UNINSTALL_BAT


def test_stop_bat_quiet_no_listener(tmp_path):
    bat = tmp_path / "停止 NMR.bat"
    bat.write_bytes(mr.bat_bytes(mr.STOP_BAT))
    proc = _run_bat(bat, ["/quiet", "65534"], tmp_path)
    assert proc.returncode == 0, proc.stdout[-500:]

def test_check_layout_reports_missing_files(tmp_path):
    missing = mr.check_layout(tmp_path)
    assert "app/nmrflow/app.py" in missing
    assert "NMR-Predictor.exe" in missing