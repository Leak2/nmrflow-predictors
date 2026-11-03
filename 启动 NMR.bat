@echo off
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
