@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist "NMR-Predictor.exe" (
    echo [ERROR] NMR-Predictor.exe not found next to this script.
    echo.
    pause
    exit /b 1
)
"NMR-Predictor.exe" --install
