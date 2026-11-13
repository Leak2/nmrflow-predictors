@echo off
chcp 65001 >nul
setlocal
set "PYTHONUTF8=1"
cd /d "%~dp0"

echo Running nmrflow eval-kit self-check...
echo.

"%~dp0runtime\python\python.exe" "%~dp0app\tools\check_eval_kit.py"
set "RC=%ERRORLEVEL%"

echo.
if "%RC%"=="0" (
    echo [OK] self-check passed.
) else (
    echo [FAIL] self-check failed, exit code %RC%.
)
echo.
pause
exit /b %RC%
