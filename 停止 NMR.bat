@echo off
chcp 65001 >nul
setlocal EnableExtensions

set "PORT=8901"
set "QUIET="

:parse
if "%~1"=="" goto parsed
if /i "%~1"=="/quiet" (
    set "QUIET=1"
    shift
    goto parse
)
set "PORT=%~1"
shift
goto parse
:parsed

set "FOUND="
set "FAILED="
echo Stopping nmrflow (port %PORT%)...
echo Note: this ends whatever process is LISTENING on that port, not only nmrflow.
echo.

for /f "tokens=5" %%a in ('netstat -ano ^| findstr /C:":%PORT% " ^| findstr /C:"LISTENING"') do (
    set "FOUND=1"
    echo Stopping PID %%a...
    taskkill /PID %%a /T /F >nul 2>&1
    if errorlevel 1 (
        echo [ERROR] taskkill failed for PID %%a.
        set "FAILED=1"
    ) else (
        echo Stopped PID %%a.
    )
)

if not defined FOUND (
    echo No listener on port %PORT%.
    set "RC=0"
    goto done
)

ping -n 2 127.0.0.1 >nul
set "STILL="
for /f "tokens=5" %%a in ('netstat -ano ^| findstr /C:":%PORT% " ^| findstr /C:"LISTENING"') do set "STILL=1"
if defined STILL (
    echo [ERROR] port %PORT% is still LISTENING.
    set "RC=1"
) else (
    echo Port %PORT% is free.
    if defined FAILED (
        echo [WARN] at least one taskkill failed.
        set "RC=1"
    ) else (
        set "RC=0"
    )
)

:done
echo.
if defined QUIET (
    endlocal & exit /b %RC%
)
pause
endlocal & exit /b %RC%
