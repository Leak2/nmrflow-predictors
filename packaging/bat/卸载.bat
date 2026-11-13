@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "ROOT=%~dp0"
if "%ROOT:~-1%"=="\" set "ROOT=%ROOT:~0,-1%"
set "EXE=%ROOT%\NMR-Predictor.exe"
set "APPDIR=%ROOT%\app\nmrflow"
set "PORT={{PORT}}"

echo ============================================
echo  nmrflow 卸载向导
echo ============================================
echo  程序目录:
echo    %ROOT%
echo.

rem ---- 安全闸：必须确认这正是 nmrflow 的安装目录 ----
rem 少了这一步，用户把 卸载.bat 单独拷到别的文件夹（或把 zip 内容直接解压到
rem D:\NMR\ 而没带顶层目录）再双击选 2，会 rmdir 掉一个与程序无关的目录。
if not exist "%EXE%" goto not_install
if not exist "%APPDIR%\app.py" goto not_install
if "%ROOT:~1,1%"==":" if "%ROOT:~2,1%"=="" goto refuse_root
if "%ROOT:~2,1%"==":" if "%ROOT:~3,1%"=="" goto refuse_root
goto menu

:not_install
echo [ERROR] 这个目录看起来不是 nmrflow 的安装目录。
echo         缺少 NMR-Predictor.exe 或 app\nmrflow\app.py。
echo         本向导只会删除它自己所在的程序目录，找不到就什么都不做。
echo.
pause
endlocal
exit /b 1

:refuse_root
echo [ERROR] 程序目录位于磁盘根目录（%ROOT%\），拒绝执行删除。
echo         请先把整个 nmrflow 文件夹放到一个子目录里。
echo.
pause
endlocal
exit /b 1

:menu
echo   [1] 只删除桌面快捷方式，保留程序与 results
echo   [2] 删除快捷方式 + 整个程序目录（可选保留 results）
echo   [3] 取消
echo.
set "CHOICE="
if not "%~1"=="" (
    set "CHOICE=%~1"
) else (
    set /p "CHOICE=请输入 1/2/3: "
)
set "CHOICE=%CHOICE: =%"
set "CHOICE=%CHOICE:~0,1%"

if "%CHOICE%"=="1" goto shortcut
if "%CHOICE%"=="2" goto full
if "%CHOICE%"=="3" goto cancel
echo.
echo [ERROR] 无效选择，什么也没做。
goto end

:shortcut
call :remove_shortcut
echo.
echo 已删除桌面快捷方式；程序目录保留：
echo   %ROOT%
echo 需要彻底删除时，再次运行本向导并选 2。
goto end

:full
echo.
echo 将删除整个程序目录：
echo   %ROOT%
echo.
set "CONFIRM="
if not "%~2"=="" (
    set "CONFIRM=%~2"
) else (
    set /p "CONFIRM=确认删除请输入 DELETE（其它任意输入取消）: "
)
if /i not "%CONFIRM%"=="DELETE" (
    echo.
    echo 已取消，什么也没做。
    goto end
)
set "KEEP="
if not "%~3"=="" (
    set "KEEP=%~3"
) else (
    set /p "KEEP=删除前把 results 文件夹移到桌面？[Y/N]: "
)
set "KEEP=%KEEP: =%"
call :stop_service
if errorlevel 1 goto end
call :remove_shortcut
if /i "%KEEP:~0,1%"=="Y" call :keep_results
echo.
echo 正在安排删除程序目录...
call :schedule_delete
echo 已安排删除；窗口关闭后约 2 秒完成。
echo 若只剩一个被占用的空目录，请手动删除它。
goto end

:remove_shortcut
if exist "%EXE%" (
    echo [1/3] 删除桌面快捷方式...
    "%EXE%" --uninstall
) else (
    echo [WARN] 未找到 NMR-Predictor.exe，跳过快捷方式删除。
)
exit /b 0

:stop_service
echo [2/3] 检查并停止 %PORT% 服务...
set "FOUND="
set "FAILED="
for /f "tokens=5" %%a in ('netstat -ano ^| findstr /C:":%PORT% " ^| findstr /C:"LISTENING"') do (
    set "FOUND=1"
    echo   停止 PID %%a
    taskkill /PID %%a /T /F >nul 2>&1
    if errorlevel 1 set "FAILED=1"
)
if not defined FOUND (
    echo   %PORT% 未在监听，跳过。
    exit /b 0
)
ping -n 2 127.0.0.1 >nul
set "STILL="
for /f "tokens=5" %%a in ('netstat -ano ^| findstr /C:":%PORT% " ^| findstr /C:"LISTENING"') do set "STILL=1"
if defined STILL (
    echo [ERROR] %PORT% 仍在监听，已中止删除。
    echo         请先手动停止服务，再重新运行卸载。
    pause
    exit /b 1
)
if defined FAILED echo [WARN] 至少一个 PID 结束失败；端口已释放，继续。
exit /b 0
:keep_results
if not exist "%ROOT%\results" (
    echo   没有 results 文件夹，跳过。
    exit /b 0
)
set "DEST=%USERPROFILE%\Desktop\nmrflow-results-backup"
if exist "%DEST%" set "DEST=%DEST%-%RANDOM%"
if not exist "%USERPROFILE%\Desktop" (
    set "DEST=%USERPROFILE%\nmrflow-results-backup-%RANDOM%"
)
echo   移动 results 到:
echo     %DEST%
move "%ROOT%\results" "%DEST%" >nul 2>nul
if errorlevel 1 (
    echo   [WARN] 移动失败，results 会随目录一起删除。
) else (
    echo   已保留 results。
)
exit /b 0

:schedule_delete
set "CLEAN=%TEMP%\nmrflow-uninstall-%RANDOM%-%RANDOM%.bat"
> "%CLEAN%" echo @echo off
>>"%CLEAN%" echo cd /d C:\
>>"%CLEAN%" echo ping -n 3 127.0.0.1 ^>nul
>>"%CLEAN%" echo for /l %%%%i in (1,1,10) do (
>>"%CLEAN%" echo     rmdir /s /q "%ROOT%" ^>nul 2^>nul
>>"%CLEAN%" echo     if not exist "%ROOT%" goto done
>>"%CLEAN%" echo     ping -n 2 127.0.0.1 ^>nul
>>"%CLEAN%" echo )
>>"%CLEAN%" echo :done
>>"%CLEAN%" echo del "%%~f0"
cd /d C:\
start "" /min /d "C:\" cmd /c ""%CLEAN%""
exit /b 0

:cancel
echo.
echo 已取消，什么也没做。
goto end

:end
echo.
pause
endlocal
exit /b 0
