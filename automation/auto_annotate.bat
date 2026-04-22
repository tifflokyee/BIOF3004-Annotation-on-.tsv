@echo off
setlocal

rem Native Windows launcher for the local annotation pipeline.
rem It does not download files, request administrator rights, or modify system settings.

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
set "PY_SCRIPT=%SCRIPT_DIR%auto_annotate_generated.py"
set SPLICEAI_PYTHON=python

if not exist "%PY_SCRIPT%" (
    echo ERROR: Cannot find "%PY_SCRIPT%".
    echo Keep auto_annotate.bat and auto_annotate_generated.py in the same folder.
    exit /b 1
)

where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python was not found on PATH.
    echo Install Python or run this from a terminal where python is available.
    exit /b 1
)

pushd "%PROJECT_ROOT%" >nul
python "%PY_SCRIPT%" %*
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%