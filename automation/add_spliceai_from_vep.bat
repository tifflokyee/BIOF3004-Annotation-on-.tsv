@echo off
setlocal

rem Windows launcher for adding real SpliceAI scores from VEP output.
rem No downloads, no administrator rights, and no system changes are performed.

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
set "PY_SCRIPT=%SCRIPT_DIR%add_spliceai_from_vep.py"

if not exist "%PY_SCRIPT%" (
    echo ERROR: Cannot find "%PY_SCRIPT%".
    exit /b 1
)

where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python was not found on PATH.
    exit /b 1
)

pushd "%PROJECT_ROOT%" >nul
python "%PY_SCRIPT%" %*
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%
