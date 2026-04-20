@echo off
setlocal

rem Windows launcher for creating the small VCF to upload to VEP + SpliceAI.

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
set "PY_SCRIPT=%SCRIPT_DIR%create_small_vcf_for_vep.py"

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
