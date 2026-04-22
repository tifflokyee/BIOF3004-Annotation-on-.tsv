@echo off
setlocal

rem Windows launcher for local SpliceAI model scoring with hg19 sequence windows from UCSC.

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
set "PY_SCRIPT=%SCRIPT_DIR%run_local_spliceai_ucsc.py"

if not exist "%PY_SCRIPT%" (
    echo ERROR: Cannot find "%PY_SCRIPT%".
    exit /b 1
)

if "%OS%"=="Windows_NT" (
    ver | find "5.1." >nul
    if not errorlevel 1 (
        echo ERROR: Local SpliceAI cannot run on Windows XP.
        echo It requires a modern Python/TensorFlow environment. Run this step on Windows 10/11 or Linux.
        exit /b 1
    )
)

if "%SPLICEAI_PYTHON%"=="" set "SPLICEAI_PYTHON=py -3.13"

pushd "%PROJECT_ROOT%" >nul
%SPLICEAI_PYTHON% "%PY_SCRIPT%" %*
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%
