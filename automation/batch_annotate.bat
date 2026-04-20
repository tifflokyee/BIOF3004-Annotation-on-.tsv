@echo off
REM Batch annotation launcher for Windows
REM This script provides an easy way to run batch annotation from the command line
REM
REM Usage:
REM   batch_annotate.bat
REM   batch_annotate.bat --sample-id sample_001
REM   batch_annotate.bat --list
REM   batch_annotate.bat --reset

setlocal enabledelayedexpansion

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."
set "PYTHON_SCRIPT=%SCRIPT_DIR%batch_annotate.py"
set "SAMPLE_LIST=%PROJECT_ROOT%\sample_list.tsv"

if not exist "%PYTHON_SCRIPT%" (
    echo ERROR: Cannot find "%PYTHON_SCRIPT%".
    echo This script must be in the automation folder.
    exit /b 1
)

if not exist "%SAMPLE_LIST%" (
    echo ERROR: Cannot find "%SAMPLE_LIST%".
    echo Please create sample_list.tsv in the project root folder.
    exit /b 1
)

where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python was not found on PATH.
    echo Install Python or run this from a terminal where python is available.
    exit /b 1
)

pushd "%PROJECT_ROOT%" >nul
python "%PYTHON_SCRIPT%" "%SAMPLE_LIST%" %*
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%
