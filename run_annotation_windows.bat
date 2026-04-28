@echo off
setlocal

echo =====================================================
echo   BIOF Annotation Pipeline (Windows)
echo =====================================================
echo.

:: Run from repository root (location of this .bat file)
set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"

:: Activate Conda environment (best in Anaconda Prompt)
call conda activate biof_annotation
if errorlevel 1 (
    if defined CONDA_EXE (
        set "CONDA_BAT=%CONDA_EXE:conda.exe=conda.bat%"
        if exist "%CONDA_BAT%" call "%CONDA_BAT%" activate biof_annotation
    )
)
if errorlevel 1 (
    echo ERROR: Could not activate conda environment "biof_annotation".
    echo Please open Anaconda Prompt, or run: conda init cmd.exe
    pause
    exit /b 1
)

echo Conda environment activated: biof_annotation
echo.

echo Running: python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample" --skip-gnomad
python3 automation/auto_annotate_generated.py --all-tsv --input-dir "sample" --skip-gnomad

if errorlevel 1 (
    echo.
    echo Pipeline failed.
    pause
    exit /b 1
)

echo.
echo Pipeline finished successfully!
pause
exit /b 0
