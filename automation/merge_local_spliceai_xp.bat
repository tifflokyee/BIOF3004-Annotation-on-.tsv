@echo off
setlocal

rem Windows XP-friendly launcher for merging an existing local SpliceAI VCF into a TSV.
rem Usage:
rem   merge_local_spliceai_xp.bat input.tsv spliceai_output.vcf output.tsv

set "INPUT_TSV=%~1"
set "SPLICEAI_VCF=%~2"
set "OUTPUT_TSV=%~3"
set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."

if "%INPUT_TSV%"=="" (
    echo Usage: merge_local_spliceai_xp.bat input.tsv spliceai_output.vcf output.tsv
    exit /b 1
)

if "%SPLICEAI_VCF%"=="" (
    echo Usage: merge_local_spliceai_xp.bat input.tsv spliceai_output.vcf output.tsv
    exit /b 1
)

if "%OUTPUT_TSV%"=="" (
    echo Usage: merge_local_spliceai_xp.bat input.tsv spliceai_output.vcf output.tsv
    exit /b 1
)

python --version >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python was not found on PATH.
    exit /b 1
)

pushd "%PROJECT_ROOT%" >nul
python "%SCRIPT_DIR%annotate_spliceai_step.py" "%INPUT_TSV%" -o "%OUTPUT_TSV%" --local-spliceai-vcf "%SPLICEAI_VCF%"
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%
