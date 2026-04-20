@echo off
setlocal

rem Stepwise Windows annotation runner.
rem Usage:
rem   run_stepwise_annotation.bat input.tsv [vep_output.vcf]
rem
rem If the VEP SpliceAI file is omitted or missing, this script stops after
rem creating the TSV that should be used to build the small VCF for VEP.

set "INPUT=%~1"
set "VEP_SPLICEAI=%~2"
set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."

if "%INPUT%"=="" (
    echo Usage: run_stepwise_annotation.bat input.tsv [vep_output.vcf]
    exit /b 1
)

where python >nul 2>nul
if errorlevel 1 (
    echo ERROR: Python was not found on PATH.
    exit /b 1
)

if not exist "result" mkdir "result"
pushd "%PROJECT_ROOT%" >nul

echo.
echo [1/6] PanelApp
python "%SCRIPT_DIR%annotate_panelapp_step.py" "%INPUT%" -o "result\step1_panelapp.tsv"
if errorlevel 1 goto fail

echo.
echo [2/6] REVEL
python "%SCRIPT_DIR%annotate_revel_step.py" "result\step1_panelapp.tsv" -o "result\step2_revel.tsv"
if errorlevel 1 goto fail

echo.
echo [3/6] AlphaMissense
python "%SCRIPT_DIR%annotate_alphamissense_step.py" "result\step2_revel.tsv" -o "result\step3_alphamissense.tsv"
if errorlevel 1 goto fail

echo.
echo [4/6] ClinVar
python "%SCRIPT_DIR%annotate_clinvar_step.py" "result\step3_alphamissense.tsv" -o "result\step4_clinvar.tsv"
if errorlevel 1 goto fail

echo.
echo [5/6] gnomAD
python "%SCRIPT_DIR%annotate_gnomad_step.py" "result\step4_clinvar.tsv" -o "result\step5_gnomad.tsv"
if errorlevel 1 goto fail

if "%VEP_SPLICEAI%"=="" goto no_vep
if not exist "%VEP_SPLICEAI%" goto no_vep

echo.
echo [6/6] SpliceAI from VEP
python "%SCRIPT_DIR%annotate_spliceai_step.py" "result\step5_gnomad.tsv" -o "result\step6_final_with_spliceai.tsv" --vep-spliceai "%VEP_SPLICEAI%"
if errorlevel 1 goto fail
echo.
echo Done. Final output: result\step6_final_with_spliceai.tsv
popd >nul
exit /b 0

:no_vep
echo.
echo SpliceAI VEP output was not provided or was not found.
echo Current output before SpliceAI: result\step5_gnomad.tsv
echo Next:
echo   1. Run: automation\create_small_vcf_for_vep.bat result\step5_gnomad.tsv -o small_clinical_variants_for_vep.vcf
echo   2. Upload small_clinical_variants_for_vep.vcf to GRCh37 Ensembl VEP with SpliceAI enabled.
echo   3. Download the VEP result.
echo   4. Run: python automation\annotate_spliceai_step.py result\step5_gnomad.tsv -o result\step6_final_with_spliceai.tsv --vep-spliceai downloaded_vep_file.vcf
popd >nul
exit /b 0

:fail
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%
