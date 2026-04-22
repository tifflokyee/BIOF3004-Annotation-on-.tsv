@echo off
setlocal

rem Stepwise Windows annotation runner.
rem Usage:
rem   run_stepwise_annotation.bat input.tsv

set "INPUT=%~1"
set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%.."

if "%INPUT%"=="" (
    echo Usage: run_stepwise_annotation.bat input.tsv
    exit /b 1
)

python --version >nul 2>nul
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

echo.
echo [6/6] Local SpliceAI
if "%OS%"=="Windows_NT" (
    ver | find "5.1." >nul
    if not errorlevel 1 goto spliceai_xp
)
call "%SCRIPT_DIR%run_local_spliceai_ucsc.bat" "result\step5_gnomad.tsv" -o "result\step6_final_with_local_spliceai.tsv" --spliceai-vcf "result\local_spliceai_model_output.vcf"
if errorlevel 1 goto fail
goto spliceai_done

:spliceai_xp
if not exist "result\local_spliceai_model_output.vcf" (
    echo ERROR: Windows XP cannot run local SpliceAI model inference.
    echo Expected existing local SpliceAI VCF: result\local_spliceai_model_output.vcf
    echo Creating VCF to run through local SpliceAI on another computer...
    python "%SCRIPT_DIR%create_small_vcf_for_vep.py" "result\step5_gnomad.tsv" -o "result\local_spliceai_input_for_external_run.vcf"
    echo.
    echo Next:
    echo   1. Copy result\local_spliceai_input_for_external_run.vcf to Windows 10/11 or Linux.
    echo   2. Run local SpliceAI there and save the output as local_spliceai_model_output.vcf.
    echo   3. Copy local_spliceai_model_output.vcf back into this result folder.
    echo   4. Rerun this automation on Windows XP.
    popd >nul
    exit /b 1
)
echo Windows XP detected; merging existing local SpliceAI VCF.
python "%SCRIPT_DIR%annotate_spliceai_step.py" "result\step5_gnomad.tsv" -o "result\step6_final_with_local_spliceai.tsv" --local-spliceai-vcf "result\local_spliceai_model_output.vcf"
if errorlevel 1 goto fail

:spliceai_done
echo.
echo Done. Final output: result\step6_final_with_local_spliceai.tsv
popd >nul
exit /b 0

:fail
set "EXIT_CODE=%ERRORLEVEL%"
popd >nul
exit /b %EXIT_CODE%
