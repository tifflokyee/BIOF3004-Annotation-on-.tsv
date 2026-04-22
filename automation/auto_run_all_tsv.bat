@echo off
setlocal enabledelayedexpansion

echo =====================================================
echo   Full Annotation Pipeline for ALL TSV files
echo   (PanelApp → REVEL → AlphaMissense → ClinVar → gnomAD → SpliceAI)
echo =====================================================

:: === CONFIGURATION ===
set INPUT_FOLDER=input
set AUTOMATION=..\automation\run_stepwise_annotation.bat

:: Check if input folder exists
if not exist "%INPUT_FOLDER%" (
    echo ERROR: Folder "%INPUT_FOLDER%" not found!
    echo Please create an "input" folder and put your .tsv files inside it.
    pause
    exit /b 1
)

for %%f in (%INPUT_FOLDER%\*.tsv) do (
    echo.
    echo =====================================================
    echo Processing file: %%~nxf
    echo =====================================================

    :: Run the full stepwise annotation (this will do all 6 steps including gnomAD)
    call %AUTOMATION% "%%f"

    echo.
    echo Finished processing %%~nxf
    echo Output should be in result\step6_final_with_local_spliceai.tsv (or similar)
    echo.
)

echo =====================================================
echo All files have been processed!
echo Check the "result" folder for outputs.
echo =====================================================
pause