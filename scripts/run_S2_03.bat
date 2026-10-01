@echo off
rem ==========================================================================
rem run_S2_03.bat -- chained short-process runner for S2_03 bulk QC (B.6/B.7)
rem Self-locating; checks exit code per cohort; prints the failing step's log.
rem Run from anywhere:  scripts\run_S2_03.bat
rem ==========================================================================
setlocal EnableExtensions
set "RSCRIPT=C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
cd /d "%~dp0\.."

if not exist "%RSCRIPT%" ( echo FATAL: Rscript not found at %RSCRIPT% & exit /b 1 )
if not exist scripts\S2_03_bulk_qc.R ( echo FATAL: scripts\S2_03_bulk_qc.R missing & exit /b 1 )

for %%C in (GSE7305 GSE6364 GSE25628 GSE7307 GSE51981 GSE11691 GSE120103 collect) do (
  echo === S2_03 step: %%C ===
  "%RSCRIPT%" scripts\S2_03_bulk_qc.R %%C > logs\S2_03_%%C.txt 2>&1
  if errorlevel 1 (
    echo STEP FAILED: %%C -- full log follows:
    type logs\S2_03_%%C.txt
    exit /b 1
  )
)
echo ALL S2_03 STEPS DONE -- upload logs\S2_03_*.txt and results\S2_bulk_outliers.csv
endlocal
