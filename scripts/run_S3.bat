@echo off
rem ==========================================================================
rem run_S3.bat -- chained short-process runner for stage S3 (B.6/B.7)
rem ORDER (dependencies): smoke -> preflight -> [FILL CONFIG from S3_00] ->
rem   bulk contrasts -> ssGSEA -> pseudobulk -> wilcoxon(sens) -> z-mean -> collect
rem The bat runs straight through; if S3_00 prints config you must fill,
 rem scripts stop FATAL at the require_config gate (by design, R3/R5).
rem Run from anywhere:  scripts\run_S3.bat
rem ==========================================================================
setlocal EnableExtensions
set "RSCRIPT=C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
cd /d "%~dp0\.."

if not exist "%RSCRIPT%" ( echo FATAL: Rscript not found at %RSCRIPT% & exit /b 1 )
if not exist scripts\S3_common.R ( echo FATAL: scripts\S3_common.R missing & exit /b 1 )
if not exist logs ( mkdir logs )

rem ---- 0. mandatory smoke (B.4) ------------------------------------------
echo === S3 step: smoke ===
"%RSCRIPT%" scripts\S3_smoke.R > logs\S3_smoke.txt 2>&1
if errorlevel 1 ( echo STEP FAILED: smoke -- log follows: & type logs\S3_smoke.txt & exit /b 1 )

rem ---- 1. read-only census -------------------------------------------------
echo === S3 step: preflight ===
"%RSCRIPT%" scripts\S3_00_preflight.R > logs\S3_00_preflight.txt 2>&1
if errorlevel 1 ( echo STEP FAILED: preflight & exit /b 1 )
echo PREFLIGHT DONE -- fill scripts\S3_common.R config block from logs\S3_00_preflight.txt, then re-run this bat.

rem ---- 2. bulk contrasts (one process per cohort) --------------------------
for %%C in (GSE7305 GSE6364 GSE25628 GSE7307 GSE51981 GSE11691 GSE120103 collect) do (
  echo === S3 step: bulk contrasts %%C ===
  "%RSCRIPT%" scripts\S3_01_bulk_contrasts.R %%C > logs\S3_01_%%C.txt 2>&1
  if errorlevel 1 ( echo STEP FAILED: bulk %%C -- log follows: & type logs\S3_01_%%C.txt & exit /b 1 )
)

rem ---- 3-6. scores + scRNA + collect ---------------------------------------
for %%S in (S3_02_ssgsea_bulk S3_03_scrna_pseudobulk S3_04_scrna_wilcoxon_sens S3_05_scrna_zmean S3_06_collect) do (
  echo === S3 step: %%S ===
  "%RSCRIPT%" scripts\%%S.R > logs\%%S.txt 2>&1
  if errorlevel 1 ( echo STEP FAILED: %%S -- log follows: & type logs\%%S.txt & exit /b 1 )
)
echo ALL S3 STEPS DONE -- upload logs\S3_*.txt, results\S3_bulk_contrasts_summary.csv, results\S3_knock_summary.md
endlocal
