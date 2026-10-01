@echo off
REM ============================================================================
REM run_S3_01_repair.bat -- Amendment A3: regenerate the 19 locked bulk
REM contrast CSVs at SYMBOL level (S3_01 v2.0 + S3_common v3.0.2; ledger #18).
REM Requires the lead's written A3 approval. Logs under logs\.
REM Overwrites results\S3_bulk_contrasts\*.csv -- script-generated outputs,
REM never hand-edited (R4). Re-running is idempotent under the same inputs.
REM ============================================================================
setlocal
set R="C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
cd /d "%~dp0\.."
if not exist logs mkdir logs

for %%C in (GSE7305 GSE6364 GSE25628 GSE7307 GSE51981 GSE11691 GSE120103) do (
  echo [S3_01 repair] %%C --^> logs\S3_01_%%C.log
  %R% scripts\S3_01_bulk_contrasts.R %%C > logs\S3_01_%%C.log 2>&1
  if errorlevel 1 type logs\S3_01_%%C.log & goto :fail
)

echo [S3_01 repair] collect --^> logs\S3_01_collect.log
%R% scripts\S3_01_bulk_contrasts.R collect > logs\S3_01_collect.log 2>&1
if errorlevel 1 type logs\S3_01_collect.log & goto :fail

echo [S3_01 repair] S3_06 collect refresh --^> logs\S3_06_collect.log
%R% scripts\S3_06_collect.R > logs\S3_06_collect.log 2>&1
if errorlevel 1 type logs\S3_06_collect.log & goto :fail

echo [S3_01 repair] DONE. Now delete scripts\S4_input_registry.tsv (A3) and run run_S4.bat
exit /b 0

:fail
echo [S3_01 repair] STEP FAILED -- log printed above and saved under logs\.
pause
exit /b 1
