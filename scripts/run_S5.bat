@echo off
REM ============================================================================
REM run_S5.bat -- stages S4_05 + S5 runner (Amendment A4, 2026-09-22).
REM Each step runs as its OWN process and writes its own log under logs\
REM (same v1.1 logging discipline as run_S4.bat -- no silent stops, R6).
REM Prerequisite: metafor installed once (same as S4):
REM   Rscript -e "utils::install.packages('metafor')"
REM If a step fails: open its log under logs\, fix the cause, re-run this file;
REM completed steps re-verify their inputs (registry) and are fast.
REM ============================================================================
setlocal
set R="C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
cd /d "%~dp0\.."
if not exist logs mkdir logs

echo [S5] 0/4 smoke test (8 legs, known-answer math) --^> logs\S5_smoke.log
%R% scripts\S5_smoke.R > logs\S5_smoke.log 2>&1
if errorlevel 1 type logs\S5_smoke.log & goto :fail

echo [S5] 1/4 axis localization map (S4_05, pure measurement) --^> logs\S4_05.log
%R% scripts\S4_05_axis_localization.R > logs\S4_05.log 2>&1
if errorlevel 1 type logs\S4_05.log & goto :fail

echo [S5] 2/4 pathway scores (z-mean axes + locked SenMayo reuse) --^> logs\S5_01.log
%R% scripts\S5_01_pathway_scores.R > logs\S5_01.log 2>&1
if errorlevel 1 type logs\S5_01.log & goto :fail

echo [S5] 3/4 pathway REML meta + GSE11691 arm + sensitivity --^> logs\S5_02.log
%R% scripts\S5_02_pathway_meta.R > logs\S5_02.log 2>&1
if errorlevel 1 type logs\S5_02.log & goto :fail

echo [S5] ALL DONE -- see results\S5_meta\S5_sensitivity_summary.md
echo          and results\S4_localization\S4_05_axis_summary.md
exit /b 0

:fail
echo [S5] STEP FAILED -- full log printed above and saved under logs\.
echo        Fix the cause, then re-run this file (done steps re-verify fast).
pause
exit /b 1
