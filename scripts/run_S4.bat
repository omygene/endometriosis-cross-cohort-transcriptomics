@echo off
REM ============================================================================
REM run_S4.bat -- stage S4 runner. Each step runs as its OWN process and
REM writes its own log under logs\ (v1.1 fix: v1.0 logged nowhere, which
REM hid the S4_03 crash -- measured 2026-09-22, R6).
REM Prerequisite: package metafor installed once:
REM   Rscript -e "utils::install.packages('metafor')"
REM If a step fails: open its log under logs\, fix the cause, re-run this file;
REM completed steps re-verify their inputs (registry/audit) and are fast.
REM ============================================================================
setlocal
set R="C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
cd /d "%~dp0\.."
if not exist logs mkdir logs

echo [S4] 0/5 smoke test (8 legs, known-answer math) --^> logs\S4_smoke.log
%R% scripts\S4_smoke.R > logs\S4_smoke.log 2>&1
if errorlevel 1 type logs\S4_smoke.log & goto :fail

echo [S4] 1/5 primary REML meta (6 cohorts) + Tier-1 core + GSE11691 arm --^> logs\S4_01.log
%R% scripts\S4_01_meta_bulk.R > logs\S4_01.log 2>&1
if errorlevel 1 type logs\S4_01.log & goto :fail

echo [S4] 2/5 sensitivity battery (sans-G7307, LOCO x6, FE-vs-RE, QC swap x2) --^> logs\S4_02.log
%R% scripts\S4_02_sensitivity.R > logs\S4_02.log 2>&1
if errorlevel 1 type logs\S4_02.log & goto :fail

echo [S4] 3/5 scRNA projection of the Tier-1 core --^> logs\S4_03.log
%R% scripts\S4_03_meta_scoring.R > logs\S4_03.log 2>&1
if errorlevel 1 type logs\S4_03.log & goto :fail

echo [S4] 4/5 immunoediting synthesis + knock summary --^> logs\S4_04.log
%R% scripts\S4_04_collect.R > logs\S4_04.log 2>&1
if errorlevel 1 type logs\S4_04.log & goto :fail

echo [S4] ALL DONE -- see results\S4_meta\S4_knock_summary.md
exit /b 0

:fail
echo [S4] STEP FAILED -- full log printed above and saved under logs\.
echo        Fix the cause, then re-run this file (done steps re-verify fast).
pause
exit /b 1
