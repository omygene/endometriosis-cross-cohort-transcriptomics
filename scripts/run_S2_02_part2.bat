@echo off
REM ============================================================
REM S2_02 part 2 chain v2 -- four short-lived processes
REM (Amendment S2_02-B.7, tooling update B.7a)
REM Self-diagnosing: checks prerequisites, prints the failing
REM log ON SCREEN if any step fails.
REM Run from anywhere; it cd's to the project root itself.
REM If a step HANGS (CPU 00, frozen log): kill Rscript.exe in
REM Task Manager, then re-run ONLY that step's Rscript line.
REM ============================================================

REM -- always work from the project root (parent of scripts\)
cd /d "%~dp0\.."
echo Working directory: %CD%

REM -- check 1: Rscript reachable?
where Rscript >nul 2>&1
if errorlevel 1 (
  echo [FAIL] Rscript.exe is not on PATH.
  echo Fix: use the full path, e.g.
  echo   "C:\Program Files\R\R-4.x.x\bin\Rscript.exe" scripts\S2_02_p2a_harmony.R
  goto :eof
)

REM -- check 2: the four step scripts present?
for %%F in (S2_02_p2a_harmony.R S2_02_p2b_cluster.R S2_02_p2c_umap.R S2_02_p2d_finish.R) do (
  if not exist "scripts\%%F" (
    echo [FAIL] missing file: scripts\%%F
    echo Copy all four p2a-p2d scripts into scripts\ then re-run this bat.
    goto :eof
  )
)

REM -- check 3: input checkpoint from PART 1 present?
if not exist "data\processed\tmp_S2_02\checkpoint_pre_harmony.rds" (
  echo [FAIL] missing checkpoint: data\processed\tmp_S2_02\checkpoint_pre_harmony.rds
  echo PART 1 must complete before this chain. Re-run S2_02_integrate_part1.R
  goto :eof
)

echo === STEP 2a Harmony ===
Rscript scripts\S2_02_p2a_harmony.R > logs\S2_02_p2a.txt 2>&1
if errorlevel 1 (
  echo STEP 2a FAILED -- full log below:
  echo ----------------------------------------
  type logs\S2_02_p2a.txt
  echo ----------------------------------------
  goto :eof
)
echo 2a OK.

echo === STEP 2b Clustering ===
Rscript scripts\S2_02_p2b_cluster.R > logs\S2_02_p2b.txt 2>&1
if errorlevel 1 (
  echo STEP 2b FAILED -- full log below:
  echo ----------------------------------------
  type logs\S2_02_p2b.txt
  echo ----------------------------------------
  goto :eof
)
echo 2b OK.

echo === STEP 2c UMAP ===
Rscript scripts\S2_02_p2c_umap.R > logs\S2_02_p2c.txt 2>&1
if errorlevel 1 (
  echo STEP 2c FAILED -- full log below:
  echo ----------------------------------------
  type logs\S2_02_p2c.txt
  echo ----------------------------------------
  goto :eof
)
echo 2c OK.

echo === STEP 2d Finish ===
Rscript scripts\S2_02_p2d_finish.R > logs\S2_02_p2d.txt 2>&1
if errorlevel 1 (
  echo STEP 2d FAILED -- full log below:
  echo ----------------------------------------
  type logs\S2_02_p2d.txt
  echo ----------------------------------------
  goto :eof
)
echo 2d OK.

echo ALL PART2 STEPS DONE
