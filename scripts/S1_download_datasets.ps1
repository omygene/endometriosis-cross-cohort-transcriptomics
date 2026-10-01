# ============================================================================
# S1_download_datasets.ps1 — Benign Immunoediting in Endometriosis (Windows)
# Run in PowerShell on the lead researcher's Windows machine.
# Downloads every pending dataset, verifies sizes, computes MD5 checksums.
# Charter refs: S1 acquisition; R2 traceability; QC gate S1.
#
# Usage:
#   1. Open PowerShell
#   2. cd to your project folder (must contain scripts\S1_GSE7307_subset.txt)
#   3. powershell -ExecutionPolicy Bypass -File .\scripts\S1_download_datasets.ps1
# ============================================================================

$ErrorActionPreference = "Stop"
$ROOT = (Get-Location).Path
$RAW  = Join-Path $ROOT "data\raw"
$LOGS = Join-Path $ROOT "logs"
New-Item -ItemType Directory -Force -Path $RAW, $LOGS | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$LOG = Join-Path $LOGS "S1_download_$stamp.log"

function Log($msg) {
    $line = "$(Get-Date -Format 'HH:mm:ss')  $msg"
    Write-Host $line
    Add-Content -Path $LOG -Value $line
}

function Get-Dataset($url, $dest) {
    Log "downloading $(Split-Path $dest -Leaf)"
    for ($try = 1; $try -le 5; $try++) {
        # curl.exe ships with Windows 10/11; -C - resumes partial files
        & curl.exe -sS -L -C - --retry 3 --max-time 3600 -o "$dest" "$url"
        if ($LASTEXITCODE -eq 0) { return }
        Log "  retry $try"; Start-Sleep -Seconds 20
    }
    throw "FAILED: $dest"
}

# ---- 1. GSE51981 (bulk primary, n=148) --------------------------------------
Get-Dataset "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE51981&format=file" "$RAW\GSE51981_RAW.tar"
Log ("GSE51981 size: {0} (expect 746424320)" -f (Get-Item "$RAW\GSE51981_RAW.tar").Length)

# ---- 2. GSE179640 (single-cell validation arm) ------------------------------
Get-Dataset "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/suppl/GSE179640_RAW.tar" "$RAW\GSE179640_RAW.tar"
Log ("GSE179640 size: {0} (expect 795351040)" -f (Get-Item "$RAW\GSE179640_RAW.tar").Length)

# ---- 3. GSE214411 (minimal/mild gradient) -----------------------------------
Get-Dataset "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214411/suppl/GSE214411_RAW.tar" "$RAW\GSE214411_RAW.tar"
Log ("GSE214411 size: {0} (expect ~1567782912)" -f (Get-Item "$RAW\GSE214411_RAW.tar").Length)

# ---- 4. GSE7307 subset (41 CELs; frozen list) --------------------------------
$SUB = Join-Path $RAW "GSE7307_subset"
New-Item -ItemType Directory -Force -Path $SUB | Out-Null
$LIST = Join-Path $ROOT "scripts\S1_GSE7307_subset.txt"
Get-Content $LIST | Where-Object { $_ -notmatch '^\s*#' -and $_ -match '\S' } | ForEach-Object {
    $parts = $_ -split '\s+'
    $gsm = $parts[1]
    $num = $gsm -replace 'GSM',''
    $prefix = "GSM" + $num.Substring(0, $num.Length - 3) + "nnn"
    Get-Dataset "https://ftp.ncbi.nlm.nih.gov/geo/samples/$prefix/$gsm/suppl/$gsm.CEL.gz" "$SUB\$gsm.CEL.gz"
}
Log ("GSE7307 subset files: {0} (expect 41)" -f (Get-ChildItem $SUB -Filter *.CEL.gz).Count)

# ---- 5. GSE213216 (primary atlas) -------------------------------------------
# OPTION A (recommended): download the processed Seurat objects via your browser
#   from https://cedars.box.com/s/1ks3eyzlpnjbrseefw3j4k7nx6p2ut02
#   into data\raw\GSE213216_processed\
# OPTION B: full raw archive (16.8 GB):
#   Get-Dataset "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/suppl/GSE213216_RAW.tar" "$RAW\GSE213216_RAW.tar"

# ---- 6. MD5 checksums ---------------------------------------------------------
$MD5FILE = Join-Path $LOGS "S1_md5_$(Get-Date -Format 'yyyyMMdd').txt"
Get-ChildItem $RAW -Recurse -File | ForEach-Object {
    $h = Get-FileHash $_.FullName -Algorithm MD5
    "{0}  {1}" -f $h.Hash.ToLower(), $_.FullName.Substring($RAW.Length + 1)
} | Tee-Object -FilePath $MD5FILE

Log "DONE. Send back to the project: $LOG and $MD5FILE"
