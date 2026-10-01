#!/usr/bin/env bash
# ============================================================================
# S1_download_datasets.sh — Benign Immunoediting in Endometriosis
# Run on the LEAD RESEARCHER's machine. Downloads every pending dataset,
# verifies sizes + MD5 checksums, and writes a local manifest.
# Charter refs: S1 acquisition; R2 traceability; QC gate S1.
# Usage: bash S1_download_datasets.sh /path/to/endometriosis_immunoediting
# ============================================================================
set -euo pipefail
ROOT="${1:-.}"
RAW="$ROOT/data/raw"
mkdir -p "$RAW"
LOG="$ROOT/logs/S1_download_$(date +%Y%m%d_%H%M%S).log"
echo "# S1 download log — $(date)" | tee "$LOG"

dl () {  # dl <url> <dest>  — resumable, retrying downloader
  local url="$1" dest="$2"
  echo "## $(date +%H:%M:%S) downloading $(basename "$dest")" | tee -a "$LOG"
  for try in 1 2 3 4 5; do
    curl -sS --retry 3 -C - --max-time 1800 -o "$dest" "$url" && break
    echo "   retry $try" | tee -a "$LOG"; sleep 20
  done
}

# ---- 1. GSE51981 (bulk primary, n=148) — full archive -----------------------
dl "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE51981&format=file" "$RAW/GSE51981_RAW.tar"
echo "GSE51981 size: $(stat -c%s "$RAW/GSE51981_RAW.tar") (expect 746424320)" | tee -a "$LOG"

# ---- 2. GSE179640 (single-cell validation arm) ------------------------------
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/suppl/GSE179640_RAW.tar" "$RAW/GSE179640_RAW.tar"
echo "GSE179640 size: $(stat -c%s "$RAW/GSE179640_RAW.tar") (expect 795351040)" | tee -a "$LOG"

# ---- 3. GSE214411 (minimal/mild gradient) -----------------------------------
dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214411/suppl/GSE214411_RAW.tar" "$RAW/GSE214411_RAW.tar"
echo "GSE214411 size: $(stat -c%s "$RAW/GSE214411_RAW.tar") (expect ~1567782912)" | tee -a "$LOG"

# ---- 4. GSE7307 subset (41 CELs; frozen list in S1_GSE7307_subset.txt) ------
mkdir -p "$RAW/GSE7307_subset"
grep -v '^#' "$(dirname "$0")/S1_GSE7307_subset.txt" | while read -r group gsm; do
  [ -z "${gsm:-}" ] && continue
  num="${gsm#GSM}"; prefix="GSM${num%???}nnn"
  dl "https://ftp.ncbi.nlm.nih.gov/geo/samples/${prefix}/${gsm}/suppl/${gsm}.CEL.gz" "$RAW/GSE7307_subset/${gsm}.CEL.gz"
done
echo "GSE7307 subset files: $(ls "$RAW/GSE7307_subset" | wc -l) (expect 41)" | tee -a "$LOG"

# ---- 5. GSE213216 (primary atlas) -------------------------------------------
# OPTION A (recommended): processed Seurat objects from the authors' Box:
#   https://cedars.box.com/s/1ks3eyzlpnjbrseefw3j4k7nx6p2ut02
#   Download via browser into data/raw/GSE213216_processed/
# OPTION B: full raw per-sample archives from GEO (16.8 GB total):
#   dl "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213216/suppl/GSE213216_RAW.tar" "$RAW/GSE213216_RAW.tar"

# ---- 6. Checksums ------------------------------------------------------------
echo "## MD5 checksums" | tee -a "$LOG"
( cd "$RAW" && find . -type f -exec md5sum {} \; ) | tee -a "$LOG" > "$ROOT/logs/S1_md5_$(date +%Y%m%d).txt"
echo "DONE. Send back: $LOG and the md5 file." | tee -a "$LOG"
