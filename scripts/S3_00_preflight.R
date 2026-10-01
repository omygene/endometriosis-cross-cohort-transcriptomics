# ============================================================================
# S3_00_preflight.R -- read-only census for stage S3 (NO analysis, R2/R9)
# Measures every input and dependency, then prints the exact values to paste
# into scripts/S3_common.R (config block TO CONFIRM AFTER S3_00 CENSUS).
# Run from project root:
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_00_preflight.R > logs\S3_00_preflight.txt 2>&1
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
dir.create("logs", showWarnings = FALSE)
say("== S3_00 preflight -- read-only census ==")

## 1. bulk inputs
meta <- load_meta()
fr <- load_frozen()
say("metadata rows:", nrow(meta), "| frozen contrasts:", nrow(fr))
cohorts <- unique(fr$cohort)
for (ch in cohorts) {
  p <- unname(CFG$bulk_rds[ch])
  ok <- !is.null(p) && file.exists(p)
  say("  ", ch, " object:", if (ok) "PRESENT" else "MISSING", if (!is.null(p)) p else "(unregistered)")
  if (ok) {
    obj <- readRDS(p)
    cn <- if (ch == CFG$agilent) colnames(obj$E) else colnames(Biobase::exprs(obj))
    m <- map_cols(ch, cn, meta)
    say("    colnames:", length(cn), "| conditions:",
        paste(names(table(m$condition)), table(m$condition), collapse = " / "))
    if (ch == "GSE51981")
      say("    _103/_134 stems present:",
          sum(grepl("GSM1256659_103", cn, fixed = TRUE)), "/",
          sum(grepl("GSM1256659_134", cn, fixed = TRUE)))
  }
}

## 2. GSVA + annotation .db packages (ssGSEA dependency, locked criterion 3)
say("== 2. ssGSEA dependencies ==")
gv <- requireNamespace("GSVA", quietly = TRUE)
say("  GSVA installed:", gv, if (gv) paste0("version ", as.character(utils::packageVersion("GSVA"))) else "-- S3_02 will STOP; needs your written approval to add (R5)")
for (db in c("hgu133plus2.db", "hgu133a.db", "hgu133a2.db"))
  say("  ", db, ":", requireNamespace(db, quietly = TRUE))
gm <- find_senmayo_gmt()
say("  SenMayo GMT:", if (is.na(gm)) "NOT FOUND (searched reference/ and data/ for *.gmt)" else gm)

## 3. scRNA atlas metadata census (fill config from THIS section)
say("== 3. atlas census ==")
if (!file.exists(CFG$atlas)) stop("FATAL: missing ", CFG$atlas, " (R2)")
obj <- readRDS(CFG$atlas)
md <- obj@meta.data
say("  atlas cells:", nrow(md), "| metadata columns:")
say("   ", paste(colnames(md), collapse = " | "))
ds_col <- if ("dataset" %in% colnames(md)) "dataset" else NA_character_
if (is.na(ds_col)) stop("FATAL: no 'dataset' column in atlas metadata")
for (ds in sort(unique(md[[ds_col]]))) {
  sub <- md[md[[ds_col]] == ds, ]
  say("  ---", ds, "samples:", length(unique(sub$sample)), "cells:", nrow(sub))
  for (cc in colnames(sub)) {
    u <- unique(sub[[cc]])
    if (length(u) <= 13 && is.character(u))
      say("    ", cc, "=", paste(sort(u), collapse = " | "))
  }
}
say("  cells matching 'N1' anywhere in metadata:")
n1hits <- apply(md, 2, function(z) sum(grepl("N1", as.character(z), fixed = TRUE)))
n1hits <- n1hits[n1hits > 0]
if (length(n1hits)) say("   ", paste(names(n1hits), n1hits, sep = ":", collapse = "  ")) else say("    (none)")

## 4. pairs audit (read-only)
pr <- read.delim(CFG$pairs11691, comment.char = "#", stringsAsFactors = FALSE)
allg <- c(pr$eutopic_gsm, pr$ectopic_gsm)
say("== 4. GSE11691 pairs:", nrow(pr), "women | distinct GSMs:", length(unique(allg)),
    "| audit:", if (length(unique(allg)) == 18 && nrow(pr) == 9) "OK" else "VIOLATION")

say("PREFLIGHT DONE -- fill scripts/S3_common.R from sections 2-3, then run S3 scripts (R3)")
