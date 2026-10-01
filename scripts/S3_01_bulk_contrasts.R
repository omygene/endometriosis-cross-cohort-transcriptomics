# ============================================================================
# S3_01_bulk_contrasts.R -- bulk primary contrasts per frozen registry
# v2.0 (2026-09-22): SYMBOL-level DE. Measured failure ledger #18: this script
# previously fitted probe-level DE (raw platform rownames), leaving the
# Affy cohorts in one probe namespace and GSE120103 (Agilent) in another --
# S4's cross-cohort meta silently lost GSE120103 (LOCO_B-G120103 identical
# to base meta, dM = 0.0). Fix: load_bulk_symbol_matrix() (S3_common v3.0.2,
# the frozen probe->symbol collapse proven in S3_02). Contrast specs, frozen
# criteria, and outputs are otherwise UNCHANGED (R5); re-run regenerates all
# contrast CSVs from the locked scripts (R4-safe).
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_01_bulk_contrasts.R <COHORT|collect>
# Locked criteria: BH-FDR < 0.05 gate; |log2FC| >= 0.5 interpretive tier
# (reported, NO exclusion); sign concordance reported in S3_06.
# One short-lived process per cohort (B.6/B.7); collect verifies completeness.
# Writes ONLY: results/S3_bulk_contrasts/<contrast_id>.csv (+ summary at collect)
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(limma); library(Biobase) })
source("scripts/S3_common.R")
dir.create(CFG$out_contrasts, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/S3_01_bulk_contrasts.R <COHORT|collect>")
cohort <- args[1]
fr <- load_frozen()
meta <- load_meta()

if (cohort == "collect") {
  say("== S3_01 collect ==")
  expect <- fr$contrast_id
  files <- file.path(CFG$out_contrasts, paste0(expect, ".csv"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("FATAL: missing contrast outputs: ", paste(missing, collapse = ", "))
  rows <- lapply(files, read.csv, stringsAsFactors = FALSE)
  summ <- data.frame(contrast_id = expect, level = fr$level,
                     n_test = vapply(rows, function(z) z$n_test[1], integer(1)),
                     n_ref = vapply(rows, function(z) z$n_ref[1], integer(1)),
                     n_FDR = vapply(rows, function(z) sum(z$adj.P.Val < CFG$fdr_gate), integer(1)),
                     n_FDR_tier = vapply(rows, function(z) sum(z$adj.P.Val < CFG$fdr_gate & abs(z$logFC) >= CFG$lfc_tier), integer(1)),
                     stringsAsFactors = FALSE)
  write.csv(summ, file.path(CFG$out_contrasts, "S3_bulk_contrasts_summary.csv"), row.names = FALSE)
  say("  contrasts verified:", length(expect), "| summary written")
  print(summ)
  say("S3_01 COLLECT DONE")
  quit(save = "no", status = 0)
}

if (!cohort %in% unique(fr$cohort)) stop("Unknown cohort: ", cohort)
specs <- fr[fr$cohort == cohort, ]
mat <- load_bulk_symbol_matrix(cohort)   # v2.0: symbol-level (ledger #18)
m <- map_cols(cohort, colnames(mat), meta)

for (i in seq_len(nrow(specs))) {
  spec <- specs[i, ]
  say("contrast:", spec$contrast_id, "(", spec$level, ")")
  tt <- run_frozen_contrast(mat, m, spec)
  out <- file.path(CFG$out_contrasts, paste0(spec$contrast_id, ".csv"))
  write.csv(tt, out, row.names = FALSE)
  say("  wrote", out, "| n genes:", nrow(tt),
      "| FDR<", CFG$fdr_gate, ":", sum(tt$adj.P.Val < CFG$fdr_gate),
      "| +tier:", sum(tt$adj.P.Val < CFG$fdr_gate & abs(tt$logFC) >= CFG$lfc_tier))
}
print(gc())
say("COHORT", cohort, "DONE")
