# ============================================================================
# S2_01b_qc_GSE183837.R — surgical v5: GSE183837 ONLY (9 samples: 3 Ctrl + 6 RIF)
# Run from project root:
#   Rscript scripts\S2_01b_qc_GSE183837.R > logs\S2_01b_run.txt 2>&1
# Runtime: ~30-60 min (9 samples, scDblFinder each).
#
# WHY SCOPED TO GSE183837 (documented decision, 2026-09-14):
# The propagation diagnostic (logs/S2_prop_diag.txt) proved the v4.1 pipeline
# CORRECT on BOTH input paths: on E09_EcPA / E10_EcO (H5) and EMS1 (trio) the
# recomputed removal+QC values equal the v4 CSVs EXACTLY (378 / 10749 / 6496),
# doublet rates match (3.85 / 10.34 / 11.63), cell order is preserved, and zero
# doublet barcodes survive the subset. The v4 GSE179640 "propagation BLOCKER"
# is therefore LIFTED (false alarm rooted in a mislabelled reference
# comparison; see logs/S2_prop_diag.txt + analysis_log entry 2026-09-14).
# GSE179640 and GSE214411 outputs are VERIFIED FINAL and must NOT be
# regenerated: a re-run would deterministically reproduce identical files
# (demonstrated) at hours of cost and nonzero overwrite risk. This script
# touches ONLY GSE183837 and NEVER opens results/S2_qc_GSE179640.csv,
# results/S2_qc_GSE214411.csv, or their rds/figures.
#
# ROOT CAUSE FIXED HERE (GSE183837 silent skip in v3/v4/v4.1):
# The tar extracts FLAT with DOT-named trios:
#   GSM5572238_Ctrl-1.barcodes.tsv.gz / .genes.tsv.gz / .matrix.mtx.gz
# while v4.1 staging only matched UNDERSCORE-named trios
# (sample_barcodes.tsv.gz) — the regex never matched, zero dirs were staged,
# and the loud "no 10x trio dirs" banner fired every run. Staging below
# accepts BOTH separators ([._]) and keeps the genes->features rename.
#
# Pre-registered thresholds: logs/S2_preregistered_design.md Section 2
# (GSE183837: 200-6000 genes, mt <= 20%, Hb <= 5%) — UNCHANGED (R5).
# qc_one_sample() below is VERBATIM IDENTICAL to S2_01 v4.1 (same
# construction, same set.seed(42) policy) so the S2_05 determinism check
# remains valid for GSE183837.
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(scDblFinder)
  library(SingleCellExperiment)
})
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

THR <- list(min.genes = 200, max.genes = 6000, max.mt = 20, max.hb = 5)  # GSE183837 (pre-registered)

# ---- qc_one_sample: VERBATIM from S2_01_qc_singlecell.R v4.1 ----------------
qc_one_sample <- function(counts, sample_id, thr) {
  seu <- CreateSeuratObject(counts, project = sample_id,
                            min.cells = 3, min.features = 0)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  hb <- grep("^HB[AB]", rownames(seu), value = TRUE)
  seu[["percent.hb"]] <- if (length(hb)) PercentageFeatureSet(seu, features = hb) else 0
  n0 <- ncol(seu)
  dbl <- NA_real_
  sce <- tryCatch({
    s <- as.SingleCellExperiment(seu)
    set.seed(42)
    scDblFinder(s, samples = NULL)
  }, error = function(e) { cat("scDblFinder failed on", sample_id, "- fallback flag\n"); NULL })
  if (!is.null(sce)) {
    keep  <- as.logical(sce$scDblFinder.class == "singlet")
    dbl   <- round(mean(!keep) * 100, 2)
    n_pre <- ncol(seu)
    seu   <- subset(seu, cells = colnames(seu)[keep])
    cat("  scDblFinder:", sum(!keep), "doublets removed (", dbl, "%), cells ",
        n_pre, " -> ", ncol(seu), "\n", sep = "")
    if (ncol(seu) != sum(keep) || ncol(seu) == n_pre)
      cat("  !! DOUBLET REMOVAL INCONSISTENT for ", sample_id,
          " — STOP and report before proceeding\n")
  }
  seu <- subset(seu, subset = nFeature_RNA >= thr$min.genes &
                          nFeature_RNA <= thr$max.genes &
                          percent.mt <= thr$max.mt & percent.hb <= thr$max.hb)
  data.frame(sample = sample_id, cells_before = n0, cells_after = ncol(seu),
             pct_kept = round(100 * ncol(seu) / n0, 1),
             dbl_pct = round(dbl, 2),
             median_genes = median(seu$nFeature_RNA),
             median_mt = round(median(seu$percent.mt), 2),
             stringsAsFactors = FALSE) -> stats
  list(obj = seu, stats = stats)
}

# ---- GSE183837: extract if needed, stage DOT- or UNDERSCORE-named trios -----
gse     <- "GSE183837"
tarfile <- file.path("data/raw", paste0(gse, "_RAW.tar"))
exdir   <- file.path("data/raw", gse)
if ((!dir.exists(exdir) || !length(list.files(exdir))) && file.exists(tarfile))
  untar(tarfile, exdir = exdir, tar = "internal")
if (!dir.exists(exdir))
  stop(gse, " folder missing and no tar found — restore data/raw/GSE183837_RAW.tar (MD5 ced92ae3a4c43b8c030a143faa4a280b) and re-run.")

# v5 staging: sample id = filename up to "_"/"." + barcodes|features|genes|matrix
# (handles both GSM..._Ctrl-1.barcodes.tsv.gz and GSM..._EMS1_barcodes.tsv.gz)
tri_files <- list.files(exdir, pattern = "[._](barcodes|features|genes|matrix)\\.(tsv|mtx)\\.gz$")
samples   <- unique(sub("[._](barcodes|features|genes|matrix)\\..*$", "", tri_files))
cat("top-level trio files:", length(tri_files), "| unique samples detected:", length(samples), "\n")
for (s in samples) {                       # move top-level trios into per-sample dirs (no-op on re-run)
  d <- file.path(exdir, s); dir.create(d, showWarnings = FALSE)
  for (ext in c("barcodes.tsv.gz", "features.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz")) {
    for (sep in c("_", ".")) {
      src <- file.path(exdir, paste0(s, sep, ext))
      # 10x v2 'genes.tsv.gz' is renamed to 'features.tsv.gz' for Read10X consistency
      if (file.exists(src)) file.rename(src, file.path(d, sub("^genes", "features", ext)))
    }
  }
}
dirs <- list.dirs(exdir, recursive = FALSE, full.names = TRUE)
dirs <- dirs[file.exists(file.path(dirs, "matrix.mtx.gz"))]
if (!length(dirs))
  stop(gse, " STILL has no staged trio dirs after v5 staging — paste `dir data\\raw\\GSE183837 /b` and report.")
if (length(dirs) != 9)
  cat("!! EXPECTED 9 samples (3 Ctrl + 6 RIF), found ", length(dirs), " — report before proceeding\n")

# condition mapping print (Amendment S2-B/B2: to be verified against GEO before freezing)
cond_map <- data.frame(sample = basename(dirs),
  condition = ifelse(grepl("RIF", basename(dirs), ignore.case = TRUE), "RIF", "Control"))
cat("\nGSE183837 condition mapping (name-based; verify vs GEO metadata before freezing):\n")
print(cond_map)

# ---- QC all 9 samples --------------------------------------------------------
cat("\n##### ", gse, " ##### (", length(dirs), " samples )\n")
objs <- list(); all_stats <- list()
for (d in dirs) {
  sid <- basename(d); cat("loading", sid, "\n")
  r <- qc_one_sample(Read10X(d), sid, THR)
  objs[[sid]] <- r$obj; all_stats[[sid]] <- r$stats
}
st <- do.call(rbind, all_stats)
write.csv(st, "results/S2_qc_GSE183837.csv", row.names = FALSE)
merged <- if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
p <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3, pt.size = 0, group.by = "orig.ident") + NoLegend()
ggsave("figures/S2_qc_GSE183837.png", p, width = 14, height = 6)
saveRDS(merged, "data/processed/GSE183837_filtered.rds")
print(st)
cat("\nS2_01b done. Wrote results/S2_qc_GSE183837.csv, figures/S2_qc_GSE183837.png,",
    "data/processed/GSE183837_filtered.rds — nothing else was touched.\n")
