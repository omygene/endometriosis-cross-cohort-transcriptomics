# ============================================================================
# S2_01_diagnose_propagation.R — ONE-SHOT diagnostic (run BEFORE any v5 re-run)
# Why: v4 run log shows an impossibility for GSE179640 — per-sample removal
# lines printed drops (e.g. C01 4133 -> 3855) yet the final QC table is
# value-identical to the pre-removal era (all 33 samples: cells_after,
# median_genes, median_mt exactly equal). Arithmetically that requires EVERY
# called doublet to have failed QC anyway (e.g. E10_EcO: all 1322 doublets in
# the 2039 QC-failures, p ~ 0) — so propagation is broken on the H5 path while
# the SAME function propagated correctly on the GSE214411 trio path.
# This script pinpoints the exact step, on real data, no edits to any result.
# Also: extracts + lists GSE183837 (skipped in v4 with an empty folder).
# Run: Rscript scripts\S2_01_diagnose_propagation.R > logs\S2_prop_diag.txt 2>&1
# Runtime: ~10-20 min.
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat); library(scDblFinder); library(SingleCellExperiment)
})
cat("Seurat", as.character(packageVersion("Seurat")),
    "| SeuratObject", as.character(packageVersion("SeuratObject")),
    "| scDblFinder", as.character(packageVersion("scDblFinder")), "\n")

diag_one <- function(counts, sid, mt_thr) {
  cat("\n=====", sid, "=====\n")
  seu <- CreateSeuratObject(counts, project = sid, min.cells = 3, min.features = 0)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  hb <- grep("^HB[AB]", rownames(seu), value = TRUE)
  seu[["percent.hb"]] <- if (length(hb)) PercentageFeatureSet(seu, features = hb) else 0
  cat("n0:", ncol(seu), "| unique cell names:", length(unique(colnames(seu))), "\n")
  sce <- as.SingleCellExperiment(seu)
  cat("ncol(sce):", ncol(sce), "| same cell order seu~sce:",
      identical(colnames(seu), colnames(sce)), "\n")
  set.seed(42)
  sce <- scDblFinder(sce, samples = NULL)
  is_dbl  <- as.logical(sce$scDblFinder.class == "doublet")
  cat("class table:", paste(names(table(sce$scDblFinder.class)),
      as.integer(table(sce$scDblFinder.class)), collapse = " / "), "\n")
  qc_pass <- seu$nFeature_RNA >= 200 & seu$nFeature_RNA <= 6000 &
             seu$percent.mt <= mt_thr & seu$percent.hb <= 5
  cat("QC-pass:", sum(qc_pass), "| doublets:", sum(is_dbl),
      "| doublets passing QC:", sum(is_dbl & qc_pass), "\n")
  keep <- !is_dbl
  cat("keep length:", length(keep), "sum:", sum(keep), "\n")
  seu2 <- subset(seu, cells = colnames(seu)[keep])
  cat("after removal:", ncol(seu2), "(expect", sum(keep), ")",
      "| doublet barcodes still inside:", sum(colnames(seu)[is_dbl] %in% colnames(seu2)), "\n")
  qc_expr <- quote(nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                   percent.mt <= mt_thr & percent.hb <= 5)
  seu3 <- subset(seu2, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                            percent.mt <= mt_thr & percent.hb <= 5)
  seu4 <- subset(seu,  subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                            percent.mt <= mt_thr & percent.hb <= 5)
  cat("QC after removal:", ncol(seu3), "| QC without removal:", ncol(seu4),
      "| identical sets:", identical(sort(colnames(seu3)), sort(colnames(seu4))), "\n")
}

# two GSE179640 H5 samples: smallest (E09_EcPA) and the worst doublet case (E10_EcO)
diag_one(Read10X_h5("data/raw/GSE179640/GSM6595259_E09_EcPA_filtered_feature_bc_matrix.h5"),
         "E09_EcPA", 20)
diag_one(Read10X_h5("data/raw/GSE179640/GSM6595261_E10_EcO_filtered_feature_bc_matrix.h5"),
         "E10_EcO", 20)
# one GSE214411 trio as the working-control contrast (this path propagated fine in v4)
diag_one(Read10X("data/raw/GSE214411/GSM6605431_EMS1"), "EMS1_trio", 10)

# ---- GSE183837: extract (if still empty) + reveal structure ------------------
ex <- "data/raw/GSE183837"
if (!dir.exists(ex) || !length(list.files(ex))) {
  cat("\nextracting GSE183837_RAW.tar ...\n")
  untar("data/raw/GSE183837_RAW.tar", exdir = ex, tar = "internal")
}
cat("\nGSE183837 top-level (first 40 entries):\n")
print(head(list.files(ex), 40))
cat("GSE183837 immediate subdirs:\n")
print(head(list.dirs(ex, recursive = FALSE, full.names = FALSE), 20))
cat("\ndiagnostic done.\n")
