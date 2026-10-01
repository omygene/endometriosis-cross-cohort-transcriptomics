# ============================================================================
# S2_01_diagnose_doublets.R — capture the EXACT scDblFinder failure
# Context: S2_01 ran with dbl_pct = NA on every GSE214411 sample, i.e. the
# pre-registered fallback flag fired silently (tryCatch). This script runs ONE
# sample through the same pipeline WITHOUT the tryCatch so the true error is
# captured verbatim, then writes logs/S2_scDblFinder_diag.txt.
# Run from project root via Rscript. Output: logs/S2_scDblFinder_diag.txt
# ============================================================================
set.seed(42)
sink("logs/S2_scDblFinder_diag.txt", split = TRUE)
cat("scDblFinder diagnostic —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
suppressPackageStartupMessages({
  library(Seurat); library(scDblFinder); library(SingleCellExperiment)
})
cat("scDblFinder version:", as.character(packageVersion("scDblFinder")), "\n")
cat("BiocParallel bpparam: ", as.character(class(BiocParallel::bpparam())), "\n\n")

# pick the first staged GSE214411 sample directory
d <- list.dirs("data/raw/GSE214411", recursive = FALSE)[1]
cat("test sample dir:", d, "\n\n")

counts <- Read10X(d)
seu <- CreateSeuratObject(counts, project = "diag", min.cells = 3)
sce <- as.SingleCellExperiment(seu)
cat("SCE dims:", nrow(sce), "genes x", ncol(sce), "cells\n\n")

cat("--- attempt 1: samples = NULL with set.seed(42) before the call ---\n")
cat("(v2 note: the original bug was seed = 42 INSIDE the call — scDblFinder\n")
cat(" 1.26.7 has no seed argument; verified empirically 2026-09-13.)\n")
r1 <- tryCatch({
  set.seed(42)
  sce1 <- scDblFinder(sce, samples = NULL)
  paste("SUCCESS — doublet rate:",
        round(mean(sce1$scDblFinder.class == "doublet") * 100, 2), "%")
}, error = function(e) paste("FAILED:", conditionMessage(e)))
cat(r1, "\n\n")

if (grepl("FAILED", r1)) {
  cat("--- full traceback of attempt 1 ---\n")
  tryCatch({ set.seed(42); scDblFinder(sce, samples = NULL) },
           error = function(e) { cat(conditionMessage(e), "\n"); print(sys.calls()) })

  cat("\n--- attempt 2: serial evaluation (BPPARAM = SerialParam) ---\n")
  r2 <- tryCatch({
    set.seed(42)
    sce2 <- scDblFinder(sce, samples = NULL,
                        BPPARAM = BiocParallel::SerialParam())
    paste("SUCCESS — doublet rate:",
          round(mean(sce2$scDblFinder.class == "doublet") * 100, 2), "%")
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(r2, "\n")
}
sink()
cat("\nDiagnostics written to logs/S2_scDblFinder_diag.txt — send this file back.\n")
