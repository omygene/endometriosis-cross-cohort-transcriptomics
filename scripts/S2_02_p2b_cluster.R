# ============================================================================
# S2_02_p2b_cluster.R — step 2b/4: FindNeighbors + FindClusters (res 0.5, frozen).
# Amendment S2_02-B.7 (2026-09-16): part 2 split into SHORT LIVED steps, one R
# process each, chained by checkpoints. Measured basis (R9): three hangs hit
# long-lived processes (Harmony x2 in-process; UMAP phase in v2.5b), while the
# SAME heavy steps completed in fresh short processes (Harmony twice: 10.3 min
# diagnostic, 6.6 min in v2.5b). Zero parameter/gate/output change vs v2.4 —
# execution topology only. Each step re-runnable alone after a kill.
# ============================================================================
set.seed(42)
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}
req <- c("Seurat")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("MISSING PACKAGES: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages(library(Seurat))
tmpdir <- "data/processed/tmp_S2_02"
load_cp <- function(name) {
  f <- file.path(tmpdir, paste0("checkpoint_", name, ".rds"))
  if (!file.exists(f)) stop("missing checkpoint: ", f, " — run the previous step first")
  obj <- readRDS(f)
  say("loaded ", name, ": ", ncol(obj), " cells")
  if (ncol(obj) != 233459) stop("cell count mismatch — investigate (R6)")
  obj
}

merged <- load_cp("post_harmony")

# ---- 4. neighbors + clusters (res 0.5, frozen) + UMAP ------------------------
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, random.seed = 42, verbose = FALSE)

saveRDS(merged, file.path(tmpdir, "checkpoint_post_cluster.rds"))
say("checkpoint saved: post_cluster — STEP 2b DONE")
