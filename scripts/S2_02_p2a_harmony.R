# ============================================================================
# S2_02_p2a_harmony.R — step 2a/4: Harmony only (fresh process).
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
req <- c("Seurat", "harmony")
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

merged <- load_cp("pre_harmony")

t_harmony <- Sys.time()
# project.dim = FALSE (Addendum S2_02-B.5, source-verified harmony 2.0.5
# RunHarmony.R): the corrected embeddings are computed and stored BEFORE the
# project.dim branch — identical either way. The TRUE default only triggers
# Seurat::ProjectDim, which needs the scale.data layer that Amendment S2_02-B
# deliberately drops before Harmony (RAM) — the smoke test caught this
# ("Layer 'scale.data' is empty" -> non-conformable arguments). No gate or
# downstream step uses harmony loadings (neighbors/clusters/UMAP/LISI read
# embeddings only).
merged <- harmony::RunHarmony(merged, group.by.vars = c("sample", "dataset"),
                              theta = c(2, 2), reduction = "pca",
                              dims.use = 1:30, reduction.save = "harmony",
                              project.dim = FALSE)
say("Harmony done in ",
    round(as.numeric(difftime(Sys.time(), t_harmony, units = "mins")), 1),
    " min")
gc(verbose = FALSE)


saveRDS(merged, file.path(tmpdir, "checkpoint_post_harmony.rds"))
say("checkpoint saved: post_harmony — STEP 2a DONE")
