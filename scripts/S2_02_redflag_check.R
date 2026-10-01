# ============================================================================
# S2_02_redflag_check.R — the ONE focused red-flag check for S2_02 (charter).
# READ-ONLY diagnostic: reads the final deliverable, writes NOTHING to results/.
# Run:  Rscript scripts\S2_02_redflag_check.R > logs\S2_02_redflag_check.txt 2>&1
# Purpose (R9 — measurement, not theory):
#   A. Independent integrity verification of the deliverable object (R2).
#   B. Direct measurement of exclusive cluster 28 (and near-exclusive 26/27):
#      real biology (patient-specific rare state) or artifact? This measurement
#      is the basis for the lead's written decision on the exclusive-cluster
#      gate (registered: "Integration FAILS if any cluster is exclusive").
# ============================================================================
set.seed(42)
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}
suppressPackageStartupMessages(library(Seurat))

obj <- readRDS("data/processed/S2_02_integrated.rds")
md  <- obj@meta.data

# ---- A. deliverable integrity ------------------------------------------------
say("== A. deliverable integrity ==")
say("cells: ", ncol(obj), " -> registered 233459: ", ncol(obj) == 233459L)
say("clusters: ", nlevels(obj$seurat_clusters), " (expected 29)")
say("genes in light assay: ", nrow(obj))
say("lineage NA cells: ", sum(is.na(obj$lineage)))
say("condition NA cells: ", sum(is.na(md$condition)))
for (v in c("cLISI_sample", "cLISI_dataset", "iLISI_lineage")) {
  say(v, " present+finite: ",
      v %in% colnames(md) && all(is.finite(md[[v]])))
}
say("reductions: ", paste(Reductions(obj), collapse = ", "))
ci <- read.csv("results/S2_02_cell_index.csv", stringsAsFactors = FALSE)
say("cell_index rows: ", nrow(ci))
say("cell order identical to index: ", identical(ci$cell, colnames(obj)))
say("cell SET identical to index: ", setequal(ci$cell, colnames(obj)))
say("meta columns: ", paste(colnames(md), collapse = ", "))

# ---- B. exclusive/near-exclusive cluster composition --------------------------
say("")
say("== B. exclusive-cluster measurement ==")
cl <- as.character(md$seurat_clusters)
for (cc in c("28", "27", "26")) {
  sub <- md[cl == cc, ]
  say("--- cluster ", cc, " | n = ", nrow(sub),
      " | share of total = ", round(100 * nrow(sub) / nrow(md), 3), "%")
  print(table(sub$sample, sub$condition))
}
# are cluster-28 cells real stromal cells? per-cell mean of the stromal panel
# (light assay: counts layer holds log-normalized values, single layer)
g <- intersect(c("PDGFRB", "DCN", "COL1A1"), rownames(obj))
say("stromal panel genes present: ", paste(g, collapse = ", "))
mat <- LayerData(obj, layer = "counts")[g, , drop = FALSE]
pc  <- colMeans(mat)                       # per-cell stromal-panel mean
stromal_clusters <- unique(as.character(md$seurat_clusters[md$lineage == "Stromal"]))
in28   <- cl == "28"
inStr  <- cl %in% stromal_clusters & !in28
say("per-cell stromal-panel mean | cluster 28:        median ",
    round(median(pc[in28]), 4), "  IQR ", round(IQR(pc[in28]), 4))
say("per-cell stromal-panel mean | other stromal:    median ",
    round(median(pc[inStr]), 4), "  IQR ", round(IQR(pc[inStr]), 4))
say("per-cell stromal-panel mean | all non-stromal:  median ",
    round(median(pc[!in28 & !inStr]), 4), "  IQR ",
    round(IQR(pc[!in28 & !inStr]), 4))

# QC metrics if carried in meta (defensive: names differ across pipelines)
qc <- intersect(colnames(md),
                c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent_mt"))
if (length(qc)) {
  say("QC columns available: ", paste(qc, collapse = ", "))
  for (v in qc) {
    say(v, " | cluster 28 median ", round(median(md[[v]][in28]), 3),
        " | all cells median ", round(median(md[[v]]), 3))
  }
} else {
  say("no QC columns in the light object meta — QC comparison N/A (expected:",
      " QC lived in the frozen per-dataset objects)")
}
say("RED FLAG CHECK DONE")
