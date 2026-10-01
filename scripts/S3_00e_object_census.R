# ============================================================================
# S3_00e_object_census.R -- READ-ONLY diagnostic (R9 direct measurement)
# Trigger 2026-09-21 (measured by S3_00d): data/processed/S2_02_integrated.rds
# is a 28-gene x 233459-cell LINEAGE PANEL object, not a full-expression
# atlas. Before proposing ANY fix, measure what actually exists in
# data/processed/: is there a full integrated object, and do the per-dataset
# objects carry full expression + the metadata columns S3 needs?
# This census loads each .rds READ-ONLY, audits it, and frees it. Nothing
# is written.
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_00e_object_census.R > logs\S3_00e_object_census.txt 2>&1
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(Seurat) })
source("scripts/S3_common.R")

say("== S3_00e object census (read-only) ==")
gm <- find_senmayo_gmt()
sen <- unique(unlist(lapply(read_gmt(gm), `[[`, "genes")))
say("SenMayo symbols in GMT:", length(sen))

d <- "data/processed"
files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("FATAL: no .rds files in ", d, " (R2)")
say("rds files in", d, ":", length(files))
for (f in files) say("  ", basename(f), "--", round(file.size(f) / 1e6, 1), "MB")

need_md <- c("dataset", "condition", "lineage", "sample", "orig.ident")
for (f in files) {
  say("----", basename(f), "----")
  x <- tryCatch(readRDS(f), error = function(e) {
    say("  readRDS FAILED:", conditionMessage(e)); NULL })
  if (is.null(x)) next
  if (!inherits(x, "Seurat")) {
    say("  non-Seurat object (class:", paste(class(x), collapse = ","),
        ") -- bulk/normalized matrix, skipping assay audit")
    rm(x); gc(); next
  }
  say("  assays:", paste(SeuratObject::Assays(x), collapse = " | "))
  a5 <- x[[DefaultAssay(x)]]
  say("  default assay:", DefaultAssay(x), "| dim:", nrow(a5), "x", ncol(a5))
  rn <- rownames(a5)
  say("  rownames[1:8]:", paste(head(rn, 8), collapse = " | "))
  say("  SenMayo overlap:", length(intersect(sen, rn)), "of", length(sen))
  md <- x@meta.data
  say("  metadata cols:", paste(colnames(md), collapse = " | "))
  have_md <- need_md %in% colnames(md)
  say("  S3-needed cols present:", paste(need_md[have_md], collapse = ", "),
      if (any(!have_md)) paste(" | MISSING:", paste(need_md[!have_md], collapse = ", ")) else "")
  if ("dataset" %in% colnames(md)) say("  md$dataset:", paste(sort(unique(as.character(md$dataset))), collapse = " | "))
  if (all(c("condition", "lineage") %in% colnames(md))) {
    say("  condition values:", paste(sort(unique(as.character(md$condition))), collapse = " | "))
    say("  lineage values:", paste(sort(unique(as.character(md$lineage))), collapse = " | "))
  }
  say("  reductions:", if (length(x@reductions)) paste(names(x@reductions), collapse = " | ") else "none")
  rm(x, a5, md); gc()
}
say("S3_00e DONE -- nothing was written or modified (read-only)")
