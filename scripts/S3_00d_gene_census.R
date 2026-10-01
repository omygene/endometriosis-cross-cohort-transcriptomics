# ============================================================================
# S3_00d_gene_census.R -- READ-ONLY diagnostic (R9 direct measurement)
# Trigger 2026-09-21: S3_03 kept only 8-28 genes after CPM screen; S3_04
# returned exactly 28 genes in every dataset x lineage; S3_05 found only
# 1 of 124 SenMayo symbols in the RNA assay rownames. Three measurements,
# one suspected root cause: the atlas RNA assay feature space is not what
# the S3 scripts assume. This census measures DIRECTLY, changes NOTHING.
# v1.1 (2026-09-21): GetAssay is NOT exported in SeuratObject >= 5.1
# (measured halt on the lead machine) -- assay access via obj[[name]].
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_00d_gene_census.R > logs\S3_00d_gene_census.txt 2>&1
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(Seurat) })
source("scripts/S3_common.R")

say("== S3_00d gene census (read-only diagnostic, v1.1) ==")
obj <- readRDS(CFG$atlas)
say("atlas:", CFG$atlas)
say("assays:", paste(SeuratObject::Assays(obj), collapse = " | "))
say("default assay:", DefaultAssay(obj))

gm <- find_senmayo_gmt()
sen <- unique(unlist(lapply(read_gmt(gm), `[[`, "genes")))
say("SenMayo symbols in GMT:", length(sen))

audit_assay <- function(obj, assay) {
  say("---- assay:", assay, "----")
  a5 <- obj[[assay]]
  say("  dim:", nrow(a5), "x", ncol(a5))
  rn <- rownames(a5)
  say("  rownames[1:10]:", paste(head(rn, 10), collapse = " | "))
  say("  n rownames:", length(rn), "| unique:", length(unique(rn)))
  say("  Ensembl-like (^ENSG):", sum(grepl("^ENSG", rn)),
      "| GSM-like (^GSM):", sum(grepl("^GSM", rn)),
      "| MT- prefix:", sum(grepl("^MT-", rn)),
      "| with version suffix (\\.\\d+$):", sum(grepl("\\.\\d+$", rn)))
  cnt <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
                  error = function(e) NULL)
  if (!is.null(cnt)) {
    nz <- Matrix::rowSums(cnt) > 0
    say("  counts: nonzero features:", sum(nz), "of", nrow(cnt))
    if (sum(nz) <= 200) {
      say("  ALL nonzero feature names:")
      for (g in rownames(cnt)[nz]) say("    ", g)
    } else {
      say("  nonzero feature names[1:20]:", paste(head(rownames(cnt)[nz], 20), collapse = " | "))
    }
  } else say("  counts layer: not retrievable")
  say("  SenMayo overlap with rownames:", length(intersect(sen, rn)), "of", length(sen))
}

for (a in SeuratObject::Assays(obj)) audit_assay(obj, a)

say("S3_00d DONE -- nothing was written or modified (read-only)")
