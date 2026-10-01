# ============================================================================
# S3_00f_join_audit.R -- READ-ONLY diagnostic (R9 direct measurement)
# Trigger 2026-09-21 (measured by S3_00d/S3_00e): the atlas carries all S3
# annotations (dataset/sample/condition/lineage) but only a 28-gene panel;
# the per-dataset objects carry full expression but QC-only metadata. The
# S3 fix will join annotations -> expression BY CELL BARCODE. Before writing
# that fix, measure the join key DIRECTLY: do atlas colnames match the
# per-dataset objects' colnames as-is, or did Seurat prefix them at merge?
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_00f_join_audit.R > logs\S3_00f_join_audit.txt 2>&1
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(Seurat) })
source("scripts/S3_common.R")

say("== S3_00f join audit (read-only) ==")
obj <- readRDS(CFG$atlas)
md <- obj@meta.data
atlas_cn <- colnames(obj)
say("atlas cells:", length(atlas_cn), "| unique:", length(unique(atlas_cn)))
say("atlas colnames[1:5]:", paste(head(atlas_cn, 5), collapse = " | "))
say("md$sample[1:3]:", paste(head(as.character(md$sample), 3), collapse = " | "))
say("md$orig.ident[1:3]:", paste(head(as.character(md$orig.ident), 3), collapse=" | "))
rm(obj); gc()

per_ds <- c(GSE179640 = "data/processed/GSE179640_filtered.rds",
            GSE183837 = "data/processed/GSE183837_filtered.rds",
            GSE214411 = "data/processed/GSE214411_filtered.rds")
for (ds in names(per_ds)) {
  say("----", ds, "----")
  x <- readRDS(per_ds[[ds]])
  cn <- colnames(x)
  sub <- md$dataset == ds
  ac <- atlas_cn[sub]
  say("  per-dataset cells:", length(cn), "| atlas cells for dataset:", sum(sub))
  say("  per-dataset colnames[1:5]:", paste(head(cn, 5), collapse = " | "))
  say("  direct barcode match:", sum(ac %in% cn), "of", length(ac))
  ## candidate prefixed keys: <orig.ident>_<barcode> ; <orig.ident>:<barcode> ;
  ## <sample>_<barcode> ; strip trailing -1 style then retry
  oi <- as.character(md$orig.ident[sub]); smp <- as.character(md$sample[sub])
  v1 <- paste0(oi, "_", ac); say("  match as <orig.ident>_<atlas>:", sum(v1 %in% cn))
  v2 <- paste0(smp, "_", ac); say("  match as <sample>_<atlas>:", sum(v2 %in% cn))
  v3 <- sub("-[0-9]+$", "", ac); say("  match after stripping -N suffix:", sum(v3 %in% cn))
  ## per-sample barcode purity: atlas sample -> per-dataset orig.ident
  if (exists("oi")) {
    tab <- table(paste(smp, oi))
    say("  sample x orig.ident pairs:", length(tab), "| atlas sample values:", paste(sort(unique(smp)), collapse=" | "))
  }
  rm(x, cn); gc()
}
say("S3_00f DONE -- nothing was written or modified (read-only)")
