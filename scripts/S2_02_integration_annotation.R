# ============================================================================
# S2_02_integration_annotation.R — Harmony integration + marker annotation
# Pre-registered: logs/S2_preregistered_design.md (Sections 3-5).
# Inputs: data/processed/*_filtered.rds  Outputs: data/processed/integrated_S2.rds,
# results/S2_batchmixing.csv, results/S2_annotation_validation.csv, figures/S2_umap_*
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2); library(patchwork)
})
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

objs <- list()
for (gse in c("GSE213216_processed", "GSE179640", "GSE214411", "GSE183837")) {
  f <- if (gse == "GSE213216_processed") {
    # merge the 6 author objects (audit said report-only; integration uses them as-is)
    rd <- list.files("data/raw/GSE213216_processed", pattern = "\\.rds$", full.names = TRUE)
    if (length(rd)) {
      ol <- lapply(rd, readRDS)
      m <- if (length(ol) > 1) merge(ol[[1]], ol[-1]) else ol[[1]]
      m$dataset <- "GSE213216"; saveRDS(m, "data/processed/GSE213216_merged.rds"); m
    }
  } else {
    fp <- file.path("data/processed", paste0(gse, "_filtered.rds"))
    if (file.exists(fp)) { o <- readRDS(fp); o$dataset <- gse; o }
  }
  if (!is.null(f)) objs[[gse]] <- f
}
stopifnot(length(objs) >= 2)

# ---- Pre-processing (fixed: 2000 HVGs vst, 30 PCs) ---------------------------
objs <- lapply(objs, function(o) {
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  o
})
merged <- merge(objs[[1]], objs[-1])
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

# ---- Harmony (pre-registered primary method) ---------------------------------
merged <- RunHarmony(merged, group.by.vars = c("orig.ident", "dataset"),
                     theta = 2, reduction = "pca", seed = 42)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, seed.use = 42)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.6, random.seed = 42)

p1 <- DimPlot(merged, group.by = "dataset", raster = TRUE) + ggtitle("by dataset")
p2 <- DimPlot(merged, group.by = "orig.ident", raster = TRUE, label = FALSE) + ggtitle("by sample")
p3 <- DimPlot(merged, label = TRUE, raster = TRUE) + ggtitle("clusters")
ggsave("figures/S2_umap_integration.png", p1 + p2 + p3, width = 21, height = 7)

# ---- Batch-mixing metrics (pre-registered: LISI) -----------------------------
lisi_ok <- requireNamespace("lisi", quietly = TRUE)
if (lisi_ok) {
  emb <- Embeddings(merged, "umap")
  meta <- merged@meta.data[, c("orig.ident", "seurat_clusters")]
  colnames(meta) <- c("batch", "celltype")
  L <- lisi::compute_lisi(emb, meta, c("batch", "celltype"))
  bm <- data.frame(metric = c("cLISI_batch_median", "iLISI_celltype_median"),
                   value = c(median(L$batch), median(L$celltype)))
  write.csv(bm, "results/S2_batchmixing.csv", row.names = FALSE); print(bm)
} else cat("lisi not installed — install.packages('lisi'); metric REQUIRED for gate.\n")

# ---- Marker-based annotation (fixed panels, design §4) -----------------------
markers <- list(
  Stromal     = c("PDGFRB", "DCN", "COL1A1"),
  Epithelial  = c("EPCAM", "KRT8", "KRT19"),
  Ciliated    = c("FOXJ1", "PIFO"),
  Secretory   = c("PAEP", "CXCL14"),
  Endothelial = c("PECAM1", "VWF"),
  Perivascular= c("RGS5", "MCAM"),
  NK_T        = c("PTPRC", "CD3D", "NCAM1", "KLRD1", "NKG7", "GZMB"),
  Myeloid     = c("CD68", "CD14", "LYZ", "ITGAX"),
  B_plasma    = c("MS4A1", "IGHG1"),
  Mast        = c("TPSB2", "KIT")
)
avg <- AverageExpression(merged, assays = "RNA", slot = "data",
                         group.by = "seurat_clusters")$RNA
score_tab <- sapply(markers, function(g) {
  g <- intersect(g, rownames(avg))
  if (!length(g)) return(rep(NA, ncol(avg)))
  colMeans(avg[g, , drop = FALSE])
})
annot <- data.frame(cluster = colnames(avg),
                    assigned = colnames(score_tab)[max.col(score_tab, ties.method = "first")],
                    max_score = apply(score_tab, 1, max, na.rm = TRUE),
                    second_score = apply(score_tab, 1, function(x) sort(x, decreasing = TRUE)[2]),
                    stringsAsFactors = FALSE)
annot$margin <- round(annot$max_score - annot$second_score, 3)
write.csv(annot, "results/S2_annotation_validation.csv", row.names = FALSE)
print(annot)

# dotplot evidence figure
allm <- unique(unlist(markers))
p4 <- DotPlot(merged, features = allm) + RotatedAxis() + ggtitle("canonical markers")
ggsave("figures/S2_umap_markers_dotplot.png", p4, width = 16, height = 8)

saveRDS(merged, "data/processed/integrated_S2.rds")
cat("\nS2_02 done. Next: verify gate criteria in logs/S2_preregistered_design.md §5,\n")
cat("then write logs/S2_verification_report.md before S3.\n")
