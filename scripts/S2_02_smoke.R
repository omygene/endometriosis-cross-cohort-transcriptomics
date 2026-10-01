# ============================================================================
# S2_02_smoke.R — 10-15 min pre-flight for S2_02_integrate.R v2.3 (Addendum
# S2_02-B.4, 2026-09-16). Exercises EVERY fragile code path of v2.3 on the
# smallest dataset (GSE183837, 40,801 cells) so API faults surface in MINUTES,
# not after an overnight run. This is a DIAGNOSTIC script (R9 culture):
#   * writes ONLY to data/processed/tmp_S2_02_smoke/ — touches no real output
#     in results/, figures/, or data/processed/ (R4: nothing here is a result)
#   * produces NO numbers used anywhere in the project (R1)
#   * the full S2_02 run and its registration are UNCHANGED
# One documented difference by construction: a single-dataset smoke gives the
# "dataset" covariate ONE level, so Harmony here is tested with
# group.by.vars = "sample" only — the test targets the RunHarmony INTERFACE
# (object plumbing), not the two-covariate call, which is unchanged in v2.3.
# Verdict line at the end: "SMOKE OK" = every fragile path executed.
# ============================================================================
set.seed(42)
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}

req <- c("Seurat", "SeuratObject", "harmony", "RANN", "ggplot2", "dplyr", "patchwork")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("MISSING PACKAGES: ", paste(miss, collapse = ", "))
stopifnot(requireNamespace("lisi", quietly = TRUE))   # engine is registered present
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})

outdir <- "data/processed/tmp_S2_02_smoke"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Same panels as v2.3 (annotation path test needs them)
PANELS <- list(
  Stromal    = c("PDGFRB", "DCN", "COL1A1"),
  Epithelial = c("EPCAM", "KRT8", "KRT19"),
  Endothelial = c("PECAM1", "VWF"),
  Perivascular = c("RGS5", "MCAM"),
  NKT        = c("PTPRC", "CD3D", "NCAM1", "KLRD1", "NKG7", "GZMB"),
  Myeloid    = c("CD68", "CD14", "LYZ", "ITGAX"),
  B_plasma   = c("MS4A1", "IGHG1"),
  Mast       = c("TPSB2", "KIT"))
SUB_EPI <- list(Ciliated = c("FOXJ1", "PIFO"), Secretory = c("PAEP", "CXCL14"))
panel_all <- unique(c(unlist(PANELS), unlist(SUB_EPI)))

align_rows <- function(m, genes) {
  miss <- setdiff(genes, rownames(m))
  if (length(miss)) {
    pad <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                dims = c(length(miss), ncol(m)),
                                dimnames = list(miss, colnames(m)))
    m <- rbind(m, pad)
  }
  m[genes, , drop = FALSE]
}

g <- "GSE183837"
o <- readRDS(file.path("data/processed", paste0(g, "_filtered.rds")))
o$dataset <- g
o$sample  <- o$orig.ident
o <- NormalizeData(o, verbose = FALSE)
o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = 2000,
                          verbose = FALSE)
union_hvgs <- VariableFeatures(o)          # smoke stand-in for the 3-dataset union
feat_set <- c(union_hvgs, setdiff(panel_all, union_hvgs))
say("smoke base: ", g, " | cells ", ncol(o), " | feat_set ", length(feat_set))

# --- v2.1 mechanism: multi-layer extraction ---------------------------------
new_names <- paste0(o$sample, "__", colnames(o))
stopifnot(!anyDuplicated(new_names))
o <- RenameCells(o, new.names = new_names)
lyrs <- Layers(o, search = "^data\\.")
nd <- do.call(cbind, lapply(lyrs, function(ly)
  align_rows(LayerData(o, layer = ly), feat_set)))
stopifnot(identical(colnames(nd), rownames(o@meta.data)))
say("PASS B mechanism OK: ", length(lyrs), " layers -> ", nrow(nd), " x ", ncol(nd))

# --- PASS C mechanisms: assemble, global ScaleData, DietSeurat, PCA ----------
merged <- CreateSeuratObject(counts = nd, meta.data = o@meta.data,
                             min.cells = 0, min.features = 0)
LayerData(merged, layer = "data") <- nd
VariableFeatures(merged) <- union_hvgs
merged <- ScaleData(merged, features = union_hvgs, verbose = FALSE)
say("ScaleData (global) OK")

panel_present <- intersect(panel_all, rownames(merged))
pd <- LayerData(merged, layer = "data")[panel_present, , drop = FALSE]
merged <- DietSeurat(merged, counts = FALSE, data = FALSE, scale.data = TRUE)
merged <- RunPCA(merged, npcs = 50, seed.use = 42, verbose = FALSE)
say("DietSeurat + PCA OK")

# --- v2.3 mechanism: Assay5 panel replacement + checkpoint round-trip --------
merged[["RNA"]] <- CreateSeuratObject(counts = pd, min.cells = 0,
                                      min.features = 0)[["RNA"]]
stopifnot(is(merged[["RNA"]], "Assay5"))
say("Assay5 panel replacement OK")
saveRDS(merged, file.path(outdir, "checkpoint_pre_harmony.rds"))
chk <- readRDS(file.path(outdir, "checkpoint_pre_harmony.rds"))
stopifnot(ncol(chk) == ncol(merged)); rm(chk); gc(verbose = FALSE)
say("checkpoint write/read OK")

# --- RunHarmony interface (single-covariate here — see header) ---------------
t0 <- Sys.time()
merged <- harmony::RunHarmony(merged, group.by.vars = "sample",
                              theta = 2, reduction = "pca",
                              dims.use = 1:30, reduction.save = "harmony",
                              project.dim = FALSE)   # matches v2.4 exactly
say("RunHarmony OK in ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# --- neighbors / clusters / UMAP ---------------------------------------------
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, random.seed = 42, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, seed.use = 42,
                  verbose = FALSE)
say("neighbors/clusters/UMAP OK: ", length(levels(merged$seurat_clusters)),
    " clusters")

# --- LISI engine (lisi package) ----------------------------------------------
emb  <- Embeddings(merged, "harmony")[, 1:30]
meta <- merged@meta.data
lx <- lisi::compute_lisi(emb, meta, c("sample", "dataset"))
rownames(lx) <- rownames(meta)
merged <- AddMetaData(merged, lx)
say("lisi package OK")

# --- annotation input path: AverageExpression(layer = "counts") on Assay5 ----
avg <- as.matrix(AverageExpression(merged, features = panel_present,
                                   group.by = "seurat_clusters",
                                   layer = "counts", verbose = FALSE)$RNA)
stopifnot(nrow(avg) > 0, ncol(avg) == length(levels(merged$seurat_clusters)))
say("AverageExpression(layer = 'counts') OK: ", nrow(avg), " genes x ",
    ncol(avg), " clusters")

# --- figures path -------------------------------------------------------------
p <- DimPlot(merged, reduction = "umap", raster = TRUE, raster.dpi = c(300, 300),
             label = TRUE) + ggtitle("S2_02 smoke")
ggsave(file.path(outdir, "smoke_umap.png"), p, width = 8.5, height = 7,
       dpi = 300, bg = "white")
ggsave(file.path(outdir, "smoke_umap.pdf"), p, width = 8.5, height = 7,
       bg = "white", device = cairo_pdf)
say("figures OK")

say("SMOKE OK — every fragile v2.3 path executed on real data. ",
    "The overnight run now carries SCALE risk only, not API risk.")
