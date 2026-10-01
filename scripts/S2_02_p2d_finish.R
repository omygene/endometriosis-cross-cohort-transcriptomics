# ============================================================================
# S2_02_p2d_finish.R — step 2d/4: LISI gates, exclusive clusters, annotation, figures, save, gate summary.
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
req <- c("Seurat", "RANN", "ggplot2", "dplyr", "patchwork")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("MISSING PACKAGES: ", paste(miss, collapse = ", "))
have_lisi <- requireNamespace("lisi", quietly = TRUE)
say("LISI engine: ", if (have_lisi) "lisi package" else
    "internal pure-R port (Amendment S2_02-A)")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
tmpdir <- "data/processed/tmp_S2_02"
load_cp <- function(name) {
  f <- file.path(tmpdir, paste0("checkpoint_", name, ".rds"))
  if (!file.exists(f)) stop("missing checkpoint: ", f, " — run the previous step first")
  obj <- readRDS(f)
  say("loaded ", name, ": ", ncol(obj), " cells")
  if (ncol(obj) != 233459) stop("cell count mismatch — investigate (R6)")
  obj
}

# ---- 0b. LISI engines (Amendment S2_02-A) — verbatim from v1 ----------------
compute_lisi_r <- function(emb, labs, perplexity = 30) {
  n <- nrow(emb); k <- min(3 * perplexity, n - 1)
  nn  <- RANN::nn2(emb, k = k + 1)              # +1: self is dropped next
  idx <- nn$nn.idx[, -1, drop = FALSE]
  D   <- nn$nn.dists[, -1, drop = FALSE]
  logU <- log(perplexity)
  lo <- rep(1e-12, n); hi <- rep(Inf, n); sg <- rep(1, n)
  for (it in seq_len(64)) {
    P <- exp(-sweep(D^2, 1, 2 * sg^2, "/"))
    rs <- rowSums(P); rs[rs == 0] <- 1
    P  <- P / rs
    H  <- -rowSums(P * log(P + 1e-12))
    big <- H > logU                              # entropy too high -> sigma too big
    hi[big] <- sg[big]; lo[!big] <- sg[!big]
    sg <- ifelse(is.infinite(hi), sg * 2, (lo + hi) / 2)
  }
  P <- exp(-sweep(D^2, 1, 2 * sg^2, "/"))        # final calibrated weights
  P <- P / rowSums(P)
  lv <- levels(factor(labs)); M <- matrix(match(labs, lv)[idx], n, k)
  s2 <- numeric(n)
  for (l in seq_along(lv)) { fr <- rowSums(P * (M == l)); s2 <- s2 + fr^2 }
  1 / s2
}
compute_lisi_any <- function(emb, mdf, cols) {
  if (have_lisi) {
    out <- lisi::compute_lisi(emb, mdf, cols)
  } else {
    out <- as.data.frame(setNames(
      lapply(cols, function(cc) compute_lisi_r(emb, mdf[[cc]])), cols))
  }
  rownames(out) <- rownames(mdf)
  out
}

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

# condition was already computed in part 1 and lives in the checkpoint meta.
COND_COLORS <- c(Ctrl = "#4DBBD5", EuE = "#E64B35", EcP = "#00A087",
                 EcPA = "#F39B7F", EcO = "#3C5488", EOR = "#8491B4",
                 EMS = "#DC0000", N = "#00BFC4", Control = "#4DBBD5",
                 RIF = "#7E6148")


merged <- load_cp("post_umap")

# ---- 5. LISI gate (design §5): cLISI sample/dataset --------------------------
emb  <- Embeddings(merged, "harmony")[, 1:30]
meta <- merged@meta.data
lx   <- compute_lisi_any(emb, meta, c("sample", "dataset"))
colnames(lx) <- paste0("cLISI_", colnames(lx))
merged <- AddMetaData(merged, lx)
say("cLISI done")

# ---- 6. exclusive-cluster rule (registered: >95% one sample = exclusive) -----
tab  <- table(cluster = merged$seurat_clusters, sample = merged$sample)
frac <- prop.table(tab, margin = 1)
exc  <- data.frame(
  cluster    = rownames(frac),
  n_cells    = as.integer(rowSums(tab)),
  top_sample = colnames(frac)[apply(frac, 1, which.max)],
  top_frac   = round(apply(frac, 1, max), 4))
exc$exclusive <- exc$top_frac > 0.95
write.csv(exc, "results/S2_02_exclusive_clusters.csv", row.names = FALSE)
say("exclusive clusters (>95% one sample): ", sum(exc$exclusive), " of ",
    nrow(exc))

# ---- 7. annotation by fixed panels (design §4 verbatim) ----------------------
# AverageExpression with layer="counts": the counts layer holds the SAME
# log-normalized values v1 held in "data" (Amendment S2_02-B) -> identical
# statistic. Seurat 5 AverageExpression already averages in non-log space and
# returns log-scale values — do NOT re-log (a second log1p would corrupt the
# 1.25x ratio).
present <- intersect(panel_all, rownames(merged))
avg <- AverageExpression(merged, features = present,
                         group.by = "seurat_clusters", layer = "counts",
                         verbose = FALSE)$RNA
avg <- as.matrix(avg)
# Seurat 5.5.1 prepends "g" to numeric group names ("0" -> "g0"; source-verified
# in locked utilities.R: CreateCategoryMatrix caller, "appending `g` to ensure
# valid variable names"). Category-matrix columns follow FACTOR LEVEL ORDER
# (sparse.model.matrix ~0+factor), so map positionally back to the true cluster
# ids — hard stop if the shape disagrees (R6: fail loud, never propagate NA).
cl_levels <- levels(merged$seurat_clusters)
if (ncol(avg) != length(cl_levels)) {
  stop("AverageExpression columns (", ncol(avg),
       ") != nlevels(seurat_clusters) (", length(cl_levels),
       ") — investigate (R6)")
}
colnames(avg) <- cl_levels
panel_score <- function(gs) colMeans(avg[intersect(gs, rownames(avg)), , drop = FALSE])
S  <- sapply(PANELS, panel_score)                # clusters x panels
Ss <- sapply(SUB_EPI, panel_score)
ord  <- order(as.integer(colnames(avg)))         # cluster order by number
S <- S[ord, ]; Ss <- Ss[ord, ]
top_i  <- apply(S, 1, which.max)
top_v  <- apply(S, 1, max)
S2nd   <- apply(S, 1, function(x) sort(x, decreasing = TRUE)[2])
val <- data.frame(
  cluster      = rownames(S),
  n_cells      = as.integer(table(merged$seurat_clusters)[rownames(S)]),
  lineage      = names(PANELS)[top_i],
  top_panel_mean    = round(top_v, 4),
  runner_up_mean    = round(S2nd, 4),
  ratio        = round(top_v / pmax(S2nd, 1e-4), 3),
  ciliated_score  = round(Ss[, "Ciliated"], 4),
  secretory_score = round(Ss[, "Secretory"], 4),
  round(as.data.frame(S), 4), check.names = FALSE)
val$support <- val$ratio >= 1.25                 # registered mechanics §3
write.csv(val, "results/S2_02_annotation_validation.csv", row.names = FALSE)
pct_ok <- mean(val$support) * 100
say(sprintf("annotation support: %d/%d clusters (%.0f%%) — >=90%% rule: %s",
            sum(val$support), nrow(val), pct_ok,
            if (pct_ok >= 90) "PASS" else "FAIL — document per §4"))
merged$lineage <- factor(val$lineage[match(as.character(merged$seurat_clusters),
                                           val$cluster)],
                         levels = names(PANELS))
if (any(is.na(merged$lineage))) {
  stop("lineage mapping produced NA for ",
       sum(is.na(merged$lineage)), " cells — investigate (R6)")
}

# ---- 8. iLISI on lineage (gate §5) -------------------------------------------
dfl <- data.frame(lineage = merged$lineage); rownames(dfl) <- rownames(meta)
li <- compute_lisi_any(emb, dfl, "lineage")
colnames(li) <- "iLISI_lineage"
if (any(!is.finite(li$iLISI_lineage))) {
  stop("iLISI has ", sum(!is.finite(li$iLISI_lineage)),
       " non-finite values — investigate (R6)")
}
merged <- AddMetaData(merged, li)
lsum <- rbind(
  data.frame(metric = "cLISI_sample",  label = "overall",
             median = median(lx$cLISI_sample),  IQR = IQR(lx$cLISI_sample)),
  data.frame(metric = "cLISI_dataset", label = "overall",
             median = median(lx$cLISI_dataset), IQR = IQR(lx$cLISI_dataset)),
  data.frame(metric = "iLISI_lineage", label = "overall",
             median = median(li$iLISI_lineage), IQR = IQR(li$iLISI_lineage)),
  merged@meta.data |>
    group_by(lineage) |>
    summarise(median = median(iLISI_lineage), IQR = IQR(iLISI_lineage),
              .groups = "drop") |>
    mutate(metric = "iLISI_lineage", label = as.character(lineage),
           lineage = NULL) |>
    dplyr::select(metric, label, median, IQR))
write.csv(lsum, "results/S2_02_lisi.csv", row.names = FALSE)
print(lsum)

# ---- 9. figures (publication pattern from S2_04: PNG 300dpi + cairo_pdf) -----
umap_base <- function(...) DimPlot(merged, reduction = "umap", raster = TRUE,
                                   raster.dpi = c(300, 300), ...) +
  theme(plot.title = element_text(face = "bold"))
figs <- list(
  clusters  = umap_base(label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("S2_02 integrated — clusters (Harmony, res 0.5)") + NoLegend(),
  # colors via NAMED scale_color_manual (DimPlot's `cols` maps by level ORDER,
  # not name — alphabetical levels would scramble COND_COLORS)
  dataset   = umap_base(group.by = "dataset") +
    scale_color_manual(values = c(GSE179640 = "#E64B35", GSE183837 = "#4DBBD5",
                                  GSE214411 = "#00A087")) +
    ggtitle("S2_02 integrated — by dataset"),
  condition = umap_base(group.by = "condition") +
    scale_color_manual(values = COND_COLORS) +
    ggtitle("S2_02 integrated — by condition"),
  lineage   = umap_base(group.by = "lineage") +
    ggtitle("S2_02 integrated — major lineage (fixed panels)"))
for (nm in names(figs)) {
  ggsave(file.path("figures", paste0("S2_02_umap_", nm, ".png")),
         figs[[nm]], width = 8.5, height = 7, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0("S2_02_umap_", nm, ".pdf")),
         figs[[nm]], width = 8.5, height = 7, bg = "white", device = cairo_pdf)
  say("figure: ", nm)
}

# ---- 10. save + gate summary --------------------------------------------------
saveRDS(merged, "data/processed/S2_02_integrated.rds")
tmp_ok <- file.remove(list.files(tmpdir, full.names = TRUE))
say("tmp slim files + checkpoints removed: ", all(tmp_ok))
say("")
say("===== S2_02 GATE SUMMARY (design §5) =====")
say("cells: ", ncol(merged), " (registered: 233,459)")
say("clusters: ", nrow(exc), " | exclusive: ", sum(exc$exclusive))
say("cLISI sample median: ", round(median(lx$cLISI_sample), 3),
    " | cLISI dataset median: ", round(median(lx$cLISI_dataset), 3))
say("iLISI lineage median: ", round(median(li$iLISI_lineage), 3))
say(sprintf("annotation >=90%% rule: %s",
            if (pct_ok >= 90) "PASS" else "FAIL"))
say("S2_02 done.")
say("STEP 2d DONE — S2_02 COMPLETE")
