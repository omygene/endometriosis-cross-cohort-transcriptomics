# ============================================================================
# S2_02_integrate.R — v2.1 (Amendment S2_02-B, 2026-09-16): RAM-isolated staged
# execution of the SAME pre-registered integration math.
# v2.1: multi-layer-aware PASS B (the frozen S2 rds carry one data layer per
# sample — extract per layer + zero-fill + cbind, mirroring StitchMatrix).
# v2.2: post-PCA scale.data drop via direct assay replacement — DietSeurat
# with all-flags-FALSE aborts ('x' cannot be empty; keeping zero layers is
# not allowed in Seurat 5.5.1).
# v2.3 (Addendum S2_02-B.3): the v2.2 run HUNG at Harmony 9/10 (measured on
# the lead's machine: Rscript.exe 0% CPU, 56 KB active working set, log frozen
# >1 h, no prompt returned, Enter did not unblock) — a hung process is
# indistinguishable from a slow one because Windows block-buffers redirected
# output and the script had no checkpoints. v2.3 changes OBSERVABILITY ONLY:
#   1) say(): timestamped + flush.console() logging after every phase line;
#   2) checkpoint rds files (pre-Harmony, post-UMAP) in tmp_S2_02 — a future
#      failure loses minutes, not hours (checkpoints auto-deleted on success);
#   3) Harmony wrapped in a timer so a hang localizes to the exact iteration;
#   4) the panel assay is built via the CreateSeuratObject path so it stays
#      Assay5 (v2.2's CreateAssayObject produced a v3 Assay — unproven for the
#      AverageExpression(layer=) call at the annotation step).
# v2.4 (Addendum S2_02-B.5): RunHarmony gains project.dim = FALSE — its TRUE
# default runs Seurat::ProjectDim AFTER convergence and needs the dropped
# scale.data layer (caught by the smoke test); the corrected embeddings are
# identical either way (source-verified) and no downstream step uses loadings.
# No parameter, gate, or output changed in v2.1/v2.2/v2.3 — execution only.
# Source-verified identity note (locked Seurat 5.5.1, preprocessing5.R):
# ScaleData.StdAssay defaults to by.layer = FALSE → it stitches all data
# layers into ONE matrix and scales GLOBALLY — exactly what PASS C does. The
# v1 math = global scaling; per-layer scaling was never in play.
#
# Why v2 exists (measured, not assumed):
#   v1 (MD5 46bcef38102a014374ee911dc03e567d) was OOM-killed by Windows at
#   Harmony iteration 7/10. Run log evidence: Vcells max used = 42,962.9 MB
#   (gc telemetry) vs 16 GB physical RAM; log ends mid-Harmony with NO Error
#   line -> external OS kill. Root cause: v1 kept all three full Seurat
#   objects (37,139 genes, counts+data layers) co-resident with the dense
#   scale.data (233,459 x 3,193 = ~6 GB) through ScaleData/PCA/Harmony.
#
# What v2 changes: execution topology ONLY (three memory-isolated passes).
#   PASS A - per dataset, sequential: readRDS -> NormalizeData ->
#            FindVariableFeatures(vst, 2000) -> keep ONLY the HVG vector ->
#            rm + gc. (v1 kept all three objects alive.)
#   PASS B - per dataset, sequential: readRDS -> deterministic cell renaming
#            (sample__barcode, hard uniqueness check; replaces v1's SILENT
#            merge() renaming) -> NormalizeData -> extract the normalized
#            data matrix restricted to (union HVGs + annotation panel genes)
#            as a SPARSE matrix + meta.data -> save slim rds -> rm + gc.
#   PASS C - assemble the three sparse matrices into ONE genes x cells matrix
#            (zero-filled rows for genes absent from a dataset = exactly what
#            merge() did in v1) -> build minimal object -> ScaleData(features
#            = union_hvgs) on ALL 233,459 cells (per-gene mean/sd across ALL
#            cells = numerically identical to v1's merged ScaleData, because
#            NormalizeData is per-cell and dataset-independent) -> drop sparse
#            layers BEFORE RunPCA -> RunPCA -> drop scale.data immediately
#            after -> Harmony / neighbors / clusters / UMAP / LISI /
#            annotation / figures = v1 code unchanged.
#
# What is FROZEN (identical to the registration and to v1): HVGs 2000 vst per
# dataset pre-merge; union-HVG rule; npcs=50, RunPCA seed.use=42; Harmony
# theta=c(2,2), covariates sample+dataset, dims 1:30; neighbors/clusters on
# harmony 1:30, resolution 0.5, random.seed=42; RunUMAP seed.use=42; LISI
# gates; exclusive-cluster >95% rule; annotation panels + 1.25x / >=90% rules;
# all output paths; set.seed(42) global.
#
# Deliverable object note: S2_02_integrated.rds is now LIGHT — its RNA assay
# holds ONLY the ~40 annotation-panel genes, and the counts layer contains
# their LOG-NORMALIZED values (identical numbers to v1's data layer; stored
# in counts because the data/scale.data layers are dropped for RAM). It is
# NOT raw counts. Genome-wide matrices remain in the frozen per-dataset S2
# rds; S3+ joins them by results/S2_02_cell_index.csv (cell -> dataset,
# sample). AverageExpression below uses layer="counts" — the IDENTICAL
# statistic to v1's layer="data" on identical values.
#
# Peak RAM estimate: PASS A ~8 GB, PASS B ~9 GB, PASS C ~12 GB transient
# (vs 43 GB measured in v1). Runtime ~2-3 h. Close heavy apps anyway.
# If PASS C is still killed -> STOP and report (R6); the sketch fallback is
# a METHOD change and requires its own Amendment (S2_02-C) — NOT automatic.
# ============================================================================
set.seed(42)

# v2.3: timestamped, flushed logging — Windows block-buffers redirected stdout,
# which made a HUNG process indistinguishable from a silent one on 2026-09-16.
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}

# ---- 0. package pre-flight (registered §6): fail in seconds, not at 2 AM ----
req <- c("Seurat", "SeuratObject", "harmony", "RANN", "ggplot2", "dplyr", "patchwork")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss))
  stop("MISSING PACKAGES: ", paste(miss, collapse = ", "),
       "\nRun: renv::install(c(",
       paste(sprintf('"%s"', miss), collapse = ", "),
       ")) then re-run S2_02. Nothing else was executed.")
# Amendment S2_02-A (2026-09-14): lisi is archived off CRAN (verified) and its
# GitHub source needs Rtools (C++). If installed we use the package; otherwise
# the faithful pure-R port below. The engine used is printed and logged.
# RESOLVED 2026-09-14: lisi 1.0 present on the machine; package engine runs.
have_lisi <- requireNamespace("lisi", quietly = TRUE)
say("LISI engine: ", if (have_lisi) "lisi package" else
    "internal pure-R port (Amendment S2_02-A)")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

gses <- c("GSE179640", "GSE214411", "GSE183837")
tmpdir <- "data/processed/tmp_S2_02"
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

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

# ---- 0c. annotation panels (design §4 verbatim — needed in PASS B) ----------
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

# ---- PASS A. per-dataset HVGs (design §3: BEFORE merge) — memory-isolated ---
hvgs <- list()
for (g in gses) {
  f <- file.path("data/processed", paste0(g, "_filtered.rds"))
  if (!file.exists(f)) stop("missing ", f, " — S2 must be complete first")
  o <- readRDS(f)
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
  hvgs[[g]] <- VariableFeatures(o)
  say(sprintf("PASS A %s: cells %d | HVGs %d", g, ncol(o), length(hvgs[[g]])))
  rm(o); gc(verbose = FALSE)
}
union_hvgs <- unique(unlist(hvgs))          # deterministic: dataset order kept
writeLines(c(sprintf("per-dataset HVG counts: %s",
                     paste(sprintf("%s=%d", gses, lengths(hvgs)), collapse = ", ")),
             sprintf("union HVG count: %d", length(union_hvgs))),
           "results/S2_02_hvg_union.txt")
say(sprintf("union HVGs: %d", length(union_hvgs)))
rm(hvgs); gc(verbose = FALSE)

# ---- PASS B. slim extraction: normalized values for union HVGs + panels -----
# feat_set = union HVGs (ScaleData input, identical math) + panel genes (for
# the §4 annotation AFTER scale.data is dropped). ScaleData still sees exactly
# union_hvgs — panel rows are passengers only.
feat_set <- c(union_hvgs, setdiff(panel_all, union_hvgs))

# Zero-fill rows for genes absent from a layer/dataset = merge()/StitchMatrix
# semantics (rowmap zero-fill).
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

for (g in gses) {
  f <- file.path("data/processed", paste0(g, "_filtered.rds"))
  o <- readRDS(f)
  o$dataset <- g
  o$sample  <- o$orig.ident
  # Deterministic cell naming (replaces v1's SILENT merge() renaming, which
  # triggered "Some cell names are duplicated... Renaming"). Hard check per R6.
  new_names <- paste0(o$sample, "__", colnames(o))
  if (anyDuplicated(new_names))
    stop("duplicate cell names after sample__barcode prefixing in ", g,
         " — investigate before proceeding (R6)")
  o <- RenameCells(o, new.names = new_names)
  o <- NormalizeData(o, verbose = FALSE)
  # v2.1 fix (2026-09-16): the frozen S2 rds are MULTI-LAYER — one data layer
  # per sample (S2_01 merged without JoinLayers; 33+13+9 layers). Calling
  # LayerData(layer="data") on >1 layers returns ONE layer with a warning and
  # the narrower per-layer feature set broke subsetting ("subscript out of
  # bounds"). Fix: extract each layer explicitly, zero-fill, cbind — the same
  # mechanism ScaleData.StdAssay itself uses (StitchMatrix, by.layer=FALSE
  # default in locked Seurat 5.5.1 — source-verified).
  lyrs <- Layers(o, search = "^data\\.")   # anchored: never matches scale.data
  nd <- do.call(cbind, lapply(lyrs, function(ly)
    align_rows(LayerData(o, layer = ly), feat_set)))
  stopifnot(identical(colnames(nd), rownames(o@meta.data)))
  say(sprintf("PASS B %s: %d data layers | extracted %d features x %d cells (sparse)",
              g, length(lyrs), nrow(nd), ncol(nd)))
  saveRDS(list(norm = nd, meta = o@meta.data),
          file.path(tmpdir, paste0(g, "_slim.rds")))
  rm(o, nd); gc(verbose = FALSE)
}

# ---- PASS C. assembly -> ScaleData -> PCA -> Harmony (design §3 verbatim) ---
slims <- lapply(gses, function(g)
  readRDS(file.path(tmpdir, paste0(g, "_slim.rds"))))
names(slims) <- gses

# Rows are already aligned to feat_set in PASS B; verify, then cbind.
stopifnot(all(vapply(slims, function(s) identical(rownames(s$norm), feat_set),
                     logical(1))))
norm_all <- do.call(cbind, lapply(slims, function(s) s$norm))
meta_all <- do.call(rbind, lapply(slims, function(s) s$meta))
rm(slims); gc(verbose = FALSE)
stopifnot(nrow(meta_all) == ncol(norm_all))
if (anyDuplicated(colnames(norm_all)))
  stop("duplicate cell names in assembled matrix — investigate (R6)")

# Traceability: cell index for joining genome-wide data in later stages.
cell_index <- data.frame(cell    = colnames(norm_all),
                         dataset = meta_all$dataset,
                         sample  = meta_all$sample,
                         stringsAsFactors = FALSE)
write.csv(cell_index, "results/S2_02_cell_index.csv", row.names = FALSE)

# ---- condition labels (same frozen maps as v1 / S2_04 / S2_05) --------------
parse_cond <- function(sid, gse) {
  if (gse == "GSE179640") {
    if (grepl("_Ctrl", sid)) return("Ctrl")
    if (grepl("_EuE",  sid)) return("EuE")
    if (grepl("_EcPA", sid)) return("EcPA")
    if (grepl("_EcP",  sid)) return("EcP")
    if (grepl("_EcO",  sid)) return("EcO")
    if (grepl("EOR",   sid)) return("EOR")
    return("other")
  }
  if (gse == "GSE214411") return(if (grepl("EMS", sid)) "EMS" else "N")
  if (grepl("RIF", sid, ignore.case = TRUE)) return("RIF")
  "Control"
}
COND_COLORS <- c(Ctrl = "#4DBBD5", EuE = "#E64B35", EcP = "#00A087",
                 EcPA = "#F39B7F", EcO = "#3C5488", EOR = "#8491B4",
                 EMS = "#DC0000", N = "#00BFC4", Control = "#4DBBD5",
                 RIF = "#7E6148")

# min.cells=0 / min.features=0: mandatory — default min.features=200 would
# drop every cell (the assembly carries ~3.2k features, but the rule stands).
merged <- CreateSeuratObject(counts = norm_all, meta.data = meta_all,
                             min.cells = 0, min.features = 0)
LayerData(merged, layer = "data") <- norm_all   # already log-normalized
rm(norm_all); gc(verbose = FALSE)
merged$condition <- mapply(parse_cond, merged$sample, merged$dataset)
VariableFeatures(merged) <- union_hvgs
say(sprintf("assembled: %d cells x %d features", ncol(merged), length(feat_set)))

# ScaleData on ALL 233,459 cells: per-gene mean/sd are GLOBAL here, numerically
# identical to v1's ScaleData on the merged object (NormalizeData is per-cell,
# merge zero-fills absent genes — both reproduced exactly above).
merged <- ScaleData(merged, features = union_hvgs, verbose = FALSE)
say("ScaleData done (global)"); gc(verbose = FALSE)

# Rescue panel values (already log-normalized) BEFORE dropping sparse layers.
panel_present <- intersect(panel_all, rownames(merged))
pd <- LayerData(merged, layer = "data")[panel_present, , drop = FALSE]

# Drop the sparse counts/data layers (multi-GB) before RunPCA's dense work.
merged <- DietSeurat(merged, counts = FALSE, data = FALSE, scale.data = TRUE)
gc(verbose = FALSE)

merged <- RunPCA(merged, npcs = 50, seed.use = 42, verbose = FALSE)
say("PCA done")

# Drop scale.data (~6 GB dense) immediately after PCA — nothing downstream
# needs it (Harmony/neighbors/UMAP read the pca/harmony embeddings only).
# v2.2 fix: DietSeurat(counts=FALSE, data=FALSE, scale.data=FALSE) aborts in
# locked Seurat 5.5.1 (".PropagateList: 'x' cannot be empty" — keeping zero
# layers is not allowed; objects.R DietSeurat, source-verified). Replacing the
# whole assay discards scale.data cleanly.
# v2.3: build the replacement via the CreateSeuratObject path so the panel
# assay stays Assay5 (v2.2's bare CreateAssayObject produced a v3 Assay —
# warned "Assay RNA changing from Assay5 to Assay" — unproven for the
# AverageExpression(layer=) call at annotation).
# Light RNA assay back: panel genes only, counts layer = their LOG-NORMALIZED
# values (identical numbers to v1's data layer — see header note).
merged[["RNA"]] <- CreateSeuratObject(counts = pd, min.cells = 0,
                                      min.features = 0)[["RNA"]]
rm(pd); gc(verbose = FALSE)

# CHECKPOINT 1 (v2.3): light object (panel assay + pca + metadata) — if
# anything after this point fails, resume material exists; a hang localizes
# instead of costing the full run.
saveRDS(merged, file.path(tmpdir, "checkpoint_pre_harmony.rds"))
say("checkpoint saved: pre_harmony (",
    file.size(file.path(tmpdir, "checkpoint_pre_harmony.rds")) %/% 1e6, " MB)")

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

# ---- 4. neighbors + clusters (res 0.5, frozen) + UMAP ------------------------
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, random.seed = 42, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, seed.use = 42,
                  verbose = FALSE)
say("clustering + UMAP done: ", length(levels(merged$seurat_clusters)), " clusters")

# CHECKPOINT 2 (v2.3): post-UMAP light object (embeddings + clusters + meta).
saveRDS(merged, file.path(tmpdir, "checkpoint_post_umap.rds"))
say("checkpoint saved: post_umap")

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

# ---- 8. iLISI on lineage (gate §5) -------------------------------------------
dfl <- data.frame(lineage = merged$lineage); rownames(dfl) <- rownames(meta)
li <- compute_lisi_any(emb, dfl, "lineage")
colnames(li) <- "iLISI_lineage"
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
