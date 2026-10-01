# ============================================================================
# S2_02_integrate_part1.R — v2.5a (Amendment S2_02-B.6, 2026-09-16)
# PART 1 of 2: everything UP TO the pre-Harmony checkpoint (PASS A, PASS B,
# assembly, global ScaleData, DietSeurat, PCA, Assay5 panel replacement,
# checkpoint save) — then STOPS gracefully.
#
# Why the split (measured, R9): harmony 2.0.5 HUNG twice inside the long-lived
# process that had churned ~29 GB (ScaleData/PCA + paging), yet the EXACT
# registered Harmony call COMPLETED in 10.3 min on all 233,459 cells when run
# in a FRESH process from the checkpoint (S2_02_harmony_diag.txt, 2026-09-16).
# Conclusion: the hang is process-context-dependent (post-heavy-phase Windows
# memory state), not a harmony fault. The fix is execution topology ONLY:
# Harmony and everything downstream run in a fresh process (part 2).
# Zero parameter/gate/output change vs v2.4. Run part 2 only after this prints
# "PART 1 DONE".
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

say("PART 1 DONE — checkpoint at data/processed/tmp_S2_02/checkpoint_pre_harmony.rds")
say("Now run: Rscript scripts\\S2_02_integrate_part2.R > logs\\S2_02_part2.txt 2>&1")
