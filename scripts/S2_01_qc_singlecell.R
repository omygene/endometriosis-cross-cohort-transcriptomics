# ============================================================================
# S2_01_qc_singlecell.R — per-dataset single-cell QC  [v4]
# Pre-registered thresholds: logs/S2_preregistered_design.md (Section 2).
# Run from project root. Inputs: data/raw/*  Outputs: results/S2_qc_*,
# figures/S2_*, data/processed/*_filtered.rds
# v4 (2026-09-13): three defects found in the v3 run log (see analysis_log.md):
#   1. dbl_pct stayed NA: `dbl <<-` assigned to the global env, not the local
#      variable -> plain `<-` now, with per-sample removal printed + self-check.
#   2. v3 QC table was byte-identical to pre-v3 (doublet removal did not
#      propagate to the object) -> explicit subset(seu, cells = colnames[keep]).
#   3. GSE214411/GSE183837 silently no-op'd on re-run (staging moved trio files
#      into per-sample dirs on first run; re-runs found nothing) and the
#      EOR02/EOR04 load-exclusion CSV was not regenerated -> staging now
#      idempotent: already-staged dirs are always detected, exclusions are
#      re-recorded every run, empty extraction = LOUD banner.
# v4.1 (2026-09-13): untar guard now also fires when the dataset folder exists
#      but is EMPTY (lead's GSE183837 case — empty dir blocked extraction of
#      the assembled, MD5-verified tar).
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(ggplot2); library(scDblFinder)
  library(SingleCellExperiment)
})
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

# ---- Pre-registered thresholds (DO NOT tune after seeing data) --------------
THRESH <- list(
  GSE179640 = list(min.genes = 200, max.genes = 6000, max.mt = 20, max.hb = 5),
  GSE214411 = list(min.genes = 200, max.genes = 6000, max.mt = 10, max.hb = 5),
  GSE183837 = list(min.genes = 200, max.genes = 6000, max.mt = 20, max.hb = 5)
)

qc_one_sample <- function(counts, sample_id, thr) {
  seu <- CreateSeuratObject(counts, project = sample_id,
                            min.cells = 3, min.features = 0)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  hb <- grep("^HB[AB]", rownames(seu), value = TRUE)
  seu[["percent.hb"]] <- if (length(hb)) PercentageFeatureSet(seu, features = hb) else 0
  n0 <- ncol(seu)
  # scDblFinder (pre-registered doublet caller)
  # v3 fix: scDblFinder 1.26.7 has NO seed argument (verified empirically on the
  # lead's machine: "unused argument (seed = 42)"). Determinism is enforced by
  # set.seed(42) immediately before the call — same seed per sample, so each
  # sample's result is independent of execution order.
  # v4 fix (verified against v3 run log): (a) `dbl <<-` assigned to the GLOBAL
  # env, never to the local `dbl` -> dbl_pct stayed NA; plain `<-` is correct
  # inside the function frame. (b) v3's QC table was byte-identical to the
  # pre-v3 table for all 33 samples -> the singlet subset did not propagate.
  # v4 restructures: explicit keep-vector, canonical subset(cells=...), and a
  # printed per-sample removal line + INCONSISTENCY self-check (silent failure
  # of this step is now impossible to miss).
  dbl <- NA_real_
  sce <- tryCatch({
    s <- as.SingleCellExperiment(seu)
    set.seed(42)
    scDblFinder(s, samples = NULL)
  }, error = function(e) { cat("scDblFinder failed on", sample_id, "- fallback flag\n"); NULL })
  if (!is.null(sce)) {
    keep  <- as.logical(sce$scDblFinder.class == "singlet")
    dbl   <- round(mean(!keep) * 100, 2)
    n_pre <- ncol(seu)
    seu   <- subset(seu, cells = colnames(seu)[keep])
    cat("  scDblFinder:", sum(!keep), "doublets removed (", dbl, "%), cells ",
        n_pre, " -> ", ncol(seu), "\n", sep = "")
    if (ncol(seu) != sum(keep) || ncol(seu) == n_pre)
      cat("  !! DOUBLET REMOVAL INCONSISTENT for ", sample_id,
          " — STOP and report before proceeding\n")
  }
  seu <- subset(seu, subset = nFeature_RNA >= thr$min.genes &
                          nFeature_RNA <= thr$max.genes &
                          percent.mt <= thr$max.mt & percent.hb <= thr$max.hb)
  data.frame(sample = sample_id, cells_before = n0, cells_after = ncol(seu),
             pct_kept = round(100 * ncol(seu) / n0, 1),
             dbl_pct = round(dbl, 2),
             median_genes = median(seu$nFeature_RNA),
             median_mt = round(median(seu$percent.mt), 2),
             stringsAsFactors = FALSE) -> stats
  list(obj = seu, stats = stats)
}

run_dataset <- function(gse, files, thr) {
  cat("\n##### ", gse, " #####\n")
  objs <- list(); all_stats <- list()
  for (f in files) {
    sid <- sub("\\.(h5|mtx.*)$", "", basename(f))
    cat("loading", sid, "\n")
    counts <- if (grepl("\\.h5$", f)) Read10X_h5(f) else Read10X(dirname(f))
    r <- qc_one_sample(counts, sid, thr)
    objs[[sid]] <- r$obj; all_stats[[sid]] <- r$stats
  }
  st <- do.call(rbind, all_stats)
  write.csv(st, file.path("results", paste0("S2_qc_", gse, ".csv")), row.names = FALSE)
  merged <- if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
  p <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3, pt.size = 0, group.by = "orig.ident") + NoLegend()
  ggsave(file.path("figures", paste0("S2_qc_", gse, ".png")), p, width = 14, height = 6)
  saveRDS(merged, file.path("data/processed", paste0(gse, "_filtered.rds")))
  print(st)
  invisible(merged)
}

# ---- GSE179640 (33 H5 + 2 mtx-trio samples; 24 bulk files recorded only) ----
# Amendment S2-A (2026-09-12, logs/S2_preregistered_design.md): GEO official
# filelist = 63 files / 59 samples = 35 single-cell + 24 bulk (bulk NOT analysed).
tar640  <- "data/raw/GSE179640_RAW.tar"
exd640  <- "data/raw/GSE179640"
if (!dir.exists(exd640) && file.exists(tar640))
  untar(tar640, exdir = exd640, tar = "internal")   # internal reader: bsdtar failed on some archives
if (dir.exists(exd640)) {
  cat("\n##### GSE179640 #####\n")
  n_bulk <- length(list.files(exd640, pattern = "bulk\\.featurecounts\\.txt\\.gz$"))
  cat("bulk featurecounts files present (recorded, NOT analysed per Amendment S2-A):", n_bulk, "\n")

  # stage the two 10x mtx-trio samples (EOR02, EOR04)
  # v4 fix: staging is now idempotent. v1-v3 moved trio files into per-sample
  # dirs on the FIRST run, so on re-runs list.files at top level found nothing
  # and the pre-registered EOR02/EOR04 load-time exclusions were NOT re-recorded
  # (results/S2_GSE179640_load_exclusions.csv silently missing after v3 run).
  exclusions <- data.frame(sample = character(), reason = character(), stringsAsFactors = FALSE)
  trio_ids <- unique(sub("_(barcodes|features|matrix)\\..*$", "",
                         list.files(exd640, pattern = "_(barcodes|features|matrix)\\.(tsv|mtx)\\.gz$")))
  for (s in trio_ids) {                      # move top-level trios into dirs (no-op on re-run)
    d <- file.path(exd640, s); dir.create(d, showWarnings = FALSE)
    for (ext in c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")) {
      src <- file.path(exd640, paste0(s, "_", ext))
      if (file.exists(src)) file.rename(src, file.path(d, ext))
    }
  }
  trio_dirs <- c()
  for (d in list.dirs(exd640, recursive = FALSE, full.names = TRUE)) {
    if (!file.exists(file.path(d, "matrix.mtx.gz"))) next
    fsz <- file.size(file.path(d, "features.tsv.gz"))
    if (is.na(fsz) || fsz < 1000) {          # pre-registered load-time exclusion (Amendment S2-A #3)
      cat("LOAD-TIME EXCLUSION:", basename(d), "- features.tsv.gz", fsz, "bytes (<1000)\n")
      exclusions <- rbind(exclusions, data.frame(sample = basename(d),
        reason = paste0("features.tsv.gz size ", fsz, " B < 1000 B (GEO upload anomaly)")))
      next
    }
    trio_dirs <- c(trio_dirs, d)
  }

  h5 <- list.files(exd640, pattern = "filtered_feature_bc_matrix\\.h5$", full.names = TRUE)
  objs <- list(); all_stats <- list()
  for (f in h5) {
    sid <- sub("\\.h5$", "", basename(f)); cat("loading", sid, "\n")
    r <- qc_one_sample(Read10X_h5(f), sid, THRESH$GSE179640)
    objs[[sid]] <- r$obj; all_stats[[sid]] <- r$stats
  }
  for (d in trio_dirs) {
    sid <- basename(d); cat("loading", sid, "(mtx trio)\n")
    counts <- tryCatch(Read10X(d), error = function(e) NULL)
    if (is.null(counts)) {                   # pre-registered load-time exclusion (Amendment S2-A #3)
      cat("LOAD-TIME EXCLUSION:", sid, "- Read10X failed\n")
      exclusions <- rbind(exclusions, data.frame(sample = sid, reason = "Read10X failed at load"))
      next
    }
    r <- qc_one_sample(counts, sid, THRESH$GSE179640)
    objs[[sid]] <- r$obj; all_stats[[sid]] <- r$stats
  }
  if (nrow(exclusions))
    write.csv(exclusions, "results/S2_GSE179640_load_exclusions.csv", row.names = FALSE)
  if (length(objs)) {
    st <- do.call(rbind, all_stats)
    write.csv(st, "results/S2_qc_GSE179640.csv", row.names = FALSE)
    merged <- if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
    p <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3, pt.size = 0, group.by = "orig.ident") + NoLegend()
    ggsave("figures/S2_qc_GSE179640.png", p, width = 14, height = 6)
    saveRDS(merged, "data/processed/GSE179640_filtered.rds")
    print(st)
  }
}

# ---- GSE214411 + GSE183837 (10x mtx trios after untar) ----------------------
for (gse in c("GSE214411", "GSE183837")) {
  tarfile <- file.path("data/raw", paste0(gse, "_RAW.tar"))
  exdir   <- file.path("data/raw", gse)
  # v4.1: untar also when the folder exists but is EMPTY — the lead's v3 run
  # had an empty data/raw/GSE183837 dir (origin unknown) which defeated the
  # !dir.exists() guard, so the assembled+MD5-verified tar was never extracted
  # and the dataset was silently skipped.
  if ((!dir.exists(exdir) || !length(list.files(exdir))) && file.exists(tarfile))
    untar(tarfile, exdir = exdir, tar = "internal")
  if (!dir.exists(exdir)) {                       # v2: LOUD skip, never silent
    cat("\n##### ", gse, " SKIPPED — no local data (", tarfile, " missing ).",
        "If portal-resident: reassemble split parts, verify MD5 vs manifest, then re-run. #####\n")
    next
  }
  # v4 fix: staging moved trio files into per-sample dirs on the first run, so
  # list.files at top level finds nothing on re-runs -> silent no-op (this is
  # why the v3 run regenerated NOTHING for GSE214411/GSE183837). Now: stage any
  # top-level files (first run), then ALWAYS detect already-staged dirs.
  samples <- unique(sub("_(barcodes|features|genes|matrix)\\..*$", "",
                        list.files(exdir, pattern = "\\.(tsv|mtx)\\.gz$")))
  # stage one directory per sample (10x trio naming convention for Read10X)
  for (s in samples) {
    d <- file.path(exdir, s); dir.create(d, showWarnings = FALSE)
    for (ext in c("barcodes.tsv.gz", "features.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz")) {
      src <- file.path(exdir, paste0(s, "_", ext))
      if (file.exists(src)) file.rename(src, file.path(d, ext))
    }
    if (file.exists(file.path(d, "genes.tsv.gz")) && !file.exists(file.path(d, "features.tsv.gz")))
      file.rename(file.path(d, "genes.tsv.gz"), file.path(d, "features.tsv.gz"))
  }
  dirs <- list.dirs(exdir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "matrix.mtx.gz"))]
  if (!length(dirs)) {                       # v4: LOUD, never silent
    cat("\n##### ", gse, " SKIPPED — ", exdir, " contains no 10x trio dirs.",
        "If the tar was just assembled: check it extracted here and file naming is",
        "*_barcodes.tsv.gz / *_features.tsv.gz / *_matrix.mtx.gz. Report the",
        "`dir` listing before re-running. #####\n")
    next
  }
  {
    cat("\n##### ", gse, " ##### (", length(dirs), " samples )\n")
    objs <- list(); all_stats <- list()
    for (d in dirs) {
      sid <- basename(d); cat("loading", sid, "\n")
      r <- qc_one_sample(Read10X(d), sid, THRESH[[gse]])
      objs[[sid]] <- r$obj; all_stats[[sid]] <- r$stats
    }
    st <- do.call(rbind, all_stats)
    write.csv(st, file.path("results", paste0("S2_qc_", gse, ".csv")), row.names = FALSE)
    merged <- if (length(objs) > 1) merge(objs[[1]], objs[-1]) else objs[[1]]
    p <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3, pt.size = 0, group.by = "orig.ident") + NoLegend()
    ggsave(file.path("figures", paste0("S2_qc_", gse, ".png")), p, width = 14, height = 6)
    saveRDS(merged, file.path("data/processed", paste0(gse, "_filtered.rds")))
    print(st)
  }
}

# ---- GSE213216 audit (author-processed rds; report-only per design §2) ------
cat("\n##### GSE213216 processed-object audit (report-only) #####\n")
expected6 <- c("aux.seurat.shared", "EnEpi", "EnS", "endothelial",
               "epithelial.annotated", "mesenchymal.annotated")
rds_files <- list.files("data/raw/GSE213216_processed", pattern = "\\.rds$", full.names = TRUE)
cat("rds files found:", length(rds_files), "(expected 6 objects)\n")
if (!length(rds_files))
  cat("NOTE: folder empty/missing — place the 6 author-processed rds here. On Windows the\n",
      "file 'aux.seurat.shared.rds' hits the RESERVED DEVICE NAME 'aux': it cannot be opened\n",
      "by R. RENAME it (e.g. to 'auxiliary.seurat.shared.rds') — renaming is allowed and the\n",
      "content MD5 (2218a6817475f69dfdeeefa7451ada0ed, manifest) is unchanged by a rename.\n")
audit <- data.frame()
for (f in rds_files) {
  info <- tryCatch({
    o <- readRDS(f)
    data.frame(object = basename(f), class = class(o)[1],
               cells = tryCatch(ncol(o), error = function(e) NA),
               genes = tryCatch(nrow(o), error = function(e) NA),
               status = "OK", stringsAsFactors = FALSE)
  }, error = function(e)
    data.frame(object = basename(f), class = NA, cells = NA, genes = NA,
               status = paste("READ_FAILED:", conditionMessage(e)), stringsAsFactors = FALSE))
  print(info); audit <- rbind(audit, info)
}
write.csv(audit, "results/S2_qc_GSE213216_audit.csv", row.names = FALSE)
cat("\nS2_01 done. QC tables in results/, figures in figures/, filtered objects in data/processed/.\n")
