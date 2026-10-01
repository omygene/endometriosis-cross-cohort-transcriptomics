# ============================================================================
# S2_03_bulk_qc.R  v3.1  -- bulk microarray QC + normalization
# Pre-registered: logs/S2_preregistered_design.md Section 6, as amended by
# logs/S2_03_amendment_A.md + Addendum S2_03-A.1 (2026-09-18):
#   - AQM runs on the RMA ExpressionSet (do.logtransform = FALSE); outlier calls
#     extracted programmatically from outlierDetection slots (locked-binary-
#     verified: statistic / threshold / which / description); rule unchanged:
#     flagged if >= 2 of 3 AQ metrics fail.
#   - AQM also for GSE7307-subset and GSE11691; GSE120103 = neqc only.
#   - GSE120103: genes$Status marshalled from ControlType (Addendum A.1).
# v3.1 fixes measured in the smoke run (logs/S2_03_smoke.txt, 2026-09-18):
#   F-A  outlierDetection@description is a length-2 vector (metric text +
#        "data-driven"/"fixed") -> collapse before grepl (R >= 4.2 rejects
#        length>1 conditions).
#   F-B  MA-plot module description is literally "<i>D<sub>a</sub></i>"
#        (threshold 0.15, fixed) -> identification pattern widened.
#   F-C  Agilent: genes$Status absent in these FE files -> neqc routed to
#        normexp.fit.detection.p and stopped "Detection p values not found
#        in the data." -> Status marshalled from ControlType with fail-loud
#        guards (coding 0/-1 measured by scripts/S2_03_agilent_diag.R).
# Execution model (B.6/B.7): ONE cohort per short-lived process, chained by
# scripts\run_S2_03.bat.
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S2_03_bulk_qc.R <COHORT|collect>
# Writes ONLY new S2_03 files (R4): data/processed/bulk_*.rds,
# results/S2_aqm_*/, results/S2_bulk_outliers_<GSE>.csv,
# results/S2_bulk_outliers.csv (collect), figures/S2_bulk_*_pca.png
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/S2_03_bulk_qc.R <COHORT|collect>")
cohort <- args[1]

suppressPackageStartupMessages({
  library(oligo); library(limma)
  library(arrayQualityMetrics); library(ggplot2)
})
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

CFG <- list(
  GSE7305   = list(tar = "data/raw/GSE7305_RAW.tar",   platform = "GPL570",  members = 20L),
  GSE6364   = list(tar = "data/raw/GSE6364_RAW.tar",   platform = "GPL570",  members = 37L),
  GSE25628  = list(tar = "data/raw/GSE25628_RAW.tar",  platform = "GPL570",  members = 22L),
  GSE7307   = list(tar = NULL, dir0 = "data/raw/GSE7307_subset", platform = "GPL570", members = 41L),
  GSE51981  = list(tar = "data/raw/GSE51981_RAW.tar",  platform = "GPL570",  members = 148L),
  GSE11691  = list(tar = "data/raw/GSE11691_RAW.tar",  platform = "GPL96",   members = 18L),
  GSE120103 = list(tar = "data/raw/GSE120103_RAW.tar", platform = "Agilent", members = 36L)
)
AQM_COHORTS <- c("GSE7305", "GSE6364", "GSE25628", "GSE7307", "GSE51981", "GSE11691")

## ---- collect mode: merge per-cohort outlier tables --------------------------
if (cohort == "collect") {
  say("== S2_03 collect: merging per-cohort outlier tables ==")
  files <- file.path("results", paste0("S2_bulk_outliers_", AQM_COHORTS, ".csv"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("FATAL: missing per-cohort outlier CSVs: ", paste(missing, collapse = ", "))
  all_rows <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  write.csv(all_rows, "results/S2_bulk_outliers.csv", row.names = FALSE)
  say("  wrote results/S2_bulk_outliers.csv --", nrow(all_rows), "arrays across", length(AQM_COHORTS), "cohorts")
  fl <- all_rows[all_rows$flagged, ]
  if (nrow(fl) > 0) {
    say("  FLAGGED arrays (>=2 of 3 AQ metrics):")
    print(fl)
  } else {
    say("  zero arrays flagged by the >=2-of-3 rule")
  }
  say("COLLECT DONE")
  quit(save = "no", status = 0)
}

## ---- cohort mode -------------------------------------------------------------
if (!cohort %in% names(CFG)) stop("Unknown cohort: ", cohort,
                                  " -- expected one of: ", paste(c(names(CFG), "collect"), collapse = ", "))
cfg <- CFG[[cohort]]
say("== S2_03 cohort:", cohort, "(", cfg$platform, ") ==")

## 1. locate raw files (fail-loud; no silent skips)
if (cohort == "GSE7307") {
  exdir <- cfg$dir0
  if (!dir.exists(exdir)) stop("FATAL: ", exdir, " missing")
  data_files <- list.files(exdir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
} else {
  exdir <- file.path("data/raw", cohort)
  if (!dir.exists(exdir)) {
    if (!file.exists(cfg$tar)) stop("FATAL: ", cfg$tar, " missing and ", exdir, " absent")
    say("  untar", cfg$tar, "->", exdir)
    untar(cfg$tar, exdir = exdir, tar = "internal")
  }
  pat <- if (cfg$platform == "Agilent") "\\.txt\\.gz$" else "\\.CEL(\\.gz)?$"
  data_files <- list.files(exdir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
}
if (length(data_files) != cfg$members)
  stop("FATAL: ", cohort, " has ", length(data_files), " data files, expected ", cfg$members,
       " -- check extraction state before re-running")
say("  data files:", length(data_files), "(expected", cfg$members, ") OK")

## 2. normalization
has_group <- FALSE
if (cfg$platform == "Agilent") {
  say("  read.maimages (source=agilent, green.only=TRUE) ...")
  rg <- limma::read.maimages(data_files, source = "agilent", green.only = TRUE)
  colnames(rg) <- basename(data_files)
  say("  RGList dim:", paste(dim(rg), collapse = " x "),
      "| genes$Status present:", "Status" %in% colnames(rg$genes))
  ## Addendum S2_03-A.1: marshal genes$Status from ControlType (fix F-C).
  ## Locked-binary evidence (limma 3.68.5): neqc -> nec routes on genes$Status;
  ## without it nec calls normexp.fit.detection.p which stops "Detection p
  ## values not found in the data." (measured in logs/S2_03_smoke.txt leg 5).
  ## ControlType coding 0/-1 measured by scripts/S2_03_agilent_diag.R first.
  if (!"Status" %in% colnames(rg$genes)) {
    ct <- rg$genes$ControlType
    if (is.null(ct))
      stop("FATAL: ", cohort, " genes$ControlType absent -- cannot marshal Status; escalate (R9)")
    tab_ct <- table(ct, useNA = "always")
    say("  ControlType table:", paste(paste(names(tab_ct), tab_ct, sep = "="), collapse = " / "))
    n_reg <- sum(ct == 0, na.rm = TRUE); n_neg <- sum(ct == -1, na.rm = TRUE)
    if (n_reg == 0 || n_neg == 0)
      stop("FATAL: ", cohort, " ControlType coding not 0/-1 (code 0: ", n_reg,
           ", code -1: ", n_neg, ") -- do NOT marshal; escalate with this log (R9)")
    rg$genes$Status <- ifelse(ct == 0, "regular", ifelse(ct == -1, "negative", "other"))
    say("  Status marshalled from ControlType:", n_reg, "regular /", n_neg, "negative /",
        sum(rg$genes$Status == "other"), "other")
  }
  eset <- limma::neqc(rg)
  rm(rg)
  saveRDS(eset, file.path("data/processed", "bulk_GSE120103_neqc.rds"))
  mat <- eset$E
} else {
  say("  oligo::read.celfiles ...")
  raw <- oligo::read.celfiles(data_files)
  print(gc())
  say("  oligo::rma ...")
  eset <- oligo::rma(raw)
  rm(raw); print(gc())
  if (cohort == "GSE7307") {
    lab <- read.delim("scripts/S1_GSE7307_subset.txt", comment.char = "#", header = FALSE,
                      col.names = c("group", "gsm"), stringsAsFactors = FALSE)
    gsm <- regmatches(colnames(eset), regexpr("GSM[0-9]+", colnames(eset), ignore.case = TRUE))
    idx <- match(gsm, lab$gsm)
    if (any(is.na(idx)))
      stop("FATAL: GSE7307 group mapping produced NA for: ",
           paste(colnames(eset)[is.na(idx)], collapse = ", "))
    eset$group <- lab$group[idx]
    has_group <- TRUE
    say("  GSE7307 groups attached:", paste(names(table(eset$group)), table(eset$group), collapse = " / "))
  }
  saveRDS(eset, file.path("data/processed", paste0("bulk_", cohort, "_rma.rds")))
  mat <- Biobase::exprs(eset)
}
say("  normalized matrix:", nrow(mat), "x", ncol(mat))

## 3. arrayQualityMetrics + registered >=2-of-3 outlier rule
if (cohort %in% AQM_COHORTS) {
  aqm_dir <- file.path("results", paste0("S2_aqm_", cohort))
  say("  arrayQualityMetrics ->", aqm_dir)
  aqm_err <- NULL
  aqm_ret <- tryCatch(
    if (has_group)
      arrayQualityMetrics(eset, outdir = aqm_dir, force = TRUE, do.logtransform = FALSE, intgroup = "group")
    else
      arrayQualityMetrics(eset, outdir = aqm_dir, force = TRUE, do.logtransform = FALSE),
    error = function(e) { aqm_err <<- e; NULL })
  if (!is.null(aqm_err))
    stop("FATAL: arrayQualityMetrics failed for ", cohort, ": ", conditionMessage(aqm_err),
         " -- escalate with this log (R9)")
  if (!file.exists(file.path(aqm_dir, "index.html")))
    stop("FATAL: AQM report index.html missing for ", cohort, " -- escalate with this log (R9)")

  mods <- NULL
  if (is.list(aqm_ret) && "modules" %in% names(aqm_ret)) mods <- aqm_ret$modules
  if (is.null(mods) && is.list(aqm_ret) && length(aqm_ret) > 0 &&
      all(vapply(aqm_ret, function(z) is(z, "aqmReportModule"), logical(1)))) mods <- aqm_ret
  if (is.null(mods)) {
    say("  return value carried no modules; rebuilding via prepdata (documented fallback)")
    pd <- arrayQualityMetrics::prepdata(eset, do.logtransform = FALSE)
    mods <- list(boxplot = arrayQualityMetrics::aqm.boxplot(pd),
                 maplot  = arrayQualityMetrics::aqm.maplot(pd),
                 heatmap = arrayQualityMetrics::aqm.heatmap(pd))
  }

  ods <- list()
  for (m in mods) {
    od <- tryCatch(slot(m, "outliers"), error = function(e) NULL)
    if (is.null(od) || !is(od, "outlierDetection")) next
    ## F-A: description is a length-2 vector (metric text + "data-driven"/"fixed",
    ## measured in logs/S2_03_smoke.txt) -> collapse to one string before grepl.
    d <- paste(slot(od, "description"), collapse = " ")
    ## F-B: the MA-plot module's description is literally "<i>D<sub>a</sub></i>"
    ## (Hoeffding's D, threshold 0.15 fixed) -- "Hoeffding" never appears in it.
    key <- if (grepl("Kolmogorov", d)) "boxplot"
      else if (grepl("Hoeffding", d) || grepl("D<sub>a</sub>", d, fixed = TRUE) ||
               grepl("D_a", d, fixed = TRUE)) "maplot"
      else if (grepl("sum of distances", d)) "heatmap"
      else NA_character_
    if (!is.na(key) && !(key %in% names(ods))) ods[[key]] <- od
  }
  if (!all(c("boxplot", "maplot", "heatmap") %in% names(ods)))
    stop("FATAL: could not identify all 3 AQ outlier metrics for ", cohort,
         " -- found: ", paste(names(ods), collapse = ", "), ". Escalate with this log (R9).")

  n <- ncol(mat)
  tab <- data.frame(gse = cohort, sample = colnames(mat), stringsAsFactors = FALSE)
  for (key in c("boxplot", "maplot", "heatmap")) {
    od  <- ods[[key]]
    st  <- slot(od, "statistic"); thr <- slot(od, "threshold"); wh <- slot(od, "which")
    if (length(st) != n)
      stop("FATAL: ", cohort, " ", key, " statistic length ", length(st), " != ", n, " arrays")
    tab[[paste0("stat_", key)]] <- as.numeric(st)
    tab[[paste0("thr_",  key)]] <- if (length(thr) == 1) as.numeric(thr) else NA_real_
    out <- rep(FALSE, n); out[wh] <- TRUE
    tab[[paste0("out_",  key)]] <- out
  }
  tab$n_fail  <- tab$out_boxplot + tab$out_maplot + tab$out_heatmap
  tab$flagged <- tab$n_fail >= 2
  out_csv <- file.path("results", paste0("S2_bulk_outliers_", cohort, ".csv"))
  write.csv(tab, out_csv, row.names = FALSE)
  say("  wrote", out_csv, "--", nrow(tab), "arrays | per-method outliers:",
      sum(tab$out_boxplot), "boxplot /", sum(tab$out_maplot), "maplot /",
      sum(tab$out_heatmap), "heatmap | FLAGGED >=2of3:", sum(tab$flagged))
} else {
  say("  AQM not registered for", cohort, "(design SS6) -- skipped by design")
}

## 4. PCA figure (every cohort)
pc  <- prcomp(t(mat), scale. = TRUE)
pct <- round(100 * summary(pc)$importance[2, 1:2], 1)
df  <- data.frame(pc$x[, 1:2], sample = colnames(mat), stringsAsFactors = FALSE)
p <- ggplot(df, aes(PC1, PC2, label = sample))
if (has_group) { df$group <- eset$group; p <- ggplot(df, aes(PC1, PC2, label = sample, color = group)) }
p <- p + geom_point(size = 2) + geom_text(size = 2, vjust = -0.6) +
  labs(title = paste(cohort, "PCA"),
       x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))
ggsave(file.path("figures", paste0("S2_bulk_", cohort, "_pca.png")), p, width = 9, height = 7, dpi = 150)

print(gc())
say("COHORT", cohort, "DONE")
