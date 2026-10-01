# ============================================================================
# S2_05_rescue_audit.R — Amendment S2-B: rescue audit + stratified QC table
# Registered BEFORE the S2_01 v4 run completed (logs/S2_preregistered_design.md
# Amendment S2-B, 2026-09-13). Oversight only — NO threshold is touched here.
# B1: characterize excluded cells (doublets + QC-fail) biologically per sample;
#     ALERT if excluded cells are immune-enriched (rules fixed in Amendment S2-B).
#     Also a determinism check: recomputed kept-count must equal the v4 CSV's
#     cells_after (independently verifies scDblFinder reproducibility).
# B2: stratified QC table per dataset x condition with fixed FLAG rules.
# Run AFTER S2_01 v4 completes. Runtime ~1.5-2.5 h (re-runs scDblFinder).
# Outputs: results/S2_rescue_audit_per_sample.csv
#          results/S2_rescue_audit_summary.csv
#          results/S2_qc_stratified.csv
#          figures/S2_rescue_audit_<GSE>.png (marker log2FC heatmap)
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(scDblFinder)
  library(SingleCellExperiment); library(ggplot2)
})

THRESH <- list(
  GSE179640 = list(min.genes = 200, max.genes = 6000, max.mt = 20, max.hb = 5),
  GSE214411 = list(min.genes = 200, max.genes = 6000, max.mt = 10, max.hb = 5),
  GSE183837 = list(min.genes = 200, max.genes = 6000, max.mt = 20, max.hb = 5)
)

# ---- Fixed panels (Amendment S2-B; selection frozen here) -------------------
IMMUNE <- c("PTPRC","CD3D","CD3E","CD4","CD8A","NCAM1","KLRD1","NKG7",
            "MS4A1","CD14","CD68","LYZ","FCGR3A","ITGAX","KIT","TPSB2")
REFERENCE <- c("EPCAM","KRT8","PECAM1","VWF","DCN","COL1A1")
CORE_PATHWAYS <- list(
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = c("NFKB1","NFKB2","RELA","RELB","TNF","TNFAIP3",
    "NFKBIA","NFKBIE","ICAM1","CXCL1","CXCL2","CXCL3","IL6","IL1B","CCL2","CCL5",
    "JUN","FOS","BIRC3","TRAF1"),
  HALLMARK_INTERFERON_GAMMA_RESPONSE = c("IFNG","STAT1","IRF1","CXCL9","CXCL10",
    "CXCL11","GBP1","GBP2","IDO1","HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1",
    "CD74","B2M","TAP1","TAP2","PSMB9","CIITA"),
  HALLMARK_INFLAMMATORY_RESPONSE = c("IL1A","IL1B","IL6","CXCL8","CCL2","CCL3",
    "CCL4","CCL5","TNF","NFKB1","PTGS2","S100A8","S100A9","TLR2","TLR4","NLRP3",
    "IL18","OSM"),
  HALLMARK_APOPTOSIS = c("BAX","BCL2","BCL2L1","BID","CASP3","CASP8","CASP9",
    "FAS","FADD","CYCS","APAF1","MCL1","XIAP","BIRC5","BBC3","PMAIP1"),
  HALLMARK_HYPOXIA = c("HIF1A","EPAS1","VEGFA","SLC2A1","LDHA","PGK1","ENO1",
    "BNIP3","BNIP3L","PDK1","CA9","ALDOA","GAPDH","HK2","PFKL")
)
IMMUNE_PATHWAYS <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     "HALLMARK_INFLAMMATORY_RESPONSE")

# prefer full MSigDB GMT if present (deterministic file), else core panels
gmt <- list.files("data/raw/gene_sets", pattern = "\\.gmt$", full.names = TRUE)
PATHWAYS <- CORE_PATHWAYS
if (length(gmt)) {
  gs <- list()
  for (ln in readLines(gmt[1])) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(f) > 2) gs[[f[1]]] <- f[-c(1, 2)]
  }
  if (all(names(CORE_PATHWAYS) %in% names(gs))) {
    PATHWAYS <- gs[names(CORE_PATHWAYS)]
    cat("pathway source: full GMT", gmt[1], "\n")
  } else cat("pathway source: core panels (GMT lacked the 5 Hallmark sets)\n")
} else cat("pathway source: core panels (no GMT found)\n")

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
  # GSE183837: name-based, to be verified against GEO metadata before freezing
  if (grepl("RIF", sid, ignore.case = TRUE)) return("RIF")
  if (grepl("con|ctrl|control|HC", sid, ignore.case = TRUE)) return("Control")
  "unmapped"
}

log2fc <- function(score, keep) {
  e <- score[!keep]; k <- score[keep]
  if (length(e) < 10 || length(k) < 10) return(NA_real_)
  log2((mean(expm1(e)) + 1) / (mean(expm1(k)) + 1))
}

audit_one <- function(counts, sid, thr, gse, csv_after) {
  seu <- CreateSeuratObject(counts, project = sid, min.cells = 3, min.features = 0)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  hb <- grep("^HB[AB]", rownames(seu), value = TRUE)
  seu[["percent.hb"]] <- if (length(hb)) PercentageFeatureSet(seu, features = hb) else 0
  # identical construction to S2_01 v4 -> identical scDblFinder calls
  sce <- as.SingleCellExperiment(seu)
  set.seed(42)
  sce <- scDblFinder(sce, samples = NULL)
  is_dbl  <- as.logical(sce$scDblFinder.class == "doublet")
  qc_fail <- !(seu$nFeature_RNA >= thr$min.genes & seu$nFeature_RNA <= thr$max.genes &
               seu$percent.mt <= thr$max.mt & seu$percent.hb <= thr$max.hb)
  keep <- (!is_dbl) & (!qc_fail)
  det <- if (!is.na(csv_after) && sum(keep) != csv_after)
    "INCONSISTENCY — kept-count != v4 cells_after" else "OK"
  if (det != "OK") cat("  !!", det, "for", sid, "(audit", sum(keep), "vs v4", csv_after, ")\n")

  seu <- NormalizeData(seu, verbose = FALSE)
  nd <- tryCatch(GetAssayData(seu, layer = "data"),
                 error = function(e) GetAssayData(seu, slot = "data"))
  row <- list(dataset = gse, sample = sid, condition = parse_cond(sid, gse),
              n_before = ncol(seu), n_kept = sum(keep), n_excluded = sum(!keep),
              n_dbl_only = sum(is_dbl & !qc_fail), n_qc_only = sum(!is_dbl & qc_fail),
              n_both = sum(is_dbl & qc_fail),
              dbl_pct_recomputed = round(mean(is_dbl) * 100, 2),
              determinism = det)
  ptprc <- tryCatch(GetAssayData(seu, layer = "counts")["PTPRC", ] > 0,
                    error = function(e) rep(NA, ncol(seu)))
  row$ptprc_frac_excluded <- if (sum(!keep) > 0) round(mean(ptprc[!keep], na.rm = TRUE), 4) else NA
  row$ptprc_frac_kept     <- round(mean(ptprc[keep], na.rm = TRUE), 4)
  for (mk in c(IMMUNE, REFERENCE)) {
    row[[paste0("l2fc_", mk)]] <- if (mk %in% rownames(nd))
      round(log2fc(as.numeric(nd[mk, ]), keep), 3) else NA_real_
  }
  for (pw in names(PATHWAYS)) {
    g <- intersect(PATHWAYS[[pw]], rownames(nd))
    row[[paste0("l2fc_", pw)]] <- if (length(g) >= 3)
      round(log2fc(colMeans(as.matrix(nd[g, , drop = FALSE])), keep), 3) else NA_real_
  }
  cat("  audited", sid, "- excluded", sum(!keep), "of", ncol(seu),
      sprintf("(%.1f%%)", 100 * mean(!keep)), "|", det, "\n")
  as.data.frame(row, check.names = FALSE)
}

run_audit <- function(gse) {
  csv <- file.path("results", paste0("S2_qc_", gse, ".csv"))
  if (!file.exists(csv)) { cat("skip", gse, "- no v4 QC CSV (run S2_01 v4 first)\n"); return(NULL) }
  after <- setNames(read.csv(csv)$cells_after,
                    sub("\\.h5$", "", basename(read.csv(csv)$sample)))
  cat("\n##### rescue audit:", gse, "#####\n")
  rows <- list()
  if (gse == "GSE179640") {
    for (f in list.files("data/raw/GSE179640",
                         pattern = "filtered_feature_bc_matrix\\.h5$", full.names = TRUE)) {
      sid <- sub("\\.h5$", "", basename(f))
      rows[[sid]] <- audit_one(Read10X_h5(f), sid, THRESH[[gse]], gse, after[[sid]])
    }
    cat("note: EOR02/EOR04 are load-excluded (Amendment S2-A) - non-auditable, recorded only.\n")
  } else {
    exdir <- file.path("data/raw", gse)
    dirs <- list.dirs(exdir, recursive = FALSE, full.names = TRUE)
    dirs <- dirs[file.exists(file.path(dirs, "matrix.mtx.gz"))]
    if (!length(dirs)) { cat("skip", gse, "- no staged trio dirs\n"); return(NULL) }
    for (d in dirs) {
      sid <- basename(d)
      rows[[sid]] <- audit_one(Read10X(d), sid, THRESH[[gse]], gse, after[[sid]])
    }
  }
  do.call(rbind, rows)
}

per_sample <- list()
for (gse in c("GSE179640", "GSE214411", "GSE183837")) {
  r <- tryCatch(run_audit(gse), error = function(e) { cat("ERROR in", gse, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) per_sample[[gse]] <- r
}
PS <- do.call(rbind, per_sample)
write.csv(PS, "results/S2_rescue_audit_per_sample.csv", row.names = FALSE)

# ---- B1 summary + ALERT rules (fixed in Amendment S2-B) ---------------------
summ <- list()
for (gse in unique(PS$dataset)) {
  d <- PS[PS$dataset == gse, ]
  med <- sapply(d[, grep("^l2fc_", names(d))], function(x) median(x, na.rm = TRUE))
  n_imm_mk <- sum(med[paste0("l2fc_", IMMUNE)] >= 1, na.rm = TRUE)
  n_imm_pw <- sum(med[paste0("l2fc_", IMMUNE_PATHWAYS)] >= 1, na.rm = TRUE)
  big <- d[d$n_excluded >= 50, ]
  ptprc_rule <- if (nrow(big))
    mean(big$ptprc_frac_excluded >= 2 * big$ptprc_frac_kept, na.rm = TRUE) > 0.5 else FALSE
  alert <- (n_imm_mk >= 2) || (n_imm_pw >= 1) || isTRUE(ptprc_rule)
  reason <- paste(c(
    if (n_imm_mk >= 2) sprintf("%d immune markers median log2FC>=1", n_imm_mk),
    if (n_imm_pw >= 1) sprintf("%d immune pathways median log2FC>=1", n_imm_pw),
    if (isTRUE(ptprc_rule)) "PTPRC+ fraction >=2x in >50% of samples"), collapse = "; ")
  summ[[gse]] <- data.frame(dataset = gse, n_samples = nrow(d),
    n_cells_before = sum(d$n_before), n_cells_excluded = sum(d$n_excluded),
    pct_excluded = round(100 * sum(d$n_excluded) / sum(d$n_before), 1),
    n_determinism_fail = sum(d$determinism != "OK"),
    immune_markers_median_l2fc_ge1 = n_imm_mk,
    immune_pathways_median_l2fc_ge1 = n_imm_pw,
    ptprc_rule_fired = isTRUE(ptprc_rule),
    ALERT = alert, alert_reason = reason)
  if (alert) cat("\n!!! RESCUE AUDIT ALERT:", gse, "-", reason,
                 "-> must be dispositioned in logs/S2_verification_report.md BEFORE S2_02 !!!\n")
  if (sum(d$determinism != "OK")) cat("\n!!! DETERMINISM FAILURE in", gse,
      "-", sum(d$determinism != "OK"), "samples: kept-count != v4 cells_after. STOP and report. !!!\n")
}
SM <- do.call(rbind, summ)
write.csv(SM, "results/S2_rescue_audit_summary.csv", row.names = FALSE)
print(SM)

# per-dataset marker heatmap (best-effort; never fatal)
for (gse in unique(PS$dataset)) {
  tryCatch({
    d <- PS[PS$dataset == gse, ]
    m <- as.matrix(d[, paste0("l2fc_", c(IMMUNE, REFERENCE))])
    rownames(m) <- sub("^GSM[0-9]+_", "", d$sample)
    colnames(m) <- sub("^l2fc_", "", colnames(m))
    library(pheatmap)
    png(file.path("figures", paste0("S2_rescue_audit_", gse, ".png")),
        width = 1400, height = max(600, 30 * nrow(m)), res = 150)
    pheatmap(t(m), cluster_rows = FALSE, cluster_cols = TRUE,
             main = paste(gse, "- excluded vs kept: marker log2FC (positive = enriched in EXCLUDED)"),
             fontsize = 8)
    dev.off()
  }, error = function(e) cat("heatmap skipped for", gse, ":", conditionMessage(e), "\n"))
}

# ---- B2 stratified QC table + FLAG rules (fixed in Amendment S2-B) ----------
strat <- list()
for (gse in c("GSE179640", "GSE214411", "GSE183837")) {
  csv <- file.path("results", paste0("S2_qc_", gse, ".csv"))
  if (!file.exists(csv)) next
  d <- read.csv(csv)
  d$condition <- sapply(d$sample, parse_cond, gse = gse)
  thr <- THRESH[[gse]]
  for (cc in unique(d$condition)) {
    x <- d[d$condition == cc, ]
    strat[[paste(gse, cc)]] <- data.frame(
      dataset = gse, condition = cc, n_samples = nrow(x),
      retention_median = median(x$pct_kept), retention_min = min(x$pct_kept),
      retention_max = max(x$pct_kept),
      median_mt_median = median(x$median_mt), median_mt_max = max(x$median_mt),
      median_genes_median = median(x$median_genes),
      dbl_pct_median = median(x$dbl_pct, na.rm = TRUE))
  }
  sp <- aggregate(pct_kept ~ condition, d, median)
  flags <- c()
  if (max(sp$pct_kept) - min(sp$pct_kept) > 20)
    flags <- c(flags, sprintf("retention spread %.1f pts > 20", max(sp$pct_kept) - min(sp$pct_kept)))
  if (max(aggregate(median_mt ~ condition, d, median)$median_mt) >= 0.75 * thr$max.mt)
    flags <- c(flags, sprintf("a condition median_mt >= 75%% of threshold (%g%%)", thr$max.mt))
  mg <- aggregate(median_genes ~ condition, d, median)$median_genes
  if (any(mg <= 220 | mg >= 5400)) flags <- c(flags, "a condition median_genes within 10% of a bound")
  if (length(flags))
    cat("\n!!! STRATIFIED QC FLAG:", gse, "-", paste(flags, collapse = "; "),
        "-> report in logs/S2_verification_report.md BEFORE S2_02 !!!\n")
}
ST <- do.call(rbind, strat)
write.csv(ST, "results/S2_qc_stratified.csv", row.names = FALSE)
print(ST)
cat("\nS2_05 done. Rescue audit + stratified QC in results/.\n")
