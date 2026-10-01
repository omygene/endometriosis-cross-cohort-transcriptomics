# ============================================================================
# S4_05_axis_localization.R -- Hypothesis Localization Diagnostic (A4, 2026-09-22)
# Pure measurement: WHERE does each pre-registered axis gene's signal live?
#   bulk arm: per locked primary contrast (S3_01 v2.0 symbol-level CSVs)
#   scRNA arm: per pseudobulk leg (S3_pb_*_A/B/ALL.csv)
# Rules locked in S4_design section 13 (R5):
#   * NO new thresholds. Only signal flag = adj.P.Val < 0.05 (the locked gate).
#   * No exclusion, no weighting, no pooling. B-G11691 reported as its own arm.
#   * Every gene x contrast/leg emits a row whether found or not (absence is a
#     measured fact, never a silently dropped row -- R9).
#   * Any missing primary contrast file -> FATAL naming it (R2). Empty pb dir
#     -> FATAL. All columns are read from the files, never assumed.
# Writes ONLY under results/S4_localization/
# Usage: Rscript scripts/S4_05_axis_localization.R
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
source("scripts/S4_common.R")

OUT <- "results/S4_localization"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## Frozen axis registry (A4, design section 13). Changes require an Amendment.
AXES <- list(
  cytotoxicity         = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GNLY"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
  senescence_dormancy  = c("CDKN1A", "CDKN2A", "GADD45A", "SERPINE1"),
  proliferation        = c("MKI67", "TOP2A", "PCNA"),
  stromal_ecm          = c("VIM", "COL1A1", "COL3A1", "FN1", "ACTA2")
)
ALL_GENES <- unique(unlist(AXES))
if (length(ALL_GENES) != 24L)
  stop("FATAL: frozen axis registry drifted (expected 24 unique genes, got ",
       length(ALL_GENES), ") -- Amendment required (R5)")
axis_of_gene <- stack(AXES)   # columns: values = gene, ind = axis
names(axis_of_gene) <- c("gene", "axis")

say("== S4_05: axis localization map (pure measurement, no exclusion) ==")
say("frozen axes:", length(AXES), "| unique genes:", length(ALL_GENES))

## ---- bulk arm: the 7 locked primary contrasts ------------------------------
fr <- load_frozen()
prim <- fr[fr$level == "primary", , drop = FALSE]
if (nrow(prim) != 7L)
  stop("FATAL: frozen registry has ", nrow(prim), " primary contrasts, expected 7 (R2)")
files <- file.path(CFG$out_contrasts, paste0(prim$contrast_id, ".csv"))
missing <- files[!file.exists(files)]
if (length(missing))
  stop("FATAL: missing primary contrast outputs (run S3_01 v2.0 repair first):\n  ",
       paste(missing, collapse = "\n  "))

sig_of <- function(adj_p) !is.na(adj_p) & adj_p < CFG$fdr_gate

bulk_rows <- list()
for (i in seq_len(nrow(prim))) {
  d <- read.csv(files[i], stringsAsFactors = FALSE)
  if (!all(c("gene", "logFC", "t", "P.Value", "adj.P.Val") %in% names(d)))
    stop("FATAL: ", basename(files[i]), " lacks locked limma columns (R2): got ",
         paste(names(d), collapse = ","))
  hit <- match(ALL_GENES, d$gene)
  present <- !is.na(hit)
  z <- d[hit[present], , drop = FALSE]
  bulk_rows[[i]] <- data.frame(
    axis = axis_of_gene$axis[match(ALL_GENES[present], axis_of_gene$gene)],
    gene = ALL_GENES[present],
    contrast_id = prim$contrast_id[i], cohort = prim$cohort[i],
    design = prim$design[i],
    logFC = z$logFC, t = z$t, P.Value = z$P.Value, adj.P.Val = z$adj.P.Val,
    sig_fdr05 = sig_of(z$adj.P.Val),
    stringsAsFactors = FALSE)
  if (any(!present))
    say("  ", prim$contrast_id[i], ": genes absent from contrast output (measured):",
        paste(ALL_GENES[!present], collapse = ", "))
}
bulk_map <- do.call(rbind, bulk_rows)
write.csv(bulk_map, file.path(OUT, "S4_05_axis_bulk_map.csv"), row.names = FALSE)
say("bulk map:", nrow(bulk_map), "rows |", sum(bulk_map$sig_fdr05),
    "significant cells (of", nrow(bulk_map), "gene x contrast)")

## ---- scRNA arm: every pseudobulk leg ---------------------------------------
pb_files <- sort(list.files(CFG4$pb_dir, pattern = "^S3_pb_.*_(A|B|ALL)\\.csv$",
                            full.names = TRUE))
if (length(pb_files) == 0L)
  stop("FATAL: no pseudobulk legs under ", CFG4$pb_dir, " (run S3_03 first, R2)")
say("scRNA legs found:", length(pb_files))

scrna_rows <- list()
for (f in pb_files) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("gene", "logFC", "t", "P.Value", "adj.P.Val",
             "dataset", "lineage", "arm", "n_test", "n_ref") %in% names(d)))
    stop("FATAL: ", basename(f), " lacks locked pseudobulk columns (R2)")
  hit <- match(ALL_GENES, d$gene)
  present <- !is.na(hit)
  z <- d[hit[present], , drop = FALSE]
  scrna_rows[[f]] <- data.frame(
    axis = axis_of_gene$axis[match(ALL_GENES[present], axis_of_gene$gene)],
    gene = ALL_GENES[present],
    dataset = z$dataset, lineage = z$lineage, arm = z$arm,
    logFC = z$logFC, t = z$t, P.Value = z$P.Value, adj.P.Val = z$adj.P.Val,
    sig_fdr05 = sig_of(z$adj.P.Val),
    n_test = z$n_test, n_ref = z$n_ref,
    stringsAsFactors = FALSE)
}
scrna_map <- do.call(rbind, scrna_rows)
write.csv(scrna_map, file.path(OUT, "S4_05_axis_scrna_map.csv"), row.names = FALSE)
say("scRNA map:", nrow(scrna_map), "rows |", sum(scrna_map$sig_fdr05),
    "significant cells | legs:", length(pb_files))

## ---- summary md (measurement statements only, no thresholds) ---------------
meta_set6 <- CFG4$meta_set
paired <- CFG4$paired_arm
arrow <- function(lfc, sig) {
  if (is.na(lfc)) return(".")
  a <- if (lfc > 0) "+" else "-"
  if (isTRUE(sig)) paste0(a, a) else a
}
L <- c()
L <- c(L, "# S4_05 -- Axis localization summary (measurement only)")
L <- c(L, "")
L <- c(L, "Locked criteria: signal flag = BH-FDR < 0.05 only (design section 13).",
       "No exclusion, no weighting, no pooling. B-G11691 is its own paired arm.",
       "'++/--' = FDR < 0.05; '+/-' = detected but not FDR-significant;",
       "'.' = gene absent from that output (measured, never dropped silently).")
for (ax in names(AXES)) {
  L <- c(L, "", paste0("## ", ax, " (", paste(AXES[[ax]], collapse = ", "), ")"))
  g <- AXES[[ax]]
  L <- c(L, "", "### bulk: gene x primary contrast (direction of TEST vs REF)")
  L <- c(L, paste0("| gene | ", paste(meta_set6, collapse = " | "), " | ",
                   paired, " |"), paste0("|", strrep("-", 6), rep("", 1), "|",
                   paste(rep("------", length(meta_set6) + 1), collapse = "|"), "|"))
  b <- bulk_map[bulk_map$axis == ax, ]
  for (gn in g) {
    cells <- vapply(c(meta_set6, paired), function(cid) {
      r <- b[b$gene == gn & b$contrast_id == cid, ]
      if (!nrow(r)) "." else arrow(r$logFC[1], r$sig_fdr05[1])
    }, character(1))
    L <- c(L, paste0("| ", gn, " | ", paste(cells, collapse = " | "), " |"))
  }
  s <- scrna_map[scrna_map$axis == ax, ]
  if (nrow(s)) {
    L <- c(L, "", "### scRNA: significant cells (FDR < 0.05) per lineage x arm")
    tab <- sort(table(paste(s$dataset, s$lineage, s$arm, sep = " | ")[s$sig_fdr05]),
                decreasing = TRUE)
    if (length(tab)) {
      L <- c(L, paste0("- ", names(tab), ": ", as.integer(tab), " of ", length(g),
                       " axis genes"))
    } else L <- c(L, "- (no significant cells in any leg)")
  }
}
writeLines(L, file.path(OUT, "S4_05_axis_summary.md"))
say("wrote", file.path(OUT, "S4_05_axis_summary.md"))
say("S4_05 DONE")
