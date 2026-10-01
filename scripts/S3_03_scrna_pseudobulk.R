# ============================================================================
# S3_03_scrna_pseudobulk.R -- PRIMARY scRNA contrasts (locked criterion 2)
# Pseudobulk per dataset x lineage (sample-level sums of RNA counts), limma
# voom. Design per dataset: ~ group [+ hormonal covariate for GSE179640 arm].
# Arms:
#   A (primary)   : all samples, N1 INCLUDED
#   B             : all samples, N1 EXCLUDED  (Amendment S2-B/B3)
#   A-HCOV        : SUSPENDED by Amendment A (2026-09-20, zero-variance covariate)
#   A-NOTX        : SUSPENDED by Amendment A (2026-09-20)
# Group map: per dataset (S3_group_map.tsv); hormonal status from the GEO
# series matrix (S3_00c extraction), joined by GSM accession = gsm_of(orig.ident)
# Gates: BH-FDR < 0.05; |log2FC| >= 0.5 interpretive tier, no exclusion.
# v3.1 (2026-09-21): expression source = per-dataset objects (measured
# S3_00d/e/f: atlas S2_02_integrated.rds is a 28-gene panel without any
# expression matrix; the 2026-09-21 first run's pseudobulk on 28 markers is
# VOID and is overwritten by this regeneration, R4). Annotations join by
# cell_key = <orig.ident>__<barcode> with a fail-loud 100% gate. The frozen
# criterion (pseudobulk per sample, voom, same gates) is UNCHANGED.
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_03_scrna_pseudobulk.R > logs\S3_03_scrna_pseudobulk.txt 2>&1
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(limma); library(Seurat) })
source("scripts/S3_common.R")
dir.create("results/S3_scrna_pseudobulk", recursive = TRUE, showWarnings = FALSE)

## Amendment A (lead decision 2026-09-20): A-HCOV / A-NOTX SUSPENDED for
## GSE179640 -- the hormonal covariate is a series-level constant (measured:
## "hormonally treated patients", zero variance; see S3_design.md section 10).
## hormonal_col is therefore NOT required. The arm code below stays, gated
## on nzchar(CFG$hormonal_col): it reactivates only if a registered external
## source is ever approved (R5); it can never run silently.
require_config(c("n1_sample", "group_col", "lineage_col"))
say("== S3_03 pseudobulk primary (per-dataset sources, v3.1) ==")
md <- atlas_annotations()
cell_grp_all <- map_cell_group(as.character(md$dataset), as.character(md[[CFG$group_col]]))
say("atlas cells:", length(cell_grp_all), "| TEST:", sum(cell_grp_all == "TEST"),
    "| REF:", sum(cell_grp_all == "REF"))
if (!(CFG$n1_sample %in% as.character(md$sample)))
  stop("FATAL: n1_sample '", CFG$n1_sample, "' not among atlas sample ids (see S3_00 census)")

## frozen gene filter for pseudobulk (pre-counts CPM screen)
cpm_screen <- function(cnt) {
  lib <- colSums(cnt)
  cpm <- sweep(cnt, 2, lib / 1e6, "/")
  keep <- rowSums(cpm > 1) >= max(3, ceiling(ncol(cnt) / 2))
  cnt[keep, , drop = FALSE]
}

fit_one <- function(cnt, grp, cov = NULL) {
  design <- if (is.null(cov)) model.matrix(~ factor(grp, levels = c("REF", "TEST"))) else
    model.matrix(~ factor(grp, levels = c("REF", "TEST")) + factor(cov))
  v <- voom(cnt, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), tt, row.names = NULL, stringsAsFactors = FALSE)
}

run_leg <- function(tag, ds, lin, counts, cell_lin, cell_smp, cell_grp, cell_horm,
                    exclude_smp = NULL, use_horm = FALSE, drop_treated = FALSE) {
  ## lin == "ALL": whole-sample reference leg (no lineage split). Latent bug
  ## measured in the 2026-09-21 log ("skip A ALL (too few cells)" every
  ## dataset, every prior version): the old selection `cell_lin == "ALL"`
  ## matched zero cells, so this leg NEVER ran. Fixed: select all cells.
  sel <- (if (lin == "ALL") rep(TRUE, length(cell_lin)) else cell_lin == lin) &
    !(cell_smp %in% exclude_smp)
  if (drop_treated) sel <- sel & cell_horm != CFG$hormonal_treated_lbl
  if (sum(sel) < CFG$min_cells_pseudobulk) { say("  skip", tag, lin, "(too few cells)"); return(invisible(NULL)) }
  smps <- sort(unique(cell_smp[sel]))
  agg <- sapply(smps, function(s) Matrix::rowSums(counts[, sel & cell_smp == s, drop = FALSE]))
  agg <- cpm_screen(as.matrix(agg))
  grp <- unname(vapply(smps, function(s) unique(cell_grp[sel & cell_smp == s]), character(1)))
  if (length(unique(grp)) < 2) { say("  skip", tag, lin, "(single group)"); return(invisible(NULL)) }
  cov <- if (use_horm) unname(vapply(smps, function(s) unique(cell_horm[sel & cell_smp == s]), character(1))) else NULL
  tt <- fit_one(agg, grp, cov)
  tt$dataset <- ds; tt$lineage <- lin; tt$arm <- tag
  tt$n_test <- sum(grp == "TEST"); tt$n_ref <- sum(grp == "REF")
  out <- file.path("results/S3_scrna_pseudobulk",
                   paste0("S3_pb_", ds, "_", lin, "_", tag, ".csv"))
  write.csv(tt, out, row.names = FALSE)
  say("  wrote", basename(out), "| genes:", nrow(tt),
      "| FDR<", CFG$fdr_gate, ":", sum(tt$adj.P.Val < CFG$fdr_gate))
}

datasets <- intersect(sort(unique(as.character(md$dataset))), names(CFG$per_dataset))
if (!length(datasets)) stop("FATAL: no atlas dataset with a registered per_dataset source (R9)")
for (ds in datasets) {
  say("dataset:", ds)
  per <- readRDS(CFG$per_dataset[[ds]])
  per <- prep_rna_norm(per)
  ann <- join_annotations(per, md, ds)
  counts <- SeuratObject::GetAssayData(per, assay = "RNA", layer = "counts")
  cell_lin <- as.character(ann$lineage)
  cell_smp <- as.character(ann$sample)
  cell_grp <- map_cell_group(rep(ds, nrow(ann)), as.character(ann$condition))
  say("lineages:", paste(sort(unique(cell_lin)), collapse = " | "))

  ## hormonal status (Amendment A gated): MEASURED ABSENT from atlas and
  ## per-dataset objects -> source is the GEO series matrix (S3_00c), joined
  ## by GSM accession = gsm_of(orig.ident). Gated: can never run silently.
  cell_horm <- rep(NA_character_, nrow(ann))
  if (nzchar(CFG$hormonal_col)) {
    hm <- load_hormonal_map()
    cell_gsm <- gsm_of(as.character(per$orig.ident))
    miss <- sort(unique(cell_gsm[is.na(match(cell_gsm, hm$gsm))]))
    if (length(miss))
      stop("FATAL: ", ds, " GSM accessions absent from hormonal map: ",
           paste(miss, collapse = ", "), " -- escalate (R9)")
    cell_horm <- hm$horm[match(cell_gsm, hm$gsm)]
    if (any(is.na(cell_horm))) stop("FATAL: unmatched cells in hormonal map (R9)")
  } else if (ds == "GSE179640") {
    say("  hormonal_col empty -- A-HCOV / A-NOTX legs skipped (fill after S3_00c census)")
  }

  for (lin in sort(unique(cell_lin))) {
    run_leg("A", ds, lin, counts, cell_lin, cell_smp, cell_grp, cell_horm)
    if (ds == "GSE214411") run_leg("B", ds, lin, counts, cell_lin, cell_smp, cell_grp, cell_horm,
                                   exclude_smp = CFG$n1_sample)
    if (ds == "GSE179640" && nzchar(CFG$hormonal_col)) {
      run_leg("A-HCOV", ds, lin, counts, cell_lin, cell_smp, cell_grp, cell_horm, use_horm = TRUE)
      run_leg("A-NOTX", ds, lin, counts, cell_lin, cell_smp, cell_grp, cell_horm, drop_treated = TRUE)
    }
  }
  run_leg("A", ds, "ALL", counts, cell_lin, cell_smp, cell_grp, cell_horm)
  if (ds == "GSE214411") run_leg("B", ds, "ALL", counts, cell_lin, cell_smp, cell_grp, cell_horm,
                                 exclude_smp = CFG$n1_sample)
  rm(per, counts, ann); gc()
}
print(gc())
say("S3_03 DONE")
