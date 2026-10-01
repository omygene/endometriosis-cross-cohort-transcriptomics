# ============================================================================
# S5_01_pathway_scores.R -- per-sample pathway scores for the 6 locked axes
# v1.0.1 (2026-09-22): registry READ fixed -- the registry is written by
#   write.csv (comma-separated) but was read back with read.delim (tab), so
#   each "file,md5" line collapsed into one column and every registered input
#   looked "vanished" (measured on lead machine, S5_02 log 14:11, ledger #20,
#   my error, R6). Fix: read.csv. The registry FILE on disk is valid as-is --
#   no deletion, no hand-edit (R4); the two fixed scripts read it correctly.
# v1.0 (2026-09-22, Amendment A4, design section 14).
# For each of the 7 bulk cohorts:
#   gene-wise z across samples WITHIN cohort -> mean z over each axis's genes.
#   SenMayo is NOT recomputed: it reuses the locked S3_02 ssGSEA scores
#   (results/S3_ssgsea_senmayo_bulk.csv) -- R4, no recomputation of frozen
#   results. Its header is verified against the schema measured on the lead
#   machine (cohort,sample,gsm,condition,senmayo_ssgsea); any drift = FATAL.
# Locked rules:
#   * structural floor: an axis needs >= 2 z-able genes (sd > 0); below that a
#     "pathway score" is a single-gene readout in disguise -> FATAL with the
#     measured count (R9: direct measurement, not theoretical inference).
#   * sample match between SenMayo CSV and the symbol matrix must be exact
#     (setequal), else FATAL -- no silent sample loss.
#   * input registry (MD5 of SenMayo CSV + 7 bulk rds + frozen TSV): created
#     loud on first run, verified on later runs, mismatch halts with the
#     amendment route (R5). Same pattern as S4.
# Writes ONLY under results/S5_pathway_scores/ (+ registry on first run)
# Usage: Rscript scripts/S5_01_pathway_scores.R
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
source("scripts/S4_common.R")

OUT <- "results/S5_pathway_scores"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

AXES <- list(
  cytotoxicity         = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GNLY"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
  senescence_dormancy  = c("CDKN1A", "CDKN2A", "GADD45A", "SERPINE1"),
  proliferation        = c("MKI67", "TOP2A", "PCNA"),
  stromal_ecm          = c("VIM", "COL1A1", "COL3A1", "FN1", "ACTA2")
)
MIN_AXIS_GENES <- 2L
SENMAYO_CSV <- file.path(CFG$out_scores, "S3_ssgsea_senmayo_bulk.csv")
SENMAYO_HDR <- c("cohort", "sample", "gsm", "condition", "senmayo_ssgsea")
REGISTRY <- "scripts/S5_input_registry.tsv"
ALL_COHORTS <- names(CFG$bulk_rds)

say("== S5_01: pathway scores (z-mean axes + locked SenMayo reuse) ==")

## ---- registry over true inputs (loud create / verify / halt) ---------------
reg_inputs <- c(SENMAYO_CSV, CFG$frozen, unname(CFG$bulk_rds))
md5_of <- function(p) unname(tools::md5sum(p))
if (file.exists(REGISTRY)) {
  prev <- read.csv(REGISTRY, stringsAsFactors = FALSE, check.names = FALSE)
  cur <- data.frame(file = reg_inputs, md5 = vapply(reg_inputs, md5_of, character(1)),
                    stringsAsFactors = FALSE)
  m <- merge(prev, cur, by = "file", suffixes = c("_prev", "_cur"), all = TRUE)
  bad <- m$md5_prev != m$md5_cur | is.na(m$md5_prev) | is.na(m$md5_cur)
  if (any(bad))
    stop("FATAL: input registry mismatch (R5). Changed/new files:\n  ",
         paste(m$file[bad], collapse = "\n  "),
         "\n  Route: document an Amendment, do not edit results by hand (R4/R5).")
  say("registry verified:", nrow(m), "inputs match", REGISTRY)
} else {
  cur <- data.frame(file = reg_inputs, md5 = vapply(reg_inputs, md5_of, character(1)),
                    stringsAsFactors = FALSE)
  if (any(!file.exists(reg_inputs)))
    stop("FATAL: missing inputs: ",
         paste(reg_inputs[!file.exists(reg_inputs)], collapse = ", "))
  write.csv(cur, REGISTRY, row.names = FALSE)
  say("registry CREATED (first run):", REGISTRY, "--", nrow(cur), "inputs")
}

## ---- SenMayo locked scores + schema guard ----------------------------------
if (!file.exists(SENMAYO_CSV)) stop("FATAL: missing ", SENMAYO_CSV, " (run S3_02)")
sen <- read.csv(SENMAYO_CSV, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(sen), SENMAYO_HDR))
  stop("FATAL: SenMayo CSV schema drift (R2). Expected: ",
       paste(SENMAYO_HDR, collapse = ","), " | got: ",
       paste(names(sen), collapse = ","),
       ". The locked S3_02 output changed -- investigate before proceeding.")
say("SenMayo schema verified:", paste(SENMAYO_HDR, collapse = ","))

## ---- per-cohort scores ------------------------------------------------------
cov_rows <- list()
for (ch in ALL_COHORTS) {
  X <- load_bulk_symbol_matrix(ch)          # v3.0.2: probe->symbol, loud
  n_s <- ncol(X)
  sd_gene <- apply(X, 1, stats::sd)
  if (any(sd_gene == 0))
    say("  ", ch, ": genes with sd == 0 (constant, cannot z):",
        sum(sd_gene == 0), "-- excluded from axis means, counted below")
  Z <- sweep(sweep(X, 1, rowMeans(X), "-"), 1, sd_gene, "/")
  Z[sd_gene == 0, ] <- NA

  sc <- matrix(NA_real_, nrow = length(AXES) + 1L, ncol = n_s,
               dimnames = list(c(names(AXES), "SenMayo"), colnames(X)))
  for (ax in names(AXES)) {
    g <- AXES[[ax]]
    hit <- g %in% rownames(Z)
    zable <- hit                      # same length as g; filled only where present
    zable[hit] <- sd_gene[g[hit]] > 0 # index by NAME: absent genes never recycled
    n_use <- sum(zable)
    cov_rows[[paste(ch, ax, sep = "|")]] <- data.frame(
      cohort = ch, pathway = ax, n_in_set = length(g),
      n_present = sum(hit), n_used = n_use,
      n_sdzero_dropped = sum(hit) - n_use, stringsAsFactors = FALSE)
    if (n_use < MIN_AXIS_GENES)
      stop("FATAL: ", ch, " / ", ax, " has only ", n_use,
           " z-able genes (< ", MIN_AXIS_GENES, " structural floor). ",
           "A pathway score over <2 genes is a single-gene readout (R9).")
    sc[ax, ] <- colMeans(Z[g[zable], , drop = FALSE], na.rm = TRUE)
    say("  ", ch, ax, "coverage:", n_use, "of", length(g), "genes")
  }
  ## SenMayo: exact sample match against this cohort's symbol matrix columns
  se <- sen[sen$cohort == ch, , drop = FALSE]
  if (!setequal(se$sample, colnames(X)))
    stop("FATAL: ", ch, " SenMayo sample mismatch (R2/R9). In matrix not CSV: ",
         paste(setdiff(colnames(X), se$sample), collapse = ","),
         " | in CSV not matrix: ", paste(setdiff(se$sample, colnames(X)), collapse = ","))
  sc["SenMayo", ] <- se$senmayo_ssgsea[match(colnames(X), se$sample)]
  cov_rows[[paste(ch, "SenMayo", sep = "|")]] <- data.frame(
    cohort = ch, pathway = "SenMayo", n_in_set = NA_integer_, n_present = NA_integer_,
    n_used = nrow(se), n_sdzero_dropped = NA_integer_, stringsAsFactors = FALSE)

  out <- data.frame(pathway = rownames(sc), sc, check.names = FALSE,
                    row.names = NULL, stringsAsFactors = FALSE)
  write.csv(out, file.path(OUT, paste0(ch, ".csv")), row.names = FALSE)
  say(ch, "scores written:", nrow(sc), "pathways x", n_s, "samples")
}
cov <- do.call(rbind, cov_rows)
write.csv(cov, file.path(OUT, "S5_01_coverage_audit.csv"), row.names = FALSE)
say("coverage audit:", file.path(OUT, "S5_01_coverage_audit.csv"))
say("S5_01 DONE")
