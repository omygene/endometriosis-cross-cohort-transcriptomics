# ============================================================================
# S4_common.R -- shared helpers + FROZEN config for stage S4 (cross-cohort
# random-effects meta-analysis of the locked S3 bulk contrasts).
# v1.1 (2026-09-22): + audit_cohorts() fail-loud gate.
# v1.1.1 (2026-09-22): audit overlap metric fixed (sapply->lapply; measured
# smoke-9 false FAIL: intersect of symbols with a count number == empty). Measured failure ledger
# #18: with probe-level S3 outputs, GSE120103 silently contributed ZERO rows
# to every meta (LOCO_B-G120103 == base meta exactly, dM = 0.0). From v1.1
# any cohort contributing zero rows, or whose genes overlap < 5% with the
# union of the other cohorts (namespace mismatch), stops the run naming the
# offender (R9). Also requires S3_01 >= v2.0 symbol-level outputs.
# Design locked 2026-09-21 by lead written approval ("ابدا التنفيذ");
# amendment A2 documented in S4_design: referee-Q&A section removed, published
# references added. Every quantitative rule maps 1:1 to S4_design sections;
# constants here are FROZEN (R5) -- any change requires a documented amendment.
#
# PRIMARY META SET -- verified 2026-09-21 against scripts/S3_contrasts_frozen.tsv
# itself (the 6 primary disease-vs-control contrasts):
#   B-G7305  B-G6364  B-G25628  B-G7307  B-G51981  B-G120103
# GSE11691 (B-G11691, paired ectopic-vs-eutopic, 9 women) is a SEPARATE
# directional concordance arm -- NEVER pooled with disease-vs-control
# (estimand separation, design section 3; Ramasamy 2008, PLoS Med 5:e184).
#
# Methods backbone (S4_design section 11):
#   REML random effects  [DerSimonian & Laird 1986; Viechtbauer 2005, 2010]
#   heterogeneity Q, I2, tau2  [Cochran 1954; Higgins & Thompson 2002]
#   SE(logFC) = |logFC/t| from limma topTable  [Ritchie 2015; Law 2014]
#   BH-FDR within each meta set  [Benjamini & Hochberg 1995]
# ============================================================================

CFG4 <- list(
  seed      = 42,
  fdr_gate  = 0.05,      # locked BH-FDR gate (unchanged from S3, R5)
  lfc_tier  = 0.5,       # interpretive tier only, NO exclusion (R5)
  k_min     = 4L,        # presence rule: gene must survive QC in >= 4 cohorts
  min_concord_k = 4L,    # Tier-1 sign concordance rule (design section 4)
  in_dir    = "results/S3_bulk_contrasts",
  pb_dir    = "results/S3_scrna_pseudobulk",
  out_meta  = "results/S4_meta",
  registry  = "scripts/S4_input_registry.tsv",
  meta_set  = c("B-G7305", "B-G6364", "B-G25628", "B-G7307",
                "B-G51981", "B-G120103"),
  paired_arm = "B-G11691",
  qc_swap   = c("B-G51981-S1", "B-G51981-S2"),
  sans_drop = "B-G7307"
)

say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }

require_metafor <- function() {
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("FATAL: package 'metafor' is required for S4.\n",
         "Install once on this machine:\n",
         "  utils::install.packages(\"metafor\")\n",
         "then re-run (R3). Reference: Viechtbauer 2010, J Stat Softw 36(3):1-48.",
         call. = FALSE)
  invisible(TRUE)
}

CONTRAST_COLS <- c("contrast_id", "gene", "logFC", "AveExpr", "t", "P.Value",
                   "adj.P.Val", "B", "n_test", "n_ref")

## Read one locked S3 bulk contrast CSV and validate its frozen structure.
load_one_contrast <- function(id) {
  p <- file.path(CFG4$in_dir, paste0(id, ".csv"))
  if (!file.exists(p))
    stop("FATAL: missing locked S3 input ", p, " -- run S3_01 collect first (R2)")
  d <- read.csv(p, stringsAsFactors = FALSE)
  miss <- setdiff(CONTRAST_COLS, colnames(d))
  if (length(miss))
    stop("FATAL: ", p, " missing columns: ", paste(miss, collapse = ", "),
         " -- S3 output structure changed? escalate (R9)")
  d$contrast_id <- id
  d
}

## Derive per-gene per-contrast SE from the locked limma topTable:
##   SE(logFC) = |logFC / t|   [Ritchie 2015, NAR 43:e47; Law 2014, Genome Biol 15:R29]
## Rows with t == 0 / NA / non-finite, or non-positive SE, are DROPPED with a
## loud count (R9) -- never repaired, never imputed.
se_from_contrast <- function(d, id) {
  bad <- is.na(d$t) | d$t == 0 | is.na(d$logFC)
  if (any(bad)) say("  ", id, ": rows dropped (t==0/NA or logFC NA):",
                    sum(bad), "of", nrow(d))
  d <- d[!bad, , drop = FALSE]
  d$se <- abs(d$logFC / d$t)
  keep <- is.finite(d$se) & d$se > 0
  if (any(!keep)) say("  ", id, ": rows dropped (non-finite/zero SE):", sum(!keep))
  d[keep, c("contrast_id", "gene", "logFC", "se"), drop = FALSE]
}

## Build the gene-wise study list from a long table (rbind of
## se_from_contrast outputs): one data.frame per gene, columns
## contrast_id / gene / logFC / se.
build_gene_list <- function(long) split(long, long$gene)

## Gene-wise meta-analysis. method = "REML" (primary) or "FE" (sensitivity).
## Genes with k < k_min are EXCLUDED from fitting and reported by the caller
## via the presence rule (design section 5) -- this function counts them.
meta_genes <- function(gene_list, method = "REML", k_min = CFG4$k_min) {
  require_metafor()
  n_in <- length(gene_list)
  ks <- vapply(gene_list, nrow, integer(1))
  keep <- ks >= k_min
  gene_list <- gene_list[keep]
  say("  genes entering", method, "meta (k >=", k_min, "):",
      length(gene_list), "of", n_in, "| excluded low-presence:", sum(!keep))
  if (!length(gene_list))
    stop("FATAL: no gene reached k >= ", k_min, " -- inputs inconsistent with S3 (R9)")
  res <- lapply(names(gene_list), function(g) {
    z <- gene_list[[g]]
    fit <- tryCatch(metafor::rma(yi = z$logFC, sei = z$se, method = method),
                    error = function(e) e)
    if (inherits(fit, "error"))
      return(data.frame(gene = g, k = nrow(z), M = NA_real_, SE = NA_real_,
                        p = NA_real_, tau2 = NA_real_, I2 = NA_real_,
                        Q = NA_real_, QEp = NA_real_,
                        k_pos = sum(z$logFC > 0), k_neg = sum(z$logFC < 0),
                        stringsAsFactors = FALSE))
    data.frame(gene = g, k = nrow(z), M = unname(fit$b), SE = fit$se,
               p = fit$pval, tau2 = fit$tau2, I2 = fit$I2, Q = fit$QE,
               QEp = fit$QEp, k_pos = sum(z$logFC > 0), k_neg = sum(z$logFC < 0),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

## Locked tier assignment (design section 4). BH-FDR is computed WITHIN the
## meta set being assigned (per-set adjustment, [Benjamini & Hochberg 1995]).
assign_tiers <- function(m) {
  m$FDR <- stats::p.adjust(m$p, method = "BH")
  m$concord_k <- pmax(m$k_pos, m$k_neg)
  m$tier <- ifelse(!is.na(m$FDR) & m$FDR < CFG4$fdr_gate &
                     abs(m$M) >= CFG4$lfc_tier &
                     m$concord_k >= CFG4$min_concord_k, "Tier1",
                   ifelse(!is.na(m$FDR) & m$FDR < CFG4$fdr_gate, "Tier2", "Tier3"))
  m
}

## Jaccard similarity of two gene sets (sensitivity battery, design section 6).
jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  length(intersect(a, b)) / u
}

## Input hash registry (design section 10). Created on FIRST run (loud),
## verified on every later run; any change halts and points to the amendment
## route (R5). ids = the contrast CSVs S4 consumes.
verify_or_create_registry <- function(ids) {
  paths <- file.path(CFG4$in_dir, paste0(ids, ".csv"))
  missing <- paths[!file.exists(paths)]
  if (length(missing))
    stop("FATAL: locked S3 inputs absent:\n  ", paste(missing, collapse = "\n  "),
         "\n-- S4 reads S3 outputs only (R4)")
  cur <- data.frame(file = paths,
                    md5 = vapply(paths, function(x) unname(tools::md5sum(x)), character(1)),
                    stringsAsFactors = FALSE)
  if (!file.exists(CFG4$registry)) {
    write.table(cur, CFG4$registry, sep = "\t", quote = FALSE, row.names = FALSE)
    say("  input registry CREATED (first run):", CFG4$registry,
        "-- now FROZEN; any later input change halts here (R5)")
    return(invisible(TRUE))
  }
  prev <- read.delim(CFG4$registry, stringsAsFactors = FALSE)
  if (!identical(sort(as.character(prev$file)), sort(cur$file)))
    stop("FATAL: registry file set changed (R5).\nprev: ",
         paste(sort(as.character(prev$file)), collapse = ", "),
         "\ncur : ", paste(sort(cur$file), collapse = ", "))
  mrg <- merge(prev, cur, by = "file", suffixes = c("_prev", "_cur"))
  bad <- mrg$file[mrg$md5_prev != mrg$md5_cur]
  if (length(bad))
    stop("FATAL: locked S3 input(s) changed after registry freeze:\n  ",
         paste(bad, collapse = "\n  "),
         "\n-- if a DOCUMENTED amendment (R5) approved this change, delete ",
         CFG4$registry, " and re-run; otherwise restore the locked S3 outputs (R4)")
  say("  input registry verified:", nrow(cur), "files, all MD5 match (frozen)")
  invisible(TRUE)
}

## ---- fail-loud cohort audit (v1.1; measured ledger #18) --------------------
## Every meta run MUST verify that EACH cohort actually contributed rows and
## shares the gene namespace. Silent zero-contribution (the GSE120103 probe/
## symbol namespace split) poisoned the first S4 run without any error.
audit_cohorts <- function(long, ids) {
  say("  cohort contribution audit:")
  sets <- lapply(ids, function(id) unique(long$gene[long$contrast_id == id]))
  names(sets) <- ids
  n_rows <- vapply(ids, function(id) sum(long$contrast_id == id), integer(1))
  for (id in ids)
    say("    ", id, ": rows =", n_rows[id], "| unique genes =", length(sets[[id]]))
  zero <- ids[n_rows == 0]
  if (length(zero))
    stop("FATAL: cohort(s) contributed ZERO rows to the meta: ",
         paste(zero, collapse = ", "),
         " -- gene namespace mismatch or empty S3 output. S3_01 must be >= v2.0 ",
         "(symbol-level). Re-run S3_01, then reset the S4 registry under a ",
         "documented amendment (R5/R9).")
  ## v1.1.1 fix (measured smoke-9 false FAIL 2026-09-22): sapply() here
  ## returned COUNTS, and intersect() of gene symbols with a number silently
  ## yielded empty -> every cohort looked disjoint. lapply keeps gene SETS.
  union_others <- lapply(ids, function(id)
    unique(unlist(sets[setdiff(ids, id)])))
  names(union_others) <- ids
  for (id in ids) {
    ov <- length(intersect(sets[[id]], union_others[[id]]))
    frac <- ov / max(1, length(sets[[id]]))
    say("    ", id, ": gene overlap with other cohorts =", ov, "of",
        length(sets[[id]]), sprintf("(%.1f%%)", 100 * frac))
    if (frac < 0.05)
      stop("FATAL: ", id, " shares < 5% of its genes with the other cohorts ",
           "(", sprintf("%.1f%%", 100 * frac), ") -- namespace mismatch (probe vs ",
           "symbol). Re-run S3_01 >= v2.0 and reset the S4 registry (R5/R9).")
  }
  invisible(TRUE)
}
