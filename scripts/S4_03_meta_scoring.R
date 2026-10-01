# ============================================================================
# S4_03_meta_scoring.R -- scRNA arm
# v1.1 (2026-09-22): fail-loud guard when no leg holds >= 10 core genes
# (v1.0 died on NULL$FDR with no useful message -- measured stop 2026-09-22). (S4_design section 7). NO gene-level meta
# is attempted on scRNA (k = 3 cohorts < k_min = 4; not defensible). Instead:
#   (a) projection of the locked bulk Tier-1 core onto every S3 scRNA
#       pseudobulk leg (S3_pb_*_A.csv / *_ALL.csv): sign-match rate vs the
#       bulk meta direction + one-sample Wilcoxon signed-rank of core logFCs;
#   (b) per-lineage readout: where does the core live.
# Requires S4_01 output. Writes results/S4_meta/S4_scrna_projection.csv
# Usage: Rscript scripts/S4_03_meta_scoring.R
# ============================================================================
set.seed(42)
source("scripts/S4_common.R")

t1_file <- file.path(CFG4$out_meta, "S4_meta_Tier1_core.csv")
if (!file.exists(t1_file))
  stop("FATAL: run S4_01 first -- missing ", t1_file, " (R2)")
t1 <- read.csv(t1_file, stringsAsFactors = FALSE)
t1_sign <- sign(t1$M); names(t1_sign) <- t1$gene
say("== S4_03: scRNA projection of bulk Tier-1 core | core n =", nrow(t1))

pb_files <- sort(list.files(CFG4$pb_dir,
                            pattern = "^S3_pb_.*_(A|ALL)\\.csv$",
                            full.names = TRUE))
if (!length(pb_files))
  stop("FATAL: no S3 pseudobulk legs under ", CFG4$pb_dir, " -- run S3_03 first (R2)")
say("pseudobulk legs found:", length(pb_files))

need_cols <- c("gene", "logFC", "dataset", "lineage", "arm")
res <- list()
for (f in pb_files) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  miss <- setdiff(need_cols, colnames(d))
  if (length(miss))
    stop("FATAL: ", f, " missing columns: ", paste(miss, collapse = ", "), " (R9)")
  sub <- d[d$gene %in% t1$gene, , drop = FALSE]
  tag <- sub("^S3_pb_", "", basename(f))
  if (nrow(sub) < 10) {
    say("  skip", tag, "-- only", nrow(sub), "core genes present (reported)")
    next
  }
  sm <- mean(sign(sub$logFC) == t1_sign[sub$gene])
  wt <- tryCatch(stats::wilcox.test(sub$logFC, mu = 0),
                 error = function(e) e)
  res[[tag]] <- data.frame(
    leg = tag, dataset = d$dataset[1], lineage = d$lineage[1], arm = d$arm[1],
    n_core_genes = nrow(sub),
    sign_match_rate = sm,
    mean_core_logFC = mean(sub$logFC),
    wilcox_p = if (inherits(wt, "error")) NA_real_ else wt$p.value,
    stringsAsFactors = FALSE)
}
if (!length(res))
  stop("FATAL: every pseudobulk leg had < 10 Tier-1 core genes -- the bulk ",
       "meta gene IDs do not intersect the scRNA gene IDs (namespace mismatch, ",
       "ledger #18 family). Check audit output; S3_01 must be >= v2.0 (R9).")
out <- do.call(rbind, res)
out$FDR <- stats::p.adjust(out$wilcox_p, method = "BH")
out <- out[order(out$dataset, out$lineage, out$arm), ]
write.csv(out, file.path(CFG4$out_meta, "S4_scrna_projection.csv"), row.names = FALSE)
say("wrote S4_scrna_projection.csv |", nrow(out), "legs |",
    "median sign-match:", round(stats::median(out$sign_match_rate), 3))

## per-lineage readout across datasets
lin <- aggregate(sign_match_rate ~ lineage, data = out, FUN = mean)
colnames(lin)[2] <- "mean_sign_match"
write.csv(lin, file.path(CFG4$out_meta, "S4_scrna_lineage_readout.csv"),
          row.names = FALSE)
say("wrote S4_scrna_lineage_readout.csv |", nrow(lin), "lineages")
print(lin)
say("S4_03 DONE")
