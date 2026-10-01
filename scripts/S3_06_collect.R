# ============================================================================
# S3_06_collect.R -- completeness + concordance collector for stage S3
# Descriptive only: NO new gates here. Reports the locked gate (BH-FDR < 0.05)
# counts, the |log2FC| >= 0.5 interpretive tier, and pre-registered
# concordance/sensitivity summaries:
#   bulk sign concordance across cohorts (+ sans-GSE7307 arm, E2 reservation)
#   GSE51981 sensitivity arms vs primary (S1/S2/S3)
#   GSE214411 N1 arms (A vs B) agreement
# Writes results/S3_knock_summary.md
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_06_collect.R > logs\S3_06_collect.txt 2>&1
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")

say("== S3_06 collect ==")
fr <- load_frozen()
cdir <- CFG$out_contrasts
files <- file.path(cdir, paste0(fr$contrast_id[fr$design %in% c("unpaired", "paired")], ".csv"))
missing <- files[!file.exists(files)]
if (length(missing)) stop("FATAL: missing bulk contrast outputs: ", paste(missing, collapse = ", "))
bulk <- lapply(files, read.csv, stringsAsFactors = FALSE)
names(bulk) <- fr$contrast_id[fr$design %in% c("unpaired", "paired")]
summ <- read.csv(file.path(cdir, "S3_bulk_contrasts_summary.csv"), stringsAsFactors = FALSE)

## --- cross-cohort sign concordance among primary bulk contrasts -----------
primary <- c("B-G7305", "B-G6364", "B-G25628", "B-G7307", "B-G51981", "B-G11691", "B-G120103")
sig <- lapply(bulk[primary], function(z) {
  zz <- z[order(z$gene), ]; zz$sign <- sign(zz$logFC)
  zz[zz$adj.P.Val < CFG$fdr_gate, c("gene", "sign", "logFC")]
})
sig_genes <- table(unlist(lapply(sig, `[[`, "gene")))
multi <- names(sig_genes)[sig_genes >= 2]
conc <- vapply(multi, function(g) {
  s <- vapply(sig, function(z) z$sign[z$gene == g][1], numeric(1))
  abs(sum(s, na.rm = TRUE)) == length(na.omit(s))
}, logical(1))
sans7307 <- setdiff(primary, "B-G7307")
sig7 <- lapply(bulk[sans7307], function(z) {
  zz <- z[z$adj.P.Val < CFG$fdr_gate, ]; zz$sign <- sign(zz$logFC); zz[, c("gene", "sign")]
})
g7 <- table(unlist(lapply(sig7, `[[`, "gene"))); m7 <- names(g7)[g7 >= 2]
conc7 <- if (length(m7)) mean(vapply(m7, function(g) {
  s <- vapply(sig7, function(z) z$sign[z$gene == g][1], numeric(1))
  abs(sum(s, na.rm = TRUE)) == length(na.omit(s))
}, logical(1))) else NA_real_

jacc <- function(a, b) { u <- union(a, b); if (!length(u)) return(NA_real_); length(intersect(a, b)) / length(u) }
g519 <- function(id) bulk[[id]]$gene[bulk[[id]]$adj.P.Val < CFG$fdr_gate]
arms519 <- data.frame(arm = c("S1", "S2", "S3"),
                      jaccard_vs_primary = c(jacc(g519("B-G51981"), g519("B-G51981-S1")),
                                             jacc(g519("B-G51981"), g519("B-G51981-S2")),
                                             jacc(g519("B-G51981"), g519("B-G51981-S3"))))

## --- scRNA N1 arms agreement ----------------------------------------------
pb <- list.files("results/S3_scrna_pseudobulk", pattern = "_A\\.csv$", full.names = TRUE)
pbB <- list.files("results/S3_scrna_pseudobulk", pattern = "_B\\.csv$", full.names = TRUE)
n1 <- data.frame()
if (length(pb) && length(pbB)) {
  a <- gsub("_A\\.csv$", "", basename(pb)); b <- gsub("_B\\.csv$", "", basename(pbB))
  common <- intersect(a, b)
  n1 <- data.frame(leg = common, stringsAsFactors = FALSE)
  n1$rankcorr_logFC <- vapply(common, function(k) {
    x <- read.csv(file.path("results/S3_scrna_pseudobulk", paste0(k, "_A.csv")), stringsAsFactors = FALSE)
    y <- read.csv(file.path("results/S3_scrna_pseudobulk", paste0(k, "_B.csv")), stringsAsFactors = FALSE)
    gg <- intersect(x$gene, y$gene)
    if (length(gg) < 10) return(NA_real_)
    cor(x$logFC[match(gg, x$gene)], y$logFC[match(gg, y$gene)], method = "spearman")
  }, numeric(1))
}

zm <- file.path("results/S3_scrna_zmean_tests.csv")
wx <- length(list.files("results/S3_scrna_wilcoxon_sens", pattern = "\\.csv$"))

L <- c("# S3 knock summary (descriptive; locked gate BH-FDR < 0.05, tier |log2FC| >= 0.5)",
       paste0("- bulk contrasts completed: ", nrow(summ), " / ", nrow(fr)),
       paste0("- genes FDR<0.05 in >=2 primary cohorts: ", length(multi),
              " | unanimous sign: ", sum(conc, na.rm = TRUE),
              " (", if (length(multi)) round(100 * mean(conc), 1) else NA, "%)"),
       paste0("- sans-GSE7307 concordance (E2 sensitivity): genes>=2: ", length(m7),
              " | unanimous: ", if (length(m7)) round(100 * conc7, 1) else NA, "%"),
       paste0("- GSE51981 arm Jaccard vs primary: S1=", round(arms519$jaccard_vs_primary[1], 3),
              " S2=", round(arms519$jaccard_vs_primary[2], 3),
              " S3=", round(arms519$jaccard_vs_primary[3], 3)),
       paste0("- scRNA pseudobulk legs (A): ", length(pb), " | N1 arms (B): ", length(pbB),
              if (nrow(n1)) paste0(" | median arm-A/B Spearman = ", round(median(n1$rankcorr_logFC, na.rm = TRUE), 3)) else ""),
       paste0("- Wilcoxon sensitivity legs: ", wx, " (sensitivity only, pseudoreplication caveat)"),
       paste0("- z-mean tests table present: ", file.exists(zm)),
       "", "## per-contrast locked-gate counts", "")
for (i in seq_len(nrow(summ)))
  L <- c(L, paste0("- ", summ$contrast_id[i], " [", summ$level[i], "] n_test=", summ$n_test[i],
                   " n_ref=", summ$n_ref[i], " FDR<0.05: ", summ$n_FDR[i],
                   " (+tier): ", summ$n_FDR_tier[i]))
writeLines(L, "results/S3_knock_summary.md")
say("wrote results/S3_knock_summary.md")
print(summ)
say("S3_06 COLLECT DONE")
