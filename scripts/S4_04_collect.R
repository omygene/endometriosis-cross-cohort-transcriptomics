# ============================================================================
# S4_04_collect.R -- immunoediting synthesis + knock summary (S4_design
# section 8). Pre-registered directional axes [Dunn 2002; Schreiber 2011;
# Mittal 2014; Symons 2018] tested on the locked Tier-1 core with exact
# binomial tests. Escape Score = fraction of Tier-1 immune-axis members that
# are DOWN (immune evasion readout). Every number here is MEASURED from the
# locked outputs -- nothing is fabricated (R1).
# Usage: Rscript scripts/S4_04_collect.R
# ============================================================================
set.seed(42)
source("scripts/S4_common.R")

t1_file <- file.path(CFG4$out_meta, "S4_meta_Tier1_core.csv")
if (!file.exists(t1_file))
  stop("FATAL: run S4_01 first -- missing ", t1_file, " (R2)")
t1 <- read.csv(t1_file, stringsAsFactors = FALSE)
say("== S4_04: immunoediting synthesis | Tier-1 core n =", nrow(t1))

## pre-registered axes (gene symbols; frozen here, R5)
AXES <- list(
  cytotoxicity = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GNLY"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
  senescence_dormancy = c("CDKN1A", "CDKN2A", "GADD45A", "SERPINE1"),
  proliferation = c("MKI67", "TOP2A", "PCNA"),
  stromal_ecm = c("VIM", "COL1A1", "COL3A1", "FN1", "ACTA2")
)

axis_rows <- list()
for (ax in names(AXES)) {
  hit <- intersect(AXES[[ax]], t1$gene)
  if (!length(hit)) {
    axis_rows[[ax]] <- data.frame(axis = ax, n_axis_genes = length(AXES[[ax]]),
                                  n_in_Tier1 = 0L, n_up = NA_integer_, n_down = NA_integer_,
                                  binom_p = NA_real_, stringsAsFactors = FALSE)
    next
  }
  sgn <- sign(t1$M[match(hit, t1$gene)])
  n_up <- sum(sgn > 0); n_dn <- sum(sgn < 0); n <- n_up + n_dn
  p <- if (n > 0) min(1, 2 * stats::pbinom(min(n_up, n_dn), n, 0.5)) else NA_real_
  axis_rows[[ax]] <- data.frame(axis = ax, n_axis_genes = length(AXES[[ax]]),
                                n_in_Tier1 = length(hit), n_up = n_up, n_down = n_dn,
                                binom_p = p, stringsAsFactors = FALSE)
}
syn <- do.call(rbind, axis_rows)
write.csv(syn, file.path(CFG4$out_meta, "S4_immunoediting_synthesis.csv"),
          row.names = FALSE)
say("wrote S4_immunoediting_synthesis.csv")
print(syn)

## Escape Score: among Tier-1 members of the two immune axes, fraction DOWN
imm <- c(AXES$cytotoxicity, AXES$antigen_presentation)
imm_hit <- intersect(imm, t1$gene)
esc <- if (length(imm_hit)) {
  sum(sign(t1$M[match(imm_hit, t1$gene)]) < 0) / length(imm_hit)
} else NA_real_
say("Escape Score (immune-axis Tier-1 fraction DOWN):",
    if (is.na(esc)) "no immune genes in Tier-1" else round(esc, 3))

## ---- collect: verify every expected artifact exists (fail-loud) ------------
expected <- c("S4_meta_all_genes.csv", "S4_meta_Tier1_core.csv",
              "S4_meta_low_presence.csv", "S4_gse11691_concordance.csv",
              "S4_sensitivity_summary.csv", "S4_scrna_projection.csv",
              "S4_scrna_lineage_readout.csv", "S4_immunoediting_synthesis.csv")
paths <- file.path(CFG4$out_meta, expected)
missing <- paths[!file.exists(paths)]
if (length(missing))
  stop("FATAL: S4 artifacts missing -- run the earlier stages first:\n  ",
       paste(missing, collapse = "\n  "))

conc  <- read.csv(file.path(CFG4$out_meta, "S4_gse11691_concordance.csv"),
                  stringsAsFactors = FALSE)
sens  <- read.csv(file.path(CFG4$out_meta, "S4_sensitivity_summary.csv"),
                  stringsAsFactors = FALSE)
proj  <- read.csv(file.path(CFG4$out_meta, "S4_scrna_projection.csv"),
                  stringsAsFactors = FALSE)
md5s  <- vapply(paths, function(x) unname(tools::md5sum(x)), character(1))

md <- c(
  "# S4 knock summary -- cross-cohort random-effects meta-analysis",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Design",
  "- Random-effects REML meta across the 6 locked primary disease-vs-control",
  "  contrasts; GSE11691 (paired ectopic-vs-eutopic) analysed as a separate",
  "  directional concordance arm (estimand separation).",
  "- Presence rule k >= 4; Tier-1 = BH-FDR < 0.05 AND |M| >= 0.5 AND sign",
  "  concordance >= 4 cohorts. Heterogeneity (tau2, I2, Q) reported for all.",
  "",
  "## Primary results (measured)",
  paste0("- Genes entering meta (k >= 4): ", sum(!is.na(read.csv(paths[1])$M))),
  paste0("- Tier-1 core genes: ", nrow(t1)),
  paste0("- Tier-2 (FDR-only) genes: ",
         sum(read.csv(paths[1], stringsAsFactors = FALSE)$tier == "Tier2")),
  paste0("- Tier-1 up / down: ", sum(t1$M > 0), " / ", sum(t1$M < 0)),
  paste0("- Median tau2 across Tier-1: ",
         round(stats::median(t1$tau2, na.rm = TRUE), 4)),
  "",
  "## GSE11691 concordance arm",
  paste0("- Common genes: ", conc$n_common_genes,
         " | spearman rho = ", round(conc$spearman_rho, 3),
         " | Tier-1 sign-match = ", round(conc$tier1_sign_match_rate, 3),
         " (", conc$tier1_genes_in_arm, " genes)"),
  "",
  "## Sensitivity battery",
  apply(sens, 1, function(r) paste0("- ", r["test"], ": ", r["metric"],
                                   " = ", r["value"], " | criterion ", r["criterion"],
                                   " | ", if (isTRUE(as.logical(r["pass"]))) "PASS" else "FAIL")),
  "",
  "## scRNA projection (bulk core -> scRNA pseudobulk legs)",
  paste0("- Legs analysed: ", nrow(proj),
         " | median sign-match rate: ",
         round(stats::median(proj$sign_match_rate), 3)),
  paste0("- Legs with FDR < 0.05 (Wilcoxon): ",
         sum(proj$FDR < CFG4$fdr_gate, na.rm = TRUE), " of ", nrow(proj)),
  "",
  "## Immunoediting synthesis",
  apply(syn, 1, function(r) paste0("- ", r["axis"], ": ", r["n_in_Tier1"], "/",
                                   r["n_axis_genes"], " in Tier-1 (up ", r["n_up"],
                                   " / down ", r["n_down"], ")")),
  paste0("- Escape Score: ", if (is.na(esc)) "NA" else round(esc, 3)),
  "",
  "## Artifact MD5s",
  paste0("- ", expected, ": ", md5s)
)
writeLines(md, file.path(CFG4$out_meta, "S4_knock_summary.md"))
say("wrote S4_knock_summary.md")
say("S4_04 DONE")
