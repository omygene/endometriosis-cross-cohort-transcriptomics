# ============================================================================
# S4_01_meta_bulk.R -- PRIMARY cross-cohort random-effects meta-analysis (REML)
# v1.1 (2026-09-22): + audit_cohorts() gate (ledger #18: silent phantom cohort).
# over the 6 locked primary disease-vs-control contrasts + the separate
# GSE11691 directional concordance arm.
# Locked criteria (S4_design sections 3-5; R5): k_min = 4 presence rule;
# Tier-1 = BH-FDR < 0.05 & |M| >= 0.5 & sign concordance >= 4 cohorts;
# heterogeneity (tau2, I2, Q) reported for every gene, never filtered.
# Writes ONLY under results/S4_meta/ (+ the frozen input registry on first run)
# Usage: Rscript scripts/S4_01_meta_bulk.R
# ============================================================================
set.seed(42)
source("scripts/S4_common.R")
require_metafor()
dir.create(CFG4$out_meta, recursive = TRUE, showWarnings = FALSE)

say("== S4_01: cross-cohort REML meta (6 primary cohorts) ==")

## -- frozen input registry: 6 primary + paired arm + 2 QC swap legs ---------
verify_or_create_registry(c(CFG4$meta_set, CFG4$paired_arm, CFG4$qc_swap))

say("-- primary meta set:", paste(CFG4$meta_set, collapse = ", "))
long <- do.call(rbind, lapply(CFG4$meta_set,
                              function(id) se_from_contrast(load_one_contrast(id), id)))
gene_list <- build_gene_list(long)
say("genes observed in any primary cohort:", length(gene_list))

audit_cohorts(long, CFG4$meta_set)
m <- meta_genes(gene_list, method = "REML")
m <- assign_tiers(m)
m <- m[order(m$FDR, -abs(m$M)), ]
write.csv(m, file.path(CFG4$out_meta, "S4_meta_all_genes.csv"), row.names = FALSE)
t1 <- m[m$tier == "Tier1", ]
write.csv(t1, file.path(CFG4$out_meta, "S4_meta_Tier1_core.csv"), row.names = FALSE)
say("  wrote S4_meta_all_genes.csv |", nrow(m), "genes | Tier-1:", nrow(t1),
    "| Tier-2:", sum(m$tier == "Tier2"))

## low-presence genes: observed but k < k_min -- reported, never hidden (sec 5)
ks <- vapply(gene_list, nrow, integer(1))
low <- data.frame(gene = names(ks)[ks < CFG4$k_min], k = ks[ks < CFG4$k_min],
                  stringsAsFactors = FALSE)
write.csv(low, file.path(CFG4$out_meta, "S4_meta_low_presence.csv"), row.names = FALSE)
say("  wrote S4_meta_low_presence.csv |", nrow(low), "genes with k <", CFG4$k_min)

## -- GSE11691 directional concordance arm (estimand separation, sec 3) ------
pa <- se_from_contrast(load_one_contrast(CFG4$paired_arm), CFG4$paired_arm)
common <- intersect(pa$gene, m$gene[!is.na(m$M)])
bx <- m$M[match(common, m$gene)]
by <- pa$logFC[match(common, pa$gene)]
rho <- suppressWarnings(stats::cor(bx, by, method = "spearman",
                                   use = "complete.obs"))
both <- intersect(t1$gene, pa$gene)
sign_match <- mean(sign(pa$logFC[match(both, pa$gene)]) ==
                   sign(m$M[match(both, m$gene)]))
conc <- data.frame(n_common_genes = length(common),
                   spearman_rho = rho,
                   tier1_genes_in_arm = length(both),
                   tier1_sign_match_rate = sign_match,
                   stringsAsFactors = FALSE)
write.csv(conc, file.path(CFG4$out_meta, "S4_gse11691_concordance.csv"),
          row.names = FALSE)
say("  GSE11691 arm:", length(common), "common genes | spearman rho =",
    round(rho, 3), "| Tier-1 sign-match:", round(sign_match, 3),
    "(", length(both), "of", nrow(t1), ")")

say("S4_01 DONE")
