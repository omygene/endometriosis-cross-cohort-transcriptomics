# ============================================================================
# S4_02_sensitivity.R -- robustness battery
# v1.1 (2026-09-22): + audit_cohorts() gate in every leg (ledger #18). (S4_design section 6; all legs
# published). Requires S4_01 output (S4_meta_Tier1_core.csv) as reference.
#   1 sans-GSE7307        5-cohort meta vs base Tier-1 (criterion: >= 70% kept)
#   2 LOCO x 6            leave-one-cohort-out (criterion: no removal > 50% loss)
#   3 FE vs RE            fixed-effects comparison on the same 6-cohort set
#   4 QC swap x 2         B-G51981 -> B-G51981-S1 / -S2 (Jaccard vs base)
# Writes results/S4_meta/S4_sensitivity_summary.csv (+ per-leg gene tables)
# Usage: Rscript scripts/S4_02_sensitivity.R
# ============================================================================
set.seed(42)
source("scripts/S4_common.R")
require_metafor()
dir.create(CFG4$out_meta, recursive = TRUE, showWarnings = FALSE)

base_file <- file.path(CFG4$out_meta, "S4_meta_Tier1_core.csv")
if (!file.exists(base_file))
  stop("FATAL: run S4_01 first -- missing ", base_file, " (R2)")
base <- read.csv(base_file, stringsAsFactors = FALSE)
base_t1 <- base$gene
say("== S4_02: sensitivity battery | base Tier-1 n =", length(base_t1))

run_set <- function(ids, tag, method = "REML") {
  say("-- leg", tag, "| cohorts:", paste(ids, collapse = ", "), "|", method)
  long <- do.call(rbind, lapply(ids,
                                function(id) se_from_contrast(load_one_contrast(id), id)))
  audit_cohorts(long, ids)
  m <- assign_tiers(meta_genes(build_gene_list(long), method = method))
  write.csv(m, file.path(CFG4$out_meta, paste0("S4_sens_", tag, "_all.csv")),
            row.names = FALSE)
  m
}
retention <- function(t1_new) {
  if (!length(base_t1)) return(NA_real_)
  length(intersect(base_t1, t1_new)) / length(base_t1)
}

rows <- list()

## -- leg 1: sans-GSE7307 -----------------------------------------------------
sans <- run_set(setdiff(CFG4$meta_set, CFG4$sans_drop), "sansG7307")
r1 <- retention(sans$gene[sans$tier == "Tier1"])
rows$sansG7307 <- data.frame(test = "sans-GSE7307", detail = "5-cohort meta w/o B-G7307",
                             metric = "Tier-1 retention vs base", value = r1,
                             criterion = ">= 0.70", pass = !is.na(r1) && r1 >= 0.70,
                             stringsAsFactors = FALSE)

## -- leg 2: LOCO --------------------------------------------------------------
for (drop in CFG4$meta_set) {
  ms <- run_set(setdiff(CFG4$meta_set, drop), paste0("LOCO_", drop))
  r <- retention(ms$gene[ms$tier == "Tier1"])
  rows[[paste0("LOCO_", drop)]] <- data.frame(
    test = paste0("LOCO:", drop), detail = "leave-one-cohort-out",
    metric = "Tier-1 retention vs base", value = r,
    criterion = ">= 0.50", pass = !is.na(r) && r >= 0.50, stringsAsFactors = FALSE)
}

## -- leg 3: FE vs RE on the full 6-cohort set ---------------------------------
re_df <- read.csv(file.path(CFG4$out_meta, "S4_meta_all_genes.csv"),
                  stringsAsFactors = FALSE)
fe <- run_set(CFG4$meta_set, "FE", method = "FE")
fe_t1 <- fe$gene[fe$tier == "Tier1"]
both_g <- intersect(re_df$gene, fe$gene)
rho_fe <- suppressWarnings(stats::cor(re_df$M[match(both_g, re_df$gene)],
                                      fe$M[match(both_g, fe$gene)],
                                      method = "spearman", use = "complete.obs"))
rows$FEvsRE <- data.frame(test = "FE-vs-RE", detail = "fixed-effects on same 6-cohort set",
                          metric = "spearman(M_RE, M_FE) | Tier-1 Jaccard",
                          value = paste(round(rho_fe, 3), round(jaccard(base_t1, fe_t1), 3), sep = " | "),
                          criterion = "rho > 0.95", pass = !is.na(rho_fe) && rho_fe > 0.95,
                          stringsAsFactors = FALSE)

## -- leg 4: QC swap (B-G51981 -> S1 / -S2; same control definition) -----------
for (sw in CFG4$qc_swap) {
  ids <- c(setdiff(CFG4$meta_set, "B-G51981"), sw)
  ms <- run_set(ids, paste0("swap_", sw))
  r <- retention(ms$gene[ms$tier == "Tier1"])
  rows[[paste0("swap_", sw)]] <- data.frame(
    test = paste0("QC-swap:", sw), detail = "B-G51981 replaced by QC leg",
    metric = "Tier-1 retention vs base", value = r,
    criterion = ">= 0.50 (reported)", pass = TRUE, stringsAsFactors = FALSE)
}

summ <- do.call(rbind, rows)
write.csv(summ, file.path(CFG4$out_meta, "S4_sensitivity_summary.csv"),
          row.names = FALSE)
say("== sensitivity summary ==")
print(summ)
say("S4_02 DONE")
