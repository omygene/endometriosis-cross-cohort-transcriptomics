# ============================================================================
# S5_02_pathway_meta.R -- pathway-level random-effects meta-analysis (REML)
# v1.0.1 (2026-09-22): registry read fixed (comma file read as tab, ledger #20,
#   same fix as S5_01 v1.0.1; my error, R6).
# v1.0 (2026-09-22, Amendment A4, design section 14).
# Score-level contrasts reuse run_frozen_contrast() VERBATIM from S3_common:
# the frozen design itself (GSE6364 phase covariate, frozen exclusions,
# duplicateCorrelation pairing for GSE11691) applied to the 6x-sample score
# matrix -- zero new contrast logic (referee defensibility).
# Locked rules (R5, identical to S4):
#   * k_min = 4 presence rule per pathway; SE = |logFC/t| (locked derivation).
#   * BH-FDR computed ACROSS the 6 pathway tests only.
#   * |M| >= 0.5 interpretive tier, NEVER an exclusion.
#   * tau2 / I2 / Q reported for every pathway, never filtered.
#   * B-G11691 paired arm: directional concordance only, NEVER pooled.
#   * Sensitivity: sans-G7307, LOCO x6, FE-vs-RE, QC-swap S1/S2.
# Writes ONLY under results/S5_meta/
# Usage: Rscript scripts/S5_02_pathway_meta.R
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
source("scripts/S4_common.R")
require_metafor()
dir.create(CFG4$out_meta, showWarnings = FALSE)
OUT <- "results/S5_meta"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REGISTRY <- "scripts/S5_input_registry.tsv"

PATHWAYS <- c("cytotoxicity", "antigen_presentation", "senescence_dormancy",
              "proliferation", "stromal_ecm", "SenMayo")

say("== S5_02: pathway-level REML meta over the 6 locked cohorts ==")

## ---- input verification (R2/R5) --------------------------------------------
if (!file.exists(REGISTRY)) stop("FATAL: ", REGISTRY, " absent -- run S5_01 first")
prev <- read.csv(REGISTRY, stringsAsFactors = FALSE, check.names = FALSE)
cur <- data.frame(file = prev$file, md5 = unname(tools::md5sum(prev$file)),
                  stringsAsFactors = FALSE)
if (any(!file.exists(prev$file)))
  stop("FATAL: registry input(s) vanished: ",
       paste(prev$file[!file.exists(prev$file)], collapse = ", "))
if (any(prev$md5 != cur$md5))
  stop("FATAL: registry input changed since S5_01 (R5): ",
       paste(prev$file[prev$md5 != cur$md5], collapse = ", "),
       " -- document an Amendment, never hand-edit results (R4).")
say("registry inputs verified (S5_01 state frozen)")

fr <- load_frozen()
need_ids <- c(CFG4$meta_set, CFG4$paired_arm, CFG4$qc_swap)
miss_spec <- setdiff(need_ids, fr$contrast_id)
if (length(miss_spec))
  stop("FATAL: frozen registry lacks: ", paste(miss_spec, collapse = ", "))

## ---- score contrasts through the FROZEN contrast machinery -----------------
score_contrast <- function(id) {
  spec <- fr[fr$contrast_id == id, , drop = FALSE]
  if (nrow(spec) != 1L) stop("FATAL: contrast id not unique: ", id)
  ch <- spec$cohort
  p <- file.path("results/S5_pathway_scores", paste0(ch, ".csv"))
  if (!file.exists(p)) stop("FATAL: missing ", p, " -- run S5_01 first (R2)")
  sc <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(sc$pathway, PATHWAYS))
    stop("FATAL: ", p, " pathway rows drifted (R2). Expected:\n  ",
         paste(PATHWAYS, collapse = ","), "\n  got:\n  ",
         paste(sc$pathway, collapse = ","))
  mat <- as.matrix(sc[, -1, drop = FALSE])
  rownames(mat) <- sc$pathway
  m <- map_cols(ch, colnames(mat), load_meta())
  tt <- run_frozen_contrast(mat, m, spec)   # FROZEN design, verbatim reuse
  se_from_contrast(tt, id)                  # locked SE = |logFC/t|
}

long <- do.call(rbind, lapply(need_ids, score_contrast))

## presence audit: exactly one estimate per pathway per cohort (R9)
tab <- table(long$gene, long$contrast_id)
say("presence table (pathways x contrasts):")
print(tab)
for (cid in CFG4$meta_set) {
  n <- sum(long$contrast_id == cid)
  if (n != length(PATHWAYS))
    stop("FATAL: ", cid, " contributed ", n, " pathway estimates, expected ",
         length(PATHWAYS), " -- score matrix or frozen design inconsistent (R9).")
}

## ---- primary REML meta per pathway ------------------------------------------
meta_set <- function(long, ids, method = "REML") {
  sub <- long[long$contrast_id %in% ids, , drop = FALSE]
  pl <- split(sub, sub$gene)
  res <- meta_genes(pl, method = method, k_min = CFG4$k_min)  # pathway == 'gene'
  names(res)[names(res) == "gene"] <- "pathway"
  res <- res[match(PATHWAYS, res$pathway), ]
  assign_tiers(res)   # BH-FDR across the 6 pathway tests + interpretive tiers
}
m_base <- meta_set(long, CFG4$meta_set, "REML")
write.csv(m_base, file.path(OUT, "S5_meta_pathways.csv"), row.names = FALSE)
say("base meta written | significant pathways (FDR < 0.05):",
    paste(m_base$pathway[m_base$FDR < CFG4$fdr_gate], collapse = ", "))

## ---- GSE11691 paired concordance arm (never pooled) -------------------------
pa <- long[long$contrast_id == CFG4$paired_arm, ]
pa <- pa[match(PATHWAYS, pa$gene), ]
both <- !is.na(m_base$M) & !is.na(pa$logFC)
sign_match <- mean(sign(m_base$M[both]) == sign(pa$logFC[both]))
rho <- suppressWarnings(stats::cor.test(m_base$M[both], pa$logFC[both],
                                       method = "spearman")$estimate)
conc <- data.frame(pathway = PATHWAYS, M_base = m_base$M,
                   logFC_GSE11691 = pa$logFC, same_sign = sign(m_base$M) == sign(pa$logFC),
                   stringsAsFactors = FALSE)
write.csv(conc, file.path(OUT, "S5_concordance_GSE11691.csv"), row.names = FALSE)
say("GSE11691 concordance: sign-match =", round(sign_match, 3),
    "| Spearman rho =", round(unname(rho), 3))

## ---- sensitivity battery (mirror of S4_02, design section 14) ---------------
sig_set <- function(m) m$pathway[!is.na(m$FDR) & m$FDR < CFG4$fdr_gate]
base_sig <- sig_set(m_base)
say("base significant set (", length(base_sig), "):",
    paste(base_sig, collapse = ", "))

cmp_scenario <- function(name, ids) {
  m <- meta_set(long, ids, "REML")
  m <- m[match(PATHWAYS, m$pathway), ]
  s <- sig_set(m)
  data.frame(scenario = name,
             n_sig = length(s), n_overlap = length(intersect(base_sig, s)),
             sign_match = mean(sign(m_base$M) == sign(m$M), na.rm = TRUE),
             max_abs_dM = max(abs(m_base$M - m$M), na.rm = TRUE),
             FDR = paste(sprintf("%s=%.3g", m$pathway, m$FDR), collapse = "; "),
             stringsAsFactors = FALSE)
}
scen <- list()
scen[["sans_G7307"]] <- cmp_scenario("sans-G7307",
                                     setdiff(CFG4$meta_set, CFG4$sans_drop))
for (drop in CFG4$meta_set)
  scen[[paste0("LOCO_", drop)]] <- cmp_scenario(paste0("LOCO-", drop),
                                                setdiff(CFG4$meta_set, drop))
for (sw in CFG4$qc_swap)
  scen[[paste0("swap_", sw)]] <- cmp_scenario(paste0("swap-", sw),
                                              c(setdiff(CFG4$meta_set, "B-G51981"), sw))
m_fe <- meta_set(long, CFG4$meta_set, "FE")
m_fe <- m_fe[match(PATHWAYS, m_fe$pathway), ]
fe_sig <- sig_set(m_fe)
scen[["FE_vs_RE"]] <- data.frame(scenario = "FE-vs-RE",
  n_sig = length(fe_sig), n_overlap = length(intersect(base_sig, fe_sig)),
  sign_match = mean(sign(m_base$M) == sign(m_fe$M), na.rm = TRUE),
  max_abs_dM = max(abs(m_base$M - m_fe$M), na.rm = TRUE),
  FDR = paste(sprintf("%s=%.3g", m_fe$pathway, m_fe$FDR), collapse = "; "),
  stringsAsFactors = FALSE)
sens <- do.call(rbind, scen)
write.csv(sens, file.path(OUT, "S5_sensitivity_scenarios.csv"), row.names = FALSE)

## ---- sensitivity summary md -------------------------------------------------
L <- c("# S5 sensitivity summary (pathway level)",
       "",
       paste0("Base significant set (BH-FDR < 0.05 across the ", length(PATHWAYS),
              " pathway tests): ", length(base_sig), " of ", length(PATHWAYS),
              " -- ", paste(base_sig, collapse = ", ")),
       "",
       "| scenario | n sig | overlap w/ base | sign match | max |dM| |",
       "|---|---|---|---|---|")
for (i in seq_len(nrow(sens)))
  L <- c(L, paste0("| ", gsub("_", "-", sens$scenario[i]), " | ", sens$n_sig[i],
                   " | ", sens$n_overlap[i], "/", length(base_sig), " | ",
                   sprintf("%.3f", sens$sign_match[i]), " | ",
                   sprintf("%.3f", sens$max_abs_dM[i]), " |"))
L <- c(L, "", "Per-pathway FE FDRs and scenario FDRs: S5_sensitivity_scenarios.csv",
       "", "Interpretation locks (design section 14): FDR < 0.05 is the only",
       "significance claim; |M| >= 0.5 is interpretive; tau2/I2 reported, never",
       "filtered; a scenario flipping a pathway across the FDR boundary is",
       "reported as a boundary flip under heterogeneity, not hidden.")
writeLines(L, file.path(OUT, "S5_sensitivity_summary.md"))
say("wrote", file.path(OUT, "S5_sensitivity_summary.md"))
say("S5_02 DONE")
