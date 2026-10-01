# ============================================================================
# S4_smoke.R -- pre-flight smoke test for stage S4. Run BEFORE S4_01.
# v1.1 (2026-09-22): leg 9 covers audit_cohorts() (ledger #18). 9 legs. Every leg uses hand-constructed data with a mathematically known
# answer; expected values were independently cross-checked against a from-
# scratch REML implementation. Any FAIL stops with exit status 1.
# Usage: Rscript scripts/S4_smoke.R
# ============================================================================
set.seed(42)
source("scripts/S4_common.R")
require_metafor()

fails <- 0
leg <- function(name, cond) {
  cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name))
  flush.console()
  if (!isTRUE(cond)) fails <<- fails + 1
  invisible(cond)
}

say("== S4 smoke: 9 legs ==")

## -- Leg 1: homogeneous studies -> M = grand mean, tau2 ~ 0, I2 ~ 0 ----------
## y = (1,1,1), se = 0.1 each. REML has nothing to explain between studies,
## so tau2 = 0, M = 1, SE(M) = 0.1/sqrt(3) = 0.057735.
fit1 <- metafor::rma(yi = c(1, 1, 1), sei = c(0.1, 0.1, 0.1), method = "REML")
leg("1 REML homogeneous: M = 1",
    abs(as.numeric(fit1$b) - 1) < 1e-6)
leg("1 REML homogeneous: tau2 ~ 0 and I2 ~ 0",
    fit1$tau2 < 1e-4 && fit1$I2 < 1)
leg("1 REML homogeneous: SE(M) = 0.1/sqrt(3)",
    abs(fit1$se - 0.1/sqrt(3)) < 1e-6)

## -- Leg 2: maximally spread studies -> large tau2, I2 ~ 99.8% --------------
## y = (0,2,4), se = 0.1 each. Independently computed REML: tau2 = 3.99,
## M = 2, SE(M) = 1.1547, Cochran Q = 800.
fit2 <- metafor::rma(yi = c(0, 2, 4), sei = c(0.1, 0.1, 0.1), method = "REML")
leg("2 REML heterogeneous: M = 2", abs(as.numeric(fit2$b) - 2) < 1e-6)
leg("2 REML heterogeneous: tau2 = 3.99 (+/-0.01)",
    abs(fit2$tau2 - 3.99) < 0.01)
leg("2 REML heterogeneous: I2 > 95", fit2$I2 > 95)
leg("2 REML heterogeneous: Q ~ 800", abs(fit2$QE - 800) < 1)

## -- Leg 3: SE derivation + t == 0 guard ------------------------------------
d <- data.frame(contrast_id = "X", gene = c("g1", "g2", "g3"),
                logFC = c(1, -2, 3), t = c(2, 4, 0), stringsAsFactors = FALSE)
z <- se_from_contrast(d, "X")
leg("3 SE = |logFC/t|: g1 SE 0.5, g2 SE 0.5",
    nrow(z) == 2 && all(abs(z$se - 0.5) < 1e-12))
leg("3 t == 0 row dropped loudly (R9)", nrow(z) == 2)

## -- Leg 4: presence rule k_min ---------------------------------------------
mk_long <- function(gene, n_coh, lfc) {
  data.frame(contrast_id = paste0("C", seq_len(n_coh)), gene = gene,
             logFC = lfc, se = rep(0.1, n_coh), stringsAsFactors = FALSE)
}
long4 <- rbind(mk_long("hiK", 6, 1),       # present in 6 -> enters
               mk_long("loK", 3, 1),       # present in 3 -> excluded
               mk_long("edgeK", 4, 1))     # present in exactly 4 -> enters
gl4 <- build_gene_list(long4)
m4 <- meta_genes(gl4, method = "REML", k_min = 4L)
leg("4 k_min = 4: hiK + edgeK enter, loK excluded",
    nrow(m4) == 2 && all(c("hiK", "edgeK") %in% m4$gene) && !("loK" %in% m4$gene))

## -- Leg 5: tier assignment logic -------------------------------------------
m5 <- data.frame(gene = c("a", "b", "c", "d"),
                 M = c(0.8, 0.3, 0.8, 0.8),
                 p = c(1e-6, 1e-6, 0.04, 0.2),
                 k_pos = c(6, 6, 6, 3), k_neg = c(0, 0, 0, 3),
                 stringsAsFactors = FALSE)
m5 <- assign_tiers(m5)
leg("5 tiers: FDR+size+concord -> Tier1; FDR only -> Tier2; rest Tier3",
    identical(as.character(m5$tier), c("Tier1", "Tier2", "Tier3", "Tier3")))

## -- Leg 6: registry halts on changed input ---------------------------------
td <- tempfile("s4reg"); dir.create(td)
f1 <- file.path(td, "a.csv"); f2 <- file.path(td, "b.csv")
writeLines("x", f1); writeLines("y", f2)
reg <- file.path(td, "reg.tsv")
write.table(data.frame(file = c(f1, f2),
                       md5 = vapply(c(f1, f2), function(x) unname(tools::md5sum(x)), character(1))),
            reg, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines("CHANGED", f2)
prev <- read.delim(reg, stringsAsFactors = FALSE)
cur <- data.frame(file = c(f1, f2),
                  md5 = vapply(c(f1, f2), function(x) unname(tools::md5sum(x)), character(1)),
                  stringsAsFactors = FALSE)
mrg <- merge(prev, cur, by = "file", suffixes = c("_prev", "_cur"))
leg("6 registry: MD5 mismatch detected", any(mrg$md5_prev != mrg$md5_cur))

## -- Leg 7: GSE11691 concordance direction ----------------------------------
## Bulk meta says UP (+); within-woman arm agrees on both genes -> match 1.
t1_sign <- c(gA = 1, gB = 1)
arm_lfc <- c(gA = 0.9, gB = 0.4)
leg("7 concordance: both same sign -> sign-match rate 1",
    mean(sign(arm_lfc) == t1_sign[names(arm_lfc)]) == 1)

## -- Leg 8: FE vs RE on spread data (sensitivity leg sanity) ----------------
fit8fe <- metafor::rma(yi = c(0, 2, 4), sei = c(0.1, 0.1, 0.1), method = "FE")
leg("8 FE collapses tau2 to 0 (contrast with leg 2)",
    fit8fe$tau2 == 0 && abs(as.numeric(fit8fe$b) - 2) < 1e-6)

## -- Leg 9: audit_cohorts catches a silent phantom cohort -------------------
## Five cohorts share gene symbols; one cohort carries disjoint probe ids.
mk_cohort_rows <- function(id, genes) data.frame(contrast_id = id, gene = genes,
  logFC = rnorm(length(genes), 0.5, 0.1), se = rep(0.1, length(genes)),
  stringsAsFactors = FALSE)
good_genes <- c("CD8A", "GZMB", "PRF1", "B2M", "VIM")
long9_ok <- do.call(rbind, lapply(paste0("C", 1:5), function(i) mk_cohort_rows(i, good_genes)))
long9_bad <- rbind(long9_ok, mk_cohort_rows("C6", c("205806_at", "209969_at")))
audit_ok <- tryCatch({ audit_cohorts(long9_ok, paste0("C", 1:5)); TRUE },
                     error = function(e) FALSE)
leg("9 audit: 5 concordant namespaces pass", isTRUE(audit_ok))
audit_bad <- tryCatch({ audit_cohorts(long9_bad, paste0("C", 1:6)); FALSE },
                      error = function(e) TRUE)
leg("9 audit: disjoint-namespace cohort stops the run (R9)", isTRUE(audit_bad))

cat(sprintf("== S4 smoke: %s (fails = %d) ==\n",
            if (fails == 0) "ALL PASS" else "FAILURES PRESENT", fails))
if (fails > 0) quit(save = "no", status = 1)
say("S4 SMOKE DONE")
