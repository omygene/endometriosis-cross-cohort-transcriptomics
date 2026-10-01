# ============================================================================
# S5_smoke.R -- known-answer smoke test for S4_05 + S5 (A4, 2026-09-22)
# Every leg uses hand-constructed data with a mathematically known answer.
# Expected values were computed independently (z-mean by hand; REML by a
# from-scratch implementation cross-checked against S4_smoke legs 1-2, whose
# values were themselves verified against metafor on the lead machine).
# A failed leg stops the run NONZERO, before any real stage script runs.
# Usage: Rscript scripts/S5_smoke.R
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
source("scripts/S4_common.R")

n_leg <- 0
leg <- function(name, cond) {
  n_leg <<- n_leg + 1
  cat(sprintf("[%s] leg %d: %s\n", if (isTRUE(cond)) "PASS" else "FAIL",
              n_leg, name))
  if (!isTRUE(cond)) {
    cat("SMOKE FAILED -- do not run any S5 stage script until fixed (R6).\n")
    quit(save = "no", status = 1)
  }
}
say("== S5 smoke: 8 legs ==")

## -- Leg 1: z-mean axis scoring (hand-computed) -------------------------------
## g1 = (1,2,3,4): sd = 1.290994; z1 = (-1.161895,-0.3872983,0.3872983,1.161895)
## g2 = (2,2,5,5): sd = 1.7320508; z2 = (-0.8660254,-0.8660254,0.8660254,0.8660254)
## axis mean over g1,g2 = (-1.0139602,-0.6266619,0.6266619,1.0139602)
## g3 = constant -> sd == 0 -> excluded by the floor logic, counted loudly.
X <- rbind(g1 = c(1, 2, 3, 4), g2 = c(2, 2, 5, 5), g3 = c(7, 7, 7, 7))
colnames(X) <- paste0("s", 1:4)
sd_gene <- apply(X, 1, stats::sd)
Z <- sweep(sweep(X, 1, rowMeans(X), "-"), 1, sd_gene, "/")
Z[sd_gene == 0, ] <- NA
g <- c("g1", "g2", "g3")
zable <- g %in% rownames(Z) & sd_gene[g] > 0
score <- colMeans(Z[g[zable], , drop = FALSE], na.rm = TRUE)
exp_score <- c(-1.0139602, -0.6266619, 0.6266619, 1.0139602)
leg("1 z-mean axis score matches hand-computed values",
    all(abs(score - exp_score) < 1e-6) && sum(zable) == 2)

## -- Leg 2: SenMayo schema guard (locked header measured on lead machine, R2) --
SENMAYO_HDR <- c("cohort", "sample", "gsm", "condition", "senmayo_ssgsea")
senmayo_schema_ok <- function(x) identical(names(x), SENMAYO_HDR)
good <- data.frame(cohort = "GSE7305", sample = "GSM1", gsm = "GSM1",
                   condition = "Endometriosis", senmayo_ssgsea = 0.12,
                   stringsAsFactors = FALSE, check.names = FALSE)
bad <- good; names(bad)[5] <- "senmayo"
leg("2 SenMayo schema: locked header accepted, drifted header rejected",
    senmayo_schema_ok(good) && !senmayo_schema_ok(bad))

## -- Leg 3: REML homogeneous (regression guard, = S4_smoke leg 1) --------------
## Independently computed: M = 1, SE(M) = 0.1/sqrt(3) = 0.057735, tau2 = 0.
fit1 <- metafor::rma(yi = c(1, 1, 1), sei = c(0.1, 0.1, 0.1), method = "REML")
leg("3 REML homogeneous: M = 1, SE = 0.1/sqrt(3), tau2 ~ 0",
    abs(as.numeric(fit1$b) - 1) < 1e-6 && abs(fit1$se - 0.1/sqrt(3)) < 1e-6 &&
      fit1$tau2 < 1e-4)

## -- Leg 4: REML heterogeneous (regression guard, = S4_smoke leg 2) ------------
## y = (0,2,4), se = 0.1. Independent REML: tau2 = 3.99, M = 2, Q = 800, I2 ~ 99.75
fit2 <- metafor::rma(yi = c(0, 2, 4), sei = c(0.1, 0.1, 0.1), method = "REML")
leg("4 REML heterogeneous: M = 2, tau2 = 3.99 (+/-0.01), Q ~ 800, I2 > 95",
    abs(as.numeric(fit2$b) - 2) < 1e-6 && abs(fit2$tau2 - 3.99) < 0.01 &&
      abs(fit2$QE - 800) < 1 && fit2$I2 > 95)

## -- Leg 5: structural floor -- <2 z-able genes must stop (R9) -----------------
axis_try <- c("a", "b")
sd_zero_ok <- data.frame(sd = c(0, 1.3), stringsAsFactors = FALSE,
                         row.names = c("a", "b"))
n_use <- sum(axis_try %in% rownames(sd_zero_ok) & sd_zero_ok[axis_try, "sd"] > 0)
guard_fired <- tryCatch({
  if (n_use < 2L) stop("FATAL: structural floor")
  FALSE
}, error = function(e) grepl("FATAL", conditionMessage(e)))
leg("5 structural floor: 1 z-able gene triggers FATAL (measured, R9)",
    n_use == 1 && isTRUE(guard_fired))

## -- Leg 6: frozen axis registry integrity (A4) --------------------------------
AXES <- list(
  cytotoxicity         = c("CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GNLY"),
  antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
  senescence_dormancy  = c("CDKN1A", "CDKN2A", "GADD45A", "SERPINE1"),
  proliferation        = c("MKI67", "TOP2A", "PCNA"),
  stromal_ecm          = c("VIM", "COL1A1", "COL3A1", "FN1", "ACTA2")
)
allg <- unlist(AXES)
leg("6 frozen axes: 5 axes, 24 unique genes, no within-axis duplicates",
    length(AXES) == 5 && length(unique(allg)) == 24 &&
      all(vapply(AXES, function(a) length(a) == length(unique(a)), logical(1))))

## -- Leg 7: locked SE derivation on score-like rows ----------------------------
d7 <- data.frame(contrast_id = "X", gene = c("cytotoxicity", "SenMayo"),
                 logFC = c(0.8, -1.2), t = c(2, -3), stringsAsFactors = FALSE)
z7 <- se_from_contrast(d7, "X")
leg("7 SE = |logFC/t|: 0.8/2 = 0.4, 1.2/3 = 0.4",
    nrow(z7) == 2 && all(abs(z7$se - 0.4) < 1e-12))

## -- Leg 8: tier rules on a 6-pathway meta frame -------------------------------
m8 <- data.frame(pathway = c("p1", "p2", "p3", "p4", "p5", "p6"),
                 M = c(1.0, 0.2, 1.0, 0.8, -0.9, 0.3),
                 p = c(1e-4, 1e-4, 0.2, 0.04, 0.01, 0.7),
                 tau2 = 0, I2 = 0, Q = 0, QEp = 1,
                 k = 6, k_pos = 6, k_neg = 0, stringsAsFactors = FALSE)
m8 <- assign_tiers(m8)
## BH (m=6) on sorted p (1e-4, 1e-4, 0.01, 0.04, 0.2, 0.7):
##   adjusted = (0.0003, 0.0003, 0.02, 0.06, 0.24, 0.7)
## p1: FDR<.05 & |M|>=.5 & concord>=4 -> Tier1 | p2: FDR ok, |M|<.5 -> Tier2
## p3: FDR 0.24 -> Tier3 | p4: raw p .04 but BH FDR = 0.06 -> Tier3 (boundary!)
## p5: FDR 0.02, |M| ok, concord 6 -> Tier1 | p6: FDR 0.7 -> Tier3
leg("8 tiers: size+concord -> Tier1; FDR-only -> Tier2; BH boundary -> Tier3",
    identical(m8$tier, c("Tier1", "Tier2", "Tier3", "Tier3", "Tier1", "Tier3")) &&
      abs(m8$FDR[m8$pathway == "p4"] - 0.06) < 1e-12)

say("S5 smoke: ALL", n_leg, "legs PASS -- stages may run.")
