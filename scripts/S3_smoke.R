# ============================================================================
# S3_smoke.R -- mandatory smoke test before ANY S3 full run (rule B.4)
# Synthetic legs exercise the same code paths as the full stage:
#   1 frozen-TSV parse + registry sanity
#   2 limma unpaired contrast (injected signal must be recovered)
#   3 limma paired contrast via duplicateCorrelation (injected signal)
#   4 GSE11691 pairing self-audit (integrity, no biology)
#   5 probe->symbol collapse (max-mean rule)
#   6 ssGSEA via GSVA on synthetic data (skipped loudly if GSVA absent)
#   7 z-mean computation equals manual reference value
#   8 pseudobulk aggregation equals manual colSums reference
#   9 get_rna() extraction path on a real SeuratObject (layer API; the
#     slot->layer deprecation halt 2026-09-20 was invisible to legs 1-8)
#  10 prep_rna_norm() fills an empty RNA data layer on a real SeuratObject
#     (S3_04 halt 2026-09-21: "Layer 'data' is empty")
#  11 join_annotations() round-trip on synthetic atlas+object, incl. its
#     fail-loud behavior on unmatched and duplicated keys (S3_00f key:
#     <orig.ident>__<barcode>, DOUBLE underscore)
#  12 pooled z-mean NA convention: genes absent from one dataset's block get
#     z = 0 after pooling, matching the frozen S3_05 definition
# Any FAIL stops the stage (exit != 0 through stop()).
# v2.8 (2026-09-20): + GSE51981 frozen-encoding guards (5 non-pooled rows must
#   exclude the 37 pelvic-pathology controls; B-G51981-S3 must pool controls).
# v2.8.1 (2026-09-20): guard bugfix -- empty frozen cells are the literal "NA"
#   marker; v2.8's nzchar() check false-FAILED B-G51981-S3 (measured on the
#   lead machine). Both new guards now treat "NA" as empty (loader convention).
# v2.9 (2026-09-20): + GSE120103 guard (all five rows must carry explicit
#   exclusions -- measured run halt at B-G120103 on Infertile control).
# v2.10 (2026-09-20): + leg 9 get_rna() on a real SeuratObject (slot->layer
#   defunct halt at S3_03 was invisible to the synthetic legs).
# v2.11 (2026-09-21): + leg 10 prep_rna_norm() (S3_04 halt: empty RNA data
#   layer; fix point shared by S3_04/S3_05 via S3_common).
# v2.12 (2026-09-21): + legs 11-12 (per-dataset source switch, S3_common v3.0:
#   join_annotations gate + pooled z-mean NA convention). The S2_02 atlas was
#   measured (S3_00d/e/f) to be a 28-gene panel; all scRNA stages now join
#   full expression from per-dataset objects by <orig.ident>__<barcode>.
# v2.12.1 (2026-09-21): leg 11 duplicated-key sub-test now injects the
# collision via join_annotations' orig/bc overrides (S3_common v3.0.1) --
# Seurat forbids duplicate colnames at rename time, so the old test halted
# inside colnames<- before the gate was ever exercised (measured on the
# lead machine).
# Run from project root:
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_smoke.R > logs\S3_smoke.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }
smoke_fail <- FALSE
fail <- function(msg) { say("  FAIL:", msg); smoke_fail <<- TRUE }

suppressPackageStartupMessages({ library(limma) })
source("scripts/S3_common.R")

say("== SMOKE 1/12: frozen registry parse ==")
fr <- load_frozen()
if (nrow(fr) != 19L) fail(paste("expected 19 frozen rows, got", nrow(fr)))
if (!all(c("B-G7305", "B-G51981-S2", "B-G120103", "B-G11691", "B-G25628-ECT") %in% fr$contrast_id))
  fail("key contrast ids missing")
## loader regression guard (measured miss 2026-09-19): read.delim's default
## na.strings converts the frozen TSV's literal "NA" covariate cells to real
## NAs, which crashes run_frozen_contrast on the first no-covariate contrast.
## Smoke must fail loudly if that ever regresses.
if (any(vapply(fr, function(z) any(is.na(z)), logical(1))))
  fail("frozen TSV parsed with NA values -- loader must keep literal 'NA' strings")
if (!identical(fr$covariate_col[fr$contrast_id == "B-G7305"], "NA"))
  fail("covariate_col 'NA' not preserved as character string")
## frozen-encoding guard (measured miss 2026-09-20): GSE25628 exploratory
## rows compare subsets of the cohort, so the untouched third group must be
## excluded explicitly; empty exclude_samples would FATAL the real run.
ex <- fr[fr$contrast_id %in% c("B-G25628-EXPL", "B-G25628-ECT", "B-G25628-EUT"), ]
if (nrow(ex) != 3L || any(!nzchar(ex$exclude_samples)))
  fail("GSE25628 exploratory rows must carry explicit exclude_samples")
## frozen-encoding guard (measured miss 2026-09-20, run halt at B-G51981):
## GSE51981 carries a third control group (Non-EMS pelvic pathology, 37 GSMs)
## that the no-pathology-ref rows must exclude explicitly; the pooled-control
## sensitivity row (B-G51981-S3) must instead KEEP both control labels in ref
## and carry NO exclusions.
ex5 <- fr[fr$contrast_id %in% c("B-G51981", "B-G51981-MILD", "B-G51981-SEV",
                                "B-G51981-S1", "B-G51981-S2"), ]
ex5_real <- nzchar(ex5$exclude_samples) & ex5$exclude_samples != "NA"
if (nrow(ex5) != 5L || any(!ex5_real))
  fail("GSE51981 non-pooled rows must exclude the 37 pelvic-pathology controls")
s3 <- fr[fr$contrast_id == "B-G51981-S3", ]
## empty cells in the frozen TSV are the LITERAL "NA" marker, not "" -- the
## v2.8 guard used nzchar() and false-FAILED on that marker (measured on the
## lead machine 2026-09-20); v2.8.1 treats "NA" as empty (loader convention).
s3_no_excl <- !nzchar(s3$exclude_samples) | s3$exclude_samples == "NA"
if (nrow(s3) != 1L || !grepl("Non-EMS pelvic pathology", s3$group_ref) ||
    !s3_no_excl)
  fail("B-G51981-S3 must pool both control labels and carry no exclusions")
## frozen-encoding guard (measured miss 2026-09-20, run halt at B-G120103):
## GSE120103 carries four groups (2A/2B/Fertile/Infertile); every frozen row
## compares only two, so all five rows must carry explicit exclusions.
ex12 <- fr[fr$cohort == "GSE120103", ]
ex12_real <- nzchar(ex12$exclude_samples) & ex12$exclude_samples != "NA"
if (nrow(ex12) != 5L || any(!ex12_real))
  fail("all five GSE120103 rows must carry explicit exclude_samples")
## group-map guard (lead decision 2026-09-20): labels are semicolon lists;
## each dataset's first test/ref label must round-trip through the mapper.
gm <- load_group_map()
for (i in seq_len(nrow(gm))) {
  lt <- split_semicolon(gm$label_test[i])[1]
  lr <- split_semicolon(gm$label_ref[i])[1]
  tg <- map_cell_group(gm$dataset[i], c(lt, lr))
  if (!identical(tg, c("TEST", "REF")))
    fail(paste("group map round-trip failed for", gm$dataset[i]))
}
say("  frozen rows:", nrow(fr), "| covariate_col NA-marker preserved: TRUE |",
    "G25628 exclusions + group map round-trip: TRUE")

say("== SMOKE 2/12: unpaired limma contrast ==")
## Power derivation (documented calibration, v2 2026-09-18): 6v6, delta=2,
## sd=1 -> noncentrality = 2/sqrt(1/6+1/6) = 3.46, df=10. Plain-t Monte
## Carlo (2000 reps): median FDR<0.05 recovery = 3/50 -- BH over m=2000 is
## power-limited BY CONSTRUCTION at these parameters; the v1 criterion
## (>=45/50) was statistically unrealistic. limma eBayes shares variance
## across genes (all sd=1 here) and recovers far more; measured on the lead
## machine (seed 42): 28/50. v2 asserts what the leg can prove: direction,
## nominal power, an FDR floor far below every simulated quantile, and
## purity of the FDR-significant set.
ng <- 2000; ns <- 12
X <- matrix(rnorm(ng * ns), ng, ns)
rownames(X) <- paste0("g", 1:ng)
up <- 1:50; X[up, 7:12] <- X[up, 7:12] + 2
grp <- factor(rep(c("REF", "TEST"), each = 6), levels = c("REF", "TEST"))
design <- model.matrix(~ grp)
fit <- eBayes(lmFit(X, design))
tt <- topTable(fit, coef = "grpTEST", number = Inf, sort.by = "none")
inj <- rownames(tt) %in% paste0("g", up)
fdr_set <- tt$adj.P.Val < 0.05
nom_hit <- sum(tt$P.Value[inj] < 0.05)
fdr_hit <- sum(fdr_set[inj])
purity  <- if (sum(fdr_set)) sum(fdr_set & inj) / sum(fdr_set) else 0
dir_ok  <- mean(tt$logFC[inj]) > 1
say("  nominal<0.05:", nom_hit, "/50 | FDR<0.05:", fdr_hit, "/50 | purity:",
    round(purity, 3), "| mean logFC injected:", round(mean(tt$logFC[inj]), 2))
if (!(dir_ok && nom_hit >= 40 && fdr_hit >= 15 && purity >= 0.8)) {
  fail("unpaired assertions not met (direction, nominal>=40, FDR>=15, purity>=0.8)")
} else {
  say("  unpaired OK")
}

say("== SMOKE 3/12: paired limma contrast ==")
## v1 used 6 pairs x delta=1.5 -> paired noncentrality = 1.5/(sqrt(2)/sqrt(6))
## = 2.60, df=5; plain-t Monte Carlo median FDR<0.05 recovery = 0/50 --
## powerless BY CONSTRUCTION (the 4 recovered by eBayes moderation alone
## confirm the pipeline works but cannot calibrate it). v2 raises the
## SYNTHETIC injection to a testable level: 8 pairs x delta=3 ->
## noncentrality = 3/(sqrt(2)/sqrt(8)) = 6.0, df=7; plain-t MC median FDR
## recovery = 33/50 (P5 = 24). This calibrates the TEST BENCH only; no real
## data or registered parameter is affected.
np <- 8; Xp <- matrix(rnorm(ng * np * 2), ng, np * 2)
rownames(Xp) <- paste0("g", 1:ng)
Xp[up, (np + 1):(np * 2)] <- Xp[up, (np + 1):(np * 2)] + 3
woman <- rep(paste0("W", 1:np), 2)
grp2 <- factor(rep(c("REF", "TEST"), each = np), levels = c("REF", "TEST"))
design2 <- model.matrix(~ grp2)
corfit <- duplicateCorrelation(Xp, design2, block = woman)
say("  synthetic paired consensus cor =", round(corfit$consensus, 3))
fit2 <- eBayes(lmFit(Xp, design2, block = woman, correlation = corfit$consensus))
tt2 <- topTable(fit2, coef = "grp2TEST", number = Inf, sort.by = "none")
inj2 <- rownames(tt2) %in% paste0("g", up)
fdr_set2 <- tt2$adj.P.Val < 0.05
nom_hit2 <- sum(tt2$P.Value[inj2] < 0.05)
fdr_hit2 <- sum(fdr_set2[inj2])
purity2  <- if (sum(fdr_set2)) sum(fdr_set2 & inj2) / sum(fdr_set2) else 0
dir_ok2  <- mean(tt2$logFC[inj2]) > 2
say("  nominal<0.05:", nom_hit2, "/50 | FDR<0.05:", fdr_hit2, "/50 | purity:",
    round(purity2, 3), "| mean logFC injected:", round(mean(tt2$logFC[inj2]), 2))
if (!(dir_ok2 && nom_hit2 >= 45 && fdr_hit2 >= 20 && purity2 >= 0.8)) {
  fail("paired assertions not met (direction, nominal>=45, FDR>=20, purity>=0.8)")
} else {
  say("  paired OK")
}

say("== SMOKE 4/12: GSE11691 pairing audit ==")
pr <- read.delim(CFG$pairs11691, comment.char = "#", stringsAsFactors = FALSE)
allg <- c(pr$eutopic_gsm, pr$ectopic_gsm)
if (!(nrow(pr) == 9 && length(unique(allg)) == 18)) fail("pairing integrity violated") else say("  9 women / 18 distinct GSMs OK")

say("== SMOKE 5/12: max-mean symbol collapse (shipped function) ==")
m2 <- matrix(c(1, 5, 2, 2, 1, 1), 3, 2)
rownames(m2) <- c("p1", "p2", "p3")
syms <- c("A", "A", "B")
col2 <- collapse_maxmean(m2, syms)
if (!(identical(sort(rownames(col2)), c("A", "B")) &&
      identical(col2["A", ], m2["p2", ]))) fail("collapse wrong") else say("  collapse keeps p2 as A, p3 as B, rownames=symbols OK")

say("== SMOKE 6/12: ssGSEA via GSVA (synthetic) ==")
if (requireNamespace("GSVA", quietly = TRUE)) {
  Xs <- matrix(rnorm(100 * 8), 100, 8)
  Xs[1:10, 5:8] <- Xs[1:10, 5:8] + 3
  rownames(Xs) <- paste0("g", 1:100)
  gv_ver <- utils::packageVersion("GSVA")
  sc <- if (gv_ver >= "1.50.0") {
    as.numeric(GSVA::gsva(GSVA::ssgseaParam(Xs, list(S = paste0("g", 1:10))))[1, ])
  } else {
    as.numeric(GSVA::gsva(Xs, list(S = paste0("g", 1:10)), method = "ssgsea")[1, ])
  }
  if (!(mean(sc[5:8]) > mean(sc[1:4]))) fail("ssGSEA did not rank enriched group higher") else say("  ssGSEA direction OK (GSVA", as.character(gv_ver), ")")
} else {
  say("  SKIP-LOUD: GSVA absent -- S3_02 cannot run until approved+installed (R5)")
}

say("== SMOKE 7/12: z-mean equals manual reference ==")
set.seed(1); M <- matrix(rpois(200, 5) + 1, 20, 10)
z1 <- scale(t(log1p(M)))   # cells x genes, gene-wise z across cells
ref <- vapply(1:5, function(g) (log1p(M)[g, 3] - mean(log1p(M)[g, ])) / sd(log1p(M)[g, ]), numeric(1))
if (!isTRUE(all.equal(as.numeric(z1[3, 1:5]), ref))) fail("z-mean mismatch") else say("  z-mean OK")

say("== SMOKE 8/12: pseudobulk aggregation ==")
set.seed(2); C0 <- matrix(rpois(50 * 6, 10), 50, 6)
agg <- rowSums(C0[, 1:3])
if (!isTRUE(all.equal(agg, Matrix::rowSums(Matrix::Matrix(C0)[, 1:3])))) fail("aggregation mismatch") else say("  aggregation OK")

say("== SMOKE 9/12: get_rna extraction path ==")
## Measured miss 2026-09-20: SeuratObject >= 5.1 makes the `slot` argument
## defunct; only a REAL SeuratObject exposes get_rna() API drift.
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  fail("SeuratObject not installed -- scRNA legs cannot run (R5)")
} else {
  C9 <- C0
  rownames(C9) <- paste0("g", 1:50); colnames(C9) <- paste0("c", 1:6)
  so9 <- suppressWarnings(SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(C9, sparse = TRUE)))
  gr <- get_rna(so9, "counts")
  if (!isTRUE(all.equal(as.matrix(gr), C9, check.attributes = FALSE)) ||
      !identical(rownames(gr), rownames(C9))) {
    fail("get_rna extraction mismatch")
  } else say("  get_rna counts extraction OK (SeuratObject",
             as.character(utils::packageVersion("SeuratObject")), ")")
}

say("== SMOKE 10/12: prep_rna_norm fills empty data layer ==")
## Measured miss 2026-09-21 (S3_04 halt): RNA assay stores counts only;
## FindMarkers died on "Layer 'data' is empty". Mirror that object shape and
## assert LogNormalize-equivalence + idempotence.
if (!requireNamespace("Seurat", quietly = TRUE)) {
  fail("Seurat not installed -- prep_rna_norm leg cannot run (R5)")
} else {
  C10 <- C0
  rownames(C10) <- paste0("g", 1:50); colnames(C10) <- paste0("c", 1:6)
  so10 <- suppressWarnings(SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(C10, sparse = TRUE)))
  so10 <- prep_rna_norm(so10)
  d10 <- SeuratObject::GetAssayData(so10, assay = "RNA", layer = "data")
  ref10 <- log1p(sweep(as.matrix(C10), 2, colSums(C10) / 1e4, "/"))
  if (nrow(d10) != nrow(C10) ||
      !isTRUE(all.equal(as.matrix(d10), ref10, tolerance = 1e-6,
                        check.attributes = FALSE))) {
    fail("prep_rna_norm did not produce LogNormalize-equivalent log1p data")
  } else {
    d10b <- as.matrix(SeuratObject::GetAssayData(prep_rna_norm(so10),
                                                 assay = "RNA", layer = "data"))
    if (!isTRUE(all.equal(d10b, as.matrix(d10), check.attributes = FALSE)))
      fail("prep_rna_norm not idempotent -- second call changed data")
    else say("  prep_rna_norm OK -- empty data layer filled (log1p CP10K), idempotent")
  }
}

say("== SMOKE 11/12: join_annotations round-trip + fail-loud gates ==")
## S3_00f measured key: atlas rowname == paste0(orig.ident, "__", barcode)
## (DOUBLE underscore). The join gate must match annotations in ANY cell
## order and must halt on unmatched or duplicated keys.
{
  fk_md <- data.frame(
    dataset   = rep("D1", 4),
    condition = c("EMS", "N", "EMS", "N"),
    lineage   = c("Stromal", "Stromal", "NKT", "NKT"),
    sample    = c("S1", "S2", "S1", "S2"),
    orig.ident= c("S1", "S2", "S1", "S2"),
    stringsAsFactors = FALSE)
  rownames(fk_md) <- c("S1__bc1", "S2__bc2", "S1__bc3", "S2__bc4")
  C11 <- matrix(1:12, 3, 4)
  colnames(C11) <- c("bc3", "bc1", "bc4", "bc2")   # shuffled vs md order
  rownames(C11) <- paste0("g", 1:3)
  so11 <- suppressWarnings(SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(C11, sparse = TRUE)))
  so11$orig.ident <- c("S1", "S1", "S2", "S2")
  ann11 <- join_annotations(so11, fk_md, "D1")
  if (!identical(as.character(ann11$condition[ann11$cell_key == "S1__bc1"]), "EMS") ||
      !identical(as.character(ann11$condition[ann11$cell_key == "S2__bc4"]), "N") ||
      !identical(as.character(ann11$lineage[ann11$cell_key == "S1__bc3"]), "NKT"))
    fail("join_annotations matched wrong rows (order or key bug)")
  ## fail-loud: one barcode with no atlas annotation must halt
  so11b <- so11
  colnames(so11b)[1] <- "bcX"
  so11b$orig.ident[1] <- "S9"
  e11 <- tryCatch({ join_annotations(so11b, fk_md, "D1"); NULL },
                  error = function(e) conditionMessage(e))
  if (is.null(e11) || !grepl("FATAL", e11))
    fail("join_annotations did not halt on unmatched key")
  ## fail-loud: duplicated (orig.ident, barcode) must halt. NOTE: we inject
  ## the collision through the orig/bc overrides -- renaming Seurat cells to
  ## a duplicate colnames is forbidden by Seurat itself before the gate
  ## would ever see it (measured 2026-09-21, smoke v2.12 halt).
  e11b <- tryCatch({
    join_annotations(so11, fk_md, "D1",
                     orig = c("S1", "S1", "S2", "S2"),
                     bc   = c("bc3", "bc3", "bc4", "bc2"))
    NULL }, error = function(e) conditionMessage(e))
  if (is.null(e11b) || !grepl("FATAL", e11b))
    fail("join_annotations did not halt on duplicated key")
  if (!smoke_fail) say("  join_annotations OK -- round-trip + both fail-loud gates")
}

say("== SMOKE 12/12: pooled z-mean NA convention ==")
## Frozen S3_05 def: gene-wise z across ALL pooled cells; genes absent from
## one dataset's block contribute NA -> z = 0 after pooling. Build two
## blocks sharing 2 of 3 genes and check the pooled z against a manual
## reference computed with that exact convention.
{
  M1 <- matrix(1:12, nrow = 3)                       # gA: 1,4,7,10 | 4 cells
  M2 <- matrix(c(2, 4, 6, 8, 10, 12), nrow = 2)      # gA: 2,6,10   | 3 cells
  rownames(M1) <- c("gA", "gB", "gC"); rownames(M2) <- c("gA", "gB")
  gn <- c("gA", "gB", "gC")
  P12 <- matrix(NA_real_, 3, 7, dimnames = list(gn, NULL))
  P12[rownames(M1), 1:4] <- M1
  P12[rownames(M2), 5:7] <- M2
  z12 <- scale(t(P12)); z12[is.na(z12)] <- 0
  zm12 <- rowMeans(z12)
  ## manual reference for cell 1 (all genes present) and cell 5 (gC absent -> 0)
  ref5 <- c((2 - mean(P12[1, ], na.rm = TRUE)) / sd(P12[1, ], na.rm = TRUE),
            (4 - mean(P12[2, ], na.rm = TRUE)) / sd(P12[2, ], na.rm = TRUE), 0)
  if (!isTRUE(all.equal(as.numeric(z12[5, ]), ref5, tolerance = 1e-8)))
    fail("pooled z NA convention mismatch (frozen definition changed?)")
  else say("  pooled z-mean OK -- NA->0 convention matches frozen definition")
}

if (smoke_fail) stop("SMOKE FAILED -- do not run S3 (rule B.4)") else say("SMOKE PASS -- S3 may run after lead approval (R3)")
