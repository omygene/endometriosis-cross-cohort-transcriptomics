#!/usr/bin/env Rscript
# ============================================================================
# S6_01_validate.R -- External validation of the locked signature in GSE213216
#                     (Fonseca et al., Nat Genet 2023 atlas) -- v1.0
#
# Project : Benign Immunoediting -- endometriosis meta-analysis
# Stage   : S6 (confirmatory, pre-registered). Opened by written PI approval
#           2026-09-22. Design mapping locked from S6_00 audit measurements
#           ONLY (R2/R9) -- see results/S6_meta/S6_00_structure_audit.md.
#
# Charter compliance:
#   R1 zero fabrication      -- every number is computed from the locked rds
#   R2 verify actual files   -- registry md5 re-verified before any compute
#   R3 stage opened by written PI confirmation only
#   R4 this script NEVER writes/modifies any input or locked results file;
#        it reads S4_meta_Tier1_core.csv and the SenMayo GMT read-only
#   R5 locked thresholds     -- BH-FDR < 0.05 is the sole claim gate
#   R6 failures/NOT_TESTABLE reported explicitly, never hidden
#   R8 PI runs this script:  Rscript scripts/S6_01_validate.R  (from root)
#   R9 every halt is a measured condition, printed before quitting
#
# Locked mapping (measured by S6_00, not guessed):
#   object   : data/raw/GSE213216_processed/auxiliary.seurat.shared.rds
#              (373,851 cells x 30,354 genes; contains ALL 9 lineages incl.
#               the identical 23,226 endothelial cells of the dedicated
#               endothelial object -> single object, single normalization)
#   patient  : meta$Patient.No.
#   tissue   : meta$Major.Class  (harmonized -- see CLASS_MAP)
#   lineage  : meta$active.cluster (9 classes)
#   doublets : meta$DF.classifications == "Singlet" kept
#
# Locked endpoints:
#   E1 = 6 scores (5 locked axes + SenMayo) per patient x lineage pseudobulk
#   E2 = SERPINE1 logCPM in Endothelial cells (primary) + Smooth muscle
#        cells (mural compartment; GSE213216 has no perivascular annotation)
#   E3 = Tier-1 directional signature = mean z(121 up) - mean z(364 down)
#
# Locked contrasts (unit = PATIENT, never cell):
#   C1 primary  : (Peritoneal lesion + Endometrioma) vs Eutopic Endometrium
#   C2 secondary: Endometrioma vs Unaffected ovary
#   C3 secondary: Peritoneal lesion vs Disease-free peritoneum
#
# Locked rules: min 20 cells/pseudobulk; min 4 patients/arm else NOT_TESTABLE;
#   Wilcoxon rank-sum; effect = delta median; ONE BH family across all primary
#   tests; independent-patient sensitivity for C1 (patients present in both
#   arms are dropped from both); verdicts P1-P4/E3 printed by locked rules.
# ============================================================================

say  <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
fail <- function(...) { message(paste0(...)); quit(save = "no", status = 1) }
need <- function(p) if (!requireNamespace(p, quietly = TRUE))
  fail("FATAL: package '", p, "' not available -- install it and re-run (R8)")
need("Matrix"); need("SeuratObject")

ROOT <- getwd()   # PI runs from project root (R8)
RES  <- file.path(ROOT, "results", "S6_meta")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS <- 20L   # per patient x lineage pseudobulk
MIN_ARM   <- 4L    # patients per arm per test

say("== S6_01: external validation in GSE213216 (locked design v1.0) ==")

# --- 1. input registry re-verification (R2/R4) -------------------------------
reg_path <- file.path(ROOT, "scripts", "S6_input_registry.tsv")
if (!file.exists(reg_path))
  fail("FATAL: scripts/S6_input_registry.tsv missing -- run S6_00 first (R2)")
reg <- readLines(reg_path, warn = FALSE)
ln  <- reg[grepl("auxiliary.seurat.shared.rds", reg, fixed = TRUE)]
if (length(ln) != 1L)
  fail("FATAL: registry entry for auxiliary.seurat.shared.rds not unique/missing")
reg_md5 <- tolower(regmatches(ln, regexpr("[0-9a-fA-F]{32}", ln)))
if (!nchar(reg_md5)) fail("FATAL: could not parse registry md5 (registry corrupt?)")

rds <- list.files(file.path(ROOT, "data"),
                  pattern = "^auxiliary\\.seurat\\.shared\\.rds$",
                  recursive = TRUE, full.names = TRUE)
if (length(rds) != 1L)
  fail("FATAL: expected exactly 1 auxiliary.seurat.shared.rds under data/, found ",
       length(rds))
fresh_md5 <- tolower(unname(tools::md5sum(rds)))
if (fresh_md5 != reg_md5)
  fail("FATAL: input md5 CHANGED since S6_00 audit (registry ", reg_md5,
       " vs now ", fresh_md5, ") -- halt (R2/R4). Re-audit before any run.")
say("registry verified: input unchanged since audit | md5 ", fresh_md5)

# --- 2. locked references -----------------------------------------------------
AXES <- list(
  cytotoxicity         = c("CD8A","CD8B","GZMB","PRF1","NKG7","GNLY"),
  antigen_presentation = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2"),
  senescence_dormancy  = c("CDKN1A","CDKN2A","GADD45A","SERPINE1"),
  proliferation        = c("MKI67","TOP2A","PCNA"),
  stromal_ecm          = c("VIM","COL1A1","COL3A1","FN1","ACTA2")
)
AXIS_GENES <- unique(unlist(AXES))
if (length(AXIS_GENES) != 24L)
  fail("FATAL: locked axis registry corrupted (n=", length(AXIS_GENES), ", expected 24)")

t1_path <- file.path(ROOT, "results", "S4_meta", "S4_meta_Tier1_core.csv")
if (!file.exists(t1_path))
  fail("FATAL: locked Tier-1 core missing: ", t1_path)
t1 <- read.csv(t1_path, stringsAsFactors = FALSE)
if (!all(c("gene", "M") %in% names(t1)))
  fail("FATAL: Tier-1 core lacks gene/M columns -- wrong file? (R4)")
UP   <- t1$gene[t1$M > 0]
DOWN <- t1$gene[t1$M < 0]
if (length(UP) != 121L || length(DOWN) != 364L)
  fail("FATAL: Tier-1 core split is ", length(UP), "/", length(DOWN),
       " -- locked registry says 121/364. File integrity violated (R4).")
say("locked references: 24 axis genes | Tier-1 core 121 up / 364 down verified")

# SenMayo GMT -- reused verbatim, never recomputed (locked)
gmt <- list.files(file.path(ROOT, "data"),
                  pattern = "^SAUL_SEN_MAYO_M45803\\.gmt$",
                  recursive = TRUE, full.names = TRUE)
if (length(gmt) != 1L)
  fail("FATAL: expected exactly 1 SAUL_SEN_MAYO_M45803.gmt under data/, found ",
       length(gmt))
gl   <- strsplit(readLines(gmt, warn = FALSE), "\t")
SEN  <- unique(unlist(lapply(gl, function(x) x[-c(1, 2)])))
SEN  <- SEN[nzchar(SEN)]
if (length(SEN) != 124L)
  fail("FATAL: SenMayo GMT has ", length(SEN), " unique genes -- audit measured 124 (R2)")
say("SenMayo GMT reused verbatim (never recomputed): ", basename(gmt), " | 124 genes")

# --- 3. load object, lock mapping ---------------------------------------------
say("loading ", rds, " (takes minutes -- 373,851 cells) ...")
obj <- readRDS(rds)
if (!is(obj, "Seurat")) fail("FATAL: object class is ", class(obj)[1], ", expected Seurat")
meta <- obj[[]]
for (cl in c("Patient.No.", "Major.Class", "active.cluster", "DF.classifications"))
  if (!cl %in% names(meta))
    fail("FATAL: expected metadata column '", cl, "' missing -- audit mapping violated (R9)")

cts <- tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts"),
                error = function(e)
                  tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", slot = "counts"),
                           error = function(e2) NULL))
if (is.null(cts)) fail("FATAL: RNA counts layer unavailable (R9)")
rm(obj); invisible(gc(verbose = FALSE))

common <- intersect(colnames(cts), rownames(meta))
if (!length(common)) fail("FATAL: zero overlap between count columns and metadata rows (R9)")
if (length(common) != nrow(meta))
  say("NOTE: using ", length(common), " of ", nrow(meta), " metadata cells (count-matrix overlap)")
meta <- meta[common, , drop = FALSE]
cts  <- cts[, common, drop = FALSE]

# gene coverage -- exact-match self-audit against S6_00 measurements
g <- rownames(cts)
if (sum(AXIS_GENES %in% g) != 24L)
  fail("FATAL: axis coverage ", sum(AXIS_GENES %in% g), "/24 -- audit measured 24/24 (R9)")
if (sum(c(UP, DOWN) %in% g) != 463L)
  fail("FATAL: Tier-1 coverage ", sum(c(UP, DOWN) %in% g), "/485 -- audit measured 463 (R9)")
if (sum(SEN %in% g) != 124L)
  fail("FATAL: SenMayo coverage ", sum(SEN %in% g), "/124 -- audit measured 124 (R9)")
say("gene coverage self-audit: axis 24/24 | Tier-1 463/485 | SenMayo 124/124  OK")

# --- 4. cell filter + locked class harmonization ------------------------------
CLASS_MAP <- c("Endometriosis"                = "Peritoneal lesion",
               "Extra-ovarian endometriosis"  = "Peritoneal lesion",
               "Endometrioma"                 = "Endometrioma",
               "Eutopic Endometrium"          = "Eutopic Endometrium",
               "Unaffected ovary"             = "Unaffected ovary",
               "No endometriosis detected"    = "Disease-free peritoneum")
keep <- meta$DF.classifications == "Singlet" &
        !is.na(meta$Patient.No.) & !is.na(meta$active.cluster) &
        !is.na(meta$Major.Class) & meta$Major.Class %in% names(CLASS_MAP)
say("cells kept (singlets, complete mapping): ", sum(keep), " of ", nrow(meta),
    " | doublets dropped: ", sum(meta$DF.classifications == "Doublet", na.rm = TRUE))
meta <- meta[keep, , drop = FALSE]
cts  <- cts[, keep, drop = FALSE]
meta$tissue <- unname(CLASS_MAP[as.character(meta$Major.Class)])
if (any(is.na(meta$tissue))) fail("FATAL: unmapped Major.Class value survived (R9)")

# --- 5. pseudobulk per patient x lineage (unit = PATIENT) ---------------------
meta$grp <- paste(meta$Patient.No., meta$active.cluster, meta$tissue, sep = "||")
G <- Matrix::sparse.model.matrix(~ 0 + grp, data = data.frame(grp = meta$grp))
colnames(G) <- sub("^grp", "", colnames(G))
ncells <- Matrix::colSums(G)          # explicit Matrix:: (base colSums does not
say("pseudobulk groups (patient x lineage x tissue): ", ncol(G),  # S4-dispatch on
    " | dropped (<", MIN_CELLS, " cells): ", sum(ncells < MIN_CELLS))  # sparse x)
G <- G[, ncells >= MIN_CELLS, drop = FALSE]
if (!ncol(G)) fail("FATAL: zero pseudobulk groups survive the cell threshold (R9)")

pb  <- cts %*% G                       # genes x groups (raw counts)
rm(cts, G); invisible(gc(verbose = FALSE))
lib <- Matrix::colSums(pb)
if (any(lib <= 0)) fail("FATAL: zero-library pseudobulk group encountered (R9)")
lcpm <- log1p(t(t(as.matrix(pb)) / lib) * 1e6)
rm(pb); invisible(gc(verbose = FALSE))

parts   <- do.call(rbind, strsplit(colnames(lcpm), "||", fixed = TRUE))
pb_meta <- data.frame(patient = parts[, 1], lineage = parts[, 2],
                      tissue = parts[, 3], stringsAsFactors = FALSE)
say("pseudobulk matrix: ", nrow(lcpm), " genes x ", ncol(lcpm), " groups | patients: ",
    length(unique(pb_meta$patient)))

# --- 6. scores (gene-wise z within lineage -- mirrors locked S5_01 logic) -----
score_tab <- pb_meta
for (ax in names(AXES)) score_tab[[ax]] <- NA_real_
score_tab$SenMayo          <- NA_real_
score_tab$tier1_signature  <- NA_real_
score_tab$SERPINE1_logCPM  <- NA_real_

for (ln2 in unique(pb_meta$lineage)) {
  idx <- which(pb_meta$lineage == ln2)
  sub <- lcpm[, idx, drop = FALSE]
  z   <- t(scale(t(sub)))                       # gene-wise z across samples
  z[is.na(z)] <- 0                              # zero-variance genes contribute 0
  for (ax in names(AXES)) score_tab[[ax]][idx] <- colMeans(z[AXES[[ax]], , drop = FALSE])
  score_tab$SenMayo[idx]         <- colMeans(z[SEN, , drop = FALSE])
  score_tab$tier1_signature[idx] <- colMeans(z[UP[UP %in% rownames(z)], , drop = FALSE]) -
                                    colMeans(z[DOWN[DOWN %in% rownames(z)], , drop = FALSE])
  score_tab$SERPINE1_logCPM[idx] <- sub["SERPINE1", ]
}
rm(lcpm); invisible(gc(verbose = FALSE))

# --- 7. locked contrasts + tests ----------------------------------------------
CONTRASTS <- list(
  C1 = list(case = c("Peritoneal lesion", "Endometrioma"), ctrl = "Eutopic Endometrium"),
  C2 = list(case = "Endometrioma",                         ctrl = "Unaffected ovary"),
  C3 = list(case = "Peritoneal lesion",                    ctrl = "Disease-free peritoneum")
)
SCORES <- c(names(AXES), "SenMayo", "tier1_signature")
E2_LINEAGES <- c("Endothelial cells", "Smooth muscle cells")   # SERPINE1 niche (E2)

run_test <- function(df, val_col, case, ctrl) {
  x <- df[[val_col]][df$tissue %in% case]
  y <- df[[val_col]][df$tissue %in% ctrl]
  if (length(x) < MIN_ARM || length(y) < MIN_ARM)
    return(list(testable = FALSE, n_case = length(x), n_ctrl = length(y),
                p = NA_real_, delta = NA_real_, W = NA_real_))
  wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE))
  list(testable = TRUE, n_case = length(x), n_ctrl = length(y),
       p = wt$p.value, delta = median(x) - median(y), W = unname(wt$statistic))
}

tests <- list(); k <- 0L
for (cn in names(CONTRASTS)) {
  cc <- CONTRASTS[[cn]]
  for (ln2 in sort(unique(score_tab$lineage))) {
    df <- score_tab[score_tab$lineage == ln2, , drop = FALSE]
    for (sc in SCORES) {
      k <- k + 1L
      r <- run_test(df, sc, cc$case, cc$ctrl)
      # independent-patient sensitivity for C1: drop patients seen in both arms
      p_ind <- NA_real_
      if (cn == "C1" && r$testable) {
        both <- intersect(df$patient[df$tissue %in% cc$case],
                          df$patient[df$tissue %in% cc$ctrl])
        if (length(both)) {
          d2 <- df[!df$patient %in% both, , drop = FALSE]
          r2 <- run_test(d2, sc, cc$case, cc$ctrl)
          if (r2$testable) p_ind <- r2$p
        } else p_ind <- r$p
      }
      tests[[k]] <- data.frame(endpoint = ifelse(sc == "tier1_signature", "E3", "E1"),
        score = sc, lineage = ln2, contrast = cn,
        n_case = r$n_case, n_ctrl = r$n_ctrl, W = r$W, delta_median = r$delta,
        p = r$p, p_indep_patient = p_ind,
        status = ifelse(r$testable, "TESTED", "NOT_TESTABLE(<4/arm)"),
        stringsAsFactors = FALSE)
    }
    if (ln2 %in% E2_LINEAGES) {   # E2: SERPINE1 niche, logCPM scale
      k <- k + 1L
      r <- run_test(df, "SERPINE1_logCPM", cc$case, cc$ctrl)
      tests[[k]] <- data.frame(endpoint = "E2", score = "SERPINE1_logCPM",
        lineage = ln2, contrast = cn, n_case = r$n_case, n_ctrl = r$n_ctrl,
        W = r$W, delta_median = r$delta, p = r$p, p_indep_patient = NA_real_,
        status = ifelse(r$testable, "TESTED", "NOT_TESTABLE(<4/arm)"),
        stringsAsFactors = FALSE)
    }
  }
}
T <- do.call(rbind, tests)
if (!any(T$endpoint == "E1" & T$contrast == "C1" & T$status == "TESTED"))
  fail("FATAL: no testable E1/C1 test in any lineage -- stage cannot proceed (R9)")
T$FDR <- NA_real_
ok <- T$status == "TESTED"
T$FDR[ok] <- p.adjust(T$p[ok], method = "BH")     # ONE family, all primary tests
say("tests: ", nrow(T), " total | tested ", sum(ok),
    " | NOT_TESTABLE ", sum(!ok), " (documented, never hidden -- R6)")

# --- 8. verdicts by LOCKED rules (printed by the machine, not interpreted) ----
vrow <- function(score, lineage, contrast) {
  r <- T[T$score == score & T$lineage == lineage & T$contrast == contrast, ]
  if (!nrow(r)) return(NULL); r[1, ]
}
verdict <- function(r, dir) {
  if (is.null(r) || r$status != "TESTED") return("UNTESTABLE")
  if (!is.na(r$FDR) && r$FDR < 0.05) {
    if ((dir == "up"   && r$delta_median > 0) ||
        (dir == "down" && r$delta_median < 0)) return("CONFIRMED")
    return("REFUTED_DIRECTION")     # significant but opposite -- reported honestly (R1/R6)
  }
  "NOT_CONFIRMED"
}
P1 <- verdict(vrow("SERPINE1_logCPM", "Endothelial cells",  "C1"), "up")
P2 <- verdict(vrow("proliferation",   "Mesenchymal cells",  "C1"), "down")
E3 <- verdict(vrow("tier1_signature", "Mesenchymal cells",  "C1"), "up")
p3a <- vrow("senescence_dormancy", "Mesenchymal cells", "C1")
p3b <- vrow("SenMayo",             "Mesenchymal cells", "C1")
P3  <- paste0("senescence_dormancy delta=",
              ifelse(is.null(p3a), "NA", signif(p3a$delta_median, 3)),
              " (raw p=", ifelse(is.null(p3a), "NA", signif(p3a$p, 3)), ")",
              " | SenMayo delta=",
              ifelse(is.null(p3b), "NA", signif(p3b$delta_median, 3)),
              " (raw p=", ifelse(is.null(p3b), "NA", signif(p3b$p, 3)), ")",
              " -- directional-only endpoint by design")
p4a <- vrow("cytotoxicity",         "Mesenchymal cells", "C1")
p4b <- vrow("antigen_presentation", "Mesenchymal cells", "C1")
P4  <- paste0("cytotoxicity delta=",
              ifelse(is.null(p4a), "NA", signif(p4a$delta_median, 3)),
              " FDR=", ifelse(is.null(p4a), "NA", signif(p4a$FDR, 3)),
              " | antigen_presentation delta=",
              ifelse(is.null(p4b), "NA", signif(p4b$delta_median, 3)),
              " FDR=", ifelse(is.null(p4b), "NA", signif(p4b$FDR, 3)),
              " -- reported as-is, no direction claimed by design")

# --- 9. outputs (new files only; nothing locked is touched -- R4) -------------
f_scores <- file.path(RES, "S6_01_scores_long.csv")
f_tests  <- file.path(RES, "S6_01_tests.csv")
f_sum    <- file.path(RES, "S6_01_summary.md")
f_man    <- file.path(RES, "S6_01_manifest.tsv")
write.csv(score_tab, f_scores, row.names = FALSE)
write.csv(T,         f_tests,  row.names = FALSE)

pxt <- as.data.frame.matrix(table(unique(score_tab[, c("patient", "tissue")])$patient,
      unique(score_tab[, c("patient", "tissue")])$tissue))
hits <- T[ok & !is.na(T$FDR) & T$FDR < 0.05, ]
hits <- hits[order(hits$FDR), ]

md <- c(
  "# S6_01 -- GSE213216 external validation (locked design v1.0)",
  "",
  paste0("Run: ", Sys.time(), " | input md5 verified vs registry: ", fresh_md5),
  paste0("Unit = patient. Pseudobulks (>=", MIN_CELLS, " cells): ", nrow(score_tab),
         " | patients: ", length(unique(score_tab$patient))),
  paste0("Tests: ", nrow(T), " | tested ", sum(ok), " | NOT_TESTABLE ", sum(!ok),
         " | one BH family across all primary tests (FDR<0.05 sole gate)"),
  "",
  "## Patients per tissue (pseudobulk level)",
  "",
  paste0("- ", colnames(pxt), ": ", colSums(pxt), " patients"),
  "",
  "```", capture.output(print(pxt)), "```",
  "",
  "## Pre-registered verdicts (locked rules, machine-printed)",
  "",
  paste0("- **P1 (SERPINE1 up, Endothelial, C1): ", P1, "**"),
  paste0("- **P2 (proliferation down, Mesenchymal, C1): ", P2, "**"),
  paste0("- **E3 (Tier-1 signature up, Mesenchymal, C1): ", E3, "**"),
  paste0("- P3 (directional only): ", P3),
  paste0("- P4 (no direction claimed): ", P4),
  "",
  "## FDR < 0.05 hits (all, unfiltered)",
  "",
  if (nrow(hits)) paste0("- ", hits$endpoint, " | ", hits$score, " | ", hits$lineage,
                         " | ", hits$contrast, " | delta=", signif(hits$delta_median, 3),
                         " | FDR=", signif(hits$FDR, 3),
                         " | n=", hits$n_case, "v", hits$n_ctrl) else "- none",
  "",
  "Negative and NOT_TESTABLE results are reported in S6_01_tests.csv with the",
  "same weight as positives (R1/R6). C1 sensitivity with patients shared across",
  "arms removed: column p_indep_patient.",
  "",
  "S6_01 DONE -- send S6_01_summary.md + S6_01_tests.csv back to the session."
)
writeLines(md, f_sum)

man <- data.frame(file = basename(c(f_scores, f_tests, f_sum)),
                  md5  = unname(tools::md5sum(c(f_scores, f_tests, f_sum))))
write.table(man, f_man, sep = "\t", row.names = FALSE, quote = FALSE)

say("wrote ", f_scores)
say("wrote ", f_tests)
say("wrote ", f_sum)
say("wrote ", f_man)
say("verdicts -> P1: ", P1, " | P2: ", P2, " | E3: ", E3)
say("S6_01 DONE -- send S6_01_summary.md + S6_01_tests.csv back to the session")
