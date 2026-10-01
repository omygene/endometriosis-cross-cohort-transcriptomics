# ============================================================================
# S3_02_ssgsea_bulk.R -- SenMayo ssGSEA scores for all 7 bulk cohorts
# Locked criterion (3): ssGSEA for bulk. Vehicle = GSVA (Bioconductor), the
# citable standard; API verified against GSVA source (ssgseaParam + gsva for
# GSVA >= 1.50; legacy gsva(..., method="ssgsea") branch kept as fallback).
# Stops FATAL if GSVA or a required annotation .db package is absent --
# adding packages needs your written approval (R5); do NOT substitute a
# hand-rolled scorer without an Amendment.
# Probe->symbol collapse (frozen): per symbol keep the probe set with the
# highest mean expression. Output: scores per array (long table).
# v1.1 (2026-09-20, measured halt at GSE120103): the neqc object carries no
#   symbol column (read.maimages default keeps Row/Col/ControlType/ProbeName/
#   SystematicName only), and the raw FE txt itself has NO GeneName column
#   (measured on GSM3393491 from GEO: 26 columns, FEATURES header). Symbols are
#   recovered from genes$SystematicName via org.Hs.eg.db -- a hard dependency
#   of the hgu133*.db packages already used, so NO new package (R5-safe).
#   Measured prefix distribution (43,376 regular probes): NM_ 29,028 / other
#   GenBank+THC+TCONS 5,556 / bare probe-id 4,171 / ENST 2,298 / NR_ 1,850 /
#   XR_ 238 / XM_ 177 / NP_ 58; 0 conflicting accessions per ProbeName.
#   Pre-measured coverage via the GEO GPL6480 annotation (same accession
#   basis): 76.2% of probes carry a symbol -> 19,565 unique symbols, and
#   123/124 SenMayo genes are covered (only CST4 absent from the platform).
#   v1.2 (2026-09-20): drop the over-strict ProbeName/rownames identity gate --
#   the EList contract (E rows == genes rows) is the actual alignment proof.
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_02_ssgsea_bulk.R > logs\S3_02_ssgsea_bulk.txt 2>&1
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
dir.create(CFG$out_scores, showWarnings = FALSE)

say("== S3_02 ssGSEA (SenMayo) ==")
if (!requireNamespace("GSVA", quietly = TRUE))
  stop("FATAL: GSVA not installed. ssGSEA (locked criterion 3) requires GSVA; ",
       "written approval needed to add it (R5), then renv::snapshot().")

gm <- find_senmayo_gmt()
if (is.na(gm)) stop("FATAL: SenMayo GMT not found under reference/ or data/ (see S3_00)")
sets <- read_gmt(gm)
sen <- unique(unlist(lapply(sets, `[[`, "genes")))
say("GMT:", gm, "| sets:", length(sets), "| unique symbols:", length(sen))

DB <- c(GSE7305 = "hgu133plus2.db", GSE6364 = "hgu133plus2.db",
        GSE25628 = "hgu133a2.db", GSE7307 = "hgu133plus2.db",
        GSE51981 = "hgu133plus2.db", GSE11691 = "hgu133a.db")

## SystematicName -> SYMBOL for the Agilent cohort via org.Hs.eg.db
## (hard dependency of the hgu133*.db packages already used -- no new package,
## R5-safe). Keytype split measured on GSM3393491 (2026-09-20): ENST* ->
## ENSEMBLTRANS; all other accessions -> ACCNUM (covers RefSeq + GenBank);
## bare probe ids (A_..) and THC/TCONS accessions are unannotated by design
## and become NA (dropped by the frozen max-mean collapse).
agilent_symbols <- function(obj, mat) {
  g <- obj$genes
  if (!all(c("ProbeName", "SystematicName") %in% colnames(g)))
    stop("FATAL: ", CFG$agilent, " genes lacks ProbeName/SystematicName; available: ",
         paste(colnames(g), collapse = ", "), " (R9)")
  ## EList class contract (limma documentation): E and genes are row-aligned
  ## BY CONSTRUCTION, and neqc returns an EList -- so genes$SystematicName is
  ## row-aligned with the matrix regardless of what rownames(mat) holds.
  ## (v1.1 bug, measured on the lead machine 2026-09-20: an over-strict
  ## identical() gate on ProbeName vs rownames false-stopped the run.)
  if (nrow(g) != nrow(mat))
    stop("FATAL: ", CFG$agilent, " genes rows ", nrow(g), " != matrix rows ",
         nrow(mat), " -- EList invariant broken (R9)")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("FATAL: org.Hs.eg.db required for ", CFG$agilent,
         " accession->symbol mapping (R5 approval + renv::snapshot)")
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  obj_db <- get("org.Hs.eg.db", envir = as.environment("package:org.Hs.eg.db"))
  sys <- as.character(g$SystematicName)
  is_enst <- startsWith(sys, "ENST")
  sym <- rep(NA_character_, length(sys))
  if (any(!is_enst))
    sym[!is_enst] <- unname(AnnotationDbi::mapIds(
      obj_db, keys = sys[!is_enst], column = "SYMBOL", keytype = "ACCNUM",
      multiVals = "first"))
  if (any(is_enst))
    sym[is_enst] <- unname(AnnotationDbi::mapIds(
      obj_db, keys = sys[is_enst], column = "SYMBOL", keytype = "ENSEMBLTRANS",
      multiVals = "first"))
  say("  ", CFG$agilent, " accession->symbol:", sum(!is.na(sym)), "of",
      length(sym), "probes mapped (bare probe ids / THC / TCONS unannotated by design)")
  say("  ", CFG$agilent, " trace -- first probe:", g$ProbeName[1], "->",
      g$SystematicName[1], "->", sym[1])
  sym
}

symbols_for <- function(cohort, obj, mat) {
  if (cohort == CFG$agilent) {
    g <- obj$genes
    cand <- intersect(c("GeneName", "Symbol", "GENE_NAME", "gene_symbol"), colnames(g))
    if (length(cand)) {
      sym <- g[[cand[1]]]
      if (is.null(sym) || all(is.na(sym))) stop("FATAL: ", cohort, " symbol column '",
                                                cand[1], "' empty")
      sym
    } else {
      ## measured 2026-09-20: genes carries only Row/Col/ControlType/ProbeName/
      ## SystematicName -- map SystematicName -> SYMBOL via org.Hs.eg.db.
      agilent_symbols(obj, mat)
    }
  } else {
    db <- unname(DB[cohort])
    if (is.null(db) || !nzchar(db))
      stop("FATAL: no .db package registered for ", cohort)
    if (!requireNamespace(db, quietly = TRUE))
      stop("FATAL: ", db, " required for ", cohort, " probe->symbol mapping; ",
           "written approval needed to add it (R5)")
    suppressPackageStartupMessages(library(db, character.only = TRUE))
    obj_db <- get(db, envir = as.environment(paste0("package:", db)))
    AnnotationDbi::mapIds(obj_db, keys = rownames(mat), column = "SYMBOL",
                          keytype = "PROBEID", multiVals = "first")
  }
}

score_one <- function(Xs) {
  gv_ver <- utils::packageVersion("GSVA")
  if (gv_ver >= "1.50.0") {
    param <- GSVA::ssgseaParam(as.matrix(Xs), list(SenMayo = sen))
    as.numeric(GSVA::gsva(param)[1, ])
  } else {
    as.numeric(GSVA::gsva(as.matrix(Xs), list(SenMayo = sen), method = "ssgsea")[1, ])
  }
}

meta <- load_meta()
all_scores <- list()
for (ch in names(CFG$bulk_rds)) {
  obj <- readRDS(unname(CFG$bulk_rds[ch]))
  mat <- if (ch == CFG$agilent) obj$E else Biobase::exprs(obj)
  sym <- symbols_for(ch, obj, mat)
  Xs <- collapse_maxmean(mat, sym)
  say(ch, ": probe matrix", nrow(mat), "-> symbol matrix", nrow(Xs))
  ## measured guard (v1.1): ssGSEA is meaningless if the gene set vanished
  ## from the symbol matrix. 10 is a structural floor, not a scientific
  ## threshold -- the Affy cohorts map ~21k symbols and GSE120103 is
  ## pre-measured at 123/124 SenMayo genes via the GEO platform annotation.
  hit <- sum(rownames(Xs) %in% sen)
  say("  ", ch, "SenMayo coverage:", hit, "of", length(sen))
  if (hit < 10L)
    stop("FATAL: SenMayo coverage collapsed for ", ch, " (R9)")
  sc <- score_one(Xs)
  m <- map_cols(ch, colnames(mat), meta)
  all_scores[[ch]] <- data.frame(cohort = ch, sample = colnames(mat),
                                 gsm = m$gsm, condition = m$condition,
                                 senmayo_ssgsea = sc, stringsAsFactors = FALSE)
}
out <- do.call(rbind, all_scores)
write.csv(out, file.path(CFG$out_scores, "S3_ssgsea_senmayo_bulk.csv"), row.names = FALSE)
say("wrote results/S3_ssgsea_senmayo_bulk.csv --", nrow(out), "arrays x", ncol(out), "cols")
say("S3_02 DONE")
