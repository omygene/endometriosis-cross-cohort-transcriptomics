# ============================================================================
# S3_common.R -- shared helpers + FROZEN config for stage S3 (primary contrasts)
# Governed by locked criteria approved 2026-09-18:
#   (1) gate BH-FDR < 0.05; |log2FC| >= 0.5 interpretive tier, NO exclusion;
#       sign concordance reported, never a filter.
#   (2) scRNA: pseudobulk per sample PRIMARY; Wilcoxon per lineage SENSITIVITY.
#   (3) scores: ssGSEA (bulk) + z-mean (scRNA).
#   (7) powercfg unchanged (no run here exceeds 2 h).
# Sourced by all S3_* scripts; every script runs with cwd = project root.
# v2.8 (2026-09-20): get_rna() uses layer= (slot defunct in SeuratObject 5.1,
# measured halt at S3_03); smoke leg 9 now covers this path.
# v2.9 (2026-09-21): + prep_rna_norm() (measured halt at S3_04: S2_02 atlas
# keeps counts only, "Layer 'data' is empty"); smoke leg 10 covers it.
# v3.0.2 (2026-09-22): + load_bulk_symbol_matrix() -- probe->symbol collapse for
# bulk DE (measured failure ledger #18: S3_01 fitted probe-level DE and never
# collapsed, so S4's cross-cohort meta silently lost GSE120103 (Agilent probe
# namespace disjoint from the Affy cohorts' shared namespace; LOCO_B-G120103
# identical to base meta, dM = 0.0, measured from uploaded S4 outputs). The
# mapping code below is ported VERBATIM from S3_02_ssgsea_bulk.R (proven path
# on this machine, renv-locked annotation packages) -- single frozen rule.
# v3.0.1 (2026-09-21): join_annotations() gains optional orig/bc overrides
# (smoke leg 11 must test the duplicated-key gate without renaming Seurat
# cells -- Seurat forbids duplicate colnames first; measured in smoke).
# v3.0 (2026-09-21): + per-dataset expression sources (measured by S3_00d/
# S3_00e/S3_00f: S2_02_integrated.rds is a 28-gene LINEAGE PANEL object, not
# a full-expression atlas; full expression lives in the per-dataset objects,
# cell-for-cell identical 1:1, joined by cell_key = <orig.ident>__<barcode>).
# New: CFG$per_dataset, atlas_annotations(), join_annotations(). The frozen
# criteria are unchanged -- same tests, same definitions, corrected source.
# ============================================================================

CFG <- list(
  seed = 42,
  meta       = "data/metadata/S1_bulk_conditions.tsv",
  frozen     = "scripts/S3_contrasts_frozen.tsv",
  pairs11691 = "scripts/S3_GSE11691_pairs.tsv",
  fdr_gate   = 0.05,
  lfc_tier   = 0.5,
  min_cells_pseudobulk = 10L,
  bulk_rds = c(
    GSE7305   = "data/processed/bulk_GSE7305_rma.rds",
    GSE6364   = "data/processed/bulk_GSE6364_rma.rds",
    GSE25628  = "data/processed/bulk_GSE25628_rma.rds",
    GSE7307   = "data/processed/bulk_GSE7307_rma.rds",
    GSE51981  = "data/processed/bulk_GSE51981_rma.rds",
    GSE11691  = "data/processed/bulk_GSE11691_rma.rds",
    GSE120103 = "data/processed/bulk_GSE120103_neqc.rds"
  ),
  agilent = "GSE120103",
  atlas   = "data/processed/S2_02_integrated.rds",  # ANNOTATIONS ONLY (S3_00d:
  # 28-gene panel, no expression matrix). Never use its assay for expression.
  per_dataset = c(   # full-expression sources (S3_00e: SenMayo 124/120/123)
    GSE179640 = "data/processed/GSE179640_filtered.rds",
    GSE183837 = "data/processed/GSE183837_filtered.rds",
    GSE214411 = "data/processed/GSE214411_filtered.rds"
  ),
  out_contrasts = "results/S3_bulk_contrasts",
  out_scores    = "results",
  ## --- TO CONFIRM AFTER S3_00 CENSUS -------------------------------------
  ## Fill each with exactly one value reported by S3_00 / S3_00b, then re-run
  ## the dependent scripts. Empty values stop them FATAL (fail-loud, no guess).
  ## group map is PER-DATASET (scripts/S3_group_map.tsv): the atlas comparisons
  ## are not all EMS-vs-N (GSE183837 = RIF vs Ctrl per the approved master
  ## plan), so a single global label pair is structurally wrong.
  n1_sample    = "GSM6605437_N1",   # sample id of frozen N1 borderline group (GSE214411)
  group_col    = "condition",   # atlas metadata column carrying group labels
  ## hormonal path (MEASURED 2026-09-19): absent from atlas AND from the
  ## per-dataset object (S3_00b: 5 QC columns only). Source = official GEO
  ## series matrix, extracted by S3_00c into hormonal_src (one row per GSM).
  ## Join key (v2.4, measured): the GSM accession = gsm_of(orig.ident).
  ## GEO Sample_title does NOT match orig.ident (audit A v1: 0 of 33).
  hormonal_obj      = "data/processed/GSE179640_filtered.rds",  # census audit
  hormonal_src      = "scripts/S3_gse179640_geo_characteristics.tsv",
  hormonal_join_col = "orig.ident",
  hormonal_col      = "",   # characteristic column name inside hormonal_src
  hormonal_treated_lbl = "",  # value meaning "treated" (from S3_00c census)
  lineage_col  = "lineage"  # atlas metadata column with the 8 S2_02 lineage patterns
)

say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }

require_config <- function(keys) {
  bad <- keys[vapply(keys, function(k) !nzchar(CFG[[k]]), logical(1))]
  if (length(bad))
    stop("FATAL: config not confirmed after S3_00 census: ", paste(bad, collapse = ", "),
         " -- fill scripts/S3_common.R from logs/S3_00_preflight.txt, then re-run (R3/R5)")
}

load_meta <- function() {
  if (!file.exists(CFG$meta)) stop("FATAL: missing ", CFG$meta)
  read.delim(CFG$meta, comment.char = "#", stringsAsFactors = FALSE)
}

load_frozen <- function() {
  if (!file.exists(CFG$frozen)) stop("FATAL: missing ", CFG$frozen)
  ## na.strings = character(): the frozen TSV uses the LITERAL string "NA" in
  ## covariate_col for contrasts without a covariate; read.delim's default
  ## would silently convert it to a real NA and crash run_frozen_contrast's
  ## covariate check (measured on the lead machine, S3_01 log 2026-09-19).
  fr <- read.delim(CFG$frozen, comment.char = "#", stringsAsFactors = FALSE,
                   colClasses = "character", na.strings = character())
  need <- c("contrast_id", "cohort", "level", "group_test", "group_ref",
            "design", "covariate_col", "exclude_samples", "note")
  if (!all(need %in% colnames(fr))) stop("FATAL: frozen TSV missing columns: ",
                                         paste(setdiff(need, colnames(fr)), collapse = ", "))
  fr
}

gsm_of <- function(x) regmatches(x, regexpr("GSM[0-9]+", x))

## map normalized object colnames -> metadata rows (one row per colname;
## GEO samples with two CEL files map both files to the same condition)
map_cols <- function(cohort, colnms, meta) {
  md <- meta[meta$gse == cohort, ]
  if (!nrow(md)) stop("FATAL: no metadata rows for ", cohort)
  m <- data.frame(colname = colnms, stringsAsFactors = FALSE)
  m$gsm <- gsm_of(colnms)
  idx <- match(m$gsm, md$gsm)
  if (any(is.na(idx)))
    stop("FATAL: ", cohort, " colnames without GSM match in metadata: ",
         paste(colnms[is.na(idx)], collapse = ", "))
  m$condition <- md$condition[idx]
  m$detail    <- md$detail[idx]
  m
}

load_bulk_matrix <- function(cohort) {
  p <- unname(CFG$bulk_rds[cohort])
  if (is.null(p)) stop("FATAL: no bulk object registered for ", cohort)
  if (!file.exists(p)) stop("FATAL: missing normalized object ", p,
                            " -- run S2_03 first (R2)")
  obj <- readRDS(p)
  mat <- if (cohort == CFG$agilent) obj$E else Biobase::exprs(obj)
  say("  ", cohort, " normalized matrix:", nrow(mat), "x", ncol(mat), "(from", p, ")")
  mat
}

split_semicolon <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

## ---- per-dataset group map (scripts/S3_group_map.tsv, frozen) --------------
## Returns per-cell TEST/REF tags from (dataset, group_col value). Any
## (dataset, condition) pair in the atlas not covered by the map stops the
## script FATAL, printing the uncovered pairs -- no silent absorption (R9).
load_group_map <- function() {
  p <- "scripts/S3_group_map.tsv"
  if (!file.exists(p)) stop("FATAL: missing ", p)
  gm <- read.delim(p, comment.char = "#", stringsAsFactors = FALSE,
                   colClasses = "character")
  need <- c("dataset", "label_test", "label_ref")
  if (!all(need %in% colnames(gm)))
    stop("FATAL: ", p, " missing columns: ", paste(setdiff(need, colnames(gm)), collapse = ", "))
  gm <- gm[nzchar(gm$dataset), , drop = FALSE]
  if (any(!nzchar(gm$label_test)) || any(!nzchar(gm$label_ref)))
    stop("FATAL: ", p, " has empty label cells -- fill every row from the census (R3)")
  gm
}

map_cell_group <- function(ds, cond) {
  gm <- load_group_map()
  out <- rep(NA_character_, length(ds))
  for (i in seq_len(nrow(gm))) {
    hit <- ds == gm$dataset[i]
    ## semicolon-separated label lists per arm (lead decision 2026-09-20:
    ## GSE179640 test = EuE;EcO;EcP;EcPA;EOR -- all patient tissues kept, no
    ## information deleted; precedent GSE25628 merged-patient contrast)
    lt <- split_semicolon(gm$label_test[i])
    lr <- split_semicolon(gm$label_ref[i])
    both <- intersect(lt, lr)
    if (length(both))
      stop("FATAL: label(s) in BOTH test and ref for ", gm$dataset[i], ": ",
           paste(both, collapse = ", "))
    hit_t <- hit & cond %in% lt
    hit_r <- hit & cond %in% lr
    if (any(hit_t & !is.na(out))) stop("FATAL: overlapping label_test in group map row ", i)
    if (any(hit_r & !is.na(out))) stop("FATAL: overlapping label_ref in group map row ", i)
    out[hit_t] <- "TEST"
    out[hit_r] <- "REF"
  }
  bad <- sort(unique(paste(ds, cond, sep = "::")[is.na(out)]))
  if (length(bad))
    stop("FATAL: (dataset, condition) pairs not covered by scripts/S3_group_map.tsv:\n  ",
         paste(bad, collapse = "\n  "),
         "\n  -- copy the exact values from logs/S3_00_preflight.txt section 3 (R2/R3)")
  out
}

## ---- hormonal map for GSE179640 (GEO series matrix, GSM-accession join) ---
## Reads CFG$hormonal_col from the GEO-extracted characteristics TSV (written
## by S3_00c) and returns GSM accession -> value, with a per-GSM consistency
## check (fail-loud). The caller joins atlas cells by
## gsm_of(md[[CFG$hormonal_join_col]]) and fail-louds on any unmatched cell.
load_hormonal_map <- function() {
  if (!nzchar(CFG$hormonal_col))
    stop("FATAL: hormonal_col empty -- run S3_00c census and fill scripts/S3_common.R (R3)")
  if (!file.exists(CFG$hormonal_src)) stop("FATAL: missing ", CFG$hormonal_src,
    " -- run scripts/S3_00c_geo_hormonal_census.R first (R2)")
  raw <- readLines(CFG$hormonal_src, warn = FALSE)
  tab <- read.delim(textConnection(raw[!startsWith(raw, "#")]),
                    stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("geo", CFG$hormonal_col)
  miss <- setdiff(need, colnames(tab))
  if (length(miss))
    stop("FATAL: ", CFG$hormonal_src, " missing columns: ", paste(miss, collapse = ", "),
         " | available: ", paste(colnames(tab), collapse = " | "))
  per <- tapply(as.character(tab[[CFG$hormonal_col]]), tab$geo,
                function(z) length(unique(z)))
  if (any(per > 1))
    stop("FATAL: hormonal label varies within GSM accession: ",
         paste(names(per)[per > 1], collapse = ", "), " -- escalate (R9)")
  map <- data.frame(
    gsm = names(per),
    horm = vapply(split(as.character(tab[[CFG$hormonal_col]]), tab$geo),
                  function(z) z[1], character(1)),
    stringsAsFactors = FALSE)
  say("  hormonal map:", nrow(map), "GSM accessions; labels:",
      paste(sort(unique(map$horm)), collapse = " | "))
  map
}

## Seurat 5 multi-layer safe assay extraction (JoinLayers route on demand)
## v2.8 (measured 2026-09-20 on the lead machine): SeuratObject >= 5.1 made the
## `slot` argument DEFUNCT (lifecycle::deprecate_stop, not a warning) -- both
## calls must use `layer =`. The JoinLayers fallback stays: on multi-layer
## objects GetAssayData(layer=...) still errors, which is the designed trigger.
get_rna <- function(obj, slot = c("counts", "data")) {
  slot <- match.arg(slot)
  tryCatch(SeuratObject::GetAssayData(obj, assay = "RNA", layer = slot),
           error = function(e) {
             say("  GetAssayData(", slot, ") failed (", conditionMessage(e), ") -- JoinLayers route")
             obj2 <- SeuratObject::JoinLayers(obj, assay = "RNA")
             SeuratObject::GetAssayData(obj2, assay = "RNA", layer = slot)
           })
}

## ---- RNA assay preparation: joined layers + guaranteed log-normalized data --
## MEASURED 2026-09-21 (S3_04 halt): expression objects here carry raw counts
## in the RNA assay and leave `data` EMPTY (the S2_02 atlas by integration
## design; the per-dataset objects likewise). Single fix point: JoinLayers
## once, then NormalizeData (LogNormalize log1p -- exactly the frozen z-mean
## normalization) ONLY when the data layer is missing/empty. Traced, never
## silent; idempotent.
prep_rna_norm <- function(obj) {
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  have_data <- tryCatch({
    d <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "data")
    length(rownames(d)) > 0
  }, error = function(e) FALSE)
  if (have_data) {
    say("  RNA data layer present -- normalization skipped (idempotent)")
  } else {
    say("  RNA data layer EMPTY (measured design) -- NormalizeData log1p (LogNormalize)")
    obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }
  obj
}

## ---- atlas annotations + per-dataset expression join (S3_00d/e/f measured) --
## The atlas (S2_02_integrated.rds) holds ALL S3 annotations (dataset, sample,
## condition, lineage) for 233,459 cells but only a 28-gene panel. The full
## expression lives in the per-dataset objects, cell-for-cell 1:1 (S3_00f:
## 109398/40801/83260 exact). The measured cell key is
##   atlas rowname == paste0(orig.ident, "__", per_dataset_colname)
## (DOUBLE underscore -- S3_00f tried single "_" and stripped-suffix forms:
## all zero matches). join_annotations() is the single join point, fail-loud.
atlas_annotations <- function() {
  if (!file.exists(CFG$atlas)) stop("FATAL: missing ", CFG$atlas, " (R2)")
  obj <- readRDS(CFG$atlas)
  md <- obj@meta.data
  rm(obj); gc()
  md
}

cell_key <- function(orig, bc) paste0(as.character(orig), "__", as.character(bc))

## orig/bc are optional overrides (default: read from the object). They exist
## so the smoke test can exercise the duplicated-key gate WITHOUT renaming
## Seurat cells -- duplicate colnames are forbidden by Seurat itself long
## before this gate would see them (measured 2026-09-21 in smoke leg 11).
join_annotations <- function(per_obj, md, ds, orig = NULL, bc = NULL) {
  sub <- md[as.character(md$dataset) == ds, , drop = FALSE]
  if (!nrow(sub)) stop("FATAL: dataset ", ds, " absent from atlas annotations")
  if (any(duplicated(rownames(sub))))
    stop("FATAL: duplicated atlas cell keys for ", ds, " (join unsafe, R9)")
  if (is.null(orig)) orig <- per_obj$orig.ident
  if (is.null(bc))   bc   <- colnames(per_obj)
  key <- cell_key(orig, bc)
  if (any(duplicated(key)))
    stop("FATAL: duplicated (orig.ident, barcode) keys inside ", ds,
         " object -- barcode collision within a sample (R9)")
  idx <- match(key, rownames(sub))
  n_miss <- sum(is.na(idx))
  if (n_miss > 0)
    stop("FATAL: ", ds, " cells without atlas annotation: ", n_miss, " of ",
         length(key), " | key examples: ",
         paste(head(key[is.na(idx)], 5), collapse = " | "),
         " -- measured key format changed? re-run S3_00f (R2/R9)")
  say("  join", ds, ":", length(key), "cells matched to atlas annotations (100%)")
  ann <- sub[idx, c("condition", "lineage", "sample"), drop = FALSE]
  ann$cell_key <- key
  ann
}

## core limma contrast from one frozen-TSV row; returns full topTable
run_frozen_contrast <- function(mat, m, spec) {
  test_lbls <- split_semicolon(spec$group_test)
  ref_lbls  <- split_semicolon(spec$group_ref)
  excl      <- split_semicolon(spec$exclude_samples)
  keep <- !(m$colname %in% excl | m$gsm %in% excl)
  mm <- m[keep, , drop = FALSE]
  X  <- mat[, keep, drop = FALSE]
  outside <- setdiff(mm$condition, c(test_lbls, ref_lbls))
  if (length(outside))
    stop("FATAL: ", spec$contrast_id, " conditions outside frozen groups: ",
         paste(outside, collapse = ", "))
  is_test <- mm$condition %in% test_lbls
  is_ref  <- mm$condition %in% ref_lbls
  if (sum(is_test) < 3 || sum(is_ref) < 3)
    stop("FATAL: ", spec$contrast_id, " too few samples: test=", sum(is_test),
         " ref=", sum(is_ref))
  grp <- factor(ifelse(is_test, "TEST", "REF"), levels = c("REF", "TEST"))

  if (spec$design == "paired") {
    pr <- read.delim(CFG$pairs11691, comment.char = "#", stringsAsFactors = FALSE,
                     na.strings = character())
    mm$woman <- pr$woman[match(mm$gsm, pr$eutopic_gsm)]
    mm$woman[is.na(mm$woman)] <- pr$woman[match(mm$gsm[is.na(mm$woman)], pr$ectopic_gsm)]
    if (any(is.na(mm$woman))) stop("FATAL: ", spec$contrast_id,
                                   " samples without woman pairing")
    tb <- table(mm$woman, mm$condition)
    if (nrow(tb) != 9 || any(tb != 1))
      stop("FATAL: ", spec$contrast_id, " pairing audit failed (expect 9 women x 2): ",
           paste(capture.output(print(tb)), collapse = " "))
    design <- model.matrix(~ grp)
    corfit <- limma::duplicateCorrelation(X, design, block = mm$woman)
    say("  ", spec$contrast_id, " paired consensus cor =", round(corfit$consensus, 3))
    fit <- limma::lmFit(X, design, block = mm$woman, correlation = corfit$consensus)
  } else {
    design <- model.matrix(~ grp)
    if (!is.null(spec$covariate_col) && nzchar(spec$covariate_col) && spec$covariate_col != "NA") {
      cc <- mm[[spec$covariate_col]]
      if (is.null(cc)) stop("FATAL: ", spec$contrast_id, " covariate column '",
                            spec$covariate_col, "' absent")
      design <- model.matrix(~ grp + factor(cc))
      say("  ", spec$contrast_id, " covariate levels:", paste(sort(unique(as.character(cc))), collapse = " | "))
    }
    fit <- limma::lmFit(X, design)
  }
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = "grpTEST", number = Inf, sort.by = "none")
  data.frame(contrast_id = spec$contrast_id, gene = rownames(tt), tt,
             n_test = sum(is_test), n_ref = sum(is_ref),
             row.names = NULL, stringsAsFactors = FALSE)
}

## find the SenMayo GMT under reference/ or data/ (census + reuse)
find_senmayo_gmt <- function() {
  cands <- c(list.files("reference", pattern = "\\.gmt$", recursive = TRUE, full.names = TRUE),
             list.files("data", pattern = "\\.gmt$", recursive = TRUE, full.names = TRUE))
  if (length(cands) == 0) return(NA_character_)
  hits <- cands[grepl("senmayo", cands, ignore.case = TRUE)]
  if (length(hits) == 1) return(hits)
  if (length(hits) == 0 && length(cands) == 1) return(cands)
  NA_character_
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sets <- lapply(lines, function(ln) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    list(name = f[1], genes = f[-c(1:2)])
  })
  names(sets) <- vapply(sets, `[[`, character(1), "name")
  sets
}

## probe -> symbol collapse (frozen rule): per symbol keep the probe with the
## highest mean expression; rownames BECOME the symbol (ssGSEA matches sets
## by rownames)
collapse_maxmean <- function(mat, symbols) {
  ok <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[ok, , drop = FALSE]; symbols <- symbols[ok]
  me <- rowMeans(mat)
  idx <- tapply(seq_len(nrow(mat)), symbols, function(ii) ii[which.max(me[ii])])
  out <- mat[unname(idx), , drop = FALSE]
  rownames(out) <- names(idx)
  out
}

## ---- probe -> symbol collapse for bulk DE (v3.0.2; ported from S3_02) ------
## Frozen rule (S3_design): per symbol keep the probe with the highest mean
## expression. The annotation DB registry and both mapping paths (Affy .db +
## Agilent SystematicName->SYMBOL) are IDENTICAL to S3_02_ssgsea_bulk.R,
## which has run successfully on this machine (R1: no new untested code).
## agilent_symbols/symbols_for below are VERBATIM copies from S3_02.
BULK_DB <- c(GSE7305 = "hgu133plus2.db", GSE6364 = "hgu133plus2.db",
             GSE25628 = "hgu133a2.db", GSE7307 = "hgu133plus2.db",
             GSE51981 = "hgu133plus2.db", GSE11691 = "hgu133a.db")

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
    db <- unname(BULK_DB[cohort])
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

## Load one bulk cohort as a SYMBOL-level matrix (probe->symbol, max-mean).
## v3.0.2: this is the single entry point for bulk DE; probe-level matrices
## must not reach limma (measured: silent cross-cohort namespace split).
load_bulk_symbol_matrix <- function(cohort) {
  p <- unname(CFG$bulk_rds[cohort])
  if (is.null(p)) stop("FATAL: no bulk object registered for ", cohort)
  if (!file.exists(p)) stop("FATAL: missing normalized object ", p,
                            " -- run S2_03 first (R2)")
  obj <- readRDS(p)
  mat <- if (cohort == CFG$agilent) obj$E else Biobase::exprs(obj)
  sym <- symbols_for(cohort, obj, mat)
  Xs <- collapse_maxmean(mat, sym)
  say("  ", cohort, " probe->symbol:", nrow(mat), "probes ->", nrow(Xs),
      "symbols (from", p, ")")
  Xs
}
