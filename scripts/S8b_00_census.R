#!/usr/bin/env Rscript
## ===========================================================================
## S8b_00_census.R -- GEO census for the S8b positive-control question (R7)
## ---------------------------------------------------------------------------
## Question (PI, locked): do known endometriosis drugs (progestins / danazol /
## GnRH analogues) reverse the Tier-1 lesion program in available GEO cohorts?
##
## This script is a CENSUS, not an analysis (R3/R7): it measures what exists
## BEFORE any hypothesis is designed. Nothing here touches expression matrices.
##
## Charter: R1 zero fabrication | R2 verify actual files | R3 locked scope
## R4 never edit results files | R5 locked thresholds | R6 failures explicit
## R7 this design | R8 PI runs on her machine | R9 every verdict from output.
##
## Data source: NCBI GEO via E-utilities (esearch + esummary), no API key.
## Packages: httr + jsonlite (already in renv from S8). Base R otherwise.
##
## Output (results/S8b_census/):
##   S8b_00_raw_hits.csv       -- every GSE returned by each query block
##   S8b_00_census_summary.md  -- machine-printed verdicts (R9)
##
## v1.0 2026-09-23 -- designed in session, locked by PI before any run.
## ===========================================================================

suppressPackageStartupMessages({library(httr); library(jsonlite)})

OUT <- "results/S8b_census"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

EUTILS <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL   <- "endometriosis_immunoediting_s8b"
EMAIL  <- "pi@example.org"   # <-- PI: set your email (NCBI polite-pool)

## ---- locked query blocks (S8b design doc section 3) ------------------------
## Each block = one drug class. "endometriosis" AND drug term, GSE only.
queries <- list(
  progestin_dienogest    = '(endometriosis[All Fields]) AND (dienogest[All Fields]) AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
  progestin_medroxy      = '(endometriosis[All Fields]) AND ("medroxyprogesterone"[All Fields] OR "MPA"[All Fields]) AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
  progestin_gestrinone   = '(endometriosis[All Fields]) AND (gestrinone[All Fields]) AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
  danazol                = '(endometriosis[All Fields]) AND (danazol[All Fields]) AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
  gnrh                   = '(endometriosis[All Fields]) AND ("gonadotropin releasing hormone"[All Fields] OR "GnRH"[All Fields] OR "leuprolide"[All Fields] OR "leuprorelin"[All Fields] OR "triptorelin"[All Fields] OR "goserelin"[All Fields] OR "nafarelin"[All Fields]) AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])',
  broad_drug_catchall    = '(endometriosis[All Fields]) AND ("drug therapy"[All Fields] OR "treatment"[All Fields] OR "hormone"[All Fields]) AND ("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])'
)

get_json <- function(url) {
  for (attempt in 1:4) {
    r <- tryCatch(GET(url, user_agent(TOOL), query = list(tool = TOOL, email = EMAIL)),
                  error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200)
      return(fromJSON(content(r, as = "text", encoding = "UTF-8")))
    Sys.sleep(2 * attempt)
  }
  stop("FATAL: eutils unreachable after 4 attempts (R6): ", url)
}

esearch_count_ids <- function(term) {
  u <- paste0(EUTILS, "/esearch.fcgi?db=gds&retmax=10000&retmode=json&term=",
              URLencode(term, reserved = TRUE))
  j <- get_json(u)
  list(count = as.integer(j$esearchresult$count),
       ids   = j$esearchresult$idlist)
}

esummary <- function(ids) {
  if (length(ids) == 0) return(data.frame())
  u <- paste0(EUTILS, "/esummary.fcgi?db=gds&retmode=json&id=",
              paste(ids, collapse = ","))
  j <- get_json(u)
  res <- j$result
  rows <- lapply(setdiff(names(res), "uids"), function(uid) {
    e <- res[[uid]]
    data.frame(
      uid        = uid,
      gse        = ifelse(is.null(e$accession), NA, e$accession),
      title      = ifelse(is.null(e$title), NA, e$title),
      summary    = ifelse(is.null(e$summary), NA, e$summary),
      taxon      = ifelse(is.null(e$taxon), NA, e$taxon),
      n_samples  = ifelse(is.null(e$n_samples), NA, e$n_samples),
      gds_type   = ifelse(is.null(e$gdstype), NA, e$gdstype),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

say("== S8b_00 census: measuring GEO availability before any design (R3/R7) ==")

all_rows <- list()
per_block <- data.frame(block = character(), n_raw = integer(), stringsAsFactors = FALSE)

for (nm in names(queries)) {
  say("-- block: ", nm)
  h <- esearch_count_ids(queries[[nm]])
  per_block <- rbind(per_block, data.frame(block = nm, n_raw = h$count, stringsAsFactors = FALSE))
  d <- esummary(h$ids)
  if (nrow(d)) {
    d$block <- nm
    all_rows[[nm]] <- d
  }
  Sys.sleep(0.4)  # NCBI: <=3 req/sec without key
}

raw <- if (length(all_rows)) do.call(rbind, all_rows) else data.frame()
if (nrow(raw)) {
  raw <- raw[!duplicated(raw$gse), ]   # a GSE may hit several blocks
  rownames(raw) <- NULL
}
write.csv(raw, file.path(OUT, "S8b_00_raw_hits.csv"), row.names = FALSE)

## ---- census verdicts (R9: printed from measured rows only) ------------------
## A GSE counts as a positive-control candidate ONLY if ALL hold:
##  (a) human, (b) names a locked drug class, (c) patient/lesion material,
##  (d) is NOT an in-vitro / animal / engineered model (those are exploratory
##      at best and are listed separately, never pooled as a patient control).
txt  <- tolower(paste(raw$title, raw$summary, sep = " | "))
is_human   <- grepl("homo sapiens", tolower(raw$taxon))
drug_rx    <- grepl("dienogest|medroxyprogesterone|danazol|gestrinone|leuprolide|leuprorelin|triptorelin|goserelin|nafarelin|gnrh agonist|progesterone treatment|progestin", txt)
patient_rx <- grepl("patient|women with|biops|in vivo|pre.?treatment|post.?treatment", txt)
invitro_rx <- grepl("in vitro|cell line|primary cell|stromal cell|epithelial cell|decidualiz|organoid|ipsc|knockdown|sirna|mouse|rat |rattus|rabbit|mus musculus|xenograft|granulosa|menstrual fluid|blood|plasma|exosomal|extracellular vesicle", txt)

cand <- raw[is_human & drug_rx & patient_rx & !invitro_rx, ]
vitro_drug <- raw[is_human & drug_rx & invitro_rx, ]

say("== S8b_00 CENSUS VERDICTS ==")
cat("BLOCK COUNTS (raw hits per drug class):\n")
print(per_block, row.names = FALSE)
cat("\nTOTAL UNIQUE GSEs AFTER DEDUP:", nrow(raw), "\n")
cat("danazol/gestrinone datasets measured:", sum(per_block$n_raw[per_block$block %in% c("danazol","progestin_gestrinone")]), "\n")
cat("\nHUMAN IN-VIVO DRUG COHORT CANDIDATES (a+b+c+not d):\n")
if (nrow(cand)) {
  print(cand[, c("gse","title","n_samples")], row.names = FALSE)
} else {
  cat("  NONE MEASURED\n")
}
cat("\nHUMAN IN-VITRO DRUG-TREATMENT SERIES (exploratory-only pool):\n")
if (nrow(vitro_drug)) {
  print(vitro_drug[, c("gse","title","n_samples")], row.names = FALSE)
} else {
  cat("  NONE MEASURED\n")
}

human_in_vivo <- nrow(cand) > 0
cat("CENSUS CONCLUSION: ",
    ifelse(human_in_vivo,
           "positive-control test in patient cohorts is DESIGNABLE -- proceed to S8b design doc gate.",
           "positive-control test in patient cohorts is NOT designable in GEO (R1: nothing to pool). Only in-vitro drug-treatment series (e.g. GSE75423 dienogest) remain; any such test is EXPLORATORY by design and is the LAST in-silico attempt (stop rule, S8b design doc section 5)."),
    "\n", sep = "")
say("== S8b_00 census DONE -- outputs under ", OUT, " ==")
