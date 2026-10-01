# ============================================================================
# S3_00c_geo_hormonal_census.R -- read-only census: extract GSE179640 sample
# characteristics from the OFFICIAL GEO series matrix and audit the join key
# vs the per-dataset object and the atlas.
# PARSER FACTS (measured 2026-09-19 from the official file):
#   - Fields are double-quoted: "GSM6102532", "condition: Control" -> quotes
#     MUST be stripped (bugs in v2.3/v2.4 kept them -> 0/33 joins).
#   - Lines are KEY \t VALUE \t VALUE ... with NO standalone "=" token ->
#     dropping a fixed first token would lose the first sample of every line
#     (bug in v2.3/v2.4: 58 instead of 59 samples).
#   - Series summary states the atlas is "of the disease in hormonally
#     treated patients" -> treatment is a SERIES-LEVEL CONSTANT, which is why
#     no per-sample treatment characteristic exists. Printed verbatim below.
# WRITES (reviewed artifact): scripts/S3_gse179640_geo_characteristics.tsv
#   (one row per GSM; feeds load_hormonal_map() in S3_common.R).
# Usage (from project root):
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_00c_geo_hormonal_census.R data\raw\GSE179640\GSE179640_series_matrix.txt.gz > logs\S3_00c_geo_hormonal_census.txt 2>&1
# If the matrix is absent at the given path it is downloaded ONCE from the
# GEO FTP (standard series-matrix URL) to that path, then reused offline.
# ============================================================================
set.seed(42)
source("scripts/S3_common.R")
say("== S3_00c GEO hormonal census (read-only, parser v3) ==")

args <- commandArgs(trailingOnly = TRUE)
sm_path <- if (length(args) >= 1) args[1] else
  file.path("data", "raw", "GSE179640", "GSE179640_series_matrix.txt.gz")
geo_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/",
                  "matrix/GSE179640_series_matrix.txt.gz")

if (!file.exists(sm_path)) {
  say("series matrix not found at", sm_path, "-- downloading once from GEO FTP")
  dir.create(dirname(sm_path), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(geo_url, destfile = sm_path, mode = "wb")
  if (!file.exists(sm_path)) stop("FATAL: download failed: ", geo_url)
}
say("series matrix:", normalizePath(sm_path))

lines <- readLines(gzfile(sm_path), warn = FALSE)

## robust tokeniser: KEY [\t(=)\t] VALUE ... ; strip surrounding double quotes
unquote <- function(v) {
  v <- sub("^\"", "", v)  # one leading double quote
  sub("\\\"$", "", v)     # one trailing double quote
}
parse_bang_line <- function(ln) {
  z <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (length(z) < 2) stop("FATAL: unexpected series-matrix line format: ",
                          substr(ln, 1, 80), " -- escalate (R9)")
  key <- z[1]
  vals <- if (z[2] == "=") z[-c(1:2)] else z[-1]
  list(key = key, vals = unquote(vals))
}

bl <- lines[startsWith(lines, "!Sample_")]
if (!length(bl)) stop("FATAL: no !Sample_ lines -- is this a GEO series matrix?")

titles <- accs <- NULL
chars <- list()
for (ln in bl) {
  p <- parse_bang_line(ln)
  if (p$key == "!Sample_title") titles <- p$vals
  else if (p$key == "!Sample_geo_accession") accs <- p$vals
  else if (p$key == "!Sample_characteristics_ch1") {
    kk <- trimws(sub(":.*$", "", p$vals))
    vv <- trimws(sub("^[^:]*:[[:space:]]*", "", p$vals))
    if (length(unique(kk)) != 1)
      stop("FATAL: mixed characteristic keys on one line: ",
           paste(unique(kk), collapse = " / "), " -- escalate (R9)")
    k <- kk[1]
    if (!nzchar(k) || grepl("!", k)) stop("FATAL: unparsable characteristic key -- escalate (R9)")
    if (!is.null(chars[[k]]))
      stop("FATAL: duplicated characteristic key '", k, "' -- escalate (R9)")
    chars[[k]] <- vv
  }
}
n <- length(accs)
if (is.null(titles) || n < 1) stop("FATAL: could not parse Sample_title / Sample_geo_accession")
if (length(titles) != n) stop("FATAL: title/accession length mismatch (", length(titles), " vs ", n, ")")

## internal cross-check: !Series_sample_id must list exactly the same count
sid <- lines[startsWith(lines, "!Series_sample_id")]
if (length(sid)) {
  p <- parse_bang_line(sid[1])
  ## value is ONE space-separated quoted field, not one token per sample
  n_sid <- if (length(p$vals) == 1) length(strsplit(p$vals[1], "\\s+")[[1]]) else length(p$vals)
  say("cross-check: !Series_sample_id lists", n_sid, "samples | parsed", n, "samples")
  if (n_sid != n) stop("FATAL: sample-count cross-check failed (", n_sid, " != ", n,
                       ") -- parser bug, escalate with the log (R9)")
} else say("cross-check: no !Series_sample_id line (not fatal)")

bad_len <- names(chars)[vapply(chars, length, integer(1)) != n]
if (length(bad_len)) stop("FATAL: characteristic length mismatch: ",
                          paste(bad_len, collapse = ", "))
say("parsed", n, "GSM samples |", length(chars), "characteristics:")
say("  ", paste(names(chars), collapse = " | "))
say("  accessions head:", paste(head(accs, 3), collapse = " "), "... tail:",
    paste(tail(accs, 3), collapse = " "))

## ---- series-level hormonal statement (verbatim, from the same file) ------
terms <- "hormon|treat|medicat|pill|gnrh|progest|therap|drug|oestrogen|estrogen|contracept"
st <- lines[startsWith(lines, "!Series_title")]
ss <- lines[startsWith(lines, "!Series_summary")]
if (length(st)) { p <- parse_bang_line(st[1]); say("series title:", p$vals[1]) }
say("-- series summary, treatment-relevant sentences (verbatim) --")
found_series_treat <- FALSE
if (length(ss)) {
  p <- parse_bang_line(ss[1])
  sents <- unlist(strsplit(p$vals[1], "(?<=[.!?])\\s+", perl = TRUE))
  for (s in sents) if (grepl(terms, s, ignore.case = TRUE)) {
    found_series_treat <- TRUE
    say("  *", s)
  }
  if (!found_series_treat) say("  (no treatment-like sentence in the series summary)")
} else say("  (no !Series_summary line)")

df <- data.frame(geo = accs, title = titles, stringsAsFactors = FALSE,
                 check.names = FALSE)
for (k in names(chars)) df[[k]] <- chars[[k]]

treat_cols <- grep("treat|horm|therap|med|drug|pill|gnrh|progest", names(chars),
                   ignore.case = TRUE, value = TRUE)
cand <- if (length(treat_cols)) treat_cols else names(chars)
if (!length(treat_cols)) say("  NO treatment/hormone-named characteristic in the matrix")
say("-- candidate columns and their value sets --")
for (cc in cand) {
  u <- sort(unique(as.character(df[[cc]])))
  say("  ", cc, "=", if (length(u) <= 20) paste(u, collapse = " | ")
      else paste0("(", length(u), " distinct values)"))
}
say("-- per-GSM consistency (must be exactly 1 value per accession) --")
for (cc in cand) {
  per <- tapply(as.character(df[[cc]]), df$geo, function(z) length(unique(z)))
  say("  ", cc, ": GSMs with >1 distinct value =", sum(per > 1), "of", length(per))
}

out_tsv <- "scripts/S3_gse179640_geo_characteristics.tsv"
con <- file(out_tsv, open = "wt")
writeLines(c(
  "# GSE179640 GEO sample characteristics from the official GEO series matrix",
  paste0("# source: ", geo_url),
  "# extracted by S3_00c_geo_hormonal_census.R (read-only census, parser v3); reviewed by lead",
  "# one row per GSM; join key: GSM accession = gsm_of(orig.ident)",
  "# MEASURED 2026-09-19: no per-sample treatment characteristic exists; the series",
  "# summary states the cohort is 'hormonally treated patients' (series-level constant)",
  "# reproducible: re-run S3_00c to regenerate identically"), con)
write.table(df, con, sep = "\t", row.names = FALSE, quote = TRUE)
close(con)
say("wrote", out_tsv, "--", nrow(df), "rows x", ncol(df), "cols")

## ---- join audits (measured, not assumed; R9) ------------------------------
fmd <- readRDS(CFG$hormonal_obj)@meta.data
origs <- sort(unique(as.character(fmd$orig.ident)))
orig_gsm <- gsm_of(origs)
if (any(!nzchar(orig_gsm)))
  stop("FATAL: orig.ident values without a GSM accession: ",
       paste(origs[!nzchar(orig_gsm)], collapse = ", "), " -- escalate (R9)")
m1 <- intersect(orig_gsm, accs)
say("join audit A (orig.ident GSM prefix <-> GEO accession): matched", length(m1),
    "of", length(orig_gsm), "orig.ident values")
if (length(m1) != length(orig_gsm)) {
  say("  orig.ident whose GSM is absent from the series matrix:")
  for (u in origs[!(orig_gsm %in% accs)]) say("   ", u)
}
amd <- readRDS(CFG$atlas)@meta.data
has_orig <- "orig.ident" %in% colnames(amd)
say("join audit B (atlas): 'orig.ident' column present =", has_orig)
if (has_orig) {
  a_orig <- sort(unique(as.character(amd$orig.ident[amd$dataset == "GSE179640"])))
  a_gsm <- gsm_of(a_orig)
  say("  atlas GSE179640 orig.ident values:", length(a_orig),
      "| with GSM prefix:", sum(nzchar(a_gsm)), "of", length(a_orig),
      "| covered by GEO accession:", sum(a_gsm %in% accs), "of", length(a_orig))
  if (any(!(a_gsm %in% accs)))
    say("  uncovered atlas orig.ident:",
        paste(a_orig[!(a_gsm %in% accs)], collapse = ", "))
}
asmp <- sort(unique(amd$sample[amd$dataset == "GSE179640"]))
asmp_gsm <- gsm_of(asmp)
say("join audit C (diagnostic): atlas 'sample' GSM-prefix coverage in GEO:",
    sum(nzchar(asmp_gsm) > 0 & asmp_gsm %in% accs), "of", length(asmp))

if (length(m1) != length(orig_gsm))
  stop("FATAL: orig.ident GSM <-> GEO accession join is not complete -- review audits A-C above and escalate with the log (R9)")
if (!has_orig)
  stop("FATAL: atlas metadata lacks 'orig.ident' -- escalate with the log (R9)")
say("JOIN KEY CONFIRMED: GSM accession = gsm_of(orig.ident)")

## ---- measured verdict for the LEAD DECISION (R3/R5) -----------------------
if (!length(treat_cols)) {
  say("VERDICT (measured): no per-sample treatment/hormonal characteristic exists")
  say("  in the official GEO series matrix.")
  if (found_series_treat)
    say("  AND the series summary states the cohort is 'hormonally treated patients'")
  say("  -> the covariate has ZERO VARIANCE (series-level constant). A-HCOV / A-NOTX")
  say("  are structurally impossible on GSE179640 as registered. LEAD DECISION:")
  say("  Amendment to suspend both arms with this documented rationale (R6), or")
  say("  register an external supplementary source (R5). Do NOT fill hormonal_col")
  say("  with any of the 5 existing columns -- none is a treatment variable.")
}
say("CENSUS DONE -- fill the GSE179640 row of scripts/S3_group_map.tsv from the")
say("  atlas condition values in logs/S3_00_preflight.txt (R3).")
