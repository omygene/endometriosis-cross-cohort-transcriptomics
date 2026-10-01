# ============================================================================
# S2_03_meta_build.R  v2 -- independently rebuild + verify S1_bulk_conditions.tsv
#
# Changelog v1 -> v2 (measured failure on the lead machine, 2026-09-18):
#   v1 stopped fail-loud at step 3 with "FATAL: unmapped samples in GSE7307"
#   because the GSE7307 series matrix contains 677 samples (cell lines and
#   other tissues of the Roth body-index series) while only the FROZEN 41-GSM
#   endometrium subset (scripts/S1_GSE7307_subset.txt) is in study scope.
#   v2 restricts GSE7307 to the frozen subset BEFORE applying mapping rules
#   (hard stop if any frozen GSM is absent from the series or if the subset
#   is not exactly 41 rows). Out-of-scope series members are excluded by
#   design, not "unmapped". All other logic, rules, and fail-loud behavior
#   are unchanged. The final cell-by-cell comparison against the shipped
#   TSV (322 rows) remains the closing gate.
#
# Purpose (governance): the shipped label file data/metadata/S1_bulk_conditions.tsv
# (322 rows, MD5 f8323c5f0e900a5bd2311b913b68762b) was derived from GEO OFFICIAL
# series-matrix metadata. This script re-downloads the 7 series matrices on the
# lead machine, re-derives every label with the frozen rules below, and compares
# the result cell-by-cell with the shipped file. ANY mismatch -> stop (R1/R2/R9).
#
# Extra hard checks:
#   - GSE7307: the 41 frozen GSMs (scripts/S1_GSE7307_subset.txt) must all be
#     present and the GEO-derived labels must agree with the frozen list 41/41.
#   - GSE25628: title classes must be Ectopic 8 / Eutopic 8 / Normal 6 and donor
#     classes patient 16 / healthy 6 (the E3 table, logs/S2_entry_checks_fix.txt).
#   - GSE120103: the 5 documented design-mismatch exclusions must be exactly
#     GSM3393522..GSM3393526 and all must carry Group 2B labels.
#
# Network: one small download per series from ftp.ncbi.nlm.nih.gov (https).
# Writes ONLY: data/metadata/geo_series_matrix/*.txt.gz (archived inputs, MD5-logged).
# R10 RAM estimate: < 0.5 GB (metadata rows only; GSE7307 series text ~170 MB
#   gz, read line-streamed). Runtime: ~2-4 min (GSE7307 download is the long pole).
#
# Run (from C:\endometriosis_immunoediting):
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S2_03_meta_build.R > logs\S2_03_meta_build.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }
options(timeout = 300)

SERIES <- c("GSE7305", "GSE6364", "GSE25628", "GSE7307", "GSE51981", "GSE11691", "GSE120103")
mdir <- "data/metadata/geo_series_matrix"
dir.create(mdir, recursive = TRUE, showWarnings = FALSE)

## ---- 0. load the frozen GSE7307 subset (scope definition) --------------------
fr <- read.delim("scripts/S1_GSE7307_subset.txt", comment.char = "#", header = FALSE,
                 col.names = c("s1_group", "gsm"), stringsAsFactors = FALSE)
if (nrow(fr) != 41) stop("FATAL: frozen GSE7307 subset must be exactly 41 rows, found ", nrow(fr))
if (!all(fr$s1_group %in% c("endometriosis", "normal_endometrium")))
  stop("FATAL: frozen GSE7307 subset has unexpected group labels")
say("== 0/4: frozen GSE7307 subset loaded:", nrow(fr), "GSMs (",
    sum(fr$s1_group == "endometriosis"), "disease /", sum(fr$s1_group == "normal_endometrium"), "normal) ==")

## ---- 1. download -------------------------------------------------------------
say("== 1/4: download 7 series matrices (GEO official) ==")
for (g in SERIES) {
  n <- as.integer(sub("GSE", "", g))
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%dnnn/%s/matrix/%s_series_matrix.txt.gz",
                 floor(n / 1000), g, g)
  dest <- file.path(mdir, paste0(g, "_series_matrix.txt.gz"))
  if (!file.exists(dest)) {
    say("  GET", url)
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE),
             error = function(e) stop("FATAL: download failed for ", g, ": ", conditionMessage(e)))
  }
  if (!file.exists(dest) || file.size(dest) < 10000)
    stop("FATAL: ", dest, " missing or implausibly small")
  say("  ", g, "->", basename(dest), "| bytes:", file.size(dest), "| md5:", unname(tools::md5sum(dest)))
}

## ---- 2. parse ----------------------------------------------------------------
say("== 2/4: parse sample metadata rows ==")
parse_matrix <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  ln  <- readLines(con, warn = FALSE)
  keep <- grep("^!Sample_(geo_accession|title|source_name_ch1|characteristics_ch1)\t", ln, value = TRUE)
  if (length(keep) < 3) stop("FATAL: ", path, " yielded ", length(keep), " sample rows (parse failure?)")
  sp   <- strsplit(keep, "\t", fixed = TRUE)
  keys <- vapply(sp, function(x) x[1], "")
  get1 <- function(key) {
    v <- sp[keys == key]
    if (!length(v)) return(NULL)
    gsub('^"|"$', "", v[[1]][-1])
  }
  acc <- get1("!Sample_geo_accession"); ttl <- get1("!Sample_title"); src <- get1("!Sample_source_name_ch1")
  ch  <- get1("!Sample_characteristics_ch1")   # first characteristics row only
  if (is.null(acc) || is.null(ttl)) stop("FATAL: ", path, " lacks accession/title rows")
  data.frame(gsm = acc, title = ttl,
             source = if (is.null(src)) "" else src,
             ch1 = if (is.null(ch)) "" else ch,
             stringsAsFactors = FALSE)
}
M <- lapply(SERIES, function(g) parse_matrix(file.path(mdir, paste0(g, "_series_matrix.txt.gz"))))
names(M) <- SERIES
for (g in SERIES) say("  ", g, ":", nrow(M[[g]]), "samples parsed")

## ---- 2b. scope filter: GSE7307 = frozen 41-GSM subset only (v2 fix) -----------
say("== 2b/4: restrict GSE7307 to the frozen 41-GSM subset (scope by design) ==")
in_series <- fr$gsm %in% M$GSE7307$gsm
if (!all(in_series))
  stop("FATAL: frozen GSMs absent from GSE7307 series: ",
       paste(fr$gsm[!in_series], collapse = ", "))
M$GSE7307 <- M$GSE7307[M$GSE7307$gsm %in% fr$gsm, ]
if (nrow(M$GSE7307) != 41) stop("FATAL: GSE7307 in-scope rows != 41: ", nrow(M$GSE7307))
say("  GSE7307 scoped to", nrow(M$GSE7307), "frozen GSMs; series members outside study scope excluded")

## ---- 3. frozen mapping rules (fail-loud: unmapped -> list + stop) ------------
say("== 3/4: apply frozen mapping rules ==")
map_cohort <- function(g, d) {
  cond <- rep(NA_character_, nrow(d))
  detail <- rep("", nrow(d))   # empty string by default (matches shipped TSV)
  if (g == "GSE7305") {
    cond[grepl("-Disease", d$title)] <- "Endometriosis"
    cond[grepl("-Normal", d$title)]  <- "Normal endometrium"
  } else if (g == "GSE6364") {
    cond[grepl(" Endometriosis$", d$ch1)] <- "Endometriosis eutopic"
    cond[grepl(" Normal$", d$ch1)]        <- "Normal eutopic"
    detail <- sub(" Phase (Endometriosis|Normal)$", "", d$ch1)
  } else if (g == "GSE25628") {
    cond[grepl("^Ectopic", d$title)] <- "Ectopic lesion"
    cond[grepl("^Eutopic", d$title)] <- "Eutopic (patient)"
    cond[grepl("^Normal", d$title)]  <- "Normal endometrium"
    detail <- sub("^donor: ", "", d$ch1)
  } else if (g == "GSE7307") {
    cond[grepl(" Disease$", d$title)] <- "Endometriosis"
    cond[grepl(" Normal$", d$title)]  <- "Normal endometrium"
  } else if (g == "GSE51981") {
    m <- c("Endometriosis_Moderate/Severe"                    = "Endometriosis moderate/severe",
           "Endometriosis_Minimal/Mild"                       = "Endometriosis minimal/mild",
           "Endometriosis_Severity Not Available"             = "Endometriosis (severity NA)",
           "Non-Endometriosis_Uterine Pelvic Pathology"       = "Non-EMS pelvic pathology",
           "Non-Endometriosis_No Uterine Pelvic Pathology"    = "Non-EMS no pelvic pathology")
    hit <- m[d$source]
    cond <- unname(ifelse(is.na(hit), NA_character_, hit))
    detail <- sub("^tissue: ", "", sub(" Endometrial tissue$", "", d$ch1))
  } else if (g == "GSE11691") {
    cond[d$source == "eutopic endometrium"]   <- "Eutopic (patient)"
    cond[grepl("^endometriosis", d$source)]   <- "Ectopic lesion"
    detail <- ifelse(cond == "Ectopic lesion", paste0("paired: ", d$source), "paired")
  } else if (g == "GSE120103") {
    grp <- sub(".*\\((Group [12][AB])\\).*", "\\1", d$source)
    m <- c("Group 1A" = "Fertile control", "Group 1B" = "Infertile control",
           "Group 2A" = "Stage IV ovarian EMS (2A)", "Group 2B" = "Stage IV ovarian EMS (2B)")
    hit <- m[grp]
    cond <- unname(ifelse(is.na(hit), NA_character_, hit))
    detail <- grp
  }
  if (any(is.na(cond)))
    stop("FATAL: unmapped samples in ", g, ":\n",
         paste(utils::capture.output(print(d[is.na(cond), ])), collapse = "\n"),
         "\nRules must be amended by the session before proceeding (R1).")
  data.frame(gse = g, gsm = d$gsm, condition = cond, detail = detail, stringsAsFactors = FALSE)
}
built <- do.call(rbind, Map(map_cohort, SERIES, M))
rownames(built) <- NULL
say("  mapped rows:", nrow(built))

## ---- 4. hard cross-checks + cell-by-cell comparison with shipped file --------
say("== 4/4: cross-checks + comparison with shipped S1_bulk_conditions.tsv ==")
b7 <- built[built$gse == "GSE7307", ]
idx <- match(fr$gsm, b7$gsm)
if (any(is.na(idx))) stop("FATAL: frozen GSMs missing from GSE7307 series: ",
                          paste(fr$gsm[is.na(idx)], collapse = ", "))
agree <- (b7$condition[idx] == "Endometriosis" & fr$s1_group == "endometriosis") |
         (b7$condition[idx] == "Normal endometrium" & fr$s1_group == "normal_endometrium")
if (!all(agree)) stop("FATAL: GSE7307 frozen-list disagreement at: ",
                      paste(fr$gsm[!agree], collapse = ", "))
say("  GSE7307 frozen-list cross-check: 41/41 agree [OK]")
t2 <- table(sub("^(\\w+).*", "\\1", M$GSE25628$title))
if (!identical(as.integer(t2[c("Ectopic", "Eutopic", "Normal")]), c(8L, 8L, 6L)))
  stop("FATAL: GSE25628 title classes diverge from E3 (8/8/6): ", paste(t2, collapse = "/"))
d2 <- table(built$detail[built$gse == "GSE25628"])
if (!identical(as.integer(d2[c("healthy donor", "patient")]), c(6L, 16L)))
  stop("FATAL: GSE25628 donor classes diverge from E3 (6 healthy / 16 patient)")
say("  GSE25628 E3 cross-check: Ectopic 8 / Eutopic 8 / Normal 6; patient 16 / healthy 6 [OK]")
b120_2b <- built$gsm[built$gse == "GSE120103" & built$detail == "Group 2B"]
if (!identical(sort(b120_2b), sort(c(paste0("GSM", 3393518:3393521), paste0("GSM", 3393522:3393526)))))
  stop("FATAL: GSE120103 Group 2B membership unexpected")
if (!all(paste0("GSM", 3393522:3393526) %in% b120_2b))
  stop("FATAL: documented GE1-v5 exclusions are not all Group 2B")
say("  GSE120103: 5 documented exclusions confirmed inside Group 2B [OK]")

ship <- read.delim("data/metadata/S1_bulk_conditions.tsv", comment.char = "#",
                   stringsAsFactors = FALSE)
built <- built[order(built$gse, built$gsm), ]; rownames(built) <- NULL
ship  <- ship[order(ship$gse, ship$gsm), ];   rownames(ship) <- NULL
if (nrow(ship) != nrow(built)) stop("FATAL: row count differs (shipped ", nrow(ship), " vs built ", nrow(built), ")")
same <- mapply(function(a, b) identical(a, b), built, ship)
if (!all(same)) {
  for (cl in names(built)[!same]) {
    w <- which(built[[cl]] != ship[[cl]])
    say("  MISMATCH column", cl, "at rows:", paste(head(w, 10), collapse = ", "))
  }
  stop("FATAL: shipped label file and fresh GEO rebuild differ -- escalate (R1/R2)")
}
say("  cell-by-cell comparison: all", nrow(built), "rows x 4 columns identical [OK]")
print(table(built$gse, built$condition))
say("META BUILD OK -- S1_bulk_conditions.tsv independently verified on this machine")
