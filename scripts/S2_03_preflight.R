# ============================================================================
# S2_03_preflight.R -- READ-ONLY measurement before S2_03 bulk QC (charter R9/B.4)
# Project: Benign Immunoediting -- endometriosis
# Rule: writes NOTHING to data/, results/ or figures/. Console + log only.
# Reference values source: tables/dataset_manifest.tsv (file, not chat; R1/R2).
# Run from project root:
#   cd /d C:\endometriosis_immunoediting
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S2_03_preflight.R > logs\S2_03_preflight.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }

issues <- character(0)
note_issue <- function(x) { issues <<- c(issues, x); say("  ISSUE: ", x) }

say("S2_03 preflight starting -- read-only measurement")
say("R: ", R.version.string)

## ---- 1. Package census -------------------------------------------------------
say("== 1. Package census (renv library) ==")
pkgs <- c("oligo", "limma", "arrayQualityMetrics", "affy", "ggplot2",
          "pd.hg.u133.plus.2", "pd.hg.u133a", "hgu133plus2cdf", "hgu133acdf")
pkg_ok <- logical(length(pkgs))
for (i in seq_along(pkgs)) {
  p  <- pkgs[i]
  ok <- requireNamespace(p, quietly = TRUE)
  pkg_ok[i] <- ok
  v <- if (ok) as.character(utils::packageVersion(p)) else "MISSING"
  say(sprintf("  %-22s %s", p, v))
}
core_missing <- pkgs[1:5][!pkg_ok[1:5]]
if (length(core_missing)) note_issue(paste("core packages missing:", paste(core_missing, collapse = ", ")))
say(sprintf("  AQM-route inputs: pd.* (oligo) present = %s/%s | affy CDF present = %s/%s",
            pkg_ok[6], pkg_ok[7], pkg_ok[8], pkg_ok[9]))

## ---- 2. Raw input census + integrity vs manifest -----------------------------
say("== 2. Raw input census vs dataset_manifest.tsv ==")
expect <- data.frame(
  tar     = file.path("data/raw", c("GSE51981_RAW.tar", "GSE7305_RAW.tar", "GSE25628_RAW.tar",
                                    "GSE6364_RAW.tar", "GSE11691_RAW.tar", "GSE120103_RAW.tar")),
  bytes   = c(746424320, 91064320, 43048960, 190904320, 62914560, 109291520),
  md5     = c("c3bc305129ce43ca97d93acad05effa4", "212fd3ae15bf7e1c1207c420a9da1586",
              "4b161cf132d486b8eaaf071b4774906e", "959690329da456c1f42a645c97dd3d09",
              "0c0a95ce7bef627ac778ad68147dd951", "9482b104a86f051ce86a1bca597f69b6"),
  members = c(148L, 20L, 22L, 37L, 18L, 36L),
  pat     = c("\\.CEL(\\.gz)?$", "\\.CEL(\\.gz)?$", "\\.CEL(\\.gz)?$",
              "\\.CEL(\\.gz)?$", "\\.CEL(\\.gz)?$", "\\.txt\\.gz$"),
  stringsAsFactors = FALSE)
for (k in seq_len(nrow(expect))) {
  f   <- expect$tar[k]
  tag <- basename(f)
  if (!file.exists(f)) { note_issue(paste(tag, "MISSING on this machine")); next }
  sz <- file.size(f)
  ok_sz <- isTRUE(!is.na(sz) && sz == expect$bytes[k])
  say(sprintf("  %s: %s bytes (expect %s) %s", tag,
              format(sz, big.mark = ","), format(expect$bytes[k], big.mark = ","),
              ifelse(ok_sz, "size-OK", "SIZE-MISMATCH")))
  if (!ok_sz) note_issue(paste(tag, "size mismatch"))
  m <- unname(tools::md5sum(f))
  say(sprintf("  %s MD5 %s %s", tag, m, ifelse(m == expect$md5[k], "MD5-OK", "MD5-MISMATCH")))
  if (m != expect$md5[k]) note_issue(paste(tag, "MD5 mismatch"))
  memb <- tryCatch(length(grep(expect$pat[k], untar(f, list = TRUE), ignore.case = TRUE)),
                   error = function(e) { note_issue(paste(tag, "untar list failed:", conditionMessage(e))); NA_integer_ })
  ok_m <- !is.na(memb) && memb == expect$members[k]
  say(sprintf("  %s members matching pattern: %s (expect %d) %s",
              tag, as.character(memb), expect$members[k], ifelse(ok_m, "count-OK", "COUNT-CHECK")))
  if (!ok_m && !is.na(memb)) note_issue(paste(tag, "member count", memb, "!=", expect$members[k]))
}

## ---- 3. Extraction-dir state (collision check) --------------------------------
say("== 3. Extraction-dir state ==")
dirs <- c("GSE51981", "GSE7305", "GSE25628", "GSE6364", "GSE11691", "GSE120103")
for (g in dirs) {
  d <- file.path("data/raw", g)
  if (!dir.exists(d)) { say(sprintf("  %s: no dir (clean untar target)", g)); next }
  fl  <- list.files(d)
  dat <- grep("\\.(CEL|txt)(\\.gz)?$", fl, ignore.case = TRUE, value = TRUE)
  prt <- grep("part_", fl, value = TRUE)
  say(sprintf("  %s: dir EXISTS - %d files (%d data-like, %d part files)",
              g, length(fl), length(dat), length(prt)))
  if (length(dat) == 0L) note_issue(paste0(g, ": dir exists with no data files -> v2 logic would SKIP untar (name collision)"))
  if (length(dat) > 0L)  say(sprintf("    note: %d data files already extracted", length(dat)))
}

## ---- 4. GSE7307 subset census + GSM mapping test ------------------------------
say("== 4. GSE7307 subset + frozen-list mapping ==")
sub_dir  <- "data/raw/GSE7307_subset"
sub_file <- "scripts/S1_GSE7307_subset.txt"
if (!file.exists(sub_file)) note_issue("scripts/S1_GSE7307_subset.txt missing from project tree")
if (!dir.exists(sub_dir)) {
  note_issue("data/raw/GSE7307_subset missing")
} else {
  cels <- list.files(sub_dir, pattern = "\\.CEL\\.gz$", ignore.case = TRUE)
  say(sprintf("  subset dir: %d CEL.gz files (expect 41)", length(cels)))
  if (length(cels) != 41L) note_issue(paste("GSE7307 subset has", length(cels), "files, expected 41"))
  if (file.exists(sub_file) && length(cels) > 0L) {
    frozen <- read.delim(sub_file, comment.char = "#", header = FALSE,
                         col.names = c("group", "gsm"), stringsAsFactors = FALSE)
    say(sprintf("  frozen list: %d GSMs (%d endometriosis + %d normal_endometrium)",
                nrow(frozen), sum(frozen$group == "endometriosis"), sum(frozen$group == "normal_endometrium")))
    bn      <- basename(cels)
    v2style <- sub("\\.CEL\\.gz$", "", bn, ignore.case = TRUE)
    rxstyle <- regmatches(bn, regexpr("GSM[0-9]+", bn, ignore.case = TRUE))
    m_v2 <- sum(v2style %in% frozen$gsm)
    m_rx <- sum(rxstyle %in% frozen$gsm)
    say(sprintf("  mapping test: v2-style (strip extension) = %d/%d matched ; GSM-regex = %d/%d matched",
                m_v2, length(bn), m_rx, length(bn)))
    if (m_rx == length(bn) && setequal(rxstyle, frozen$gsm)) {
      say("  GSM-regex mapping = EXACT 41/41 set match with frozen list")
    } else {
      note_issue("GSE7307 GSM-regex mapping incomplete vs frozen list")
    }
    if (m_v2 != length(bn)) note_issue("v2-style label matching would leave NA groups (F3 confirmed by measurement)")
    say("  example filenames: ", paste(utils::head(bn, 3), collapse = " | "))
  }
}

## ---- 5. Disk space on C: ------------------------------------------------------
say("== 5. Disk space (C:) ==")
tryCatch({
  di <- shell("dir C:\\", intern = TRUE)
  say("  ", tail(di, 1))
}, error = function(e) say("  disk check skipped: ", conditionMessage(e)))

## ---- 6. Summary -----------------------------------------------------------------
say("== PREFLIGHT SUMMARY ==")
if (length(issues) == 0L) {
  say("ALL CLEAN -- ready for Amendment S2_03-A route decision and smoke test")
} else {
  for (x in issues) say("  OPEN ISSUE: ", x)
}
say("Preflight done. Upload logs/S2_03_preflight.txt for verification (R2).")
