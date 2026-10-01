# ============================================================================
# S2_00b_entry_checks_fix.R — corrective re-run of S2_00 items
# Documents and fixes, per Amendment S2-A (logs/S2_preregistered_design.md):
#   E1-GSE179640: reconcile tar member NAMES vs official GEO filelist (63 files)
#   E2: re-list GSE7305 CELs with R's internal tar reader (bsdtar failed, status 1)
#   E3: re-run GSE25628 metadata parse (requires GEOquery -> run AFTER S0 lock)
# Prerequisite: copy logs/S2_GSE179640_GEO_filelist.txt (shared storage) into
# your local logs/ folder before running.
# Run from project root. Output: logs/S2_entry_checks_fix.txt
# ============================================================================
set.seed(42)
dir.create("logs", showWarnings = FALSE)
sink("logs/S2_entry_checks_fix.txt", split = TRUE)
cat("S2 entry checks — corrective re-run —", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- E1-GSE179640: member names vs official GEO filelist --------------------
cat("== E1-GSE179640: tar members vs GEO official filelist ==\n")
members <- tryCatch(basename(untar("data/raw/GSE179640_RAW.tar", list = TRUE,
                                   tar = "internal")),
                    error = function(e) { cat("TAR READ FAILED:", conditionMessage(e), "\n"); character(0) })
cat("local tar members:", length(members), "\n")

geo_file <- "logs/S2_GSE179640_GEO_filelist.txt"
if (file.exists(geo_file)) {
  # header line starts with '#', so skip it explicitly (v2 fix: comment.char
  # had eaten the header and silently produced 0 official files)
  fl <- read.delim(geo_file, skip = 1, header = FALSE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE,
                   col.names = c("kind", "Name", "Time", "Size", "Type"))
  geo_names <- fl$Name[fl$kind == "File"]
  cat("GEO official files:", length(geo_names), "\n")
  missing_local <- setdiff(geo_names, members)
  extra_local   <- setdiff(members, geo_names)
  cat("in GEO but NOT in local tar:", length(missing_local), "\n")
  if (length(missing_local)) cat(paste("  -", missing_local, collapse = "\n"), "\n")
  cat("in local tar but NOT in GEO:", length(extra_local), "\n")
  if (length(extra_local)) cat(paste("  -", extra_local, collapse = "\n"), "\n")
  verdict <- if (length(members) == 63 && length(missing_local) == 0 &&
                 length(extra_local) == 0) "PASS" else "CHECK"
  cat("E1-GSE179640 verdict:", verdict, "(criterion: 63/63 exact name match)\n")
} else {
  cat("GEO filelist copy missing in logs/ — copy S2_GSE179640_GEO_filelist.txt from shared storage first.\n")
}

# ---- E2 redo: GSE7305 vs GSE7307 overlap — GSM-level via official GEO filelists
# (authoritative evidence saved 2026-09-12; local tar not required for this check.
#  Root cause of the earlier failure: GSE7305_RAW.tar was portal-resident, not yet
#  copied to the local machine — it IS required locally before S2_03 bulk QC.)
cat("\n== E2 redo: GSE7305 vs GSE7307 overlap (GSM-level, official GEO records) ==\n")
read_gsm <- function(path) {
  fl <- read.delim(path, skip = 1, header = FALSE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE,
                   col.names = c("kind", "Name", "Time", "Size", "Type"))
  sub("\\..*$", "", fl$Name[fl$kind == "File"])
}
if (file.exists("logs/S2_GSE7305_GEO_filelist.txt") &&
    file.exists("logs/S2_GSE7307_GEO_filelist.txt")) {
  g7305 <- read_gsm("logs/S2_GSE7305_GEO_filelist.txt")
  g7307 <- read_gsm("logs/S2_GSE7307_GEO_filelist.txt")
  subset41 <- unique(unlist(regmatches(
    readLines("scripts/S1_GSE7307_subset.txt"),
    gregexpr("GSM[0-9]+", readLines("scripts/S1_GSE7307_subset.txt")))))
  cat("GSE7305 official GSMs:", length(g7305), "(manifest: 20)\n")
  cat("GSE7307 official GSMs:", length(g7307), "\n")
  cat("frozen 41-subset all within official GSE7307:", all(subset41 %in% g7307), "\n")
  common <- intersect(g7305, subset41)
  cat("GSM accessions shared by GSE7305 and the frozen subset:", length(common), "\n")
  if (length(common)) cat(paste("  -", common, collapse = "\n"), "\n")
} else {
  cat("evidence filelists missing — copy S2_GSE7305_GEO_filelist.txt and\n")
  cat("S2_GSE7307_GEO_filelist.txt from shared storage into local logs/.\n")
}
cat("NOTE stands: distinct GSM accessions do not exclude shared donors (same tissue\n")
cat("bank, Neurocrine submissions 2007). Pre-registered action unchanged: in S4, any\n")
cat("score-stage association relying on BOTH GSE7305 and GSE7307 is re-run without\n")
cat("GSE7307 as a sensitivity analysis.\n")

# ---- E3 redo: GSE25628 metadata parse (needs GEOquery from S0) --------------
cat("\n== E3 redo: GSE25628 metadata (22 samples, GEO official) ==\n")
if (requireNamespace("GEOquery", quietly = TRUE)) {
  g <- tryCatch(GEOquery::getGEO("GSE25628", GSEMatrix = TRUE, getGPL = FALSE),
                error = function(e) NULL)
  if (!is.null(g)) {
    pd <- Biobase::pData(g[[1]])
    cols <- grep("characteristics|source_name|title", colnames(pd), value = TRUE)
    print(pd[, c("geo_accession", utils::head(cols, 4)), drop = FALSE])
  } else cat("GEOquery fetch failed — record manually from GEO page at S2 execution.\n")
} else cat("GEOquery not installed yet — run scripts/S0_environment_lock.R first, then re-run this script.\n")

sink()
cat("\nCorrective entry checks written to logs/S2_entry_checks_fix.txt\n")
