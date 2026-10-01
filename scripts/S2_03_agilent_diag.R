# ============================================================================
# S2_03_agilent_diag.R  v1  -- R9 diagnostic for the measured Agilent neqc failure
#
# MEASURED failure (logs/S2_03_smoke.txt, leg 5, 2026-09-18):
#   read.maimages OK, genes$Status absent -> limma::neqc ERROR:
#   "Detection p values not found in the data."
#
# Locked-binary evidence (limma 3.68.5, reference/locked_sources/limma):
#   neqc -> nec; nec routes on genes$Status:
#     present -> normexp.fit.control   (matches tolower(Status) against
#                                       "regular"/"negative"; hard stops:
#                                       "No regular probes found" /
#                                       "Fewer than two negative control probes found")
#     absent  -> normexp.fit.detection.p -> requires a detection p-value
#                                       column these FE files do not carry
#                                       -> the exact measured stop.
#
# HYPOTHESIS UNDER MEASUREMENT (not assumed): Agilent ControlType coding is
#   0 = regular probe, -1 = negative control. This script MEASURES
#   table(ControlType) on 2 real files; the marshalling is applied ONLY if
#   both codes are present; otherwise it stops before touching anything.
#
# Writes ONLY under data/processed/tmp_S2_03_diag/ (R4: no result files).
# R10 RAM estimate: 2 Agilent FE files (~60k probes x 2 arrays) -> < 0.5 GB.
#
# Run (from C:\endometriosis_immunoediting):
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S2_03_agilent_diag.R > logs\S2_03_agilent_diag.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }
suppressPackageStartupMessages(library(limma))

tmpd <- "data/processed/tmp_S2_03_diag/GSE120103"
dir.create(tmpd, recursive = TRUE, showWarnings = FALSE)
tarf <- "data/raw/GSE120103_RAW.tar"
if (!file.exists(tarf)) stop("FATAL: ", tarf, " missing")

say("== DIAG 1/3: untar 2 members + read.maimages ==")
members <- untar(tarf, list = TRUE)
data_members <- grep("\\.txt\\.gz$", members, value = TRUE, ignore.case = TRUE)
say("  txt.gz members in tar:", length(data_members), "(expected 36)")
if (length(data_members) != 36L) stop("FATAL: tar member count mismatch: ", length(data_members))
untar(tarf, files = data_members[1:2], exdir = tmpd, tar = "internal")
f <- list.files(tmpd, pattern = "\\.txt\\.gz$", full.names = TRUE,
                recursive = TRUE, ignore.case = TRUE)
if (length(f) != 2L) stop("FATAL: expected 2 extracted files, found ", length(f))
say("  files:", paste(basename(f), collapse = " | "))
rg <- read.maimages(f, source = "agilent", green.only = TRUE)
colnames(rg) <- basename(f)
say("  class:", class(rg), "| dim:", paste(dim(rg), collapse = " x "))
say("  RGList components:", paste(names(rg), collapse = ", "))
say("  genes columns:", paste(colnames(rg$genes), collapse = ", "))
say("  genes$Status present:", "Status" %in% colnames(rg$genes))

say("== DIAG 2/3: MEASURE ControlType coding ==")
ct <- rg$genes$ControlType
if (is.null(ct)) stop("FATAL: genes$ControlType absent -- hypothesis wrong; escalate with this log")
print(table(ct, useNA = "always"))
n0   <- sum(ct == 0,  na.rm = TRUE)
nm1  <- sum(ct == -1, na.rm = TRUE)
noth <- sum(is.na(ct) | !(ct %in% c(0, -1)))
say("  ControlType==0 (regular?):", n0, "| ==-1 (negative?):", nm1, "| other/NA:", noth)
if (n0 == 0 || nm1 == 0)
  stop("FATAL: ControlType coding is NOT 0/-1 as hypothesized (code 0: ", n0,
       ", code -1: ", nm1, ") -- do NOT marshal Status; escalate with this log")
say("  coding 0/-1 CONFIRMED on these files (regular:", n0, ", negative:", nm1, ")")

say("== DIAG 3/3: apply Status marshalling + neqc ==")
rg$genes$Status <- ifelse(ct == 0, "regular", ifelse(ct == -1, "negative", "other"))
print(table(rg$genes$Status, useNA = "always"))
eset <- neqc(rg)
say("  neqc OK | class:", class(eset), "| dim:", paste(dim(eset), collapse = " x "))
say("  E range:", paste(round(range(eset$E), 3), collapse = " .. "))
if ("Status" %in% colnames(eset$genes)) {
  say("  output retains Status (controls kept in object; filtering is an S2_04 decision):")
  print(table(eset$genes$Status, useNA = "always"))
}
saveRDS(eset, "data/processed/tmp_S2_03_diag/neqc_diag.rds")
print(gc())
say("DIAG OK -- Status marshalling resolves the measured neqc failure")
