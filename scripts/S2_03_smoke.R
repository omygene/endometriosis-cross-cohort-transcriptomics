# ============================================================================
# S2_03_smoke.R -- mandatory smoke test before the S2_03 full run (rule B.4)
# Exercises every fragile path of S2_03_bulk_qc.R v3 on small real data:
#   GSE7305 (20 CEL, GPL570): untar -> read.celfiles -> rma -> AQM on the RMA
#   ExpressionSet (Amendment S2_03-A route) -> module/outlier extraction.
#   GSE120103 (2 Agilent files): read.maimages -> neqc (settles the genes$Status
#   question by direct measurement, R9).
# Writes ONLY under data/processed/tmp_S2_03_smoke/ (left in place for review).
# Run from project root:
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S2_03_smoke.R > logs\S2_03_smoke.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }
suppressPackageStartupMessages({ library(oligo); library(limma); library(arrayQualityMetrics) })

tmp <- "data/processed/tmp_S2_03_smoke"
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
smoke_fail <- FALSE

say("== SMOKE 1/5: untar + count (GSE7305, expect 20 CEL) ==")
exdir <- file.path(tmp, "GSE7305")
untar("data/raw/GSE7305_RAW.tar", exdir = exdir, tar = "internal")
cels <- list.files(exdir, pattern = "\\.CEL(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
say("  CELs found:", length(cels))
if (length(cels) != 20L) { say("  FAIL: expected 20"); smoke_fail <- TRUE }

say("== SMOKE 2/5: read.celfiles + rma ==")
raw <- tryCatch(oligo::read.celfiles(cels), error = function(e) { say("  READ ERROR:", conditionMessage(e)); NULL })
if (is.null(raw)) { smoke_fail <- TRUE } else {
  say("  raw class:", class(raw)[1], "| dims:", paste(dim(raw), collapse = " x "))
  say("  sampleNames head:", paste(utils::head(Biobase::sampleNames(raw), 3), collapse = " | "))
  eset <- tryCatch(oligo::rma(raw), error = function(e) { say("  RMA ERROR:", conditionMessage(e)); NULL })
  rm(raw); print(gc())
  if (is.null(eset)) { smoke_fail <- TRUE } else {
    say("  rma class:", class(eset)[1], "| dims:", paste(dim(eset), collapse = " x "))
    say("  exprs range:", paste(signif(range(Biobase::exprs(eset)), 4), collapse = " .. "))

    say("== SMOKE 3/5: arrayQualityMetrics on the RMA ExpressionSet (Amendment S2_03-A route) ==")
    aqm_err <- NULL
    aqm_ret <- tryCatch(
      arrayQualityMetrics(eset, outdir = file.path(tmp, "aqm"), force = TRUE, do.logtransform = FALSE),
      error = function(e) { aqm_err <<- e; NULL })
    if (!is.null(aqm_err)) {
      say("  AQM ERROR:", conditionMessage(aqm_err)); smoke_fail <- TRUE
    } else if (!file.exists(file.path(tmp, "aqm", "index.html"))) {
      say("  FAIL: AQM ran without error but index.html missing"); smoke_fail <- TRUE
    } else {
      say("  AQM report written OK (index.html present)")
      if (is.null(aqm_ret)) say("  return value: NULL -> v3 will use the prepdata fallback path")
      say("  AQM returned class:", paste(class(aqm_ret), collapse = "/"))
      if (is.list(aqm_ret)) say("  return names:", paste(names(aqm_ret), collapse = ", "))
      mods <- NULL
      if (is.list(aqm_ret) && "modules" %in% names(aqm_ret)) mods <- aqm_ret$modules
      if (is.null(mods) && is.list(aqm_ret) && length(aqm_ret) > 0 &&
          all(vapply(aqm_ret, function(z) is(z, "aqmReportModule"), logical(1)))) mods <- aqm_ret
      if (is.null(mods)) say("  return value carried NO modules -> v3 will use the prepdata fallback path")
      if (!is.null(mods)) {
        for (j in seq_along(mods)) {
          m  <- mods[[j]]
          od <- tryCatch(slot(m, "outliers"), error = function(e) NULL)
          sec <- tryCatch(as.character(slot(m, "section"))[1], error = function(e) "?")
          if (is(od, "outlierDetection"))
            say(sprintf("   module[%d] section=%s | n_stat=%d thr=%s n_out=%d | desc=%.70s",
                        j, sec, length(slot(od, "statistic")),
                        paste(signif(slot(od, "threshold"), 3), collapse = ","),
                        length(slot(od, "which")), slot(od, "description")))
          else
            say(sprintf("   module[%d] section=%s | no outlierDetection slot", j, sec))
        }
      }
    }

    say("== SMOKE 4/5: API formals dump (B.4 measurement) ==")
    for (fn in c("prepdata", "aqm.boxplot", "aqm.maplot", "aqm.heatmap", "outliers", "boxplotOutliers")) {
      f <- tryCatch(get(fn, envir = asNamespace("arrayQualityMetrics")), error = function(e) NULL)
      if (is.null(f)) say("  ", fn, ": NOT FOUND") else
        say("  ", fn, "(", paste(names(formals(f)), collapse = ", "), ")")
    }
  }
}

say("== SMOKE 5/5: Agilent leg (2 files from GSE120103) ==")
mem  <- untar("data/raw/GSE120103_RAW.tar", list = TRUE)
txts <- grep("\\.txt\\.gz$", mem, value = TRUE)[1:2]
say("  members chosen:", paste(txts, collapse = " | "))
agdir <- file.path(tmp, "GSE120103"); dir.create(agdir, showWarnings = FALSE)
untar("data/raw/GSE120103_RAW.tar", files = txts, exdir = agdir, tar = "internal")
ag <- list.files(agdir, pattern = "\\.txt\\.gz$", full.names = TRUE, recursive = TRUE)
rg <- tryCatch(limma::read.maimages(ag, source = "agilent", green.only = TRUE),
               error = function(e) { say("  read.maimages ERROR:", conditionMessage(e)); NULL })
if (is.null(rg)) { smoke_fail <- TRUE } else {
  say("  RG dims:", nrow(rg$G), "x", ncol(rg$G))
  say("  genes columns:", paste(colnames(rg$genes), collapse = ", "))
  say("  genes$Status present:", "Status" %in% colnames(rg$genes))
  y <- tryCatch(limma::neqc(rg), error = function(e) { say("  NEQC ERROR:", conditionMessage(e)); NULL })
  if (is.null(y)) { smoke_fail <- TRUE; say("  NEQC FAILED -> escalate (R9)") }
  else say("  neqc OK:", nrow(y$E), "x", ncol(y$E))
}

print(gc())
say(if (smoke_fail) "SMOKE DONE WITH FAILURES -- see lines above" else "SMOKE OK")
