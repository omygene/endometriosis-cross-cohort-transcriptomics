# ============================================================================
# S2_03_GSE120103_fix.R -- completes cohort GSE120103 after the v3.1 run
# stopped on ONE file (measured, logs/S2_03_GSE120103.txt 2026-09-17):
#   GSM3393522_..._GE1-v5_95_Feb07_1_1.txt.gz -> different chip print design
#   (GE1-v5_95_Feb07) vs the other 35 files (GE1_107_Sep09); read.maimages:
#   "fewer rows than files previously read" + "embedded nul(s) found in input".
# Decision (lead, documented in analysis_log.md): exclude that one array.
# Mirrors S2_03_bulk_qc.R v3.1 Agilent leg step-for-step (same calls, same
# outputs). GSE120103 is neqc-only by design (S2_03-A) -> no AQM, no outlier CSV.
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }

suppressPackageStartupMessages({ library(limma); library(ggplot2) })
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

cohort <- "GSE120103"
exdir  <- file.path("data/raw", cohort)
data_files <- list.files(exdir, pattern = "\\.txt\\.gz$", full.names = TRUE, ignore.case = TRUE)

## 1. documented exclusion of the mismatched-design array
drop_me <- data_files[grepl("GE1-v5_95_Feb07", basename(data_files))]
if (length(drop_me) != 1)
  stop("FATAL: expected exactly 1 GE1-v5 file, found ", length(drop_me), " -- re-check before running")
data_files <- setdiff(data_files, drop_me)
if (length(data_files) != 35)
  stop("FATAL: expected 35 files after exclusion, found ", length(data_files))
say("excluded (documented):", basename(drop_me))
say("data files:", length(data_files), "(36 minus 1 documented exclusion)")

## 2. read + Status marshalling (identical to v3.1 / Addendum A.1)
say("read.maimages (source=agilent, green.only=TRUE) ...")
rg <- limma::read.maimages(data_files, source = "agilent", green.only = TRUE)
colnames(rg) <- basename(data_files)
say("RGList dim:", paste(dim(rg), collapse = " x "))
if (!"Status" %in% colnames(rg$genes)) {
  ct <- rg$genes$ControlType
  if (is.null(ct)) stop("FATAL: genes$ControlType absent -- escalate (R9)")
  tab_ct <- table(ct, useNA = "always")
  say("ControlType table:", paste(paste(names(tab_ct), tab_ct, sep = "="), collapse = " / "))
  n_reg <- sum(ct == 0, na.rm = TRUE); n_neg <- sum(ct == -1, na.rm = TRUE)
  if (n_reg == 0 || n_neg == 0)
    stop("FATAL: ControlType coding not 0/-1 (0:", n_reg, " -1:", n_neg, ") -- do NOT marshal; escalate (R9)")
  rg$genes$Status <- ifelse(ct == 0, "regular", ifelse(ct == -1, "negative", "other"))
  say("Status marshalled:", n_reg, "regular /", n_neg, "negative /",
      sum(rg$genes$Status == "other"), "other")
}

## 3. neqc + save (identical outputs to v3.1)
eset <- limma::neqc(rg)
rm(rg)
saveRDS(eset, "data/processed/bulk_GSE120103_neqc.rds")
mat <- eset$E
say("normalized matrix:", nrow(mat), "x", ncol(mat),
    "| E range:", round(min(mat), 3), "..", round(max(mat), 3))

## 4. PCA figure (identical to v3.1)
pc  <- prcomp(t(mat), scale. = TRUE)
pct <- round(100 * summary(pc)$importance[2, 1:2], 1)
df  <- data.frame(pc$x[, 1:2], sample = colnames(mat))
p <- ggplot(df, aes(PC1, PC2, label = sample)) +
  geom_point(size = 2) + geom_text(size = 2, vjust = -0.6) +
  labs(title = "GSE120103 PCA (35 arrays; 1 documented exclusion)",
       x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))
ggsave("figures/S2_bulk_GSE120103_pca.png", p, width = 9, height = 7, dpi = 150)
say("wrote figures/S2_bulk_GSE120103_pca.png")
print(gc())
say("COHORT GSE120103 DONE")