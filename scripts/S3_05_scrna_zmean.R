# ============================================================================
# S3_05_scrna_zmean.R -- SenMayo z-mean scores for scRNA (locked criterion 3)
# Definition (frozen): per cell, z-mean = mean over SenMayo genes (present in
# the object) of gene-wise z-scores computed on log1p-normalized expression
# across ALL cells of the atlas. Summarized per dataset x lineage x sample,
# then a sample-level limma comparison (with N1 arms for GSE214411).
# Scores are descriptive inputs to S4 (random-effects there).
# v2.4 (2026-09-21): per-dataset expression sources (S3_00d/e/f measured:
# atlas is a 28-gene panel). The frozen definition is preserved EXACTLY:
# per dataset we take the SenMayo rows of the log1p-normalized data
# (prep_rna_norm, same normalization), POOL the three blocks cell-for-cell
# (all 233,459 cells, keyed by cell_key), then compute gene-wise z-scores
# across ALL pooled cells. Genes absent from one dataset's platform leave
# that block as NA; NA z-scores are set to 0 (the original convention),
# so every cell still gets a z-mean over the genes it can express.
# v2.4.1 (2026-09-21): model.matrix now receives data = d. LATENT BUG in the
# original S3_05 (carried verbatim into v2.4): the call lacked data=, so R
# looked for a standalone `grp` object that never existed (grp is a column
# of summ). Unreachable until today -- every earlier run died upstream (old:
# SenMayo gate at 1/124; new: data pipeline). Measured halt on the lead
# machine: object 'grp' not found after the summary CSV was written.
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_05_scrna_zmean.R > logs\S3_05_scrna_zmean.txt 2>&1
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(limma); library(Seurat) })
source("scripts/S3_common.R")
dir.create("results", showWarnings = FALSE)

require_config(c("n1_sample", "group_col", "lineage_col"))
say("== S3_05 z-mean (SenMayo; per-dataset sources, v2.4) ==")
gm <- find_senmayo_gmt()
if (is.na(gm)) stop("FATAL: SenMayo GMT not found (see S3_00)")
sen <- unique(unlist(lapply(read_gmt(gm), `[[`, "genes")))
say("SenMayo symbols in GMT:", length(sen))

md <- atlas_annotations()
datasets <- intersect(sort(unique(as.character(md$dataset))), names(CFG$per_dataset))
if (!length(datasets)) stop("FATAL: no atlas dataset with a registered per_dataset source (R9)")

blocks <- list()   # per dataset: log1p-normalized SenMayo rows (genes x cells)
keys <- character(0)
for (ds in datasets) {
  per <- readRDS(CFG$per_dataset[[ds]])
  per <- prep_rna_norm(per)
  ann <- join_annotations(per, md, ds)
  dat <- SeuratObject::GetAssayData(per, assay = "RNA", layer = "data")
  use <- intersect(sen, rownames(dat))
  say(ds, "SenMayo genes present:", length(use), "of", length(sen))
  blocks[[ds]] <- as.matrix(dat[use, , drop = FALSE])
  keys <- c(keys, ann$cell_key)
  rm(per, dat, ann); gc()
}

## frozen definition: gene-wise z across ALL pooled cells
all_genes <- sort(unique(unlist(lapply(blocks, rownames))))
P <- matrix(NA_real_, length(all_genes), length(keys),
            dimnames = list(all_genes, keys))
col0 <- 1L
for (ds in datasets) {
  b <- blocks[[ds]]
  P[rownames(b), col0:(col0 + ncol(b) - 1)] <- b
  col0 <- col0 + ncol(b)
}
rm(blocks); gc()
say("pooled SenMayo matrix:", nrow(P), "genes x", ncol(P), "cells (atlas = 233459)")
if (ncol(P) != 233459) say("  NOTE: pooled cell count differs from atlas census -- reported, not filtered")

z <- scale(t(P))          # cells x genes, gene-wise z across ALL pooled cells
z[is.na(z)] <- 0          # gene absent from a dataset's platform -> z = 0
zm <- rowMeans(z)
rm(P, z); gc()

md$senmayo_zmean <- zm[match(rownames(md), names(zm))]
if (any(is.na(md$senmayo_zmean)))
  stop("FATAL: atlas cells without z-mean after pooled join (R9)")
md$grp <- map_cell_group(as.character(md$dataset), as.character(md[[CFG$group_col]]))

summ <- aggregate(senmayo_zmean ~ dataset + lineage + sample + grp,
                  data = data.frame(dataset = md$dataset, lineage = md[, CFG$lineage_col],
                                    sample = md$sample, grp = md$grp, senmayo_zmean = md$senmayo_zmean),
                  FUN = mean)
write.csv(summ, "results/S3_scrna_zmean_summary.csv", row.names = FALSE)
say("wrote results/S3_scrna_zmean_summary.csv --", nrow(summ), "rows")

## sample-level comparison per dataset x lineage (with N1 arms)
res <- list()
for (ds in sort(unique(summ$dataset))) {
  for (lin in sort(unique(summ$lineage[summ$dataset == ds]))) {
    for (arm in if (ds == "GSE214411") c("A", "B") else "A") {
      d <- summ[summ$dataset == ds & summ$lineage == lin, ]
      if (arm == "B") d <- d[d$sample != CFG$n1_sample, ]
      if (length(unique(d$grp)) < 2 || min(table(d$grp)) < 2) next
      design <- model.matrix(~ factor(grp, levels = c("REF", "TEST")), data = d)
      fit <- eBayes(lmFit(d$senmayo_zmean, design))
      tt <- topTable(fit, coef = 2, number = 1)
      res[[paste(ds, lin, arm, sep = "_")]] <- data.frame(
        dataset = ds, lineage = lin, arm = arm, n_test = sum(d$grp == "TEST"),
        n_ref = sum(d$grp == "REF"), logFC = tt$logFC[1], P = tt$P.Value[1],
        FDR = tt$adj.P.Val[1], stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, res)
write.csv(out, "results/S3_scrna_zmean_tests.csv", row.names = FALSE)
say("wrote results/S3_scrna_zmean_tests.csv --", nrow(out), "tests")
print(out)
say("S3_05 DONE")
