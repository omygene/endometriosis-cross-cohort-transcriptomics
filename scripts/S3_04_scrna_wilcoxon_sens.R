# ============================================================================
# S3_04_scrna_wilcoxon_sens.R -- SENSITIVITY ONLY (locked criterion 2)
# Cell-level Wilcoxon per dataset x lineage on log-normalized data.
# Pseudoreplication caveat: cell-level tests inflate n; that is WHY this leg
# is sensitivity-only and the pseudobulk (S3_03) is primary. Every output row
# carries sensitivity = TRUE.
# v2.5 (2026-09-21): crash-safety + presto gate.
#   MEASURED on the lead machine 2026-09-21: without presto, FindMarkers
#   (wilcox) on GSE179640 (109,398 cells x 29,448 genes) ran >9 h with no
#   log line and the R process was killed silently (twice). Fixes:
#   (a) fail-loud gate: dataset > 50k cells requires presto -- the exact
#       engine Seurat itself recommends in its own warning message; the
#       Wilcoxon test is IDENTICAL (presto is Seurat's drop-in fast path),
#       this is a compute-engine change, NOT a criterion change;
#   (b) per-lineage progress line with timestamp before each FindMarkers;
#   (c) each CSV is written the moment its test completes (a crash never
#       loses more than the single running test).
# Usage:  "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S3_04_scrna_wilcoxon_sens.R > logs\S3_04_scrna_wilcoxon_sens.txt 2>&1
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(Seurat) })
source("scripts/S3_common.R")
dir.create("results/S3_scrna_wilcoxon_sens", recursive = TRUE, showWarnings = FALSE)

require_config(c("group_col", "lineage_col"))
say("== S3_04 Wilcoxon SENSITIVITY (not primary; per-dataset sources, v2.5) ==")
has_presto <- requireNamespace("presto", quietly = TRUE)
say("presto engine:", if (has_presto) "ACTIVE (fast wilcox, Seurat drop-in)" else
  "NOT INSTALLED -- datasets > 50k cells will stop FATAL (see gate below)")
md <- atlas_annotations()

datasets <- intersect(sort(unique(as.character(md$dataset))), names(CFG$per_dataset))
if (!length(datasets)) stop("FATAL: no atlas dataset with a registered per_dataset source (R9)")
for (ds in datasets) {
  say("dataset:", ds)
  per <- readRDS(CFG$per_dataset[[ds]])
  per <- prep_rna_norm(per)
  ann <- join_annotations(per, md, ds)
  if (!has_presto && ncol(per) > 50000)
    stop("FATAL: ", ds, " has ", ncol(per), " cells and presto is not installed.\n",
         "  Without presto, FindMarkers(wilcox) on >50k cells x full gene space was\n",
         "  MEASURED to run >9 h and die silently on this machine (2026-09-21).\n",
         "  Install presto (binary, ~2 min, no compiler):\n",
         "    Rscript -e \"install.packages('presto', repos=c('https://immunogenomics.r-universe.dev','https://cloud.r-project.org'))\"\n",
         "  The Wilcoxon test itself is UNCHANGED -- presto is the engine Seurat calls.\n",
         "  Then re-run S3_04 only (R3).")
  grp <- map_cell_group(rep(ds, nrow(ann)), as.character(ann$condition))
  cell_lin <- as.character(ann$lineage)

  for (lin in sort(unique(cell_lin))) {
    sel <- cell_lin == lin
    g <- grp[sel]
    if (length(unique(g)) < 2) { say("  skip", ds, lin, "(single group)"); next }
    say("  [", format(Sys.time(), "%H:%M:%S"), "] FindMarkers start:", ds, lin,
        "| cells:", sum(sel), "| genes:", nrow(per))
    so <- subset(per, cells = colnames(per)[sel])
    DefaultAssay(so) <- "RNA"
    Idents(so) <- factor(g, levels = c("REF", "TEST"))
    sub <- FindMarkers(so, ident.1 = "TEST", ident.2 = "REF", test.use = "wilcox",
                       logfc.threshold = 0, min.pct = 0, only.pos = FALSE, verbose = FALSE)
    sub$gene <- rownames(sub); sub$dataset <- ds; sub$lineage <- lin
    sub$sensitivity <- TRUE
    out <- file.path("results/S3_scrna_wilcoxon_sens",
                     paste0("S3_wx_", ds, "_", lin, ".csv"))
    write.csv(sub, out, row.names = FALSE)   # immediate: crash loses at most this test
    say("  wrote", basename(out), "| genes:", nrow(sub),
        "| nominal P<0.05:", sum(sub$p_val < 0.05))
    rm(so, sub); gc()
  }
  rm(per, ann); gc()
}
print(gc())
say("S3_04 DONE")
