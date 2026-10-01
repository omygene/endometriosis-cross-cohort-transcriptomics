# ============================================================================
# S2_02_journal_figures.R — cosmetic publication pass over the FINAL S2_02 rds.
# Requested in writing by the lead (2026-09-17): "عاوزة الصور جاهزة للرفع للمجلة".
#
# DISPLAY-ONLY script (R4-safe): reads data/processed/S2_02_integrated.rds,
# writes NEW files under figures/journal/, never touches any existing results
# file, never saves the object. Zero parameter/data change.
#
# Cosmetic changes vs the QC figures (all approved/registered):
#  1. No internal stage IDs in titles; no plot titles at all (journal style —
#     the caption carries the description).
#  2. Axis labels "UMAP 1"/"UMAP 2".
#  3. Condition harmonization FOR DISPLAY ONLY (lead's written decision
#     2026-09-17): Ctrl (GSE179640) / Control (GSE183837) / N (GSE214411) are
#     biologically "healthy, no disease" from 3 cohorts -> one display
#     category "Healthy control" (one color; caption states the 3 cohorts).
#     The metadata column `condition` is NOT modified or saved.
#  4. Colorblind-safe palettes (Okabe-Ito based) for condition/dataset/lineage.
#  5. Cluster labels boxed (label.box verified in locked Seurat 5.5.1 source,
#     visualization.R DimPlot signature) to reduce overlap artifacts.
# Output: figures/journal/S2_02_umap_{clusters,dataset,condition,lineage}.{png,pdf}
#         + figures/journal/S2_02_fig_integrated_4panel.{png,pdf} (A-D tags).
# Run:  Rscript scripts\S2_02_journal_figures.R > logs\S2_02_journal_figures.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}
req <- c("Seurat", "ggplot2", "patchwork", "ggrastr", "ggrepel")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("MISSING PACKAGES: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})

obj <- readRDS("data/processed/S2_02_integrated.rds")
stopifnot(ncol(obj) == 233459L)
dir.create(file.path("figures", "journal"), showWarnings = FALSE, recursive = TRUE)

# ---- display-only condition harmonization -----------------------------------
cond_map <- c(Ctrl = "Healthy control", Control = "Healthy control",
              N = "Healthy control", EuE = "EuE", EcP = "EcP", EcPA = "EcPA",
              EcO = "EcO", EOR = "EOR", EMS = "EMS", RIF = "RIF")
cond_levels <- c("Healthy control", "EuE", "EcP", "EcPA", "EcO", "EOR",
                 "EMS", "RIF")
cond_pub <- unname(cond_map[as.character(obj$condition)])
if (any(is.na(cond_pub))) {
  stop("unknown condition label(s): ",
       paste(setdiff(unique(as.character(obj$condition)), names(cond_map)),
             collapse = ", "), " — investigate (R6)")
}
obj$condition_pub <- factor(cond_pub, levels = cond_levels)
say("display conditions: ", paste(levels(obj$condition_pub), collapse = ", "))

# ---- palettes (colorblind-safe; Okabe-Ito based) -----------------------------
COND_COLORS_PUB <- c("Healthy control" = "#009E73", "EuE" = "#0072B2",
                     "EcP" = "#E69F00", "EcPA" = "#CC79A7", "EcO" = "#56B4E9",
                     "EOR" = "#7E6148", "EMS" = "#D55E00", "RIF" = "#999999")
DATASET_COLORS_PUB <- c(GSE179640 = "#D55E00", GSE183837 = "#0072B2",
                        GSE214411 = "#009E73")
LINEAGE_COLORS_PUB <- c(Stromal = "#E69F00", Epithelial = "#D55E00",
                        Endothelial = "#009E73", Perivascular = "#56B4E9",
                        NKT = "#0072B2", Myeloid = "#CC79A7",
                        B_plasma = "#999999", Mast = "#000000")

# ---- base plot ---------------------------------------------------------------
pub_base <- function(...) {
  DimPlot(obj, reduction = "umap", raster = TRUE, raster.dpi = c(300, 300),
          ...) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_classic(base_size = 13) +
    theme(axis.text = element_text(color = "black"),
          legend.title = element_text(face = "bold"))
}

figs <- list(
  clusters = pub_base(label = TRUE, repel = TRUE, label.box = TRUE,
                      label.size = 3.2) + NoLegend(),
  dataset  = pub_base(group.by = "dataset") +
    scale_color_manual(values = DATASET_COLORS_PUB) +
    labs(color = "Dataset"),
  condition = pub_base(group.by = "condition_pub") +
    scale_color_manual(values = COND_COLORS_PUB) +
    labs(color = "Condition"),
  lineage  = pub_base(group.by = "lineage") +
    scale_color_manual(values = LINEAGE_COLORS_PUB) +
    labs(color = "Lineage")
)

for (nm in names(figs)) {
  ggsave(file.path("figures", "journal", paste0("S2_02_umap_", nm, ".png")),
         figs[[nm]], width = 7.5, height = 6.5, dpi = 300, bg = "white")
  ggsave(file.path("figures", "journal", paste0("S2_02_umap_", nm, ".pdf")),
         figs[[nm]], width = 7.5, height = 6.5, bg = "white",
         device = cairo_pdf)
  say("figure: ", nm)
}

# ---- combined 4-panel figure (A-D) -------------------------------------------
combo <- (figs$clusters | figs$lineage) / (figs$condition | figs$dataset) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))
ggsave(file.path("figures", "journal", "S2_02_fig_integrated_4panel.png"),
       combo, width = 13, height = 11.5, dpi = 300, bg = "white")
ggsave(file.path("figures", "journal", "S2_02_fig_integrated_4panel.pdf"),
       combo, width = 13, height = 11.5, bg = "white", device = cairo_pdf)
say("figure: 4-panel combo")
say("JOURNAL FIGURES DONE — nothing else written, object not saved")
