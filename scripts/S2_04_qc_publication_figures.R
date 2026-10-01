# ============================================================================
# S2_04_qc_publication_figures.R — publication-grade QC figures (v4.1)
# v2: vertical-stacked violin figure per dataset (3 metrics stacked, one
#     shared x-axis, short labels colored by condition, thin boxplot overlay,
#     threshold lines).
# v3 (2026-09-14): GSE183837 branch in parse_sample (Control/RIF coloring,
#     mapping frozen vs GEO per Amendment S2-B/B2); COND_COLORS extended.
# v4 (2026-09-14) — fixes from lead's v3-run report (documented, R4: only
#     this script changed; no results/CSV touched; no thresholds touched, R5):
#   F1. CRASH "factor level [2] is duplicated" in violin_stack AND the ECDF
#       section: parse_sample built factor levels from
#       df$short[order(df$cond, df$short)] — at CELL level (millions of
#       rows) the levels vector is massively duplicated and factor()
#       refuses it. Summaries survived because they parse the per-sample
#       CSV (unique rows). Fix: levels = unique(...) — lead's own proposed
#       fix, adopted verbatim in effect. One fix covers violins + ECDF.
#   F2. Duplicated legend (points for colour + squares for fill, same
#       labels): summary panels now use shape-21 FILLED points (fill = cond)
#       so panels A–D share ONE fill scale -> patchwork collects a single
#       legend.
#   F3. GSE179640 x-labels in panels B/C/D crowded (33 samples): x text
#       size now scales down when > 20 samples.
#   F4. sprintf("C  Mito%% vs threshold (%d%%)", mt_thr, mt_thr) had an
#       extra argument (warning only) -> single argument.
#   F5. "→" in panel A title failed glyph conversion in the PDF device ->
#       PDF outputs now use device = cairo_pdf (keeps the arrow; PNG
#       unchanged).
# v4.1 (2026-09-14) — post-run figure review (documented; script only, R4):
#   F6. F2 was incomplete: panel D's geom_col still contributed its own
#       fill legend (square keys), so patchwork's guides="collect" kept TWO
#       legends (circle keys from panels A–C + square keys from D) with
#       identical labels/colors. Fix: show.legend = FALSE on geom_col —
#       bar colors unchanged; the single shape-21 fill legend now speaks
#       for all four panels. No other change; no data/threshold touched.
# Inputs : results/S2_qc_<GSE>.csv (summary level)
#          data/processed/<GSE>_filtered.rds (cell level)  [ECDF + violins]
# Outputs: figures/S2_qc_pub_<GSE>_summary.{png,pdf}
#          figures/S2_qc_pub_<GSE>_violins.{png,pdf}
#          figures/S2_qc_pub_ecdf.{png,pdf}
# Note: panel D (doublets) is active for all datasets (dbl_pct present).
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})
dir.create("figures", showWarnings = FALSE)

COND_COLORS <- c("Ctrl" = "#4DBBD5", "N" = "#4DBBD5", "Control" = "#4DBBD5",
                 "EuE" = "#E64B35", "EMS" = "#D1495B", "EcP" = "#00A087",
                 "EcO" = "#3C5488", "EcPA" = "#F39B7F", "EOR" = "#8491B4",
                 "RIF" = "#7E6148", "sample" = "#757575")
MT_THR <- c(GSE179640 = 20, GSE214411 = 10, GSE183837 = 20)   # design §2

parse_sample <- function(df, gse) {
  if (gse == "GSE179640") {
    df$short <- sub("^GSM\\d+_(.+?)_filtered.*$", "\\1", df$sample, perl = TRUE)
    df$cond <- dplyr::case_when(
      grepl("_Ctrl", df$sample) ~ "Ctrl", grepl("_EuE",  df$sample) ~ "EuE",
      grepl("_EcPA", df$sample) ~ "EcPA", grepl("_EcP",  df$sample) ~ "EcP",
      grepl("_EcO",  df$sample) ~ "EcO",  grepl("EOR",   df$sample) ~ "EOR",
      TRUE ~ "other")
  } else if (gse == "GSE214411") {
    df$short <- sub("^GSM\\d+_", "", df$sample)
    df$cond  <- ifelse(grepl("EMS", df$sample), "EMS", "N")
  } else if (gse == "GSE183837") {
    # v3: mapping frozen vs GEO (2026-09-14, Amendment S2-B/B2)
    df$short <- sub("^GSM\\d+_", "", df$sample)
    df$cond  <- ifelse(grepl("RIF", df$sample, ignore.case = TRUE), "RIF", "Control")
  } else {
    df$short <- df$sample; df$cond <- "sample"
  }
  # v4 F1: unique() — at cell level the ordered vector is duplicated millions
  # of times; factor(levels = ...) must be duplicate-free.
  df$short <- factor(df$short, levels = unique(df$short[order(df$cond, df$short)]))
  df
}

theme_pub <- theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 6.5),
        axis.text.x = element_text(size = 6.5),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 10))

summary_row <- function(gse) {
  csv <- file.path("results", paste0("S2_qc_", gse, ".csv"))
  if (!file.exists(csv)) { cat("skip", gse, "- no CSV\n"); return(NULL) }
  d <- parse_sample(read.csv(csv), gse)
  mt_thr <- MT_THR[[gse]]
  xsz <- if (nrow(d) > 20) 4.5 else 6.5          # v4 F3

  pA <- ggplot(d) +
    geom_segment(aes(x = cells_after, xend = cells_before, y = short, yend = short),
                 color = "grey80", linewidth = 0.5) +
    geom_point(aes(x = cells_before, y = short), color = "grey65", size = 1.2) +
    geom_point(aes(x = cells_after, y = short, fill = cond),                 # v4 F2
               shape = 21, size = 1.8, color = "white", stroke = 0.3) +
    geom_text(aes(x = Inf, y = short, label = paste0(round(pct_kept), "%")),
              hjust = -0.15, size = 1.9, color = "grey40") +
    scale_fill_manual(values = COND_COLORS) +
    coord_cartesian(clip = "off") +
    labs(title = "A  Retention (before → after QC)", x = "cells per sample", y = NULL) +
    theme_pub + theme(plot.margin = margin(5, 28, 5, 5))

  pB <- ggplot(d, aes(short, median_genes)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 200, ymax = 6000,
             fill = "#2E7D32", alpha = 0.08) +
    geom_hline(yintercept = c(200, 6000), color = "#2E7D32", linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(fill = cond), shape = 21, size = 1.8, color = "white", stroke = 0.3) +  # v4 F2
    scale_fill_manual(values = COND_COLORS) +
    labs(title = "B  Complexity vs pre-registered band", x = NULL, y = "median genes / cell") +
    theme_pub + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = xsz))

  pC <- ggplot(d, aes(short, median_mt)) +
    geom_hline(yintercept = mt_thr, color = "#C62828", linetype = "dashed", linewidth = 0.5) +
    geom_point(aes(fill = cond), shape = 21, size = 1.8, color = "white", stroke = 0.3) +  # v4 F2
    scale_fill_manual(values = COND_COLORS) +
    coord_cartesian(ylim = c(0, mt_thr * 1.18)) +
    labs(title = sprintf("C  Mito%% vs threshold (%d%%)", mt_thr),           # v4 F4
         x = NULL, y = "median mitochondrial %") +
    theme_pub + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = xsz))

  panels <- list(pA, pB, pC)
  if (!all(is.na(d$dbl_pct))) {
    pD <- ggplot(d, aes(short, dbl_pct, fill = cond)) +
      geom_col(width = 0.7, show.legend = FALSE) +                        # v4.1 F6
      scale_fill_manual(values = COND_COLORS) +
      labs(title = "D  Doublets (scDblFinder)", x = NULL, y = "doublet %") +
      theme_pub + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = xsz))
    panels <- c(panels, list(pD))
  } else {
    cat(gse, ": dbl_pct all NA — panel D skipped.\n")
  }
  h <- if (nrow(d) > 20) 7.5 else 4.5
  row_fig <- wrap_plots(panels, nrow = 1, guides = "collect") +
    plot_annotation(title = sprintf("%s — %d single-cell samples", gse, nrow(d)),
                    theme = theme(plot.title = element_text(face = "bold", size = 12))) &
    theme(legend.position = "bottom")
  ggsave(file.path("figures", paste0("S2_qc_pub_", gse, "_summary.png")),
         row_fig, width = 14, height = h, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0("S2_qc_pub_", gse, "_summary.pdf")),     # v4 F5
         row_fig, width = 14, height = h, bg = "white", device = cairo_pdf)
  cat("summary figure written for", gse, "\n")
  invisible(d)
}

all_d <- list()
for (gse in c("GSE179640", "GSE214411", "GSE183837")) {
  r <- summary_row(gse)
  if (!is.null(r)) { r$dataset <- gse; all_d[[gse]] <- r }
}

# ---- v2: vertical violin stack (same look lead approved, un-cramped) --------
# 3 metrics stacked vertically, ONE shared x-axis, short labels colored by
# condition, thin boxplot overlay, pre-registered thresholds as dashed lines.
violin_stack <- function(gse) {
  f <- file.path("data/processed", paste0(gse, "_filtered.rds"))
  if (!file.exists(f)) { cat("skip violins for", gse, "- no rds\n"); return(invisible(NULL)) }
  m <- readRDS(f)@meta.data
  m$sample <- m$orig.ident
  m <- parse_sample(as.data.frame(m), gse)          # v4 F1: no longer crashes
  sl <- unique(m[, c("short", "cond")])
  sl <- sl[order(match(as.character(sl$short), levels(m$short))), ]
  L <- rbind(
    data.frame(short = m$short, metric = "genes / cell",      value = m$nFeature_RNA),
    data.frame(short = m$short, metric = "UMI counts / cell", value = m$nCount_RNA),
    data.frame(short = m$short, metric = "mitochondrial %",   value = m$percent.mt))
  L$metric <- factor(L$metric, levels = c("genes / cell", "UMI counts / cell",
                                          "mitochondrial %"))
  n <- length(levels(m$short))
  sample_cols <- setNames(scales::hue_pal()(n), levels(m$short))  # same hue family as Seurat default
  thr <- data.frame(metric = factor(c("genes / cell", "genes / cell", "mitochondrial %"),
                                    levels = levels(L$metric)),
                    yint = c(200, 6000, MT_THR[[gse]]),
                    col  = c("#2E7D32", "#2E7D32", "#C62828"))
  lbl_cols <- COND_COLORS[as.character(sl$cond)]
  p <- ggplot(L, aes(short, value, fill = short)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.15) +
    geom_boxplot(width = 0.07, outlier.shape = NA, alpha = 0.5, linewidth = 0.2) +
    geom_hline(data = thr, aes(yintercept = yint, color = col),
               linetype = "dashed", linewidth = 0.4) +
    scale_color_identity() +
    facet_grid(metric ~ ., scales = "free_y") +
    scale_fill_manual(values = sample_cols, guide = "none") +
    labs(title = sprintf("%s — post-QC cell-level QC (%d samples, %s cells)",
                         gse, n, format(nrow(m), big.mark = ",")),
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 7.5, colour = lbl_cols),
          axis.text.y = element_text(size = 7.5),
          strip.text = element_text(face = "bold", size = 9),
          plot.title = element_text(face = "bold", size = 11))
  w <- max(6.5, 0.32 * n + 2)
  ggsave(file.path("figures", paste0("S2_qc_pub_", gse, "_violins.png")), p,
         width = w, height = 9.5, dpi = 300, bg = "white")
  ggsave(file.path("figures", paste0("S2_qc_pub_", gse, "_violins.pdf")), p,   # v4 F5
         width = w, height = 9.5, bg = "white", device = cairo_pdf)
  cat("violin stack written for", gse, "\n")
}
for (gse in c("GSE179640", "GSE214411", "GSE183837")) violin_stack(gse)

# ---- Cell-level ECDF panels (solve the unreadable-violin problem) ------------
rds <- list.files("data/processed", pattern = "_filtered\\.rds$", full.names = TRUE)
if (length(rds)) {
  meta_all <- list()
  for (f in rds) {
    gse <- sub("_filtered\\.rds$", "", basename(f))
    m <- readRDS(f)@meta.data
    m$sample <- m$orig.ident
    m <- parse_sample(as.data.frame(m), gse)        # v4 F1: no longer crashes
    m$dataset <- gse
    meta_all[[gse]] <- m[, c("dataset", "sample", "short", "cond",
                             "nFeature_RNA", "nCount_RNA", "percent.mt")]
  }
  M <- do.call(rbind, meta_all)
  thr <- data.frame(dataset = names(MT_THR), metric = "percent.mt",
                    x = as.numeric(MT_THR))
  Mm <- rbind(
    data.frame(M[, c("dataset","sample","cond")], metric = "nFeature_RNA", value = M$nFeature_RNA),
    data.frame(M[, c("dataset","sample","cond")], metric = "nCount_RNA",  value = M$nCount_RNA),
    data.frame(M[, c("dataset","sample","cond")], metric = "percent.mt",  value = M$percent.mt))
  pE <- ggplot(Mm, aes(value, color = cond, group = interaction(dataset, sample))) +
    stat_ecdf(linewidth = 0.35, pad = FALSE) +
    geom_vline(data = thr, aes(xintercept = x), color = "#C62828",
               linetype = "dashed", linewidth = 0.4) +
    scale_color_manual(values = COND_COLORS) +
    facet_grid(dataset ~ metric, scales = "free") +
    labs(title = "Cell-level QC distributions (ECDF) — one curve per sample",
         x = "value", y = "cumulative fraction of cells") +
    theme_pub + theme(strip.text = element_text(face = "bold", size = 8),
                      legend.position = "bottom")
  ggsave("figures/S2_qc_pub_ecdf.png", pE, width = 12, height = 3.2 * length(rds) + 1.5,
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave("figures/S2_qc_pub_ecdf.pdf", pE, width = 12, height = 3.2 * length(rds) + 1.5,
         bg = "white", limitsize = FALSE, device = cairo_pdf)                  # v4 F5
  cat("ECDF figure written\n")
}
cat("S2_04 done.\n")
