## ===========================================================================
## make_fig1_overview.R -- Fig 1 (study overview), STANDALONE (v1.0, 2026-09-29)
##   Four panels per the locked figure-architecture plan:
##     Fig1a_dataset_map  | Fig1b_locked_contrasts | Fig1c_pipeline |
##     Fig1d_analytical_branches
##   Split out of make_figures.R v2.x so Fig 1 can be re-rendered alone.
##   Presentation only (A5/R5): no new estimands, no new thresholds.
##   New numbers (GSE213216 cells / patients / pseudobulks) are READ from
##   locked files, fail-loud (R1/R2):
##     results/S2_qc_GSE213216_audit.csv
##     results/S6_meta/S6_01_scores_long.csv
##   Reviewer-facing fixes vs v2.x:
##   * contrasts title de-conflated: six cohorts = meta set; the GSE11691
##     paired arm is NOT a seventh "primary" contrast.
##   * pipeline: "(pending R3)" removed (S6 locked & complete); S8 (LINCS)
##     and S8b (exploratory drug response) boxes added; robust core named.
## Usage (project root):
##   Rscript scripts/make_fig1_overview.R
## ===========================================================================
set.seed(42)
source("scripts/S4_common.R")   # CFG4 paths + meta_set + paired_arm (frozen)

FIGDIR <- "figures"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("FATAL: ggplot2 required (in renv lock).", call. = FALSE)
suppressPackageStartupMessages(library(ggplot2))
have_grid <- requireNamespace("gridExtra", quietly = TRUE)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

## ---- locked palette (identical to make_figures.R) -------------------------
C_DOWN  <- "#1F4E79"   # strong navy
C_UP    <- "#C00000"   # strong crimson
C_NS    <- "#BFBFBF"   # light grey
C_ACC   <- "#2E7D32"   # green
C_WARN  <- "#C55A11"   # amber
C_PURP  <- "#5B2D8E"   # purple (validation atlas)
C_FAIL  <- "#8B1A1A"   # dark red
INK     <- "#262626"
GRID    <- "#E6E6E6"

theme_pub <- function(base = 8.5) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.line  = element_line(color = INK, linewidth = 0.4),
      axis.ticks = element_line(color = INK, linewidth = 0.4),
      axis.title = element_text(color = INK, face = "bold", size = base + 1),
      axis.text  = element_text(color = INK, size = base),
      plot.title    = element_text(color = INK, face = "bold", size = base + 2.5,
                                   hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(color = "#595959", size = base - 0.5,
                                   hjust = 0, margin = margin(b = 4)),
      plot.margin = margin(6, 10, 4, 6),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text  = element_text(size = base - 0.5, color = INK),
      legend.key.size = grid::unit(3.5, "mm"),
      legend.margin = margin(b = 0),
      strip.background = element_rect(fill = "#EFEFEF", color = "#BFBFBF",
                                      linewidth = 0.3),
      strip.text = element_text(face = "bold", size = base, color = INK)
    )
}

## ---- helpers (identical discipline to make_figures.R) ---------------------
need <- function(p) {
  if (!file.exists(p))
    stop("FATAL: missing locked results file ", p,
         " -- run the owning stage first (R2).", call. = FALSE)
  invisible(p)
}
read_locked <- function(p) { need(p); read.csv(p, stringsAsFactors = FALSE) }
check_cols <- function(d, cols, p) {
  miss <- setdiff(cols, names(d))
  if (length(miss))
    stop("FATAL: ", p, " lacks locked columns: ", paste(miss, collapse = ", "),
         " | got: ", paste(names(d), collapse = ","), " (R2)", call. = FALSE)
}
## v2.1 #22: plain "png"/"pdf" devices only (no cairo_pdf / type="cairo")
save_both <- function(plot, name, w, h) {
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".png")), plot, width = w,
                  height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".pdf")), plot, width = w,
                  height = h, device = "pdf", bg = "white")
  say("wrote figures/", name, ".png/.pdf")
}
combo <- function(plots, name, ws, h) {
  if (have_grid) {
    p <- gridExtra::grid.arrange(grobs = plots, ncol = length(plots),
                                 widths = grid::unit(ws, "in"),
                                 padding = grid::unit(3, "mm"))
    save_both(p, name, sum(ws) + 0.3, h)
  } else {
    say("gridExtra absent -- writing panels individually under figures/",
        name, "_panel*.png")
    for (i in seq_along(plots))
      save_both(plots[[i]], paste0(name, "_panel", LETTERS[i]), ws[i], h)
  }
}

say("== make_fig1_overview v1.0: presentation only, zero new computation (A5) ==")

## ---- Fig1a: dataset map ---------------------------------------------------
f_audit <- "results/S2_qc_GSE213216_audit.csv"
f_s6l   <- "results/S6_meta/S6_01_scores_long.csv"
audit <- read_locked(f_audit)
s6l   <- read_locked(f_s6l)
check_cols(audit, c("object", "cells"), f_audit)
check_cols(s6l, c("patient", "lineage"), f_s6l)
n_cells <- audit$cells[audit$object == "auxiliary.seurat.shared.rds"]
if (length(n_cells) != 1 || is.na(n_cells))
  stop("FATAL: GSE213216 cell count not derivable from ", f_audit, call. = FALSE)
n_pat <- length(unique(s6l$patient))
n_pb  <- nrow(s6l)

dm <- rbind(
  data.frame(x = 1:6, y = 4,
             gse = c("GSE7305", "GSE6364", "GSE25628", "GSE7307",
                     "GSE51981", "GSE120103"),
             nn  = c("10 v 10", "21 v 16", "16 v 6", "18 v 23",
                     "77 v 34", "13 v 9"),
             role = "bulk meta set (k_min = 4)"),
  data.frame(x = 1, y = 3, gse = "GSE11691", nn = "9 women",
             role = "paired arm (never pooled)"),
  data.frame(x = 1:3, y = 2,
             gse = c("GSE179640", "GSE183837", "GSE214411"),
             nn  = c("scRNA", "scRNA", "scRNA"),
             role = "scRNA localization"),
  data.frame(x = 1, y = 1, gse = "GSE213216",
             nn = sprintf("%s cells\n%d patients\n%d pseudobulks",
                          format(n_cells, big.mark = ","), n_pat, n_pb),
             role = "validation atlas"))

role_lvls <- c("bulk meta set (k_min = 4)", "paired arm (never pooled)",
               "scRNA localization", "validation atlas")
dm$role <- factor(dm$role, levels = role_lvls)
role_cols <- c("bulk meta set (k_min = 4)" = C_DOWN,
               "paired arm (never pooled)" = C_WARN,
               "scRNA localization"        = C_ACC,
               "validation atlas"          = C_PURP)

fig1a <- ggplot(dm, aes(x = x, y = y)) +
  geom_tile(aes(fill = role), height = 0.62, width = 0.88,
            color = "white", linewidth = 0.8, show.legend = TRUE) +
  geom_text(aes(y = y + 0.09, label = gse), size = 2.8, fontface = "bold",
            color = "white", show.legend = FALSE) +
  geom_text(aes(y = y - 0.12, label = nn), size = 2.2, lineheight = 0.95,
            color = "white", show.legend = FALSE) +
  scale_fill_manual(values = role_cols) +
  scale_y_continuous(breaks = 4:1,
    labels = c("Bulk meta set\n(k_min = 4)", "Paired arm\n(never pooled)",
               "scRNA\nlocalization", "Validation atlas\n(external)"),
    expand = c(0.03, 0.03)) +
  scale_x_continuous(expand = c(0.015, 0.015)) +
  coord_cartesian(xlim = c(0.4, 6.6), ylim = c(0.55, 4.45)) +
  labs(title = "Dataset map",
       subtitle = paste0("six disease-vs-control bulk cohorts; one paired ",
                         "ectopic-vs-eutopic arm; three scRNA localization ",
                         "datasets; one independent single-cell validation ",
                         "atlas (GSE213216)")) +
  theme_pub() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 6.5, color = INK, lineheight = 0.95),
        legend.position = "bottom")

## ---- Fig1b: locked contrasts ----------------------------------------------
coh <- data.frame(
  id   = c(CFG4$meta_set, CFG4$paired_arm),
  gse  = c("GSE7305", "GSE6364", "GSE25628", "GSE7307", "GSE51981",
           "GSE120103", "GSE11691"),
  nn   = c("10 v 10", "21 v 16", "16 v 6", "18 v 23", "77 v 34",
           "13 v 9", "9 women"),
  cmp  = c("Endo vs\nNormal", "Eutopic vs Normal\n(+ phase)",
           "Lesion/Eutopic\nvs Normal", "Endo vs\nNormal",
           "Endo (all stages)\nvs Control", "Stage IV vs\nFertile control",
           "Ectopic vs Eutopic\n(paired)"),
  role = c(rep("meta set (k_min = 4)", 6), "paired arm (never pooled)"))
coh$x <- seq_len(nrow(coh))
fig1b <- ggplot(coh, aes(x = x, y = 1)) +
  geom_tile(aes(fill = role), height = 0.62, width = 0.88, color = "white",
            linewidth = 0.8, show.legend = TRUE) +
  geom_text(aes(y = 1.09, label = gse, color = role), size = 3.0,
            fontface = "bold", show.legend = FALSE) +
  geom_text(aes(y = 0.90, label = nn, color = role), size = 2.5,
            show.legend = FALSE) +
  scale_color_manual(values = c("meta set (k_min = 4)" = "white",
                                "paired arm (never pooled)" = INK)) +
  scale_fill_manual(values = c("meta set (k_min = 4)" = C_DOWN,
                               "paired arm (never pooled)" = C_WARN)) +
  scale_x_continuous(breaks = coh$x, labels = coh$cmp,
                     expand = c(0.015, 0.015)) +
  coord_cartesian(ylim = c(0.60, 1.40)) +
  labs(title = "Locked contrasts: six cohorts (meta set) + one paired arm",
       subtitle = paste0("six disease-vs-control cohorts enter the ",
                         "random-effects meta; the paired ectopic-vs-eutopic ",
                         "arm is reported separately and never pooled")) +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.x = element_text(size = 6, color = INK, lineheight = 0.95),
        legend.position = "bottom")

## ---- Fig1c: locked pipeline -----------------------------------------------
flow <- data.frame(
  x = seq_len(6),
  lbl = c("S3\n19 frozen contrasts",
          "S4\ngene-level REML meta\nTier-1 + robust core",
          "S4_05 + S5\nlocalization +\npathway-level meta",
          "S6\nexternal validation\n(GSE213216)",
          "S8\nLINCS reversal\nanalysis",
          "S8b\nexploratory drug\nresponse"),
  fill = c("#DCE6F1", "#B8CCE4", "#8FAADC",
           "#F2DCB3", "#E4D3EE", "#F0E3C8"))
fig1c <- ggplot(flow, aes(x = x, y = 1)) +
  geom_tile(aes(fill = fill), height = 0.78, width = 0.80,
            color = "#7F7F7F", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 2.3, lineheight = 0.95, color = INK) +
  geom_segment(data = data.frame(x = c(1.42, 2.42, 3.42, 4.42, 5.42)),
               aes(x = x, xend = x + 0.14, y = 1, yend = 1),
               arrow = grid::arrow(length = grid::unit(1.6, "mm"),
                                   type = "closed"),
               color = "#595959", linewidth = 0.5, inherit.aes = FALSE) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  ylim(0.45, 1.55) +
  labs(title = "Locked analysis pipeline") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

## ---- Fig1d: analytical branches -------------------------------------------
br <- data.frame(
  x = seq_len(4),
  lbl = c("Pseudobulk localization\n3 scRNA datasets\n(8 lineages)",
          "External validation\nGSE213216\n(21 patients)",
          "Perturbational reversal\nLINCS L1000 connectivity\n(SERPINE1-KD, estradiol)",
          "Exploratory\ndecidualization drug\nresponse (S8b)"),
  fill = c("#B8CCE4", "#F2DCB3", "#E4D3EE", "#F0E3C8"))
fig1d <- ggplot(br, aes(x = x, y = 1)) +
  geom_tile(aes(fill = fill), height = 0.78, width = 0.80,
            color = "#7F7F7F", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 2.3, lineheight = 0.95, color = INK) +
  geom_segment(data = data.frame(x = c(1.42, 2.42, 3.42)),
               aes(x = x, xend = x + 0.14, y = 1, yend = 1),
               arrow = grid::arrow(length = grid::unit(1.6, "mm"),
                                   type = "closed"),
               color = "#595959", linewidth = 0.5, inherit.aes = FALSE) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  ylim(0.45, 1.55) +
  labs(title = "Analytical branches",
       subtitle = paste0("discovery to validation to reversal; BH-FDR < 0.05 ",
                         "is the sole claim gate; exploratory analyses carry ",
                         "no claim")) +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

## ---- save -----------------------------------------------------------------
save_both(fig1a, "Fig1a_dataset_map",          7.2, 3.0)
save_both(fig1b, "Fig1b_locked_contrasts",     7.2, 2.6)
save_both(fig1c, "Fig1c_pipeline",             7.6, 1.9)
save_both(fig1d, "Fig1d_analytical_branches",  7.2, 1.9)
combo(list(fig1a, fig1b, fig1c, fig1d), "fig1_overview",
      c(7.2, 7.2, 7.6, 7.2), 8.6)

say("== make_fig1_overview DONE ==")
