# ============================================================================
# make_figures.R -- publication figures, locked results S4 / S4_05 / S5 / S6
# v3.0 (2026-09-23): Fig 6 implemented -- S6 GSE213216 external validation,
#   drawn ONLY from locked S6_01 outputs (presentation only, R5/A5):
#   6a direction heatmap (7 scores x tested lineages, C1) with FDR/raw-p marks;
#   6b the sole FDR<0.05 hit (tier1_signature, mural compartment) at patient
#   level incl. independent-patient sensitivity p; 6c the honest negative
#   (pre-registered P1: SERPINE1 endothelium NOT confirmed). Negatives are
#   plotted with the same weight as positives (R1/R6).
# v3.1 (2026-09-23): PI visual feedback -- fig6b/c median crossbars halved
#   (linewidth 0.6 -> 0.3) so medians no longer dominate the panels.
# v3.2 (2026-09-23): Fig 7 added -- S8 reverse-confirmation in L1000 space,
#   drawn ONLY from locked S8_01 outputs (7 SERPINE1-KD lines, 10 estradiol
#   cell-line aggregates; honest negatives plotted with equal weight).
# v3.3 (2026-09-25, PI review): Fig 7 colored by reversal direction (crimson
#   = moves opposite to lesion / navy = with lesion); HA1E (strongest
#   single-line signal, hypothesis-generating only) ringed amber + subtitle
#   flag. Presentation only (A5/R5) -- no number changed, colors/labels only.
# v2.2 (2026-09-22): fix measured on the lead's REAL-data v2.1 render:
#   #24 fig4b title/subtitle clipped at the right panel edge (real counts
#       "21 of 715 ... (all in GSE214411)" made the one-line title ~90 chars
#       in a 5.6in panel). Fix: short fixed title; counts moved into a
#       two-line subtitle.
# v2.1 (2026-09-22): live-fire fixes measured on the lead's machine (R6/R9):
#   #22 CRASH: ggsave(device="cairo_pdf") -> "Unknown graphics device" -- the
#       cairo_pdf device name is not accepted by the current ggplot2 device
#       validator, so the run died on the FIRST save and every figure after
#       fig1a was never drawn. Fix: base "pdf" device (cairo-capable R uses
#       cairo automatically); PNG saved without the legacy type="cairo" arg
#       (kills the ragg fallback warning too).
#   #23 fig1a tile labels STILL overflowed their tiles in v2.0 (measured on
#       the returned PNG: multi-line labels ~33mm wide in ~20mm tiles).
#       Fix: restructure -- only GSE id + n inside the tile (max 9 chars),
#       contrast description moved to a two-line axis label below each tile.
#   Also measured and fixed in v2.1 (sandbox rehearsal render):
#     * combo plates squeezed panels -> titles/subtitles clipped; combo() now
#       passes each panel's design width to grid.arrange (no squeeze).
#     * fig2c median-I2 label could clip at the right edge (median 85.7%) ->
#       label placed left of the line.
#     * fig3b base-line label overlapped the dashes -> moved above the line.
#     * fig4 strip labels clipped ("senescence_dorma...") -> wrapped two-line
#       strip labels; bulk contrast tick labels angled to stop touching.
#     * sig_fdr05 parse guard: fail-loud unless logical/TRUE/FALSE (R9) --
#       a character column silently produced an all-NA panel in testing.
#     * fig5a in-plot legend overlapped the senescence_dormancy CI whisker ->
#       legend moved below the panel.
#     * fig5b "senescence_dormancy" point label clipped at the panel edge ->
#       long names wrapped to two lines.
#   Changes remain PRESENTATION ONLY (R5/A5): no new estimands, no new
#   thresholds; every plotted number is read verbatim from a locked results
#   file; rounding/formatting only.
# Fail-loud on any missing input (R2). Deterministic: seed fixed.
# Usage: Rscript scripts/make_figures.R
# ============================================================================
set.seed(42)
source("scripts/S4_common.R")   # CFG4 paths + meta_set + paired_arm

FIGDIR <- "figures"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("FATAL: ggplot2 required (in renv lock).", call. = FALSE)
suppressPackageStartupMessages(library(ggplot2))
have_grid <- requireNamespace("gridExtra", quietly = TRUE)

## ---- locked palette (bold, print-safe, colorblind-tolerated) ---------------
C_DOWN  <- "#1F4E79"   # strong navy
C_UP    <- "#C00000"   # strong crimson
C_NS    <- "#BFBFBF"   # light grey
C_ACC   <- "#2E7D32"   # green (pass / concordance)
C_WARN  <- "#C55A11"   # amber (reported / pending)
C_PURP  <- "#5B2D8E"   # purple (heterogeneity)
C_FAIL  <- "#8B1A1A"   # dark red (fail)
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

r2 <- function(x) sprintf("%.2f", x)
r3 <- function(x) sprintf("%.3f", x)
pct1 <- function(x) sprintf("%.1f%%", 100 * x)

## ---- helpers ---------------------------------------------------------------
need <- function(p) {
  if (!file.exists(p))
    stop("FATAL: missing locked results file ", p,
         " -- run the owning stage first (R2).", call. = FALSE)
  invisible(p)
}
read_locked <- function(p) { need(p); read.csv(p, stringsAsFactors = FALSE) }

## v2.1 #22: no cairo_pdf, no type="cairo" -- plain "png"/"pdf" devices work
## on every R build (ragg is picked automatically for PNG when installed).
save_both <- function(plot, name, w, h) {
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".png")), plot, width = w,
                  height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".pdf")), plot, width = w,
                  height = h, device = "pdf", bg = "white")
  say("wrote figures/", name, ".png/.pdf")
}

## v2.1: ws = per-panel design widths; the plate gets sum(ws) so no panel is
## ever squeezed (measured cause of clipped titles/labels in v2.0 combos).
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
CORE_COLS <- c("gene","k","M","SE","p","tau2","I2","Q","QEp","k_pos","k_neg",
               "FDR","concord_k","tier")
check_cols <- function(d, cols, p) {
  miss <- setdiff(cols, names(d))
  if (length(miss))
    stop("FATAL: ", p, " lacks locked columns: ", paste(miss, collapse = ", "),
         " | got: ", paste(names(d), collapse = ","), " (R2)", call. = FALSE)
}
## v2.1: fail-loud flag parser (R9) -- a non-logical sig column must never
## silently produce an empty/NA panel.
as_flag <- function(d, col, p) {
  if (is.logical(d[[col]])) return(d[[col]])
  v <- toupper(trimws(as.character(d[[col]])))
  if (any(!v %in% c("TRUE", "FALSE")))
    stop("FATAL: ", p, " column ", col, " is not logical/TRUE/FALSE (R2/R9)",
         call. = FALSE)
  v == "TRUE"
}

say("== make_figures v2.1: presentation only, zero new computation (A5) ==")

## ---- inputs ----------------------------------------------------------------
f_all   <- file.path(CFG4$out_meta, "S4_meta_all_genes.csv")
f_t1    <- file.path(CFG4$out_meta, "S4_meta_Tier1_core.csv")
f_sens  <- file.path(CFG4$out_meta, "S4_sensitivity_summary.csv")
f_bulk5 <- "results/S4_localization/S4_05_axis_bulk_map.csv"
f_scr5  <- "results/S4_localization/S4_05_axis_scrna_map.csv"
f_pmet  <- "results/S5_meta/S5_meta_pathways.csv"
f_pconc <- "results/S5_meta/S5_concordance_GSE11691.csv"
f_psens <- "results/S5_meta/S5_sensitivity_scenarios.csv"
M  <- read_locked(f_all);  check_cols(M,  CORE_COLS, f_all)
t1 <- read_locked(f_t1);  check_cols(t1, CORE_COLS, f_t1)
sens  <- read_locked(f_sens)
bulk5 <- read_locked(f_bulk5)
scr5  <- read_locked(f_scr5)
pm  <- read_locked(f_pmet)
pconc <- read_locked(f_pconc)
psens <- read_locked(f_psens)
check_cols(bulk5, c("axis","gene","contrast_id","logFC","sig_fdr05"), f_bulk5)
check_cols(scr5,  c("axis","gene","dataset","lineage","arm","logFC","sig_fdr05"), f_scr5)
check_cols(pm,    c("pathway","M","SE","I2","FDR"), f_pmet)
check_cols(pconc, c("pathway","M_base","logFC_GSE11691"), f_pconc)
check_cols(psens, c("scenario","FDR"), f_psens)
bulk5$sig_fdr05 <- as_flag(bulk5, "sig_fdr05", f_bulk5)
scr5$sig_fdr05  <- as_flag(scr5,  "sig_fdr05", f_scr5)

## ===========================================================================
## Fig 1 -- study overview
## ===========================================================================
## v2.1 #23: tiles carry ONLY the GSE id and the n (<= 9 chars, fits easily);
## the contrast description lives in a two-line axis label below the tile.
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
fig1a <- ggplot(coh, aes(x = x, y = 1)) +
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
  labs(title = "Seven locked primary contrasts",
       subtitle = "six disease-vs-control cohorts enter the random-effects meta; the paired ectopic-vs-eutopic arm is reported separately") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.x = element_text(size = 6, color = INK, lineheight = 0.95),
        legend.position = "bottom")

flow <- data.frame(
  x = seq_len(4),
  lbl = c("S3\n19 frozen contrasts",
          "S4\ngene-level REML meta\nTier-1 core",
          "S4_05 + S5\nlocalization +\npathway-level meta",
          "S6\nexternal validation\n(pending R3)"),
  fill = c("#DCE6F1", "#B8CCE4", "#8FAADC", "#F2DCB3"))
fig1b <- ggplot(flow, aes(x = x, y = 1)) +
  geom_tile(aes(fill = fill), height = 0.78, width = 0.84,
            color = "#7F7F7F", linewidth = 0.4) +
  geom_text(aes(label = lbl), size = 2.5, lineheight = 0.95, color = INK) +
  geom_segment(data = data.frame(x = c(1.44, 2.44, 3.44)),
               aes(x = x, xend = x + 0.12, y = 1, yend = 1),
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
save_both(fig1a, "fig1a_contrasts", 7.2, 2.6)
save_both(fig1b, "fig1b_pipeline", 7.2, 1.9)
combo(list(fig1a, fig1b), "fig1_overview", c(7.2, 7.2), 2.8)

## ===========================================================================
## Fig 2 -- Tier-1 core
## ===========================================================================
Mv <- M
Mv$neglog <- ifelse(Mv$FDR > 0, -log10(Mv$FDR),
                    max(-log10(Mv$FDR[Mv$FDR > 0]), na.rm = TRUE) + 0.5)
Mv$tier <- factor(Mv$tier, levels = c("Tier3", "Tier2", "Tier1"))
fig2a <- ggplot(Mv, aes(x = M, y = neglog, color = tier)) +
  geom_point(size = 0.6, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = c(Tier3 = C_NS, Tier2 = "#7A9CC4",
                                Tier1 = C_UP),
                     breaks = c("Tier1", "Tier2", "Tier3")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = INK, linewidth = 0.4) +
  annotate("text", x = -Inf, y = -log10(0.05), label = "  FDR = 0.05",
           hjust = 0, vjust = -0.4, size = 2.6, color = INK, fontface = "bold") +
  labs(x = "REML meta effect size M (log2 FC)", y = expression(-log[10] ~ FDR),
       title = paste0("Cross-cohort random-effects meta (",
                      format(nrow(Mv), big.mark = ","), " genes)"),
       subtitle = paste0("Tier-1 core: ", nrow(t1),
                         " genes (FDR < 0.05, |M| >= 0.5, same sign in >= 4 cohorts)")) +
  theme_pub() +
  theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "#BFBFBF", linewidth = 0.3))

t1d <- t1[order(-abs(t1$M)), ][seq_len(min(12, nrow(t1))), ]
t1d$gene <- factor(t1d$gene, levels = rev(t1d$gene))
fig2b <- ggplot(t1d, aes(x = gene, y = M)) +
  geom_hline(yintercept = 0, color = "#7F7F7F", linewidth = 0.4) +
  geom_errorbar(aes(ymin = M - 1.96 * SE, ymax = M + 1.96 * SE),
                width = 0.22, color = "#595959", linewidth = 0.5) +
  geom_point(aes(color = M > 0), size = 2.2) +
  scale_color_manual(values = c(`TRUE` = C_UP, `FALSE` = C_DOWN),
                     breaks = c(TRUE, FALSE), labels = c("up", "down")) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "M (95% CI)",
       title = "Strongest Tier-1 genes by |M|",
       subtitle = "display ranking only; full statistics in S4_meta_Tier1_core.csv") +
  theme_pub() +
  theme(legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "#BFBFBF", linewidth = 0.3))

med_i2 <- stats::median(M$I2, na.rm = TRUE)
fig2c <- ggplot(M, aes(x = I2)) +
  geom_histogram(bins = 34, fill = "#D9D9D9", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = med_i2, color = C_PURP, linewidth = 0.7) +
  annotate("text", x = med_i2, y = Inf,
           label = paste0("median I2 = ", round(med_i2, 1), "% "),
           vjust = 1.6, hjust = 1, size = 2.7, color = C_PURP, fontface = "bold") +
  labs(x = expression(I^2 ~ "(%)"), y = "genes",
       title = "Between-cohort heterogeneity",
       subtitle = "heterogeneity is reported for every gene and never filtered (locked rule)") +
  theme_pub()
save_both(fig2a, "fig2a_volcano", 4.6, 4.0)
save_both(fig2b, "fig2b_forest", 4.2, 4.0)
save_both(fig2c, "fig2c_heterogeneity", 4.0, 4.0)
combo(list(fig2a, fig2b, fig2c), "fig2_tier1_core", c(4.6, 4.2, 4.0), 4.2)

say("fig1 + fig2 written")

## ===========================================================================
## Fig 3 -- sensitivity battery
## ===========================================================================
sens$status <- ifelse(sens$pass, "PASS",
                      ifelse(grepl("reported", sens$criterion), "REPORTED", "FAIL"))
sens$dot_col <- ifelse(sens$status == "PASS", C_ACC,
                       ifelse(sens$status == "REPORTED", C_WARN, C_FAIL))
sens$val_disp <- ifelse(grepl("|", sens$value, fixed = TRUE), sens$value,
                 ifelse(grepl("retention", sens$metric),
                        pct1(as.numeric(sens$value)), r3(as.numeric(sens$value))))
sens$test <- factor(sens$test, levels = rev(unique(sens$test)))
fig3a <- ggplot(sens, aes(x = 1, y = test)) +
  geom_point(aes(color = I(dot_col)), size = 3, show.legend = FALSE) +
  geom_text(aes(x = 1.06, label = paste0(metric, "  |  ", val_disp,
                                         "  (", criterion, ")"),
                hjust = 0), size = 2.7, color = INK) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.9, 3.8), expand = c(0, 0)) +
  labs(title = "Locked sensitivity battery",
       subtitle = "PASS / FAIL against pre-registered criteria; swaps and boundary tests are reported, never hidden") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.x = element_blank())

sens_files <- list.files(CFG4$out_meta, pattern = "^S4_sens_.*_all\\.csv$",
                         full.names = TRUE)
if (length(sens_files) == 0L)
  stop("FATAL: no S4_sens_*_all.csv files -- run S4_02 first (R2).")
base_set <- t1$gene
ret <- do.call(rbind, lapply(sens_files, function(p) {
  d <- read.csv(p, stringsAsFactors = FALSE)
  tag <- sub("^S4_sens_", "", sub("_all\\.csv$", "", basename(p)))
  set <- d$gene[d$tier == "Tier1"]
  data.frame(scenario = tag, retained = length(intersect(set, base_set)),
             stringsAsFactors = FALSE)
}))
grp <- ifelse(grepl("^swap", ret$scenario), "QC swap",
              ifelse(ret$scenario == "FE", "FE vs RE",
                     ifelse(grepl("^sans", ret$scenario), "sans", "LOCO")))
ret$grp <- factor(grp, levels = c("sans", "LOCO", "FE vs RE", "QC swap"))
ret$base <- length(base_set)
ret <- ret[order(ret$grp, -ret$retained), ]
ret$scenario <- factor(ret$scenario, levels = ret$scenario)
fig3b <- ggplot(ret, aes(x = scenario, y = retained, fill = grp)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = length(base_set), linetype = "dashed",
             color = C_UP, linewidth = 0.6) +
  annotate("text", x = 1, y = length(base_set), vjust = -0.5, hjust = 0,
           size = 2.6, color = C_UP, fontface = "bold",
           label = paste0("base Tier-1 = ", length(base_set))) +
  geom_text(aes(label = retained), vjust = -0.35, size = 2.5, color = INK) +
  scale_fill_manual(values = c(sans = "#8FAADC", LOCO = C_DOWN,
                               `FE vs RE` = C_PURP, `QC swap` = C_ACC)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Tier-1 genes retained",
       title = "Tier-1 set stability across scenarios",
       subtitle = "leave-one-cohort-out, sans-G7307, fixed-effects contrast, and frozen QC-swap replications") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        legend.position = "right")
save_both(fig3a, "fig3a_battery", 6.6, 3.6)
save_both(fig3b, "fig3b_retention", 6.6, 3.6)
combo(list(fig3a, fig3b), "fig3_sensitivity", c(6.6, 6.6), 3.8)

## ===========================================================================
## Fig 4 -- S4_05 localization
## ===========================================================================
AXES_ORD <- c("cytotoxicity", "antigen_presentation", "senescence_dormancy",
              "proliferation", "stromal_ecm")
## v2.1: two-line strip labels (measured clip of "senescence_dormancy")
AXES_LBL <- c("cytotoxicity" = "cytotoxicity",
              "antigen_presentation" = "antigen\npresentation",
              "senescence_dormancy" = "senescence\ndormancy",
              "proliferation" = "proliferation",
              "stromal_ecm" = "stromal / ECM")
ord5 <- c(CFG4$meta_set, CFG4$paired_arm)
b5 <- bulk5
b5$contrast_id <- factor(b5$contrast_id, levels = ord5)
b5$axis <- factor(b5$axis, levels = AXES_ORD)
b5$fillv <- pmax(pmin(b5$logFC, 2), -2)
fig4a <- ggplot(b5, aes(x = contrast_id, y = gene)) +
  geom_tile(aes(fill = fillv), color = "white", linewidth = 0.4) +
  geom_tile(data = b5[b5$sig_fdr05, ], aes(x = contrast_id, y = gene),
            fill = NA, color = INK, linewidth = 0.55) +
  scale_fill_gradient2(low = C_DOWN, mid = "white", high = C_UP,
                       midpoint = 0, limits = c(-2, 2),
                       name = "log2 FC (clipped)") +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(axis = AXES_LBL)) +
  labs(x = NULL, y = NULL,
       title = "Axis genes across the locked primary contrasts (bulk)",
       subtitle = "tile color = log2 FC; black border = BH-FDR < 0.05; last column = paired arm (never pooled)") +
  theme_pub() +
  theme(panel.spacing.y = grid::unit(2.5, "mm"),
        panel.spacing.x = grid::unit(1.5, "mm"),
        axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1,
                                   size = 7, face = "bold"),
        strip.text.y.left = element_text(size = 7.5, lineheight = 0.9),
        legend.position = "right",
        legend.title = element_text(size = 7.5, face = "bold"),
        strip.placement = "outside")

s5 <- scr5[scr5$sig_fdr05, , drop = FALSE]
if (nrow(s5) == 0L)
  stop("FATAL: zero significant scRNA cells -- an empty figure would be ",
       "silently delivered (R6).", call. = FALSE)
s5_ds <- paste(unique(s5$dataset), collapse = ", ")
s5$leg <- gsub("GSE214411 | ", "", paste(s5$dataset, s5$lineage, s5$arm,
                                         sep = " | "), fixed = TRUE)
s5$axis <- factor(s5$axis, levels = AXES_ORD)
s5$gene <- factor(s5$gene, levels = unique(s5$gene[order(s5$axis, s5$gene)]))
fig4b <- ggplot(s5, aes(x = leg, y = gene)) +
  geom_point(aes(size = abs(logFC), color = logFC > 0), alpha = 0.9) +
  scale_color_manual(values = c(`TRUE` = C_UP, `FALSE` = C_DOWN),
                     breaks = c(TRUE, FALSE), labels = c("up", "down")) +
  scale_size_continuous(range = c(2.2, 5.2), breaks = c(1.5, 2, 3),
                        name = "|log2 FC|") +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(axis = AXES_LBL)) +
  labs(x = NULL, y = NULL,
       title = "Significant axis-gene cells in scRNA pseudobulk",
       subtitle = paste0(nrow(s5), " of ", nrow(scr5),
                         " measured cells pass BH-FDR < 0.05 (all in ",
                         s5_ds, ")\n",
                         "blank = not significant (absence is measured, not dropped)")) +
  theme_pub() +
  theme(panel.spacing.y = grid::unit(2.5, "mm"),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.8),
        strip.text.y.left = element_text(size = 7.5, lineheight = 0.9),
        legend.position = "right",
        legend.title = element_text(size = 7.5, face = "bold"),
        strip.placement = "outside")
save_both(fig4a, "fig4a_bulk_heatmap", 5.6, 6.4)
save_both(fig4b, "fig4b_scrna_dots", 5.6, 6.4)
combo(list(fig4a, fig4b), "fig4_localization", c(5.6, 5.6), 6.6)

say("fig3 + fig4 written")

## ===========================================================================
## Fig 5 -- S5 pathway meta
## ===========================================================================
pm$pathway <- factor(pm$pathway,
                     levels = pm$pathway[order(pm$M, decreasing = TRUE)])
fig5a <- ggplot(pm, aes(x = pathway, y = M)) +
  geom_hline(yintercept = 0, color = "#7F7F7F", linewidth = 0.4) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted",
             color = "#A6A6A6", linewidth = 0.4) +
  geom_errorbar(aes(ymin = M - 1.96 * SE, ymax = M + 1.96 * SE),
                width = 0.24, color = "#595959", linewidth = 0.55) +
  geom_point(aes(color = FDR < 0.05), size = 2.6) +
  geom_text(aes(y = M + 1.96 * SE, label = paste0("I2 = ", round(I2))),
            vjust = -0.5, size = 2.3, color = "#595959") +
  scale_color_manual(values = c(`TRUE` = C_UP, `FALSE` = "#8C8C8C"),
                     breaks = c(TRUE, FALSE),
                     labels = c("BH-FDR < 0.05", "n.s.")) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "pathway score log2 FC (95% CI)",
       title = "Pathway-level random-effects meta (six cohorts)",
       subtitle = "dotted lines = interpretive |M| = 0.5 (never a filter); BH-FDR across six pathway tests") +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.margin = margin(t = -4))

xr <- diff(range(pconc$M_base, na.rm = TRUE))
yr <- diff(range(pconc$logFC_GSE11691, na.rm = TRUE))
pconc$lx <- pconc$M_base
pconc$lx[pconc$M_base > 0] <- pconc$M_base[pconc$M_base > 0] + 0.04 * xr
pconc$lx[pconc$M_base <= 0] <- pconc$M_base[pconc$M_base <= 0] - 0.04 * xr
pconc$hj <- ifelse(pconc$M_base > 0, 0, 1)
## v2.1: wrap long pathway names (measured right-edge clip of senescence_dormancy)
pconc$lbl <- as.character(pconc$pathway)
w <- AXES_LBL[pconc$lbl]
pconc$lbl[!is.na(w)] <- w[!is.na(w)]
fig5b <- ggplot(pconc, aes(x = M_base, y = logFC_GSE11691)) +
  geom_hline(yintercept = 0, color = "#D9D9D9", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "#D9D9D9", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#8C8C8C", linewidth = 0.5) +
  geom_point(size = 2.8, color = C_DOWN) +
  geom_text(aes(x = lx, label = lbl, hjust = hj), size = 2.6,
            color = INK, fontface = "bold", vjust = -0.5, lineheight = 0.9) +
  scale_x_continuous(expand = expansion(mult = 0.42)) +
  scale_y_continuous(expand = expansion(mult = 0.30)) +
  labs(x = "M (six-cohort base meta)", y = "paired ectopic-vs-eutopic log2 FC",
       title = "GSE11691 paired concordance arm",
       subtitle = "directional concordance only -- the paired arm is never pooled into the meta") +
  theme_pub()

fdr_mat <- do.call(rbind, lapply(strsplit(psens$FDR, "; ", fixed = TRUE),
                                 function(x) {
  v <- strsplit(x, "=", fixed = TRUE)
  data.frame(pathway = vapply(v, `[[`, character(1), 1),
             FDR = as.numeric(vapply(v, `[[`, character(1), 2)),
             stringsAsFactors = FALSE)
}))
fdr_mat$scenario <- rep(psens$scenario, times = vapply(strsplit(psens$FDR, "; ", fixed = TRUE), length, integer(1)))
PATH_ORD <- pm$pathway[order(pm$M, decreasing = TRUE)]
SC_ORD <- c("sans-G7307", paste0("LOCO-", CFG4$meta_set),
            paste0("swap-", CFG4$qc_swap), "FE-vs-RE")
fdr_mat$pathway <- factor(fdr_mat$pathway, levels = PATH_ORD)
fdr_mat$scenario <- factor(fdr_mat$scenario, levels = SC_ORD)
if (any(is.na(fdr_mat$FDR)) || anyNA(fdr_mat$pathway) || anyNA(fdr_mat$scenario))
  stop("FATAL: scenario/pathway labels in S5_sensitivity_scenarios.csv do not ",
       "match the locked meta_set / qc_swap / pathway names (R2): got scenario[",
       paste(unique(psens$scenario), collapse = ", "), "] pathway[",
       paste(unique(fdr_mat$pathway[is.na(fdr_mat$pathway)]), collapse = ","),
       "]", call. = FALSE)
fdr_mat$nl <- pmin(-log10(fdr_mat$FDR), 3)
fdr_mat$txt <- ifelse(fdr_mat$FDR < 0.001, format(fdr_mat$FDR, digits = 1,
                                                  scientific = TRUE),
                      r3(fdr_mat$FDR))
fig5c <- ggplot(fdr_mat, aes(x = scenario, y = pathway)) +
  geom_tile(aes(fill = nl), color = "white", linewidth = 0.5) +
  geom_text(aes(label = txt, color = nl > 1.7), size = 2.4, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient(low = "white", high = C_UP, limits = c(0, 3),
                      breaks = c(0, 1.3, 3),
                      labels = c("0", "1.3", ">=3 (FDR<=0.001)"),
                      name = expression(-log[10] ~ FDR)) +
  labs(x = NULL, y = NULL,
       title = "Sensitivity scenarios at pathway level",
       subtitle = "cell text = BH-FDR (six pathway tests); tile depth = -log10(FDR)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 7.5, face = "bold"))

pm2 <- pm
pm2$pathway <- factor(pm2$pathway, levels = PATH_ORD)
fig5d <- ggplot(pm2, aes(x = pathway, y = I2)) +
  geom_col(width = 0.68, fill = C_PURP, color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(I2)), hjust = -0.25, size = 2.5, color = INK) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression(I^2 ~ "(%)"),
       title = "Between-cohort heterogeneity per pathway",
       subtitle = "reported, never filtered (locked rule)") +
  theme_pub()
save_both(fig5a, "fig5a_pathway_forest", 5.6, 3.9)
save_both(fig5b, "fig5b_concordance", 5.2, 4.4)
save_both(fig5c, "fig5c_scenarios", 7.0, 3.9)
save_both(fig5d, "fig5d_pathway_I2", 5.2, 3.6)
combo(list(fig5a, fig5b), "fig5_pathway_meta_AB", c(5.6, 5.2), 4.3)
combo(list(fig5c, fig5d), "fig5_pathway_meta_CD", c(7.0, 5.2), 4.0)

## ===========================================================================
## Fig 6 -- S6 GSE213216 external validation (locked S6_01 outputs)
## ===========================================================================
## >>> FIG6_BLOCK
f_s6t <- "results/S6_meta/S6_01_tests.csv"
f_s6s <- "results/S6_meta/S6_01_scores_long.csv"
s6t <- read_locked(f_s6t)
s6s <- read_locked(f_s6s)
check_cols(s6t, c("endpoint","score","lineage","contrast","n_case","n_ctrl",
                  "delta_median","p","p_indep_patient","status","FDR"), f_s6t)
check_cols(s6s, c("patient","lineage","tissue","tier1_signature",
                  "SERPINE1_logCPM"), f_s6s)

S6_SCORES <- c("cytotoxicity","antigen_presentation","senescence_dormancy",
               "proliferation","stromal_ecm","SenMayo","tier1_signature")
c1t <- s6t[s6t$contrast == "C1" & s6t$status == "TESTED" &
             s6t$score %in% S6_SCORES, ]
if (!nrow(c1t))
  stop("FATAL: no TESTED C1 rows in ", f_s6t, " (R2)", call. = FALSE)

LIN_SHORT <- c("Mesenchymal cells"="Mesenchymal","Endothelial cells"="Endothelial",
               "Smooth muscle cells"="Smooth muscle","T/NK cells"="T/NK",
               "B/Plasma cells"="B/Plasma","Myeloid cells"="Myeloid",
               "Epithelial cells"="Epithelial","Erythrocytes"="Erythrocytes",
               "Mast cells"="Mast")
if (any(is.na(LIN_SHORT[c1t$lineage])))
  stop("FATAL: unmapped lineage label in S6 tests (R2/R9)", call. = FALSE)
ord <- c1t[c1t$score == "tier1_signature", ]
ord <- ord$lineage[order(ord$delta_median, decreasing = TRUE)]
c1t$lin  <- factor(LIN_SHORT[c1t$lineage], levels = LIN_SHORT[ord])
c1t$scr  <- factor(gsub("_", " ", c1t$score), levels = gsub("_", " ", S6_SCORES))
c1t$mark <- ifelse(!is.na(c1t$FDR) & c1t$FDR < 0.05, "*",
                   ifelse(!is.na(c1t$p) & c1t$p < 0.05, "+", ""))
lim <- max(abs(c1t$delta_median), na.rm = TRUE)
fig6a <- ggplot(c1t, aes(x = lin, y = scr, fill = delta_median)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = mark), size = 4.5, color = INK, fontface = "bold") +
  scale_fill_gradient2(low = C_DOWN, mid = "white", high = C_UP,
                       midpoint = 0, limits = c(-lim, lim), name = NULL) +
  labs(x = NULL, y = NULL,
       title = "External validation in GSE213216: direction across lineages",
       subtitle = "C1: ectopic lesions vs eutopic endometrium, patient-level pseudobulk (Wilcoxon delta-median)\n* FDR < 0.05    + raw p < 0.05    (negative results shown with equal weight)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "right")

## patient-level panels: the validated hit (6b) and the honest negative (6c)
## v3.0 #25: 3.1in panels need SHORT titles (<= ~24 chars) -- stats and the
## pre-registered verdict live in a two-line subtitle (same lesson as #24).
CASE_T <- c("Peritoneal lesion", "Endometrioma"); CTRL_T <- "Eutopic Endometrium"
mk_s6_panel <- function(lin, ycol, ylab, ttl, verdict) {
  d <- s6s[s6s$lineage == lin & s6s$tissue %in% c(CASE_T, CTRL_T), ]
  tr <- s6t[s6t$contrast == "C1" & s6t$status == "TESTED" &
              s6t$lineage == lin & s6t$score == ycol, ]
  if (nrow(d) == 0L || nrow(tr) != 1L)
    stop("FATAL: S6 panel data mismatch for ", lin, " / ", ycol, " (R2)",
         call. = FALSE)
  d$grp <- ifelse(d$tissue %in% CASE_T, "Ectopic lesion", "Eutopic")
  d$grp <- factor(d$grp, levels = c("Eutopic", "Ectopic lesion"))
  sub1 <- paste0("FDR = ", signif(tr$FDR, 3), " | delta = ",
                 round(tr$delta_median, 2), " | n = ", tr$n_case, " v ",
                 tr$n_ctrl)
  sub2 <- verdict
  if (!is.na(tr$p_indep_patient))
    sub2 <- paste0(sub2, " | indep p = ", signif(tr$p_indep_patient, 3))
  sub <- paste0(sub1, "\n", sub2)
  ggplot(d, aes(x = grp, y = .data[[ycol]], color = grp)) +
    geom_point(position = position_jitter(width = 0.16, height = 0, seed = 42),
               size = 1.9, alpha = 0.9) +
    stat_summary(fun = median, geom = "crossbar", width = 0.45,
                 linewidth = 0.3, color = INK) +
    scale_color_manual(values = c("Eutopic" = C_DOWN, "Ectopic lesion" = C_UP),
                       guide = "none") +
    labs(x = NULL, y = ylab, title = ttl, subtitle = sub) +
    theme_pub()
}
fig6b <- mk_s6_panel("Smooth muscle cells", "tier1_signature",
                     "Tier-1 signature (z)",
                     "Validated: mural Tier-1",
                     "sole FDR<0.05 hit")
fig6c <- mk_s6_panel("Endothelial cells", "SERPINE1_logCPM",
                     "SERPINE1 logCPM",
                     "SERPINE1 - endothelium",
                     "pre-registered P1: NOT CONFIRMED")
save_both(fig6a, "fig6a_validation_heatmap", 6.8, 4.4)
save_both(fig6b, "fig6b_tier1_mural", 3.1, 4.4)
save_both(fig6c, "fig6c_serpine1_endothelium", 3.1, 4.4)
combo(list(fig6a, fig6b, fig6c), "fig6_GSE213216_validation",
      c(6.8, 3.1, 3.1), 4.6)
## <<< FIG6_BLOCK

## >>> FIG7_BLOCK
## Fig 7 -- S8 reverse-confirmation in L1000 space (presentation only, R5/A5):
##   drawn ONLY from the locked S8_01 outputs. Negative r = perturbation moves
##   Tier-1 landmark genes OPPOSITE to the lesion (reversal direction).
##   Both hypotheses NOT CONFIRMED at the locked BH-FDR gate -- plotted with
##   the same weight as a positive would be (R1/R6).
## v3.3 (2026-09-25, PI review): points colored by direction (reversal =
##   crimson, lesion direction = navy; colorblind-tolerated pairing);
##   HA1E -- the strongest single-line signal (hypothesis-generating only) --
##   ringed amber and flagged in the subtitle (requested by the lead).
f_s7k <- "results/S8_perturb/S8_01_kd_scores.csv"
f_s7e <- "results/S8_perturb/S8_01_estradiol_scores.csv"
s7k <- read_locked(f_s7k)
check_cols(s7k, c("arm","signatureid","n_overlap","r_tier1","cellline"), f_s7k)
s7e <- read_locked(f_s7e)
check_cols(s7e, c("arm","signatureid","n_overlap","r_tier1"), f_s7e)
s7k <- s7k[!is.na(s7k$r_tier1), ]
s7e <- s7e[!is.na(s7e$r_tier1), ]
if (nrow(s7k) < 4L || nrow(s7e) < 4L)
  stop("FATAL: fig7 needs >=4 scored signatures per arm; an empty or gutted ",
       "figure would be silently delivered (R6).", call. = FALSE)

mk_s7_panel <- function(d, ylab, ttl, sub, hl = NA) {
  d <- d[order(d$r_tier1), ]
  d$lab <- factor(d[[ylab]], levels = d[[ylab]])
  d$dir <- factor(ifelse(d$r_tier1 < 0, "reversal direction",
                         "lesion direction"),
                  levels = c("reversal direction", "lesion direction"))
  p <- ggplot(d, aes(x = r_tier1, y = lab)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
    geom_vline(xintercept = median(d$r_tier1), color = INK, linewidth = 0.3,
               linetype = "dashed") +
    geom_point(aes(color = dir), size = 2.2) +
    scale_color_manual(values = c("reversal direction" = C_UP,
                                  "lesion direction" = C_DOWN)) +
    labs(x = "r vs lesion Tier-1 signature (46 landmark genes)",
         y = NULL, title = ttl, subtitle = sub) +
    theme_pub()
  if (!is.na(hl) && hl %in% levels(d$lab)) {
    p <- p +
      geom_point(data = d[d$lab == hl, ], shape = 21, size = 4.2,
                 stroke = 1.3, color = C_WARN, fill = NA)
  }
  p
}
fig7a <- mk_s7_panel(s7k, "cellline", "SERPINE1 knockdown vs lesion program",
  paste0("7 cell lines | median r = ", round(median(s7k$r_tier1), 2),
         " | P1a NOT CONFIRMED (FDR = 0.313)
",
         "HA1E (ringed): strongest single-line signal r = ",
         round(s7k$r_tier1[s7k$cellline == "HA1E"][1], 2),
         " -- hypothesis-generating only"),
  hl = "HA1E")
fig7b <- mk_s7_panel(s7e, "signatureid", "Estradiol vs lesion program",
  paste0(nrow(s7e), " cell lines (aggregated) | median r = ",
         round(median(s7e$r_tier1), 2),
         " | P2a NOT CONFIRMED (FDR = 0.313)"))
save_both(fig7a, "fig7a_serpine1_kd", 4.9, 3.4)
save_both(fig7b, "fig7b_estradiol", 4.9, 3.4)
combo(list(fig7a, fig7b), "fig7_s8_reverse_confirmation", c(4.9, 4.9), 3.6)
## <<< FIG7_BLOCK

say("make_figures v3.3 DONE -- outputs under figures/")
