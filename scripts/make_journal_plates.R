## ===========================================================================
## make_journal_plates.R -- FINAL journal plates (v1.0, 2026-09-30)
##   One folder "FIGURES PLATES FOR JOURNAL SUBMISSION" in the project root,
##   containing Fig1..Fig7 and FigS1..FigS7 (PNG + PDF, 300 dpi), each panel
##   lettered, all numbers read from locked results (R1/R2), fail-loud.
##   Fig 1: 2x2 manual-style layout (dataset map | contrasts / pipeline | branches)
##   Applied fixes vs earlier drafts:
##     * panel letters (a, b, c...) on every panel
##     * Fig3 battery split: cohort-removal scenarios (robust-core gate) vs
##       robustness contrasts (reported, not gated)
##     * Fig5b annotated: gene level rho = 0.324 | 391 Tier-1 | 80.3%
##     * Fig7 subtitles take FDR from S8_01_hypotheses.csv (0.312, machine-printed)
##     * Fig6 lineage label "Smooth muscle (mural)"; FigS4 caption maps
##       Perivascular = mural (smooth muscle); Stromal = mesenchymal
##     * FigS2/S3/S4/S5 assembled from existing PNG rasters (source objects
##       not in results/) -- black-box-free plates; FigS1 needs her y-label
##       rerun first (see NOTE below)
##   NOTE on FigS1: edit make_figS1_exploratory_drugs.R line
##       y = "paired log2 FC, dienogest + cAMP vs vehicle",
##     to
##       y = "combined dienogest + dibutyryl-cAMP (decidualization) vs vehicle",
##     rerun it, THEN run this script.
## Usage (project root):
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/make_journal_plates.R
## ===========================================================================
set.seed(42)
source("scripts/S4_common.R")

OUT <- "FIGURES PLATES FOR JOURNAL SUBMISSION"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIGDIR <- "figures"

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("FATAL: ggplot2 required (in renv lock).", call. = FALSE)
suppressPackageStartupMessages({
  library(ggplot2); library(grid)
})
have_grid <- requireNamespace("gridExtra", quietly = TRUE)
have_png  <- requireNamespace("png", quietly = TRUE)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
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
as_flag <- function(d, col, p) {
  if (is.logical(d[[col]])) return(d[[col]])
  v <- toupper(trimws(as.character(d[[col]])))
  if (any(!v %in% c("TRUE", "FALSE")))
    stop("FATAL: ", p, " column ", col, " is not logical/TRUE/FALSE (R2/R9)",
         call. = FALSE)
  v == "TRUE"
}
save_plate <- function(plot, name, w, h) {
  ggplot2::ggsave(file.path(OUT, paste0(name, ".png")), plot, width = w,
                  height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(OUT, paste0(name, ".pdf")), plot, width = w,
                  height = h, device = "pdf", bg = "white")
  say("wrote '", OUT, "/", name, ".png/.pdf'")
}
tg <- function(p, letter) {
  p + labs(tag = letter) +
    theme(plot.tag = element_text(face = "bold", size = 11),
          plot.tag.position = c(0.01, 0.99))
}
CORE_COLS <- c("gene","k","M","SE","p","tau2","I2","Q","QEp","k_pos","k_neg",
               "FDR","concord_k","tier")
r2 <- function(x) sprintf("%.2f", x)
r3 <- function(x) sprintf("%.3f", x)
pct1 <- function(x) sprintf("%.1f%%", 100 * x)

## locked palette
C_DOWN  <- "#1F4E79"; C_UP <- "#C00000"; C_NS <- "#BFBFBF"
C_ACC <- "#2E7D32"; C_WARN <- "#C55A11"; C_PURP <- "#5B2D8E"; C_FAIL <- "#8B1A1A"
INK <- "#262626"

## compact journal theme (panels ~2.4-3.6in wide)
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(axis.line  = element_line(color = INK, linewidth = 0.35),
          axis.ticks = element_line(color = INK, linewidth = 0.35),
          axis.title = element_text(color = INK, face = "bold", size = base),
          axis.text  = element_text(color = INK, size = base - 1),
          plot.title = element_text(color = INK, face = "bold",
                                    size = base + 1.5, hjust = 0,
                                    margin = margin(b = 1.5)),
          plot.subtitle = element_text(color = "#595959", size = base - 1,
                                       hjust = 0, margin = margin(b = 2),
                                       lineheight = 0.95),
          plot.margin = margin(5, 7, 3, 5),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = base - 1),
          legend.key.size = grid::unit(2.6, "mm"),
          legend.margin = margin(b = 0),
          strip.background = element_rect(fill = "#EFEFEF", color = "#BFBFBF",
                                          linewidth = 0.25),
          strip.text = element_text(face = "bold", size = base - 0.5))
}

## raster panel (for figures whose source objects are not in results/)
raster_panel <- function(file, letter = NULL, title = NULL) {
  need(file)
  if (!have_png)
    stop("FATAL: package 'png' required to embed raster panel ", file,
         " -- install.packages('png') (R2).", call. = FALSE)
  img <- png::readPNG(file)
  g <- grid::rasterGrob(img, interpolate = TRUE)
  labs_g <- if (!is.null(letter)) {
    grid::textGrob(letter, x = unit(1, "mm"), y = unit(1, "mm"),
                   just = c("left", "bottom"),
                   gp = gpar(fontface = "bold", fontsize = 11))
  } else NULL
  tit_g <- if (!is.null(title)) {
    grid::textGrob(title, x = unit(1, "mm"), y = unit(0, "npc"),
                   just = c("left", "bottom"),
                   gp = gpar(fontsize = 8, fontface = "bold", col = INK))
  } else NULL
  gridExtra::arrangeGrob(
    ggplot2::ggplot() + ggplot2::annotation_custom(g) +
      ggplot2::geom_blank() + ggplot2::theme_void() +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                     size = 8, hjust = 0, margin = ggplot2::margin(b = 2))),
    labs_g, ncol = 1, heights = grid::unit(c(1, 0), c("null", "null")))
}
arr2 <- function(grobs, ncol, widths_in, heights_in, name, w_pad = 0.1, h_pad = 0.1) {
  if (!have_grid)
    stop("FATAL: gridExtra required for plate assembly (R2).", call. = FALSE)
  p <- gridExtra::arrangeGrob(grobs = grobs, ncol = ncol,
                              widths = grid::unit(widths_in, "in"),
                              heights = grid::unit(heights_in, "in"),
                              padding = grid::unit(2, "mm"))
  save_plate(p, name, sum(widths_in) + w_pad, sum(heights_in) + h_pad)
}

say("== make_journal_plates v1.0: presentation only, zero new computation (A5) ==")

## ---- locked inputs ---------------------------------------------------------
f_t1    <- file.path(CFG4$out_meta, "S4_meta_Tier1_core.csv")
f_all   <- file.path(CFG4$out_meta, "S4_meta_all_genes.csv")
f_sens  <- file.path(CFG4$out_meta, "S4_sensitivity_summary.csv")
f_bulk5 <- "results/S4_localization/S4_05_axis_bulk_map.csv"
f_scr5  <- "results/S4_localization/S4_05_axis_scrna_map.csv"
f_pmet  <- "results/S5_meta/S5_meta_pathways.csv"
f_pconc <- "results/S5_meta/S5_concordance_GSE11691.csv"
f_psens <- "results/S5_meta/S5_sensitivity_scenarios.csv"
f_core  <- "results/S8_perturb/S8_00_core_gene_coverage.csv"
f_cov   <- "results/S8_perturb/S8_00_coverage.csv"
f_hyp   <- "results/S8_perturb/S8_01_hypotheses.csv"
f_s6t   <- "results/S6_meta/S6_01_tests.csv"
f_s6s   <- "results/S6_meta/S6_01_scores_long.csv"
f_s6a   <- "results/S2_qc_GSE213216_audit.csv"

t1    <- read_locked(f_t1);   check_cols(t1, CORE_COLS, f_t1)
M     <- read_locked(f_all);  check_cols(M, CORE_COLS, f_all)
sens  <- read_locked(f_sens)
bulk5 <- read_locked(f_bulk5); check_cols(bulk5, c("axis","gene","contrast_id","logFC","sig_fdr05"), f_bulk5)
scr5  <- read_locked(f_scr5);  check_cols(scr5, c("axis","gene","dataset","lineage","arm","logFC","sig_fdr05"), f_scr5)
pm    <- read_locked(f_pmet);  check_cols(pm, c("pathway","M","SE","I2","FDR"), f_pmet)
pconc <- read_locked(f_pconc); check_cols(pconc, c("pathway","M_base","logFC_GSE11691"), f_pconc)
psens <- read_locked(f_psens); check_cols(psens, c("scenario","FDR"), f_psens)
core  <- read_locked(f_core);  check_cols(core, c("gene","M","in_L1000","is_landmark"), f_core)
cov   <- read_locked(f_cov);   check_cols(cov, c("measure","n"), f_cov)
hyp   <- read_locked(f_hyp);   check_cols(hyp, c("hypothesis","test","n","p_value","BH_FDR","verdict"), f_hyp)
s6t   <- read_locked(f_s6t)
s6s   <- read_locked(f_s6s);   check_cols(s6s, c("patient","lineage","tissue","tier1_signature","SERPINE1_logCPM"), f_s6s)
audit <- read_locked(f_s6a);   check_cols(audit, c("object","cells"), f_s6a)
bulk5$sig_fdr05 <- as_flag(bulk5, "sig_fdr05", f_bulk5)
scr5$sig_fdr05  <- as_flag(scr5,  "sig_fdr05", f_scr5)
hyp_txt <- function(h)
  sprintf("%s (FDR = %.3g)", h, hyp$BH_FDR[hyp$hypothesis == h][1])

## ===========================================================================
## FIG 1 -- study overview (2x2, letters a-d)
## ===========================================================================
## a: dataset map (compact lanes)
n_cells <- audit$cells[audit$object == "auxiliary.seurat.shared.rds"]
n_pat <- length(unique(s6s$patient)); n_pb <- nrow(s6s)
dm <- rbind(
  data.frame(x = 1:6, y = 4,
             gse = c("GSE7305","GSE6364","GSE25628","GSE7307","GSE51981","GSE120103"),
             nn = c("10v10","21v16","16v6","18v23","77v34","13v9"),
             role = "bulk meta set (k_min = 4)"),
  data.frame(x = 1, y = 3, gse = "GSE11691", nn = "9 women",
             role = "paired arm (never pooled)"),
  data.frame(x = 1:3, y = 2, gse = c("GSE179640","GSE183837","GSE214411"),
             nn = c("scRNA","scRNA","scRNA"), role = "scRNA localization"),
  data.frame(x = 1, y = 1, gse = "GSE213216",
             nn = sprintf("%s cells\n%d patients\n%d pseudobulks",
                          format(n_cells, big.mark = ","), n_pat, n_pb),
             role = "validation atlas"))
role_lvls <- c("bulk meta set (k_min = 4)","paired arm (never pooled)",
               "scRNA localization","validation atlas")
dm$role <- factor(dm$role, levels = role_lvls)
role_cols <- c("bulk meta set (k_min = 4)" = C_DOWN,
               "paired arm (never pooled)" = C_WARN,
               "scRNA localization" = C_ACC, "validation atlas" = C_PURP)
dm$w <- 0.88; dm$w[dm$role == "validation atlas"] <- 2.9
fig1a <- ggplot(dm, aes(x, y)) +
  geom_tile(aes(fill = role, width = w), height = 0.66, color = "white",
            linewidth = 0.6, show.legend = TRUE) +
  geom_text(aes(y = y + 0.11, label = gse), size = 1.9, fontface = "bold",
            color = "white", show.legend = FALSE) +
  geom_text(aes(y = y - 0.13, label = nn), size = 1.5, lineheight = 0.9,
            color = "white", show.legend = FALSE) +
  scale_fill_manual(values = role_cols) +
  scale_y_continuous(breaks = 4:1, labels = c("Bulk meta set","Paired arm",
                                              "scRNA localization","Validation atlas"),
                     expand = c(0.03, 0.03)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  coord_cartesian(xlim = c(-0.8, 6.6), ylim = c(0.55, 4.45)) +
  labs(title = "Dataset map",
       subtitle = "six bulk cohorts (meta set); one paired arm; three scRNA localization datasets; one independent validation atlas") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.y = element_text(size = 5.5),
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.key.size = unit(2, "mm"))

## b: locked contrasts
coh <- data.frame(
  gse = c("GSE7305","GSE6364","GSE25628","GSE7307","GSE51981","GSE120103","GSE11691"),
  nn  = c("10 v 10","21 v 16","16 v 6","18 v 23","77 v 34","13 v 9","9 women"),
  cmp = c("Endo vs\nNormal","Eutopic vs Normal\n(+ phase)","Lesion/Eutopic\nvs Normal",
          "Endo vs\nNormal","Endo (all stages)\nvs Control","Stage IV vs\nFertile control",
          "Ectopic vs Eutopic\n(paired)"),
  role = c(rep("meta set (k_min = 4)", 6), "paired arm (never pooled)"))
coh$x <- seq_len(nrow(coh))
fig1b <- ggplot(coh, aes(x, 1)) +
  geom_tile(aes(fill = role), height = 0.6, width = 0.86, color = "white",
            linewidth = 0.6, show.legend = TRUE) +
  geom_text(aes(y = 1.10, label = gse, color = role), size = 1.9,
            fontface = "bold", show.legend = FALSE) +
  geom_text(aes(y = 0.88, label = nn, color = role), size = 1.5,
            show.legend = FALSE) +
  scale_color_manual(values = c("meta set (k_min = 4)" = "white",
                                "paired arm (never pooled)" = INK)) +
  scale_fill_manual(values = c("meta set (k_min = 4)" = C_DOWN,
                               "paired arm (never pooled)" = C_WARN)) +
  scale_x_continuous(breaks = coh$x, labels = coh$cmp, expand = c(0.01, 0.01)) +
  coord_cartesian(ylim = c(0.55, 1.45)) +
  labs(title = "Locked contrasts: six cohorts + one paired arm",
       subtitle = "six disease-vs-control cohorts enter the random-effects meta; the paired arm is never pooled") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.x = element_text(size = 4.3, lineheight = 0.9),
        legend.position = "bottom", legend.text = element_text(size = 5),
        legend.key.size = unit(2, "mm"))

## c: pipeline (2 rows x 3)
pl <- data.frame(
  x = c(1, 2, 3, 3, 2, 1),
  y = c(2, 2, 2, 1, 1, 1),
  lbl = c("S3\n19 frozen\ncontrasts",
          "S4\ngene-level REML\nTier-1 + robust core",
          "S4_05 + S5\nlocalization +\npathway meta",
          "S6\nexternal validation\n(GSE213216)",
          "S8\nLINCS reversal\nanalysis",
          "S8b\nexploratory drug\nresponse"),
  fill = c("#DCE6F1","#B8CCE4","#8FAADC","#F2DCB3","#E4D3EE","#F0E3C8"))
pl$x <- as.numeric(pl$x); pl$y <- as.numeric(pl$y)
fig1c <- ggplot(pl, aes(x, y)) +
  geom_tile(aes(fill = fill), height = 0.72, width = 0.82,
            color = "#7F7F7F", linewidth = 0.3) +
  geom_text(aes(label = lbl), size = 1.55, lineheight = 0.9, color = INK) +
  geom_segment(data = data.frame(x = c(1.43, 2.43, 3.43, 3.43, 2.43, 1.43),
                                 y = c(2, 2, 2, 1.36, 1, 1),
                                 xend = c(1.57, 2.57, 3.57, 2.57, 1.57, 0.57),
                                 yend = c(2, 2, 1.64, 1, 1, 1)),
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed"),
               color = "#595959", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  coord_cartesian(xlim = c(0.4, 3.6), ylim = c(0.45, 2.55)) +
  labs(title = "Locked analysis pipeline") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

## d: analytical branches (vertical)
br <- data.frame(
  x = 1, y = 4:1,
  lbl = c("Pseudobulk localization\n3 scRNA datasets (8 lineages)",
          "External validation\nGSE213216 (21 patients)",
          "Perturbational reversal\nLINCS L1000 (SERPINE1-KD, estradiol)",
          "Exploratory decidualization\ndrug response (S8b)"),
  fill = c("#B8CCE4", "#F2DCB3", "#E4D3EE", "#F0E3C8"))
fig1d <- ggplot(br, aes(x, y)) +
  geom_tile(aes(fill = fill), height = 0.74, width = 0.9,
            color = "#7F7F7F", linewidth = 0.3) +
  geom_text(aes(label = lbl), size = 1.55, lineheight = 0.9, color = INK) +
  geom_segment(data = data.frame(x = 1, xend = 1, y = c(3.62, 2.62, 1.62),
                                 yend = c(3.36, 2.36, 1.36)),
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed"),
               color = "#595959", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0.4, 1.6), ylim = c(0.4, 4.6)) +
  labs(title = "Analytical branches",
       subtitle = "discovery to validation to reversal; BH-FDR < 0.05 is the sole claim gate; exploratory analyses carry no claim") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

arr2(list(tg(fig1a, "a"), tg(fig1b, "b"), tg(fig1c, "c"), tg(fig1d, "d")),
     ncol = 2, widths_in = c(3.6, 3.6), heights_in = c(3.1, 2.4),
     name = "Fig1", w_pad = 0.15, h_pad = 0.15)

## ===========================================================================
## FIG 2 -- Tier-1 core (1x3: a volcano, b top genes, c heterogeneity)
## ===========================================================================
Mv <- M
Mv$neglog <- ifelse(Mv$FDR > 0, -log10(Mv$FDR),
                    max(-log10(Mv$FDR[Mv$FDR > 0]), na.rm = TRUE) + 0.5)
Mv$tier <- factor(Mv$tier, levels = c("Tier3", "Tier2", "Tier1"))
fig2a <- ggplot(Mv, aes(x = M, y = neglog, color = tier)) +
  geom_point(size = 0.35, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = c(Tier3 = C_NS, Tier2 = "#7A9CC4", Tier1 = C_UP),
                     breaks = c("Tier1", "Tier2", "Tier3")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = INK,
             linewidth = 0.35) +
  annotate("text", x = -Inf, y = -log10(0.05), label = " FDR = 0.05",
           hjust = 0, vjust = -0.4, size = 1.8, color = INK, fontface = "bold") +
  labs(x = "REML meta effect size M (log2 FC)", y = expression(-log[10] ~ FDR),
       title = sprintf("Cross-cohort random-effects meta (%s genes)",
                       format(nrow(Mv), big.mark = ",")),
       subtitle = sprintf("Tier-1 core: %d genes (FDR < 0.05, |M| >= 0.5, same sign in >= 4 cohorts)",
                          nrow(t1))) +
  theme_pub() +
  theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "#BFBFBF",
                                         linewidth = 0.25),
        legend.text = element_text(size = 4.5),
        legend.key.size = unit(2, "mm"))

t1d <- t1[order(-abs(t1$M)), ][seq_len(min(12, nrow(t1))), ]
t1d$gene <- factor(t1d$gene, levels = rev(t1d$gene))
fig2b <- ggplot(t1d, aes(x = gene, y = M)) +
  geom_hline(yintercept = 0, color = "#7F7F7F", linewidth = 0.35) +
  geom_errorbar(aes(ymin = M - 1.96 * SE, ymax = M + 1.96 * SE),
                width = 0.22, color = "#595959", linewidth = 0.4) +
  geom_point(aes(color = M > 0), size = 1.6) +
  scale_color_manual(values = c(`TRUE` = C_UP, `FALSE` = C_DOWN),
                     breaks = c(TRUE, FALSE), labels = c("up", "down")) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "M (95% CI)",
       title = "Strongest Tier-1 genes by |M|",
       subtitle = "display ranking only; full statistics in S4_meta_Tier1_core.csv") +
  theme_pub() +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "#BFBFBF",
                                         linewidth = 0.25),
        legend.text = element_text(size = 4.5), legend.key.size = unit(2, "mm"),
        axis.text.y = element_text(size = 5, face = "italic"))

med_i2 <- stats::median(M$I2, na.rm = TRUE)
fig2c <- ggplot(M, aes(x = I2)) +
  geom_histogram(bins = 30, fill = "#D9D9D9", color = "white", linewidth = 0.25) +
  geom_vline(xintercept = med_i2, color = C_PURP, linewidth = 0.6) +
  annotate("text", x = med_i2, y = Inf,
           label = paste0("median I2 = ", round(med_i2, 1), "% "),
           vjust = 1.6, hjust = 1, size = 1.9, color = C_PURP, fontface = "bold") +
  labs(x = expression(I^2 ~ "(%)"), y = "genes",
       title = "Between-cohort heterogeneity",
       subtitle = "reported for every gene and never filtered (locked rule)") +
  theme_pub()

arr2(list(tg(fig2a, "a"), tg(fig2b, "b"), tg(fig2c, "c")),
     ncol = 3, widths_in = c(2.6, 2.3, 2.1), heights_in = c(3.0),
     name = "Fig2", w_pad = 0.15, h_pad = 0.12)

say("Fig1 + Fig2 plates written")

## ===========================================================================
## FIG 3 -- robust-core stress test (2x2: a retention, b battery, c stability, d core)
## ===========================================================================
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
ret <- ret[order(ret$grp, -ret$retained), ]
ret$scenario <- factor(ret$scenario, levels = ret$scenario)
fig3a <- ggplot(ret, aes(scenario, retained, fill = grp)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_hline(yintercept = length(base_set), linetype = "dashed", color = C_UP,
             linewidth = 0.5) +
  annotate("text", x = 1, y = length(base_set), vjust = -0.5, hjust = 0,
           size = 1.9, color = C_UP, fontface = "bold",
           label = paste0("base Tier-1 = ", length(base_set))) +
  geom_text(aes(label = retained), vjust = -0.35, size = 1.9, color = INK) +
  scale_fill_manual(values = c(sans = "#8FAADC", LOCO = C_DOWN,
                               `FE vs RE` = C_PURP, `QC swap` = C_ACC)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Tier-1 genes retained",
       title = "Tier-1 set stability across scenarios",
       subtitle = "leave-one-cohort-out, sans-G7307, fixed-effects contrast, frozen QC-swap replications") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 4.6),
        legend.position = "right", legend.text = element_text(size = 4.5),
        legend.key.size = unit(2, "mm"))

## b: sensitivity battery -- SPLIT into two blocks
sens$val_disp <- ifelse(grepl("|", sens$value, fixed = TRUE), sens$value,
                 ifelse(grepl("retention", sens$metric),
                        pct1(as.numeric(sens$value)),
                        r3(as.numeric(sens$value))))
blk1 <- sens[seq_len(7), ]; blk2 <- sens[8:nrow(sens), ]
n1 <- nrow(blk1); n2 <- nrow(blk2)
blk1$y <- (n1 + 2):4; blk2$y <- 2:(3 - n2)
fig3b <- ggplot() +
  geom_point(data = blk1, aes(x = 1, y = y),
             color = ifelse(blk1$pass, C_ACC, C_FAIL), size = 2) +
  geom_text(data = blk1,
            aes(x = 1.05, y = y,
                label = paste0(test, ": ", metric, " = ", val_disp,
                               " (", criterion, ")")),
            hjust = 0, size = 1.9, color = INK) +
  geom_point(data = blk2, aes(x = 1, y = y), color = C_WARN, size = 2) +
  geom_text(data = blk2,
            aes(x = 1.05, y = y,
                label = paste0(test, ": ", metric, " = ", val_disp,
                               " (", criterion, ")")),
            hjust = 0, size = 1.9, color = INK) +
  annotate("text", x = 0.9, y = n1 + 3.3, hjust = 0, size = 2.1,
           fontface = "bold", color = INK,
           label = "Cohort-removal scenarios (robust-core gate)") +
  annotate("text", x = 0.9, y = 3.2, hjust = 0, size = 2.1,
           fontface = "bold", color = C_WARN,
           label = "Robustness contrasts (reported, not gated)") +
  annotate("segment", x = 0.9, xend = 3.8, y = 3.5, yend = 3.5,
           linewidth = 0.3, color = "grey70") +
  annotate("text", x = 0.9, y = 0.2, hjust = 0, size = 1.7, color = "#595959",
           label = "Robust-core membership requires retention under all seven cohort-removal scenarios; per-scenario floors are not analysis failures.") +
  scale_x_continuous(limits = c(0.85, 4.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, n1 + 4), expand = c(0, 0)) +
  labs(title = "Locked sensitivity battery",
       subtitle = "PASS / FAIL against pre-registered criteria; swaps and boundary tests are reported, never hidden") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

## d: robust core summary
core <- core[order(-core$M), ]
core$rank <- seq_len(nrow(core))
core$dir <- ifelse(core$M > 0, "up-regulated", "down-regulated")
slots <- data.frame(
  gene = c("FOSB", "FOS", "EGR1", "C1QB", "IL6"),
  lx = c(2.0, 1.2, 3.8, 5.2, 11.4),
  ly = c(2.60, 2.20, 1.88, 0.70, 0.95),
  hj = c(0, 1, 0, 1, 0))
lab <- merge(core, slots, by = "gene")
fig3d <- ggplot(core, aes(rank, M)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = dir), size = 1.2) +
  geom_point(data = core[core$is_landmark == TRUE, ], size = 2.2, shape = 21,
             fill = NA, color = "black", stroke = 0.45) +
  geom_text(data = lab, aes(x = lx, y = ly, label = gene, hjust = hj),
            size = 2.0, color = INK, fontface = "italic") +
  scale_color_manual(values = c("up-regulated" = C_UP, "down-regulated" = C_DOWN)) +
  labs(title = sprintf("Robust core: %d genes (up %d, down %d) retained under all seven cohort-removal scenarios",
                       nrow(core), sum(core$M > 0), sum(core$M < 0)),
       subtitle = paste0("six leave-one-cohort-out reanalyses plus the sans-GSE7307 boundary test;\n",
                         "rings = the 8 L1000 landmark genes; up module = AP-1 / immediate-early response\n",
                         "(FOS, FOSB, EGR1, JUNB, IER2) plus IL6, SOCS3 and C1QB."),
       x = "robust-core genes ranked by pooled M", y = "pooled M (log2 FC)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.text = element_text(size = 4.5),
        legend.key.size = unit(2, "mm"))

## a: robust-core retention matrix (66 x 7 cohort-removal scenarios, all retained)
scen7 <- ret$scenario[ret$grp %in% c("sans", "LOCO")]
mat <- do.call(rbind, lapply(sens_files, function(p) {
  d <- read.csv(p, stringsAsFactors = FALSE)
  tag <- sub("^S4_sens_", "", sub("_all\\.csv$", "", basename(p)))
  if (!tag %in% scen7) return(NULL)
  data.frame(scenario = tag, gene = d$gene[d$tier == "Tier1"],
             stringsAsFactors = FALSE)
}))
mat$kept <- mat$gene %in% core$gene
mat$gene <- factor(mat$gene, levels = core$gene[order(-core$M)])
mat$scenario <- factor(mat$scenario, levels = scen7)
fig3m <- ggplot(mat, aes(scenario, gene)) +
  geom_tile(aes(fill = kept), color = "white", linewidth = 0.1) +
  scale_fill_manual(values = c("TRUE" = C_ACC, "FALSE" = C_FAIL)) +
  annotate("text", x = length(scen7) + 0.6,
           y = nrow(core) * 0.5, hjust = 0, size = 2.0, color = INK,
           label = sprintf("%d of %d robust-core genes\nretained in every scenario\n(7 of 7 required)", nrow(core), nrow(core))) +
  coord_cartesian(xlim = c(0.5, length(scen7) + 4.2), clip = "off") +
  labs(title = "Robust-core retention matrix",
       subtitle = "rows = 66 robust-core genes (ranked by M); columns = the seven cohort-removal scenarios; all retained") +
  theme_pub() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 4.3, angle = 30, hjust = 1),
        legend.position = "none")

arr2(list(tg(fig3m, "a"), tg(fig3b, "b"), tg(fig3a, "c"), tg(fig3d, "d")),
     ncol = 2, widths_in = c(3.6, 3.6), heights_in = c(3.0, 3.0),
     name = "Fig3", w_pad = 0.3, h_pad = 0.15)


say("Fig3 plate written")

## ===========================================================================
## FIG 4 -- functional-axis landscape (1x2: a bulk heatmap, b scRNA dots)
## ===========================================================================
AXES_ORD <- c("cytotoxicity", "antigen_presentation", "senescence_dormancy",
              "proliferation", "stromal_ecm")
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
  geom_tile(aes(fill = fillv), color = "white", linewidth = 0.3) +
  geom_tile(data = b5[b5$sig_fdr05, ], aes(x = contrast_id, y = gene),
            fill = NA, color = INK, linewidth = 0.45) +
  scale_fill_gradient2(low = C_DOWN, mid = "white", high = C_UP, midpoint = 0,
                       limits = c(-2, 2), name = "log2 FC (clipped)") +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(axis = AXES_LBL)) +
  labs(x = NULL, y = NULL,
       title = "Axis genes across the locked primary contrasts (bulk)",
       subtitle = "tile color = log2 FC; black border = BH-FDR < 0.05; last column = paired arm (never pooled)") +
  theme_pub() +
  theme(panel.spacing.y = unit(1.8, "mm"), panel.spacing.x = unit(1, "mm"),
        axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1, size = 5,
                                   face = "bold"),
        axis.text.y = element_text(size = 4.6),
        strip.text.y.left = element_text(size = 5.5, lineheight = 0.85),
        legend.position = "right",
        legend.title = element_text(size = 5.5, face = "bold"),
        legend.text = element_text(size = 4.5), legend.key.size = unit(2, "mm"),
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
  scale_size_continuous(range = c(1.6, 4.0), breaks = c(1.5, 2, 3),
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
  theme(panel.spacing.y = unit(1.8, "mm"),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5),
        axis.text.y = element_text(size = 4.6),
        strip.text.y.left = element_text(size = 5.5, lineheight = 0.85),
        legend.position = "right",
        legend.title = element_text(size = 5.5, face = "bold"),
        legend.text = element_text(size = 4.5), legend.key.size = unit(2, "mm"),
        strip.placement = "outside")

arr2(list(tg(fig4a, "a"), tg(fig4b, "b")),
     ncol = 2, widths_in = c(3.6, 3.6), heights_in = c(5.2),
     name = "Fig4", w_pad = 0.15, h_pad = 0.12)

## ===========================================================================
## FIG 5 -- pathway meta + paired concordance (2x2: a forest, b concordance,
##          c sensitivity matrix, d pathway I2)
## ===========================================================================
pm$pathway <- factor(pm$pathway, levels = pm$pathway[order(pm$M, decreasing = TRUE)])
fig5a <- ggplot(pm, aes(x = pathway, y = M)) +
  geom_hline(yintercept = 0, color = "#7F7F7F", linewidth = 0.35) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", color = "#A6A6A6",
             linewidth = 0.35) +
  geom_errorbar(aes(ymin = M - 1.96 * SE, ymax = M + 1.96 * SE),
                width = 0.24, color = "#595959", linewidth = 0.45) +
  geom_point(aes(color = FDR < 0.05), size = 2.0) +
  geom_text(aes(y = M + 1.96 * SE, label = paste0("I2 = ", round(I2))),
            vjust = -0.5, size = 1.7, color = "#595959") +
  scale_color_manual(values = c(`TRUE` = C_UP, `FALSE` = "#8C8C8C"),
                     breaks = c(TRUE, FALSE), labels = c("BH-FDR < 0.05", "n.s.")) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "pathway score log2 FC (95% CI)",
       title = "Pathway-level random-effects meta (six cohorts)",
       subtitle = "dotted lines = interpretive |M| = 0.5 (never a filter); BH-FDR across six pathway tests") +
  theme_pub() +
  theme(legend.position = "bottom", legend.text = element_text(size = 4.5),
        legend.key.size = unit(2, "mm"), legend.margin = margin(t = -3),
        axis.text.y = element_text(size = 5))

xr <- diff(range(pconc$M_base, na.rm = TRUE))
pconc$lx <- pconc$M_base
pconc$lx[pconc$M_base > 0] <- pconc$M_base[pconc$M_base > 0] + 0.05 * xr
pconc$lx[pconc$M_base <= 0] <- pconc$M_base[pconc$M_base <= 0] - 0.05 * xr
pconc$hj <- ifelse(pconc$M_base > 0, 0, 1)
pconc$lbl <- as.character(pconc$pathway)
w <- AXES_LBL[pconc$lbl]
pconc$lbl[!is.na(w)] <- w[!is.na(w)]
fig5b <- ggplot(pconc, aes(x = M_base, y = logFC_GSE11691)) +
  geom_hline(yintercept = 0, color = "#D9D9D9", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#D9D9D9", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#8C8C8C",
              linewidth = 0.4) +
  geom_point(size = 2.2, color = C_DOWN) +
  geom_text(aes(x = lx, label = lbl, hjust = hj), size = 2.0, color = INK,
            fontface = "bold", vjust = -0.5, lineheight = 0.9) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.6, size = 1.9,
           color = INK, fontface = "bold",
           label = "gene level: rho = 0.324 | 391 Tier-1 genes | 80.3% sign agreement\n(paired arm never pooled into the meta)") +
  scale_x_continuous(expand = expansion(mult = 0.42)) +
  scale_y_continuous(expand = expansion(mult = 0.30)) +
  labs(x = "M (six-cohort base meta)", y = "paired ectopic-vs-eutopic log2 FC",
       title = "GSE11691 paired concordance arm",
       subtitle = "directional concordance only") +
  theme_pub()

fdr_mat <- do.call(rbind, lapply(strsplit(psens$FDR, "; ", fixed = TRUE),
                                 function(x) {
  v <- strsplit(x, "=", fixed = TRUE)
  data.frame(pathway = vapply(v, `[[`, character(1), 1),
             FDR = as.numeric(vapply(v, `[[`, character(1), 2)),
             stringsAsFactors = FALSE)
}))
fdr_mat$scenario <- rep(psens$scenario,
                        times = vapply(strsplit(psens$FDR, "; ", fixed = TRUE),
                                       length, integer(1)))
PATH_ORD <- pm$pathway[order(pm$M, decreasing = TRUE)]
SC_ORD <- c("sans-G7307", paste0("LOCO-", CFG4$meta_set),
            paste0("swap-", CFG4$qc_swap), "FE-vs-RE")
fdr_mat$pathway <- factor(fdr_mat$pathway, levels = PATH_ORD)
fdr_mat$scenario <- factor(fdr_mat$scenario, levels = SC_ORD)
if (any(is.na(fdr_mat$FDR)) || anyNA(fdr_mat$pathway) || anyNA(fdr_mat$scenario))
  stop("FATAL: scenario/pathway labels in S5_sensitivity_scenarios.csv do not ",
       "match the locked meta_set / qc_swap / pathway names (R2).", call. = FALSE)
fdr_mat$nl <- pmin(-log10(fdr_mat$FDR), 3)
fdr_mat$txt <- ifelse(fdr_mat$FDR < 0.001,
                      format(fdr_mat$FDR, digits = 1, scientific = TRUE),
                      r3(fdr_mat$FDR))
fig5c <- ggplot(fdr_mat, aes(x = scenario, y = pathway)) +
  geom_tile(aes(fill = nl), color = "white", linewidth = 0.4) +
  geom_text(aes(label = txt, color = nl > 1.7), size = 1.7, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient(low = "white", high = C_UP, limits = c(0, 3),
                      breaks = c(0, 1.3, 3), labels = c("0", "1.3", ">=3"),
                      name = expression(-log[10] ~ FDR)) +
  labs(x = NULL, y = NULL,
       title = "Sensitivity scenarios at pathway level",
       subtitle = "cell text = BH-FDR (six pathway tests)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, size = 4.2),
        axis.text.y = element_text(size = 5),
        legend.position = "right",
        legend.title = element_text(size = 5.5, face = "bold"),
        legend.text = element_text(size = 4.5), legend.key.size = unit(2, "mm"))

pm2 <- pm
pm2$pathway <- factor(pm2$pathway, levels = PATH_ORD)
fig5d <- ggplot(pm2, aes(x = pathway, y = I2)) +
  geom_col(width = 0.68, fill = C_PURP, color = "white", linewidth = 0.25) +
  geom_text(aes(label = round(I2)), hjust = -0.25, size = 1.9, color = INK) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression(I^2 ~ "(%)"),
       title = "Between-cohort heterogeneity per pathway",
       subtitle = "reported, never filtered (locked rule)") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 5))

arr2(list(tg(fig5a, "a"), tg(fig5b, "b"), tg(fig5c, "c"), tg(fig5d, "d")),
     ncol = 2, widths_in = c(3.6, 3.6), heights_in = c(3.4, 3.4),
     name = "Fig5", w_pad = 0.15, h_pad = 0.15)

say("Fig4 + Fig5 plates written")

## ===========================================================================
## FIG 6 -- external validation (a design strip spanning; b direction heatmap;
##          c mural hit; d SERPINE1 non-replication)
## ===========================================================================
S6_SCORES <- c("cytotoxicity","antigen_presentation","senescence_dormancy",
               "proliferation","stromal_ecm","SenMayo","tier1_signature")
c1t <- s6t[s6t$contrast == "C1" & s6t$status == "TESTED" &
             s6t$score %in% S6_SCORES, ]
if (!nrow(c1t))
  stop("FATAL: no TESTED C1 rows in ", f_s6t, " (R2)", call. = FALSE)
## unified lineage nomenclature: mural = smooth muscle
LIN_SHORT <- c("Mesenchymal cells" = "Mesenchymal",
               "Endothelial cells" = "Endothelial",
               "Smooth muscle cells" = "Smooth muscle (mural)",
               "T/NK cells" = "T/NK", "B/Plasma cells" = "B/Plasma",
               "Myeloid cells" = "Myeloid", "Epithelial cells" = "Epithelial",
               "Erythrocytes" = "Erythrocytes", "Mast cells" = "Mast")
if (any(is.na(LIN_SHORT[c1t$lineage])))
  stop("FATAL: unmapped lineage label in S6 tests (R2/R9)", call. = FALSE)
ord <- c1t[c1t$score == "tier1_signature", ]
ord <- ord$lineage[order(ord$delta_median, decreasing = TRUE)]
c1t$lin <- factor(LIN_SHORT[c1t$lineage], levels = LIN_SHORT[ord])
c1t$scr <- factor(gsub("_", " ", c1t$score), levels = gsub("_", " ", S6_SCORES))
c1t$mark <- ifelse(!is.na(c1t$FDR) & c1t$FDR < 0.05, "*",
                   ifelse(!is.na(c1t$p) & c1t$p < 0.05, "+", ""))
lim <- max(abs(c1t$delta_median), na.rm = TRUE)
fig6b <- ggplot(c1t, aes(x = lin, y = scr, fill = delta_median)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = mark), size = 3.6, color = INK, fontface = "bold") +
  scale_fill_gradient2(low = C_DOWN, mid = "white", high = C_UP, midpoint = 0,
                       limits = c(-lim, lim), name = NULL) +
  labs(x = NULL, y = NULL,
       title = "External validation in GSE213216: direction across lineages",
       subtitle = "C1: ectopic lesions vs eutopic endometrium, patient-level pseudobulk (Wilcoxon delta-median)\n* FDR < 0.05    + raw p < 0.05    (negative results shown with equal weight)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "right", legend.text = element_text(size = 4.5),
        legend.key.size = unit(2, "mm"))

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
                 round(tr$delta_median, 2), " | n = ", tr$n_case, " v ", tr$n_ctrl)
  sub2 <- verdict
  if (!is.na(tr$p_indep_patient))
    sub2 <- paste0(sub2, " | indep p = ", signif(tr$p_indep_patient, 3))
  ggplot(d, aes(x = grp, y = .data[[ycol]], color = grp)) +
    geom_point(position = position_jitter(width = 0.16, height = 0, seed = 42),
               size = 1.4, alpha = 0.9) +
    stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.25,
                 color = INK) +
    scale_color_manual(values = c("Eutopic" = C_DOWN, "Ectopic lesion" = C_UP),
                       guide = "none") +
    labs(x = NULL, y = ylab, title = ttl,
         subtitle = paste0(sub1, "\n", sub2)) +
    theme_pub() +
    theme(axis.text.x = element_text(size = 5))
}
fig6c <- mk_s6_panel("Smooth muscle cells", "tier1_signature",
                     "Tier-1 signature (z)",
                     "Validated: mural Tier-1",
                     "sole FDR<0.05 hit")
fig6d <- mk_s6_panel("Endothelial cells", "SERPINE1_logCPM",
                     "SERPINE1 logCPM",
                     "SERPINE1 - endothelium",
                     "pre-registered P1: NOT CONFIRMED")

## a: design strip (spans the top)
n_plan <- nrow(s6t)
n_not <- sum(grepl("not", tolower(s6t$status)))
n_ok <- n_plan - n_not
dsg <- data.frame(
  x = 1:4,
  lbl = c(sprintf("%s cells", format(n_cells, big.mark = ",")),
          sprintf("%d patients", n_pat),
          sprintf("%d pseudobulks", n_pb),
          sprintf("%d planned tests: %d testable | %d not testable",
                  n_plan, n_ok, n_not)),
  fill = c(C_PURP, C_PURP, C_PURP, C_FAIL))
fig6a <- ggplot(dsg, aes(x, 1)) +
  geom_tile(aes(fill = fill), height = 0.7, width = 0.92, color = "white",
            linewidth = 0.8) +
  geom_text(aes(label = lbl), color = "white", size = 2.1, lineheight = 0.95) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  ylim(0.5, 1.5) +
  labs(title = "External-validation design (GSE213216)",
       subtitle = "independent single-cell atlas; patient-level pseudobulks; verdicts machine-printed") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

p6 <- gridExtra::arrangeGrob(
  grobs = list(tg(fig6a, "a"), tg(fig6b, "b"), tg(fig6c, "c"), tg(fig6d, "d")),
  layout_matrix = rbind(c(1, 1), c(2, 3), c(2, 4)),
  widths = grid::unit(c(4.4, 2.6), "in"),
  heights = grid::unit(c(1.05, 1.55, 1.55), "in"),
  padding = grid::unit(2, "mm"))
save_plate(p6, "Fig6", 7.0 + 0.15, 4.15 + 0.15)

## ===========================================================================
## FIG 7 -- perturbational reversal: negative result (2x2)
## ===========================================================================
s7k <- read_locked("results/S8_perturb/S8_01_kd_scores.csv")
s7e <- read_locked("results/S8_perturb/S8_01_estradiol_scores.csv")
check_cols(s7k, c("arm","signatureid","n_overlap","r_tier1","cellline"),
           "results/S8_perturb/S8_01_kd_scores.csv")
check_cols(s7e, c("arm","signatureid","n_overlap","r_tier1"),
           "results/S8_perturb/S8_01_estradiol_scores.csv")
s7k <- s7k[!is.na(s7k$r_tier1), ]
s7e <- s7e[!is.na(s7e$r_tier1), ]
if (nrow(s7k) < 4L || nrow(s7e) < 4L)
  stop("FATAL: fig7 needs >=4 scored signatures per arm (R6).", call. = FALSE)

g_cov <- function(x) cov$n[cov$measure == x][1]
fun <- data.frame(
  x = factor(c("Tier-1 genes", "L1000 landmark", "Robust core (66)",
               "Core in L1000", "Core landmark"),
             levels = c("Tier-1 genes", "L1000 landmark", "Robust core (66)",
                        "Core in L1000", "Core landmark")),
  n = c(485, g_cov("tier1_total_landmark"), 66, g_cov("robust_core_in_L1000"),
        g_cov("robust_core_landmark")),
  grp = c("Tier-1", "Tier-1", "robust core", "robust core", "robust core"))
fig7a <- ggplot(fun, aes(x, n)) +
  geom_col(aes(fill = grp), width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 2.4, fontface = "bold",
            color = INK) +
  scale_fill_manual(values = c("Tier-1" = "#8FAADC", "robust core" = C_PURP)) +
  labs(title = "LINCS coverage audit",
       subtitle = sprintf("Tier-1: %d landmark genes (%d up / %d down of 485);\nrobust core: %d of 66 in L1000, %d landmark genes",
                          g_cov("tier1_total_landmark"), g_cov("tier1_up_landmark"),
                          g_cov("tier1_down_landmark"), g_cov("robust_core_in_L1000"),
                          g_cov("robust_core_landmark")),
       x = NULL, y = "genes") +
  theme_pub() + theme(legend.position = "none",
        axis.text.x = element_text(size = 4.8, lineheight = 0.9))

mk_s7_panel <- function(d, ylab, ttl, sub, hl = NA) {
  d <- d[order(d$r_tier1), ]
  d$lab <- factor(d[[ylab]], levels = d[[ylab]])
  d$dir <- factor(ifelse(d$r_tier1 < 0, "reversal direction", "lesion direction"),
                  levels = c("reversal direction", "lesion direction"))
  p <- ggplot(d, aes(x = r_tier1, y = lab)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_vline(xintercept = median(d$r_tier1), color = INK, linewidth = 0.25,
               linetype = "dashed") +
    geom_point(aes(color = dir), size = 1.8) +
    scale_color_manual(values = c("reversal direction" = C_UP,
                                  "lesion direction" = C_DOWN)) +
    labs(x = "r vs lesion Tier-1 signature (46 landmark genes)", y = NULL,
         title = ttl, subtitle = sub) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 5),
          legend.position = "bottom", legend.text = element_text(size = 4.5),
          legend.key.size = unit(2, "mm"), legend.margin = margin(t = -3))
  if (!is.na(hl) && hl %in% levels(d$lab)) {
    p <- p + geom_point(data = d[d$lab == hl, ], shape = 21, size = 3.6,
                        stroke = 1.1, color = C_WARN, fill = NA)
  }
  p
}
fig7b <- mk_s7_panel(s7k, "cellline", "SERPINE1 knockdown vs lesion program",
  paste0("7 cell lines | median r = ", round(median(s7k$r_tier1), 2), " | ",
         hyp_txt("P1a"), "\nHA1E (ringed): strongest single-line signal r = ",
         round(s7k$r_tier1[s7k$cellline == "HA1E"][1], 2),
         " -- hypothesis-generating only"), hl = "HA1E")
fig7c <- mk_s7_panel(s7e, "signatureid", "Estradiol vs lesion program",
  paste0(nrow(s7e), " cell lines (aggregated) | median r = ",
         round(median(s7e$r_tier1), 2), " | ", hyp_txt("P2a")))

hyp$txt_p <- sprintf("p = %.3g", hyp$p_value)
hyp$txt_fdr <- sprintf("FDR = %.3g", hyp$BH_FDR)
fig7d <- ggplot(hyp, aes(x = 1, y = seq(nrow(hyp), 1))) +
  geom_text(aes(x = 0.0, label = hypothesis), hjust = 0, size = 2.2,
            fontface = "bold", color = INK) +
  geom_text(aes(x = 0.33, label = test), hjust = 0, size = 1.7, color = "#595959") +
  geom_text(aes(x = 0.85, label = txt_p), hjust = 0, size = 1.9, color = INK) +
  geom_text(aes(x = 1.02, label = txt_fdr), hjust = 0, size = 1.9, color = INK) +
  geom_text(aes(x = 1.22, label = verdict), hjust = 0, size = 1.9,
            fontface = "bold", color = C_FAIL) +
  xlim(0, 1.45) + ylim(0.4, nrow(hyp) + 0.9) +
  labs(title = "Pre-registered reversal hypotheses: none confirmed",
       subtitle = "best BH-FDR = 0.129 (P2b); reported as a negative result") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.background = element_blank())

arr2(list(tg(fig7a, "a"), tg(fig7b, "b"), tg(fig7c, "c"), tg(fig7d, "d")),
     ncol = 2, widths_in = c(3.6, 3.6), heights_in = c(3.0, 2.2),
     name = "Fig7", w_pad = 0.15, h_pad = 0.15)

say("Fig6 + Fig7 plates written")

## ===========================================================================
## SUPPLEMENTARY PLATES (rasters where source objects are not in results/)
## ===========================================================================
lab_g <- function(letter)
  grid::textGrob(letter, gp = gpar(fontface = "bold", fontsize = 11))
rp <- function(file, letter, w_in, h_in) {
  need(file)
  if (!have_png)
    stop("FATAL: package 'png' required for raster panel ", file, " (R2).",
         call. = FALSE)
  img <- png::readPNG(file)
  g <- grid::rasterGrob(img, interpolate = TRUE)
  gridExtra::arrangeGrob(
    lab_g(letter),
    gridExtra::arrangeGrob(
      ggplot2::ggplot() + ggplot2::annotation_custom(g) +
        ggplot2::geom_blank() + ggplot2::theme_void(),
      ncol = 1, widths = grid::unit(w_in, "in"),
      heights = grid::unit(h_in, "in")),
    ncol = 2, widths = grid::unit(c(0.16, w_in), "in"))
}

## FigS1 -- exploratory decidualization (requires her re-run with the new y-label)
arr2(list(rp(file.path(FIGDIR, "figS1a_drugs_ecsc.png"), "a", 3.4, 2.9),
          rp(file.path(FIGDIR, "figS1b_drugs_nesc.png"), "b", 3.4, 2.9)),
     ncol = 2, widths_in = c(3.56, 3.56), heights_in = c(2.9),
     name = "FigS1", w_pad = 0.15, h_pad = 0.12)

## FigS2 -- bulk cohort PCA / QC (7 cohorts; rasters re-used, no black boxes)
pca_files <- file.path(FIGDIR, c("GSE6364_pca.png", "GSE7305_pca.png",
                                 "GSE7307_pca.png", "GSE11691_pca.png",
                                 "GSE25628_pca.png", "GSE51981_pca.png",
                                 "GSE120103_pca.png"))
arr2(list(rp(pca_files[1], "a", 1.72, 1.5), rp(pca_files[2], "b", 1.72, 1.5),
          rp(pca_files[3], "c", 1.72, 1.5), rp(pca_files[4], "d", 1.72, 1.5),
          rp(pca_files[5], "e", 1.72, 1.5), rp(pca_files[6], "f", 1.72, 1.5),
          rp(pca_files[7], "g", 1.72, 1.5)),
     ncol = 4, widths_in = c(1.88, 1.88, 1.88, 1.88), heights_in = c(1.5, 1.5),
     name = "FigS2", w_pad = 0.15, h_pad = 0.15)

## FigS3 -- scRNA QC summaries (3 datasets)
s3_files <- file.path(FIGDIR, c("S2_qc_pub_GSE179640_summary.png",
                                "S2_qc_pub_GSE183837_summary.png",
                                "S2_qc_pub_GSE214411_summary.png"))
arr2(list(rp(s3_files[1], "a", 2.3, 1.5), rp(s3_files[2], "b", 2.3, 1.5),
          rp(s3_files[3], "c", 2.3, 1.5)),
     ncol = 3, widths_in = c(2.46, 2.46, 2.46), heights_in = c(1.5),
     name = "FigS3", w_pad = 0.15, h_pad = 0.12)

## FigS4 -- integrated embedding (nomenclature caption added on the plate)
s4_file <- file.path(FIGDIR, "S2_02_fig_integrated_4panel.png")
need(s4_file)
img4 <- png::readPNG(s4_file)
g4 <- grid::rasterGrob(img4, interpolate = TRUE)
cap4 <- grid::textGrob(
  "Lineage nomenclature: Perivascular = mural (smooth muscle); Stromal = mesenchymal.",
  gp = gpar(fontsize = 8, col = INK), x = unit(0, "npc"), just = "left")
pS4 <- gridExtra::arrangeGrob(
  gridExtra::arrangeGrob(
    ggplot2::ggplot() + ggplot2::annotation_custom(g4) + ggplot2::geom_blank() +
      ggplot2::theme_void(), ncol = 1, widths = grid::unit(4.6, "in"),
      heights = grid::unit(4.4, "in")),
  cap4, ncol = 1, heights = grid::unit(c(4.4, 0.3), c("in", "in")))
save_plate(gridExtra::arrangeGrob(lab_g("a"),
          gridExtra::arrangeGrob(pS4, ncol = 1),
          ncol = 2, widths = grid::unit(c(0.16, 4.6), "in")),
          "FigS4", 4.75 + 0.15, 4.7 + 0.15)

## FigS5 -- post-QC distributions (ECDF + violins)
s5_files <- file.path(FIGDIR, c("S2_qc_pub_ecdf.png",
                                "S2_qc_pub_GSE179640_violins.png",
                                "S2_qc_pub_GSE183837_violins.png",
                                "S2_qc_pub_GSE214411_violins.png"))
arr2(list(rp(s5_files[1], "a", 3.4, 2.6), rp(s5_files[2], "b", 3.4, 2.6),
          rp(s5_files[3], "c", 3.4, 2.6), rp(s5_files[4], "d", 3.4, 2.6)),
     ncol = 2, widths_in = c(3.56, 3.56), heights_in = c(2.6, 2.6),
     name = "FigS5", w_pad = 0.15, h_pad = 0.15)

## FigS7 -- external-validation test matrix (drawn; locked endpoint codes)
s6t$testable <- !grepl("not", tolower(s6t$status))
s6t$ep <- factor(s6t$endpoint, levels = c("E1", "E2", "E3"),
                 labels = c("six pathway axes", "SERPINE1 logCPM",
                            "Tier-1 signature"))
tt <- aggregate(testable ~ lineage + ep, data = s6t,
                FUN = function(z) c(planned = length(z), tested = sum(z)))
tt <- data.frame(lineage = tt$lineage, ep = tt$ep, planned = tt$testable[, 1],
                 tested = tt$testable[, 2])
tt$frac <- tt$tested / tt$planned
tt$txt <- sprintf("%d/%d", tt$tested, tt$planned)
figs7 <- ggplot(tt, aes(ep, lineage)) +
  geom_tile(aes(fill = frac), color = "white", linewidth = 0.7) +
  geom_text(aes(label = txt), size = 2.2, color = "white", fontface = "bold") +
  scale_fill_gradient2(low = C_FAIL, mid = "#B07C3E", high = C_DOWN,
                       midpoint = 0.5, limits = c(0, 1),
                       name = "fraction testable") +
  labs(title = sprintf("External-validation test matrix (%d planned contrasts)",
                       nrow(s6t)),
       subtitle = sprintf("%d testable | %d not testable; cell text = tested / planned\nnot-testable = lineage absent in a patient or gene set unmeasurable; measured, never dropped",
                          n_ok, n_not),
       x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 6, color = INK),
        legend.text = element_text(size = 4.5), legend.key.size = unit(2, "mm"),
        panel.background = element_blank())
save_plate(tg(figs7, "a"), "FigS7", 3.8, 2.9)

say("== make_journal_plates DONE -- see folder '", OUT, "' ==")
