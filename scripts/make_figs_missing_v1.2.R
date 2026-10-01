## ===========================================================================
## make_figs_missing.R -- panels missing from figures/ (v1.2, 2026-09-29)
##   Fig3d  robust-core summary (66 genes)
##   Fig6a  external-validation design strip
##   Fig7a + Fig7d  LINCS coverage audit and hypothesis verdicts
##   (Fig7b/c already exist as fig7_s8_reverse_confirmation.png -- kept as-is)
##   FigS7  external-validation test matrix (195 / 109 / 86)
## ALL numbers are read from locked results files (R1/R2), fail-loud.
## NOTE: the gene-level GSE11691 scatter cannot be drawn from shipped results
##   (only the summary row exists); add the annotation
##   "rho = 0.324 | 391 Tier-1 genes | 80.3% sign agreement" to Fig5b instead.
## Usage (project root):
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/make_figs_missing.R
## ===========================================================================
set.seed(42)
source("scripts/S4_common.R")
FIGDIR <- "figures"; dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("FATAL: ggplot2 required (in renv lock).", call. = FALSE)
suppressPackageStartupMessages(library(ggplot2))

say  <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")
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
save_both <- function(plot, name, w, h) {
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".png")), plot, width = w,
                  height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".pdf")), plot, width = w,
                  height = h, device = "pdf", bg = "white")
  say("wrote figures/", name, ".png/.pdf")
}
INK <- "#262626"
theme_pub <- function(base = 8.5) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(axis.title = element_text(face = "bold", size = base),
          axis.text = element_text(color = INK, size = base - 1),
          plot.title = element_text(face = "bold", size = base + 2, hjust = 0),
          plot.subtitle = element_text(color = "#595959", size = base - 1, hjust = 0),
          plot.margin = margin(6, 10, 4, 6),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = base - 1),
          strip.background = element_rect(fill = "#EFEFEF", color = "#BFBFBF", linewidth = 0.3),
          strip.text = element_text(face = "bold", size = base))
}

say("== make_figs_missing v1.2: all numbers from locked results (R2) ==")

## ---- Fig3d: robust core (66 genes) --------------------------------------
f <- "results/S8_perturb/S8_00_core_gene_coverage.csv"
core <- read_locked(f)
check_cols(core, c("gene", "M", "in_L1000", "is_landmark"), f)
core <- core[order(-core$M), ]
core$rank <- seq_len(nrow(core))
core$dir  <- ifelse(core$M > 0, "up-regulated", "down-regulated")
n_up <- sum(core$M > 0); n_dn <- sum(core$M < 0)

mod_up <- intersect(c("FOSB", "FOS", "EGR1", "JUNB", "IL6", "SOCS3", "C1QB"),
                    core$gene[core$M > 0])
lab <- core[core$gene %in% mod_up, ]
lab$ylab <- lab$M + ifelse(lab$M > 0, 0.12, -0.12)
## stagger overlapping labels (FOS/FOSB, EGR1/JUNB/SOCS3)
lab$ylab[lab$gene == "FOS"]    <- lab$ylab[lab$gene == "FOS"] + 0.22
lab$ylab[lab$gene == "JUNB"]   <- lab$ylab[lab$gene == "JUNB"] + 0.16
lab$ylab[lab$gene == "SOCS3"]  <- lab$ylab[lab$gene == "SOCS3"] - 0.20
lab$ylab[lab$gene == "C1QB"]   <- lab$ylab[lab$gene == "C1QB"] - 0.16

fig3d <- ggplot(core, aes(rank, M)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = dir), size = 1.5) +
  geom_point(data = core[core$is_landmark == TRUE, ], size = 2.8, shape = 21,
             fill = NA, color = "black", stroke = 0.5) +
  geom_text(data = lab, aes(y = ylab, label = gene), size = 2.4,
            color = INK, fontface = "italic") +
  scale_color_manual(values = c("up-regulated" = "#C00000",
                                "down-regulated" = "#1F4E79")) +
  labs(title = sprintf("Robust core: %d genes (up %d, down %d) retained under all seven cohort-removal scenarios",
                       nrow(core), n_up, n_dn),
       subtitle = paste0("Six leave-one-cohort-out reanalyses plus the sans-GSE7307 boundary test;\n",
                         "rings = the 8 L1000 landmark genes; up module = AP-1 / immediate-early response\n",
                         "(FOS, FOSB, EGR1, JUNB, IER2) plus IL6, SOCS3 and C1QB."),
       x = "robust-core genes ranked by pooled M", y = "pooled M (log2 FC)") +
  theme_pub()
save_both(fig3d, "Fig3d_robust_core_summary", 7.2, 3.2)

## ---- Fig6a: validation design strip -------------------------------------
f1 <- "results/S2_qc_GSE213216_audit.csv"
f2 <- "results/S6_meta/S6_01_scores_long.csv"
f3 <- "results/S6_meta/S6_01_tests.csv"
audit <- read_locked(f1); check_cols(audit, c("object", "cells"), f1)
s6l <- read_locked(f2);   check_cols(s6l, c("patient", "lineage"), f2)
s6t <- read_locked(f3);   check_cols(s6t, c("status"), f3)
n_cells <- audit$cells[audit$object == "auxiliary.seurat.shared.rds"]
n_pat  <- length(unique(s6l$patient))
n_pb   <- nrow(s6l)
n_plan <- nrow(s6t)
n_not  <- sum(grepl("not", tolower(s6t$status)))
n_ok   <- n_plan - n_not
dsg <- data.frame(
  x = 1:4,
  lbl = c(sprintf("%s cells", format(n_cells, big.mark = ",")),
          sprintf("%d patients", n_pat),
          sprintf("%d pseudobulks", n_pb),
          sprintf("%d planned tests\n%d testable | %d not testable", n_plan, n_ok, n_not)),
  fill = c("#5B2D8E", "#5B2D8E", "#5B2D8E", "#8B1A1A"))
fig6a <- ggplot(dsg, aes(x, 1)) +
  geom_tile(aes(fill = fill), height = 0.7, width = 0.92, color = "white", linewidth = 1) +
  geom_text(aes(label = lbl), color = "white", size = 2.8, lineheight = 0.95) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  ylim(0.5, 1.5) +
  labs(title = "External-validation design (GSE213216)",
       subtitle = "independent single-cell atlas; patient-level pseudobulks; verdicts machine-printed") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())
save_both(fig6a, "Fig6a_validation_design", 7.2, 1.6)

## ---- Fig7: LINCS connectivity-based reversal analysis --------------------
fc <- "results/S8_perturb/S8_00_coverage.csv"
fh <- "results/S8_perturb/S8_01_hypotheses.csv"
cov <- read_locked(fc); check_cols(cov, c("measure", "n"), fc)
hyp <- read_locked(fh); check_cols(hyp, c("hypothesis", "test", "n", "p_value", "BH_FDR", "verdict"), fh)
## NOTE (v1.1): Fig7b/c already exist as figures/fig7_s8_reverse_confirmation.png
## (dot plots, direction-coloured) -- superior to bar versions; NOT redrawn here.

g <- function(x) cov$n[cov$measure == x][1]
fun <- data.frame(
  x = factor(c("Tier-1 genes", "L1000 landmark", "Robust core (66)",
               "Core in L1000", "Core landmark"),
             levels = c("Tier-1 genes", "L1000 landmark", "Robust core (66)",
                        "Core in L1000", "Core landmark")),
  n = c(485, g("tier1_total_landmark"), 66, g("robust_core_in_L1000"),
        g("robust_core_landmark")),
  grp = c("Tier-1", "Tier-1", "robust core", "robust core", "robust core"))
fig7a <- ggplot(fun, aes(x, n)) +
  geom_col(aes(fill = grp), width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3, fontface = "bold", color = INK) +
  scale_fill_manual(values = c("Tier-1" = "#8FAADC", "robust core" = "#5B2D8E")) +
  labs(title = "LINCS coverage audit",
       subtitle = sprintf("Tier-1: %d landmark genes (%d up / %d down of 485); robust core: %d of 66 in L1000,\n%d landmark (EGR1, FOS, GADD45B, CCND1, LYPLA1, PLA2G4A, NNT, UGDH)",
                          g("tier1_total_landmark"), g("tier1_up_landmark"),
                          g("tier1_down_landmark"), g("robust_core_in_L1000"),
                          g("robust_core_landmark")),
       x = NULL, y = "genes") +
  theme_pub() + theme(legend.position = "none",
        axis.text.x = element_text(size = 6.5, lineheight = 0.95))
save_both(fig7a, "Fig7a_lincs_coverage", 4.6, 3.2)

hyp$txt_p   <- sprintf("p = %.3g", hyp$p_value)
hyp$txt_fdr <- sprintf("FDR = %.3g", hyp$BH_FDR)
fig7d <- ggplot(hyp, aes(x = 1, y = seq(nrow(hyp), 1))) +
  geom_text(aes(x = 0.0, label = hypothesis), hjust = 0, size = 2.9, fontface = "bold", color = INK) +
  geom_text(aes(x = 0.35, label = test), hjust = 0, size = 2.3, color = "#595959") +
  geom_text(aes(x = 0.85, label = txt_p), hjust = 0, size = 2.6, color = INK) +
  geom_text(aes(x = 1.02, label = txt_fdr), hjust = 0, size = 2.6, color = INK) +
  geom_text(aes(x = 1.22, label = verdict), hjust = 0, size = 2.6, fontface = "bold",
            color = "#8B1A1A") +
  xlim(0, 1.45) + ylim(0.4, nrow(hyp) + 0.8) +
  labs(title = "Pre-registered reversal hypotheses: none confirmed",
       subtitle = "best BH-FDR = 0.129 (P2b: estradiol - proliferation axis); reported as a negative result") +
  theme_pub() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.background = element_blank())
save_both(fig7d, "Fig7d_hypotheses", 7.2, 2.2)

## ---- FigS7: external-validation test matrix ------------------------------
## locked endpoint codes (verified vs S6_01_tests.csv):
##   E1 = six pathway axes (cytotoxicity, antigen_presentation,
##        senescence_dormancy, proliferation, stromal_ecm, SenMayo)
##   E2 = SERPINE1_logCPM (Endothelial / Smooth muscle x C1-C3)
##   E3 = tier1_signature (all lineages x C1-C3)
s6t$testable <- !grepl("not", tolower(s6t$status))
s6t$ep <- factor(s6t$endpoint, levels = c("E1", "E2", "E3"),
                 labels = c("six pathway axes", "SERPINE1 logCPM",
                            "Tier-1 signature"))
tt <- aggregate(testable ~ lineage + ep, data = s6t,
                FUN = function(z) c(planned = length(z), tested = sum(z)))
tt <- data.frame(lineage = tt$lineage, ep = tt$ep, planned = tt$testable[, 1],
                 tested = tt$testable[, 2])
tt$frac <- tt$tested / tt$planned
tt$txt  <- sprintf("%d/%d", tt$tested, tt$planned)
figs7 <- ggplot(tt, aes(ep, lineage)) +
  geom_tile(aes(fill = frac), color = "white", linewidth = 0.8) +
  geom_text(aes(label = txt), size = 2.7, color = "white", fontface = "bold") +
  scale_fill_gradient2(low = "#8B1A1A", mid = "#B07C3E", high = "#1F4E79",
                       midpoint = 0.5, limits = c(0, 1), name = "fraction testable") +
  labs(title = sprintf("External-validation test matrix (%d planned contrasts; %d testable, %d not testable)",
                       nrow(s6t), n_ok, n_not),
       subtitle = "cell text = tested / planned; not-testable = lineage absent in a patient or gene set unmeasurable; measured, never dropped",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(size = 7.5, color = INK),
        panel.background = element_blank())
save_both(figs7, "FigS7_test_matrix", 5.2, 3.4)

say("== make_figs_missing DONE ==")
