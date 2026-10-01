## ===========================================================================
## make_figS1_exploratory_drugs.R -- Supplementary Fig S1 (EXPLORATORY ONLY)
##   dienogest + cAMP vs the lesion Tier-1 program (S8b), per-gene scatter.
## v1.1 (2026-09-25): PI review -- legend restored (points = Tier-1 genes;
##   dashed = linear fit with 95% CI band); titles/subtitles shortened to stop
##   clipping in the combined plate. No numeric change (presentation only).
## ---------------------------------------------------------------------------
## Written under the same project discipline (R1/R2/R4/R5/R6/R8/R9):
##   * READ-ONLY on locked results; the only computation is a RE-DERIVATION of
##     the per-gene (M, logFC) pairs from the GEO cache that S8b_01 already
##     downloaded (no new download logic here -- the cache must exist).
##   * CERTIFICATION GATE: the re-derived r and bootstrap CI must match the
##     locked S8b summary (S8b_01_exploratory_drug.csv, md5-registry checked)
##     within 1e-4 (r) / 0.01 (CI) BEFORE any pixel is drawn; otherwise FATAL.
##   * Presentation only: no new estimands, no new thresholds (A5/R5).
##   * The figure is watermarked EXPLORATORY and carries the locked stop-rule.
## Locked inputs:
##   results/S4_meta/S4_meta_Tier1_core.csv          (registry md5 below)
##   results/S8_perturb/S8b_01_exploratory_drug.csv  (registry md5 below)
##   results/S8_perturb/cache_geo/GSE75423_series_matrix.txt.gz  (ECSC)
##   results/S8_perturb/cache_geo/GSE75425_series_matrix.txt.gz  (NESC)
##   results/S8_perturb/cache_geo/GPL13497.soft                  (annotation)
## Output: figures/figS1a_drugs_ecsc.* + figS1b_drugs_nesc.* + combined plate
## Usage (project root):
##   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts/make_figS1_exploratory_drugs.R
## ===========================================================================

set.seed(42)

say <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

FIGDIR <- "figures"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("FATAL: ggplot2 required (in renv lock).", call. = FALSE)
suppressPackageStartupMessages(library(ggplot2))
have_grid <- requireNamespace("gridExtra", quietly = TRUE)

C_DOWN  <- "#1F4E79"   # navy (consistent with make_figures palette)
C_UP    <- "#C00000"   # crimson
C_WARN  <- "#C55A11"   # amber (exploratory watermark)
INK     <- "#262626"

theme_s1 <- function(base = 8.5) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.line  = element_line(color = INK, linewidth = 0.4),
      axis.ticks = element_line(color = INK, linewidth = 0.4),
      axis.title = element_text(color = INK, face = "bold", size = base + 1),
      axis.text  = element_text(color = INK, size = base),
      plot.title    = element_text(color = C_WARN, face = "bold",
                                   size = base + 2.5, hjust = 0,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(color = "#595959", size = base - 0.5,
                                   hjust = 0, margin = margin(b = 4),
                                   lineheight = 0.95),
      plot.margin = margin(6, 10, 4, 6),
      legend.position = "top",
      legend.text  = element_text(size = 8, color = INK),
      legend.key.width = grid::unit(9, "mm"),
      legend.margin = margin(b = 0),
      plot.caption = element_text(color = C_WARN, size = base - 0.5,
                                  face = "italic", hjust = 0)
    )
}

save_both <- function(plot, name, w, h) {
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".png")), plot, width = w,
                  height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(FIGDIR, paste0(name, ".pdf")), plot, width = w,
                  height = h, device = "pdf", bg = "white")
  say("wrote figures/", name, ".png/.pdf")
}

say("== make_figS1 v1.1: EXPLORATORY supplementary figure, certified re-derivation ==")

## ---- locked inputs + registry (R2) -----------------------------------------
F_T1  <- "results/S4_meta/S4_meta_Tier1_core.csv"
F_SUM <- "results/S8_perturb/S8b_01_exploratory_drug.csv"
F_G23 <- "results/S8_perturb/cache_geo/GSE75423_series_matrix.txt.gz"
F_G25 <- "results/S8_perturb/cache_geo/GSE75425_series_matrix.txt.gz"
F_ANN <- "results/S8_perturb/cache_geo/GPL13497.soft"

REG <- c("results/S4_meta/S4_meta_Tier1_core.csv" = "3afbc9d85bec0b4409137cb73cc854bb",
         "results/S8_perturb/S8b_01_exploratory_drug.csv" = "5307042af5acfcabb5985a6ca3cb9710")
for (p in names(REG)) {
  if (!file.exists(p))
    stop("FATAL: locked input missing: ", p, " (R2)", call. = FALSE)
  if (unname(tools::md5sum(p)) != REG[[p]])
    stop("FATAL: md5 mismatch for ", p, " (locked registry)", call. = FALSE)
}
for (p in c(F_G23, F_G25, F_ANN))
  if (!file.exists(p))
    stop("FATAL: GEO cache missing: ", p,
         " -- run scripts/S8b_01_exploratory_drug.R first (R2).", call. = FALSE)
say("figS1: locked inputs verified (md5 match); GEO cache present")

t1 <- read.csv(F_T1, stringsAsFactors = FALSE)
sm <- read.csv(F_SUM, stringsAsFactors = FALSE)
need_cols <- function(d, cols, p) {
  miss <- setdiff(cols, names(d))
  if (length(miss))
    stop("FATAL: ", p, " lacks columns: ", paste(miss, collapse = ", "),
         " (R2)", call. = FALSE)
}
need_cols(t1, c("gene", "M"), F_T1)
need_cols(sm, c("dataset", "r", "ci_lo", "ci_hi"), F_SUM)

## ---- annotation parsing (identical logic to S8b_01) -------------------------
lns <- readLines(F_ANN, warn = FALSE)
h   <- grep("^ID\t", lns)[1]
if (is.na(h)) stop("FATAL: platform table header not found in GPL13497.soft",
                   call. = FALSE)
tbl <- lns[(h + 1):length(lns)]
tbl <- tbl[grepl("^A_23_P", tbl)]
if (length(tbl) < 1000L)
  stop("FATAL: platform table looks truncated (", length(tbl), " rows)",
       call. = FALSE)
col_ann <- strsplit(lns[h], "\t", fixed = TRUE)[[1]]
i_id  <- match("ID", col_ann)
i_sym <- match("GENE_SYMBOL", col_ann)
ann <- vapply(strsplit(tbl, "\t", fixed = TRUE),
              function(p) if (length(p) >= i_sym) p[i_sym] else "", character(1))
names(ann) <- vapply(strsplit(tbl, "\t", fixed = TRUE),
                     function(p) if (length(p) >= i_id) p[i_id] else "",
                     character(1))
ann <- ann[ann != ""]

## ---- re-derivation (identical logic to S8b_01) ------------------------------
read_gse <- function(f) {
  con <- gzfile(f, "rt")
  x <- read.delim(con, comment.char = "!", check.names = FALSE,
                  stringsAsFactors = FALSE)
  close(con)
  x
}
collapse_lfc <- function(x) {
  x$sym <- ann[as.character(x[[1]])]
  x <- x[!is.na(x$sym), ]
  nums <- vapply(x[, -c(1, ncol(x))], as.numeric, numeric(nrow(x)))
  meds <- tapply(seq_len(nrow(x)), x$sym,
                 function(ii) apply(nums[ii, , drop = FALSE], 2, median))
  g <- do.call(rbind, meds)
  rownames(g) <- names(meds)
  lg <- log2(g + 1)
  stopifnot(ncol(lg) == 8L)  ## 4 vehicle + 4 treated, same patient order
  rowMeans(lg[, 5:8] - lg[, 1:4])
}
boot_ci <- function(Mv, Lv, B = 2000, seed = 42) {
  set.seed(seed)
  out <- numeric(B)
  n <- length(Mv)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(Lv[idx])) > 1L) out[b] <- cor(Mv[idx], Lv[idx])
  }
  out <- out[!is.na(out)]
  quantile(out, c(0.025, 0.975))
}
prep <- function(f, label) {
  lfc <- collapse_lfc(read_gse(f))
  ov <- t1$gene[t1$gene %in% names(lfc)]
  Mv <- t1$M[match(ov, t1$gene)]
  Lv <- unname(lfc[ov])
  data.frame(dataset = label, gene = ov, M = Mv, logFC = Lv,
             stringsAsFactors = FALSE)
}
d23 <- prep(F_G23, "GSE75423_ECSC_dienogest")
d25 <- prep(F_G25, "GSE75425_NESC_dienogest")
say("figS1: per-gene pairs re-derived | ECSC n = ", nrow(d23),
    " | NESC n = ", nrow(d25), " Tier-1 genes")

## ---- certification against the locked S8b summary (R9) ----------------------
certify <- function(d) {
  row <- sm[sm$dataset == d$dataset[1], ]
  if (nrow(row) != 1L)
    stop("FATAL: dataset ", d$dataset[1],
         " absent from locked S8b summary (R2)", call. = FALSE)
  r  <- cor(d$M, d$logFC)
  ci <- boot_ci(d$M, d$logFC)
  say("  ", d$dataset[1], " | recomputed r = ", sprintf("%.6f", r),
      " | locked r = ", sprintf("%.6f", row$r),
      " | recomputed CI [", sprintf("%.4f", unname(ci[1])), ", ",
      sprintf("%.4f", unname(ci[2])), "]")
  if (abs(r - row$r) > 1e-4)
    stop("FATAL: re-derived r does not match locked S8b run for ",
         d$dataset[1], " (|delta| = ", abs(r - row$r),
         "). Refusing to draw from unverified numbers (R9).", call. = FALSE)
  if (abs(unname(ci[1]) - row$ci_lo) > 0.01 ||
      abs(unname(ci[2]) - row$ci_hi) > 0.01)
    stop("FATAL: re-derived bootstrap CI does not match locked S8b run for ",
         d$dataset[1], " (R9).", call. = FALSE)
  list(r = r, ci = ci)
}
c23 <- certify(d23)
c25 <- certify(d25)
say("figS1: CERTIFIED -- re-derivation matches locked S8b summary (R9)")

## ---- plot (presentation only) ----------------------------------------------
KEY <- c("TOP2A", "PCNA", "IL6", "SOCS3", "FOSB", "C1QB")  ## locked key genes
mk_panel <- function(d, c, ttl_arm) {
  kd <- d[d$gene %in% KEY, ]
  ggplot(d, aes(x = M, y = logFC)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.5) +
    geom_point(aes(color = "Tier-1 genes"), size = 1.8, alpha = 0.6) +
    geom_smooth(aes(linetype = "linear fit (95% CI)"), method = "lm",
                se = TRUE, color = C_WARN, linewidth = 0.6, fill = "#F2DCB3") +
    scale_color_manual(values = c("Tier-1 genes" = C_DOWN)) +
    scale_linetype_manual(values = c("linear fit (95% CI)" = "dashed")) +
    geom_text(data = kd, aes(label = gene), size = 2.4, color = INK,
              fontface = "bold", vjust = -0.9) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             size = 2.8, color = INK, fontface = "bold",
             label = paste0("r = ", round(c$r, 3),
                            "  |  bootstrap 95% CI [",
                            round(unname(c$ci[1]), 3), ", ",
                            round(unname(c$ci[2]), 3), "]")) +
    labs(x = "lesion Tier-1 effect size M (locked S4 meta)",
         y = "paired log2 FC, dienogest + cAMP vs vehicle",
         title = paste0("EXPLORATORY — ", ttl_arm),
         subtitle = paste0(nrow(d),
           " Tier-1 genes | single dataset, paired 4v4 | NO claim gate\npositive r = drug moves genes WITH the lesion program (no reversal)"),
         caption = "hypothesis-generating only — locked stop rule: in-silico chapter closed") +
    theme_s1()
}
figS1a <- mk_panel(d23, c23, "ECSC (GSE75423)")
figS1b <- mk_panel(d25, c25, "NESC (GSE75425)")
save_both(figS1a, "figS1a_drugs_ecsc", 4.6, 4.2)
save_both(figS1b, "figS1b_drugs_nesc", 4.6, 4.2)
if (have_grid) {
  p <- gridExtra::grid.arrange(figS1a, figS1b, ncol = 2,
                               widths = grid::unit(c(4.6, 4.6), "in"),
                               padding = grid::unit(3, "mm"))
  save_both(p, "figS1_drugs_exploratory", 9.5, 4.2)
} else {
  say("gridExtra absent -- individual panels only (figS1a / figS1b)")
}
say("figS1 DONE -- exploratory supplementary figure, certified (R9)")
