# ============================================================================
# S2_03_journal_figures_v2.R  v2 -- publication-quality PCA panels + plate
#
# Changelog v1 -> v2 (lead review 2026-09-18, defects 1-5 + plate request):
#   1. Labels: FLAGGED-ONLY in every cohort (the n<=30 full-label rule is
#      abolished). Template unified on the GSE51981/GSE120103 look:
#      colored points + ringed flagged outliers only.
#   2. No clipped labels: adaptive hjust/vjust by position + 16% axis
#      expansion + coord_cartesian(clip = "off") + plot margin.
#   3. GSE7305 label-crowding: resolved by (1) -- no unflagged labels remain.
#   4. GSE6364: phase kept as shape (phase is a design covariate needed at
#      S4), legend title shortened to "Phase".
#   5. Legend title "Condition" removed everywhere (caption carries it).
#   6. NEW: PLATE_all_cohorts (7 panels A-G, facet_wrap, free scales,
#      shared color legend at bottom, flagged = black ring, no in-panel
#      text labels, phase shapes omitted in the plate only).
#
# Reads (never modifies): data/processed/bulk_<GSE>_rma.rds (6 cohorts),
#   data/processed/bulk_GSE120103_neqc.rds,
#   data/metadata/S1_bulk_conditions.tsv (GEO-verified by S2_03_meta_build.R),
#   results/S2_bulk_outliers_<GSE>.csv (the 6 AQM cohorts).
# Writes ONLY to figures_journal/ (overwrites v1 PNG/PDF outputs there; the
#   old figures/ folder and all data stay untouched).
#
# Fail-loud: every plotted GSM must carry a condition; every plotted condition
#   must have a frozen Okabe-Ito color; GSE120103's labeled-but-absent set must
#   be exactly the 5 documented GE1-v5 exclusions; outlier CSV mandatory for
#   the 6 AQM cohorts.
# R10 RAM estimate: largest matrix 54675 x 148 + plate render -> < 2 GB.
#   Runtime ~2 min.
#
# Run (from C:\endometriosis_immunoediting):
#   "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe" scripts\S2_03_journal_figures_v2.R > logs\S2_03_journal_figures_v2.txt 2>&1
# ============================================================================
set.seed(42)
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n"); flush.console() }
suppressPackageStartupMessages({ library(Biobase); library(ggplot2) })

outdir <- "figures_journal"
dir.create(outdir, showWarnings = FALSE)
if (!dir.exists(outdir)) stop("FATAL: could not create ", outdir)

## ---- frozen Okabe-Ito mapping (fail-loud on unmapped condition) --------------
OI <- c("Normal endometrium"              = "#0072B2",
        "Normal eutopic"                  = "#0072B2",
        "Non-EMS no pelvic pathology"     = "#0072B2",
        "Non-EMS pelvic pathology"        = "#56B4E9",
        "Infertile control"               = "#56B4E9",
        "Fertile control"                 = "#009E73",
        "Endometriosis"                   = "#D55E00",
        "Endometriosis eutopic"           = "#D55E00",
        "Endometriosis moderate/severe"   = "#D55E00",
        "Endometriosis minimal/mild"      = "#E69F00",
        "Endometriosis (severity NA)"     = "#F0E442",
        "Eutopic (patient)"               = "#E69F00",
        "Ectopic lesion"                  = "#CC79A7",
        "Stage IV ovarian EMS (2A)"       = "#D55E00",
        "Stage IV ovarian EMS (2B)"       = "#CC79A7")

COHORTS <- c("GSE7305", "GSE6364", "GSE25628", "GSE7307", "GSE51981", "GSE11691", "GSE120103")
RDS <- c(GSE7305 = "data/processed/bulk_GSE7305_rma.rds",
         GSE6364 = "data/processed/bulk_GSE6364_rma.rds",
         GSE25628 = "data/processed/bulk_GSE25628_rma.rds",
         GSE7307 = "data/processed/bulk_GSE7307_rma.rds",
         GSE51981 = "data/processed/bulk_GSE51981_rma.rds",
         GSE11691 = "data/processed/bulk_GSE11691_rma.rds",
         GSE120103 = "data/processed/bulk_GSE120103_neqc.rds")

lab <- read.delim("data/metadata/S1_bulk_conditions.tsv", comment.char = "#",
                  stringsAsFactors = FALSE)
say("labels loaded:", nrow(lab), "rows across", length(unique(lab$gse)), "cohorts")

plate <- NULL  # accumulator for the composite plate

for (g in COHORTS) {
  say("== ", g, " ==")
  if (!file.exists(RDS[[g]])) stop("FATAL: ", RDS[[g]], " missing")
  obj <- readRDS(RDS[[g]])
  if (is(obj, "ExpressionSet")) mat <- Biobase::exprs(obj)
  else if (is.list(obj) && !is.null(obj$E)) mat <- obj$E
  else stop("FATAL: unrecognized object class in ", RDS[[g]], ": ", paste(class(obj), collapse = "/"))
  n <- ncol(mat)
  say("  matrix:", nrow(mat), "x", n)

  pc  <- prcomp(t(mat), scale. = TRUE)
  pct <- round(100 * summary(pc)$importance[2, 1:2], 1)
  df  <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], sample = colnames(mat),
                    gsm = regmatches(colnames(mat), regexpr("GSM[0-9]+", colnames(mat))),
                    lbl = regmatches(colnames(mat), regexpr("GSM[0-9]+(_[0-9]+)?", colnames(mat))),
                    stringsAsFactors = FALSE)

  lg <- lab[lab$gse == g, ]
  idx <- match(df$gsm, lg$gsm)
  if (any(is.na(idx)))
    stop("FATAL: ", g, " samples without condition labels: ",
         paste(df$gsm[is.na(idx)], collapse = ", "), " -- escalate (R1)")
  df$condition <- lg$condition[idx]
  df$phase     <- lg$detail[idx]
  absent <- setdiff(lg$gsm, df$gsm)
  if (g == "GSE120103") {
    if (!identical(sort(absent), sort(paste0("GSM", 3393522:3393526))))
      stop("FATAL: GSE120103 labeled-but-absent set != documented 5 GE1-v5 exclusions: ",
           paste(absent, collapse = ", "))
    say("  labeled-but-absent = exactly the 5 documented GE1-v5 exclusions [OK]")
  } else if (length(absent)) {
    say("  note: labeled samples not plotted in", g, ":", paste(absent, collapse = ", "))
  }
  unmapped <- setdiff(unique(df$condition), names(OI))
  if (length(unmapped)) stop("FATAL: no frozen Okabe-Ito color for: ", paste(unmapped, collapse = ", "))

  ## flagged outliers (AQM cohorts only; CSV mandatory for them)
  ocsv <- file.path("results", paste0("S2_bulk_outliers_", g, ".csv"))
  flagged <- character(0)
  if (g != "GSE120103" && !file.exists(ocsv))
    stop("FATAL: outlier CSV missing for ", g, " -- run state inconsistent (R9)")
  if (file.exists(ocsv)) {
    ot <- read.csv(ocsv, stringsAsFactors = FALSE)
    flagged <- ot$sample[ot$flagged]
  }
  df$flagged <- df$sample %in% flagged
  say("  conditions:", paste(names(table(df$condition)), table(df$condition), sep = "=", collapse = " / "),
      "| flagged:", sum(df$flagged),
      if (sum(df$flagged)) paste0("(", paste(df$gsm[df$flagged], collapse = ", "), ")") else "")

  ## accumulate plate rows (phase shapes intentionally NOT carried to plate)
  plate <- rbind(plate, data.frame(cohort = g, PC1 = df$PC1, PC2 = df$PC2,
                                   condition = df$condition, flagged = df$flagged,
                                   stringsAsFactors = FALSE))

  ttl <- paste0(g, " - PCA (n = ", n, ")")
  if (g == "GSE120103") ttl <- paste0(g, " - PCA (n = 31 of 36; 5 design-mismatch exclusions documented)")

  p <- ggplot(df, aes(PC1, PC2))
  if (g == "GSE6364") {
    df$phase <- factor(df$phase, levels = c("Proliferative", "Early Secretory", "Mid Secretory"))
    if (any(is.na(df$phase))) stop("FATAL: GSE6364 phase values outside the frozen 3-level set")
    p <- p + geom_point(aes(color = condition, shape = phase), size = 2.6) +
      scale_shape_manual(values = c("Early Secretory" = 16, "Mid Secretory" = 17,
                                    "Proliferative" = 15),
                         name = "Phase")
  } else {
    p <- p + geom_point(aes(color = condition), size = 2.2, alpha = 0.9)
  }
  p <- p + scale_color_manual(values = OI, name = NULL)

  ## flagged: black ring + bold label with adaptive placement (no clipping)
  if (any(df$flagged)) {
    fl <- df[df$flagged, , drop = FALSE]
    xr <- range(df$PC1); dxr <- diff(xr)
    yr <- range(df$PC2); dyr <- diff(yr)
    fx <- (fl$PC1 - xr[1]) / dxr
    fy <- (fl$PC2 - yr[1]) / dyr
    fl$hj <- ifelse(fx > 0.75, 1, ifelse(fx < 0.25, 0, 0.5))
    fl$nx <- ifelse(fx > 0.75, -0.015 * dxr, ifelse(fx < 0.25, 0.015 * dxr, 0))
    fl$vj <- ifelse(fy > 0.85, 1.4, -0.7)
    fl$ny <- ifelse(fy > 0.85, -0.02 * dyr, 0.02 * dyr)
    p <- p + geom_point(data = fl, shape = 1, size = 5, stroke = 1.4, color = "black") +
      geom_text(data = fl,
                aes(x = PC1 + nx, y = PC2 + ny, label = lbl, hjust = hj, vjust = vj),
                size = 3.0, fontface = "bold", color = "black", show.legend = FALSE)
  }

  p <- p +
    scale_x_continuous(expand = expansion(mult = 0.16)) +
    scale_y_continuous(expand = expansion(mult = 0.16)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(plot.margin = margin(12, 18, 10, 10)) +
    labs(title = ttl, x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))

  fn <- file.path(outdir, paste0(g, "_pca"))
  ggsave(paste0(fn, ".png"), p, width = 7.5, height = 6, dpi = 300)
  ggsave(paste0(fn, ".pdf"), p, width = 7.5, height = 6)
  say("  wrote ", fn, ".png (300 dpi) + .pdf")
}

## ---- composite plate (panels A-G, free scales, shared legend) ----------------
say("== PLATE ==")
plate$strip <- factor(
  sprintf("%s. %s (n = %d)",
          LETTERS[match(plate$cohort, COHORTS)],
          plate$cohort,
          ave(seq_len(nrow(plate)), plate$cohort, FUN = length)),
  levels = sprintf("%s. %s (n = %d)", LETTERS[seq_along(COHORTS)], COHORTS,
                   sapply(COHORTS, function(cc) sum(plate$cohort == cc))))

pp <- ggplot(plate, aes(PC1, PC2)) +
  geom_point(aes(color = condition), size = 1.6, alpha = 0.9) +
  geom_point(data = plate[plate$flagged, ], shape = 1, size = 3.4, stroke = 1.0,
             color = "black") +
  scale_color_manual(values = OI, name = NULL,
                     guide = guide_legend(ncol = 3, byrow = TRUE,
                                          override.aes = list(size = 2.6, alpha = 1))) +
  facet_wrap(~ strip, scales = "free", ncol = 2) +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 9),
        strip.background = element_rect(fill = "grey92", color = NA),
        legend.position = "bottom",
        legend.text = element_text(size = 7.5),
        legend.key.size = grid::unit(0.35, "cm"),
        panel.grid.minor = element_blank()) +
  labs(x = "PC1", y = "PC2")

ggsave(file.path(outdir, "PLATE_all_cohorts.png"), pp, width = 9.5, height = 13.5, dpi = 300)
ggsave(file.path(outdir, "PLATE_all_cohorts.pdf"), pp, width = 9.5, height = 13.5)
say("  wrote PLATE_all_cohorts.png (300 dpi) + .pdf | panels A-G,",
    sum(plate$flagged), "flagged ringed; PC variance shown in individual figures only")

print(gc())
say("JOURNAL FIGURES V2 DONE")
