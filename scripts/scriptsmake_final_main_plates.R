# scripts/make_final_main_plates.R

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

required <- c(
  "ggplot2", "dplyr", "readr", "tidyr", "forcats",
  "patchwork", "scales", "stringr", "grid", "gridExtra"
)

missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Install missing packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(forcats)
  library(patchwork)
  library(scales)
  library(stringr)
  library(grid)
  library(gridExtra)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
res_dir <- file.path(root, "results")
fig_dir <- file.path(root, "figures", "final_plates")
src_dir <- file.path(fig_dir, "figure_source_data")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)

find_file <- function(patterns, root = res_dir) {
  hits <- unlist(lapply(patterns, function(p) {
    list.files(root, pattern = p, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }))
  hits <- unique(hits[file.exists(hits)])
  if (length(hits) == 0) {
    stop("Required file not found. Patterns: ", paste(patterns, collapse = " | "))
  }
  hits[1]
}

save_plate <- function(plot, stem, width, height, dpi = 600) {
  tiff_path <- file.path(fig_dir, paste0(stem, ".tiff"))
  pdf_path  <- file.path(fig_dir, paste0(stem, ".pdf"))

  ggsave(
    tiff_path, plot, width = width, height = height, units = "in",
    dpi = dpi, compression = "lzw", bg = "white"
  )

  ggsave(
    pdf_path, plot, width = width, height = height, units = "in",
    device = grDevices::pdf, bg = "white", useDingbats = FALSE
  )
}

theme_plate <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = "grey35", hjust = 0),
      plot.caption = element_text(size = base_size - 2, colour = "grey35", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(7, 10, 7, 7)
    )
}

theme_blank_plate <- function() {
  theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

label_panel <- function(tag) {
  ggplot() +
    annotate("text", x = 0, y = 1, label = tag, fontface = "bold", size = 7, hjust = 0, vjust = 1) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

blue <- "#235784"
blue_light <- "#88A8D8"
red <- "#C00000"
red_light <- "#E77C63"
green <- "#2E7D32"
purple <- "#5B2C92"
grey <- "#9B9B9B"
dark <- "#222222"
gold <- "#C96A0A"

meta_file <- find_file(c("^S4_meta_all_genes.*\\.csv$", "S4_meta_all_genes"))
tier_file <- find_file(c("^S4_meta_Tier1_core.*\\.csv$", "S4_meta_Tier1_core"))
sens_file <- find_file(c("^S4_sensitivity_summary.*\\.csv$", "S4_sensitivity_summary"))
path_file <- find_file(c("^S5_meta_pathways.*\\.csv$", "S5_meta_pathways"))
pair_path_file <- find_file(c("^S5_concordance_GSE11691.*\\.csv$", "S5_concordance_GSE11691"))
s6_file <- find_file(c("^S6_01_tests.*\\.csv$", "S6_01_tests"))
s6_scores_file <- find_file(c("^S6_01_scores_long.*\\.csv$", "S6_01_scores_long"))
contrast_file <- find_file(c("^S3_contrasts_frozen.*\\.tsv$", "S3_contrasts_frozen"))
lincs_file <- find_file(c("^S8_01_hypotheses.*\\.csv$", "S8_01_hypotheses"))

meta <- read_csv(meta_file, show_col_types = FALSE)
tier <- read_csv(tier_file, show_col_types = FALSE)
sens <- read_csv(sens_file, show_col_types = FALSE)
paths <- read_csv(path_file, show_col_types = FALSE)
pair_paths <- read_csv(pair_path_file, show_col_types = FALSE)
s6 <- read_csv(s6_file, show_col_types = FALSE)
s6_scores <- read_csv(s6_scores_file, show_col_types = FALSE)
contrasts <- read_tsv(contrast_file, comment = "#", show_col_types = FALSE)
lincs <- read_csv(lincs_file, show_col_types = FALSE)

write_csv(meta, file.path(src_dir, "Fig2_meta_all_genes.csv"))
write_csv(tier, file.path(src_dir, "Fig2_Fig3_Tier1_core.csv"))
write_csv(sens, file.path(src_dir, "Fig3_sensitivity_summary.csv"))
write_csv(paths, file.path(src_dir, "Fig4_pathway_meta.csv"))
write_csv(pair_paths, file.path(src_dir, "Fig4_paired_concordance.csv"))
write_csv(s6, file.path(src_dir, "Fig6_validation_all_tests.csv"))
write_csv(lincs, file.path(src_dir, "FigS1_LINCS_hypotheses.csv"))

names(meta) <- tolower(names(meta))
names(tier) <- tolower(names(tier))
names(sens) <- tolower(names(sens))
names(paths) <- tolower(names(paths))
names(pair_paths) <- tolower(names(pair_paths))
names(s6) <- tolower(names(s6))
names(contrasts) <- tolower(names(contrasts))

pick_col <- function(df, candidates, required = TRUE) {
  nms <- names(df)
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Missing expected column. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

gene_col <- pick_col(meta, c("gene", "symbol", "genesymbol", "hgnc_symbol"))
m_col <- pick_col(meta, c("m", "pooled_m", "effect", "estimate", "meta_effect"))
fdr_col <- pick_col(meta, c("fdr", "padj", "bh_fdr", "qvalue"))
i2_col <- pick_col(meta, c("i2", "i_squared"))
tier_col_meta <- pick_col(meta, c("tier", "tier_class"), required = FALSE)

tier_gene_col <- pick_col(tier, c("gene", "symbol", "genesymbol", "hgnc_symbol"))
tier_m_col <- pick_col(tier, c("m", "pooled_m", "effect", "estimate", "meta_effect"))
tier_fdr_col <- pick_col(tier, c("fdr", "padj", "bh_fdr", "qvalue"), required = FALSE)
robust_col <- pick_col(tier, c("robust_core", "robust", "is_robust_core", "robustcore"), required = FALSE)

if (is.na(tier_col_meta)) {
  meta$tier_plot <- case_when(
    meta[[fdr_col]] < 0.05 & abs(meta[[m_col]]) >= 0.5 ~ "Tier1",
    meta[[fdr_col]] < 0.05 ~ "Tier2",
    TRUE ~ "Tier3"
  )
} else {
  meta$tier_plot <- as.character(meta[[tier_col_meta]])
}

meta$tier_plot <- factor(meta$tier_plot, levels = c("Tier3", "Tier2", "Tier1"))

if (is.na(robust_col)) {
  robust <- tier %>%
    filter(abs(.data[[tier_m_col]]) >= 0.5) %>%
    slice_head(n = 66)
} else {
  robust <- tier %>%
    filter(as.logical(.data[[robust_col]]))
}

if (nrow(robust) == 0) {
  stop("No robust-core rows found. Check the robust-core column in ", tier_file)
}

robust <- robust %>%
  mutate(
    direction = if_else(.data[[tier_m_col]] > 0, "Up-regulated", "Down-regulated"),
    gene_plot = fct_reorder(.data[[tier_gene_col]], .data[[tier_m_col]])
  )

sig_path <- paths %>%
  mutate(
    significant = fdr < 0.05,
    pathway_label = str_replace_all(pathway, "_", " ")
  ) %>%
  arrange(M) %>%
  mutate(pathway_label = factor(pathway_label, levels = pathway_label))

path_pair <- paths %>%
  select(pathway, M) %>%
  left_join(pair_paths %>% select(pathway, logfc_gse11691, same_sign), by = "pathway") %>%
  mutate(
    pathway_label = str_replace_all(pathway, "_", " "),
    concordance = if_else(same_sign, "Same direction", "Opposite direction")
  ) %>%
  arrange(M) %>%
  mutate(pathway_label = factor(pathway_label, levels = pathway_label))

s6_tested <- s6 %>%
  filter(status == "TESTED")

mural <- s6_tested %>%
  filter(score == "tier1_signature", lineage %in% c("Smooth muscle cells", "Mural cells", "Perivascular cells")) %>%
  arrange(FDR) %>%
  slice(1)

serp <- s6_tested %>%
  filter(score == "SERPINE1_logCPM", lineage == "Endothelial cells") %>%
  slice(1)

if (nrow(mural) == 0) stop("Mural/smooth-muscle Tier-1 validation row not found.")
if (nrow(serp) == 0) stop("Endothelial SERPINE1 validation row not found.")

# -----------------------
# Fig. 1: design
# -----------------------

primary <- contrasts %>%
  filter(level == "primary", cohort != "GSE11691") %>%
  transmute(
    cohort,
    n = str_extract(note, "^[0-9]+\\s+vs\\s+[0-9]+"),
    comparison = case_when(
      cohort == "GSE7305" ~ "Endometriosis-associated tissue vs normal endometrium",
      cohort == "GSE6364" ~ "Eutopic endometrium: endometriosis vs control",
      cohort == "GSE25628" ~ "Patient tissue (eutopic + ectopic) vs normal endometrium",
      cohort == "GSE7307" ~ "Endometriosis-associated tissue vs normal endometrium",
      cohort == "GSE51981" ~ "Endometriosis vs non-EMS/no pelvic pathology",
      cohort == "GSE120103" ~ "Stage IV ovarian EMS vs fertile control",
      TRUE ~ note
    )
  )

fig1a <- primary %>%
  mutate(cohort = factor(cohort, levels = rev(cohort))) %>%
  ggplot(aes(y = cohort)) +
  geom_segment(aes(x = 0, xend = 1), colour = "grey88", linewidth = 7, lineend = "butt") +
  geom_point(aes(x = 0.08), size = 8, colour = blue) +
  geom_text(aes(x = 0.08, label = cohort), colour = "white", fontface = "bold", size = 3.2) +
  geom_text(aes(x = 0.22, label = comparison), hjust = 0, size = 3.1, colour = dark) +
  geom_text(aes(x = 0.97, label = n), hjust = 1, size = 3.1, colour = "grey30") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Discovery meta-analysis: six heterogeneous bulk contrasts",
    subtitle = "Random-effects meta-analysis; primary contrasts are locked before effect-size extraction",
    x = NULL, y = NULL
  ) +
  theme_plate(10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )

fig1b <- tibble(
  stage = c(
    "Bulk discovery",
    "Sensitivity",
    "Pathway meta-analysis",
    "scRNA localization",
    "External validation",
    "Perturbational analysis"
  ),
  detail = c(
    "Six primary cohorts; random-effects gene-level meta-analysis",
    "LOCO ×6; sans-GSE7307; fixed-effects; frozen QC swaps",
    "Six prespecified axes; BH correction across pathway family",
    "Three scRNA datasets; sample/donor-level pseudobulk",
    "GSE213216: patient-level pseudobulk; 195 planned tests",
    "LINCS reversal hypotheses; confirmatory BH-FDR gate"
  ),
  x = 1:6
) %>%
  ggplot(aes(x = x, y = 0)) +
  geom_segment(aes(x = 1, xend = 6), linewidth = 0.8, colour = "grey60",
               arrow = arrow(type = "closed", length = unit(2.5, "mm"))) +
  geom_point(size = 12, colour = blue_light) +
  geom_text(aes(label = stage), fontface = "bold", size = 3.2, vjust = 0.2) +
  geom_text(aes(label = detail), size = 2.5, vjust = -2.4, lineheight = 0.95) +
  scale_x_continuous(limits = c(0.4, 6.7)) +
  coord_cartesian(ylim = c(-1.2, 0.55), clip = "off") +
  labs(
    title = "Locked evidence framework",
    subtitle = "Discovery, robustness, localization, validation and perturbation stages are analytically separated"
  ) +
  theme_plate(10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )

fig1c <- tibble(
  evidence = c(
    "Primary claim",
    "Robustness support",
    "External validation",
    "Exploratory evidence"
  ),
  definition = c(
    "BH-FDR < 0.05 within the prespecified family",
    "Cohort-removal / model / QC-swap analyses reported transparently",
    "Independent atlas; patient-level pseudobulk; planned test matrix",
    "Reported separately; no confirmatory claim"
  )
) %>%
  ggplot(aes(y = fct_rev(evidence), x = 1)) +
  geom_tile(fill = c(red, blue, purple, gold), width = 0.95, height = 0.75) +
  geom_text(aes(label = evidence), colour = "white", fontface = "bold", size = 3.4) +
  geom_text(aes(x = 1.55, label = definition), hjust = 0, size = 3.05, colour = dark) +
  coord_cartesian(xlim = c(0.45, 4.5), clip = "off") +
  labs(title = "Interpretive hierarchy", x = NULL, y = NULL) +
  theme_plate(10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )

Fig1 <- (fig1a / fig1b / fig1c) +
  plot_layout(heights = c(1.25, 1.05, 0.85)) +
  plot_annotation(tag_levels = "a")

save_plate(Fig1, "Fig1", 12.5, 10.0)

# -----------------------
# Fig. 2: gene meta-analysis
# -----------------------

fig2a <- meta %>%
  mutate(
    neglog10_fdr = -log10(pmax(.data[[fdr_col]], 1e-300)),
    tier_plot = factor(tier_plot, levels = c("Tier3", "Tier2", "Tier1"))
  ) %>%
  ggplot(aes(x = .data[[m_col]], y = neglog10_fdr, colour = tier_plot)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = dark, linewidth = 0.45) +
  geom_point(alpha = 0.55, size = 0.65) +
  scale_colour_manual(values = c("Tier3" = "grey75", "Tier2" = blue_light, "Tier1" = red)) +
  labs(
    title = "Cross-cohort random-effects meta-analysis",
    subtitle = "All analysed genes; dashed line denotes BH-FDR = 0.05",
    x = expression("Pooled meta effect, M (log"[2]*" fold change)"),
    y = expression("-log"[10]*" BH-FDR"),
    colour = NULL
  ) +
  theme_plate(10) +
  theme(legend.position = c(0.88, 0.13))

top_tier <- tier %>%
  arrange(desc(abs(.data[[tier_m_col]]))) %>%
  slice_head(n = 12) %>%
  mutate(
    direction = if_else(.data[[tier_m_col]] > 0, "Up", "Down"),
    gene_plot = fct_reorder(.data[[tier_gene_col]], .data[[tier_m_col]])
  )

fig2b <- top_tier %>%
  ggplot(aes(x = .data[[tier_m_col]], y = gene_plot, colour = direction)) +
  geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.4) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Up" = red, "Down" = blue)) +
  labs(
    title = "Selected Tier-1 genes",
    subtitle = "Ranked by absolute pooled effect size",
    x = expression("Pooled M (log"[2]*" fold change)"),
    y = NULL,
    colour = NULL
  ) +
  theme_plate(10) +
  theme(legend.position = "bottom")

fig2c <- meta %>%
  filter(!is.na(.data[[i2_col]])) %>%
  ggplot(aes(x = .data[[i2_col]])) +
  geom_histogram(binwidth = 5, fill = "grey78", colour = "white") +
  geom_vline(xintercept = median(meta[[i2_col]], na.rm = TRUE), colour = purple, linewidth = 0.9) +
  annotate(
    "text",
    x = median(meta[[i2_col]], na.rm = TRUE),
    y = Inf,
    label = paste0("Median I² = ", round(median(meta[[i2_col]], na.rm = TRUE), 1), "%"),
    hjust = 1.05, vjust = 1.6, colour = purple, fontface = "bold", size = 3.2
  ) +
  labs(
    title = "Between-cohort heterogeneity",
    subtitle = "Reported for every gene; never used as a filter",
    x = expression(I^2*" (%)"),
    y = "Genes"
  ) +
  theme_plate(10)

Fig2 <- (fig2a | fig2b | fig2c) +
  plot_layout(widths = c(1.2, 0.9, 0.85)) +
  plot_annotation(tag_levels = "a")

save_plate(Fig2, "Fig2", 15.5, 5.6)

# -----------------------
# Fig. 3: robustness
# -----------------------

sens_n_col <- pick_col(sens, c("tier1_n", "n_tier1", "tier1_count", "n_retained", "retained"))
sens_label_col <- pick_col(sens, c("scenario", "scenario_id", "name", "analysis"))

sens_plot <- sens %>%
  transmute(
    scenario = .data[[sens_label_col]],
    n_tier1 = as.numeric(.data[[sens_n_col]])
  ) %>%
  mutate(
    scenario = str_replace_all(scenario, "_", " "),
    scenario = fct_reorder(scenario, n_tier1)
  )

fig3a <- sens_plot %>%
  ggplot(aes(x = scenario, y = n_tier1)) +
  geom_col(fill = blue) +
  geom_hline(yintercept = 485, linetype = "dashed", colour = red, linewidth = 0.8) +
  geom_text(aes(label = n_tier1), hjust = -0.15, size = 3.2) +
  coord_flip(clip = "off") +
  labs(
    title = "Tier-1 set stability across prespecified scenarios",
    subtitle = "Dashed line: base Tier-1 set (n = 485)",
    x = NULL,
    y = "Tier-1 genes retained"
  ) +
  theme_plate(10) +
  theme(plot.margin = margin(7, 35, 7, 7))

robust_n <- nrow(robust)
robust_up <- sum(robust[[tier_m_col]] > 0, na.rm = TRUE)
robust_down <- sum(robust[[tier_m_col]] < 0, na.rm = TRUE)

fig3b <- robust %>%
  ggplot(aes(x = seq_len(n()), y = .data[[tier_m_col]], colour = direction)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.45) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c("Up-regulated" = red, "Down-regulated" = blue)) +
  labs(
    title = paste0("Robust core: ", robust_n, " genes"),
    subtitle = paste0(robust_up, " up-regulated; ", robust_down, " down-regulated; retained across all prespecified cohort-removal scenarios"),
    x = "Robust-core genes ranked by pooled M",
    y = expression("Pooled M (log"[2]*" fold change)"),
    colour = NULL
  ) +
  theme_plate(10) +
  theme(legend.position = "top")

label_candidates <- robust %>%
  arrange(desc(abs(.data[[tier_m_col]]))) %>%
  slice_head(n = min(10, nrow(robust)))

fig3b <- fig3b +
  geom_text(
    data = label_candidates,
    aes(label = .data[[tier_gene_col]]),
    colour = dark, size = 2.7, vjust = ifelse(label_candidates[[tier_m_col]] > 0, -0.8, 1.3),
    show.legend = FALSE
  )

fig3c <- tibble(
  category = c("Base Tier-1", "Robust core"),
  n = c(485, robust_n),
  label = c("Base Tier-1", "Retained in all LOCO + sans-GSE7307")
) %>%
  ggplot(aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), vjust = -0.3, fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Base Tier-1" = blue_light, "Robust core" = purple)) +
  labs(
    title = "A restricted core withstands cohort removal",
    x = NULL, y = "Genes"
  ) +
  ylim(0, 550) +
  theme_plate(10) +
  theme(legend.position = "none")

Fig3 <- (fig3a | fig3b | fig3c) +
  plot_layout(widths = c(1.05, 1.35, 0.75)) +
  plot_annotation(tag_levels = "a")

save_plate(Fig3, "Fig3", 16.0, 5.8)

# -----------------------
# Fig. 4: pathway evidence
# -----------------------

fig4a <- sig_path %>%
  mutate(
    lo = M - 1.96 * SE,
    hi = M + 1.96 * SE,
    significance = if_else(significant, "BH-FDR < 0.05", "Not significant")
  ) %>%
  ggplot(aes(x = M, y = pathway_label, colour = significance)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "grey75", linetype = "dotted") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.55) +
  geom_point(size = 3) +
  geom_text(
    aes(label = paste0("I²=", round(I2), "%")),
    x = 0.55, hjust = 0, size = 2.7, colour = "grey35", show.legend = FALSE
  ) +
  scale_colour_manual(values = c("BH-FDR < 0.05" = red, "Not significant" = grey)) +
  coord_cartesian(xlim = c(-1.55, 0.95), clip = "off") +
  labs(
    title = "Pathway-level random-effects meta-analysis",
    subtitle = "Only the proliferation axis meets BH-FDR < 0.05; dotted lines are interpretive |M| = 0.5",
    x = expression("Pooled pathway score (log"[2]*" fold change; 95% CI)"),
    y = NULL,
    colour = NULL
  ) +
  theme_plate(10) +
  theme(legend.position = "bottom")

fig4b <- path_pair %>%
  ggplot(aes(x = M, y = logfc_gse11691, colour = concordance, label = pathway_label)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey55") +
  geom_point(size = 3.2) +
  geom_text_repel <- NULL

fig4b <- fig4b +
  geom_text(
    nudge_y = 0.08, size = 2.8, check_overlap = TRUE, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Same direction" = blue, "Opposite direction" = red)) +
  labs(
    title = "Reserved paired arm: directional comparison",
    subtitle = "GSE11691 was analysed separately and was never pooled",
    x = expression("Six-cohort pooled M"),
    y = expression("Paired GSE11691 log"[2]*" fold change"),
    colour = NULL
  ) +
  theme_plate(10) +
  theme(legend.position = "bottom")

fig4c <- sig_path %>%
  ggplot(aes(x = fct_reorder(pathway_label, I2), y = I2, fill = significant)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = paste0(round(I2), "%")), vjust = -0.25, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = red, `FALSE` = purple)) +
  labs(
    title = "Heterogeneity is substantial across all pathway tests",
    x = NULL,
    y = expression(I^2*" (%)")
  ) +
  ylim(0, 105) +
  theme_plate(10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "none"
  )

Fig4 <- (fig4a | fig4b | fig4c) +
  plot_layout(widths = c(1.25, 1, 0.8)) +
  plot_annotation(tag_levels = "a")

save_plate(Fig4, "Fig4", 16.0, 5.8)

# -----------------------
# Fig. 5: discovery localization
# -----------------------

loc_img <- c(
  file.path(root, "figures", "fig4_localization.jpg"),
  file.path(root, "figures", "main", "fig4_localization.jpg"),
  file.path(root, "fig4_localization.jpg")
)

loc_img <- loc_img[file.exists(loc_img)]

if (length(loc_img) > 0 && requireNamespace("magick", quietly = TRUE)) {
  img <- magick::image_read(loc_img[1])
  fig5a <- ggplot() +
    annotation_custom(rasterGrob(as.raster(img), interpolate = TRUE), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    labs(
      title = "Discovery-stage single-cell pseudobulk localization",
      subtitle = "Axis-gene localization is reported as a discovery-stage analysis"
    ) +
    theme_blank_plate()
} else {
  fig5a <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.58,
      label = "Insert final scRNA pseudobulk localization panel\nfrom fig4_localization.jpg",
      size = 5, colour = "grey35", hjust = 0.5
    ) +
    annotate(
      "text", x = 0.5, y = 0.40,
      label = "Required source image not found in figures/ or magick package unavailable.",
      size = 3.5, colour = red, hjust = 0.5
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_blank_plate()
}

sc_sig_file <- find_file(c("^S3_scrna_zmean_tests.*\\.csv$", "S3_scrna_zmean_tests"))
sc_sig <- read_csv(sc_sig_file, show_col_types = FALSE)
names(sc_sig) <- tolower(names(sc_sig))

sc_fdr <- pick_col(sc_sig, c("fdr", "padj", "bh_fdr"), required = FALSE)
sc_p <- pick_col(sc_sig, c("p", "p_value", "pvalue"), required = FALSE)

n_sig <- if (!is.na(sc_fdr)) sum(sc_sig[[sc_fdr]] < 0.05, na.rm = TRUE) else NA_integer_
n_total <- nrow(sc_sig)

fig5b <- tibble(
  metric = c("Measured gene-by-lineage tests", "BH-FDR significant tests"),
  n = c(n_total, n_sig)
) %>%
  ggplot(aes(x = metric, y = n, fill = metric)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.35, fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Measured gene-by-lineage tests" = blue_light, "BH-FDR significant tests" = red)) +
  labs(
    title = "Localization evidence is limited and explicitly bounded",
    subtitle = "Use full test table as Supplementary Data",
    x = NULL, y = "Tests"
  ) +
  theme_plate(10) +
  theme(legend.position = "none")

fig5c <- tibble(
  statement = c(
    "Discovery localization",
    "Patient/sample-level pseudobulk",
    "SERPINE1 finding",
    "Interpretation"
  ),
  text = c(
    "Three single-cell datasets were used for discovery-stage localization.",
    "The biological sample/donor—not individual cells—was the statistical unit.",
    "Endothelial/perivascular localization is discovery-stage only.",
    "External lineage-level SERPINE1 replication was not supported."
  )
) %>%
  ggplot(aes(y = fct_rev(statement), x = 1)) +
  geom_tile(fill = c(blue, blue_light, gold, grey), width = 0.95, height = 0.76) +
  geom_text(aes(label = statement), colour = "white", fontface = "bold", size = 3.1) +
  geom_text(aes(x = 1.55, label = text), hjust = 0, size = 3.0, colour = dark) +
  coord_cartesian(xlim = c(0.5, 5.2), clip = "off") +
  labs(title = "Pre-specified interpretation boundary", x = NULL, y = NULL) +
  theme_plate(10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank()
  )

Fig5 <- (fig5a | (fig5b / fig5c)) +
  plot_layout(widths = c(1.45, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")

save_plate(Fig5, "Fig5", 15.5, 8.0)

# -----------------------
# Fig. 6: external validation
# -----------------------

tier1_rows <- s6_tested %>%
  filter(score == "tier1_signature", contrast == "C1") %>%
  mutate(
    lineage = fct_reorder(lineage, delta_median),
    sig = FDR < 0.05
  )

fig6a <- tier1_rows %>%
  ggplot(aes(x = delta_median, y = lineage, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.45) +
  geom_point(size = 3) +
  geom_text(
    aes(label = paste0("FDR=", formatC(FDR, format = "f", digits = 3))),
    nudge_x = 0.06, hjust = 0, size = 2.6, colour = dark, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = grey), labels = c(`TRUE` = "BH-FDR < 0.05", `FALSE` = "Not significant")) +
  coord_cartesian(xlim = c(min(tier1_rows$delta_median, na.rm = TRUE) - 0.35,
                           max(tier1_rows$delta_median, na.rm = TRUE) + 0.55),
                  clip = "off") +
  labs(
    title = "Independent validation across lineages",
    subtitle = "Tier-1 directional signature; patient-level pseudobulks",
    x = expression("Median delta (lesion - eutopic)"),
    y = NULL,
    colour = NULL
  ) +
  theme_plate(10) +
  theme(legend.position = "bottom")

smooth_scores <- s6_scores %>%
  filter(
    lineage %in% c("Smooth muscle cells", "Mural cells", "Perivascular cells"),
    tissue %in% c("Eutopic Endometrium", "Endometrioma", "Peritoneal lesion")
  ) %>%
  mutate(
    group = case_when(
      tissue == "Eutopic Endometrium" ~ "Eutopic",
      TRUE ~ "Ectopic lesion"
    )
  )

fig6b <- smooth_scores %>%
  ggplot(aes(x = group, y = tier1_signature, colour = group)) +
  geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.12) +
  geom_jitter(width = 0.10, height = 0, size = 2.2, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.52, colour = dark, linewidth = 0.45) +
  scale_colour_manual(values = c("Eutopic" = blue, "Ectopic lesion" = red)) +
  annotate(
    "text", x = 1.5, y = max(smooth_scores$tier1_signature, na.rm = TRUE) + 0.15,
    label = paste0(
      "Mural/smooth-muscle Tier-1 signature\n",
      "Δ=", sprintf("%.2f", mural$delta_median),
      "; BH-FDR=", sprintf("%.3f", mural$FDR)
    ),
    size = 3.2, fontface = "bold"
  ) +
  labs(
    title = "Validated mural-compartment Tier-1 signal",
    subtitle = "C1: ectopic lesion versus eutopic endometrium",
    x = NULL, y = "Tier-1 directional signature"
  ) +
  theme_plate(10) +
  theme(legend.position = "none")

fig6c <- tibble(
  endpoint = c("Tier-1 signature\nmural/smooth muscle", "SERPINE1\nendothelium"),
  delta = c(mural$delta_median, serp$delta_median),
  FDR = c(mural$FDR, serp$FDR),
  verdict = c("Validated", "Not replicated")
) %>%
  ggplot(aes(x = endpoint, y = delta, fill = verdict)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.45) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = paste0("Δ=", sprintf("%.2f", delta), "\nFDR=", sprintf("%.3f", FDR))),
    vjust = ifelse(c(mural$delta_median, serp$delta_median) >= 0, -0.25, 1.15),
    size = 3.2
  ) +
  scale_fill_manual(values = c("Validated" = green, "Not replicated" = grey)) +
  labs(
    title = "Positive and negative endpoints are reported together",
    x = NULL, y = expression("Median delta")
  ) +
  theme_plate(10) +
  theme(legend.position = "none")

test_counts <- s6 %>%
  mutate(testable = if_else(status == "TESTED", "Tested", "Not testable")) %>%
  count(testable)

fig6d <- test_counts %>%
  ggplot(aes(x = testable, y = n, fill = testable)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.35, fontface = "bold", size = 4) +
  scale_fill_manual(values = c("Tested" = purple, "Not testable" = "grey75")) +
  labs(
    title = "Validation testability was prespecified",
    subtitle = "All planned endpoints are retained in the test matrix",
    x = NULL, y = "Planned tests"
  ) +
  theme_plate(10) +
  theme(legend.position = "none")

Fig6 <- (fig6a | fig6b | fig6c | fig6d) +
  plot_layout(widths = c(1.15, 1.1, 0.85, 0.7)) +
  plot_annotation(tag_levels = "a")

save_plate(Fig6, "Fig6", 18.0, 5.7)

# -----------------------
# Supplementary LINCS figure
# -----------------------

figs1 <- lincs %>%
  mutate(
    hypothesis = factor(hypothesis, levels = hypothesis),
    confirmed = verdict == "CONFIRMED"
  ) %>%
  ggplot(aes(x = hypothesis, y = BH_FDR, fill = confirmed)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = red, linewidth = 0.7) +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0("FDR=", sprintf("%.3f", BH_FDR))), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(`TRUE` = green, `FALSE` = grey)) +
  labs(
    title = "Pre-registered reversal hypotheses: none confirmed",
    subtitle = "Dashed line denotes BH-FDR = 0.05",
    x = NULL, y = "BH-FDR"
  ) +
  ylim(0, 1.12) +
  theme_plate(10) +
  theme(legend.position = "none")

save_plate(figs1, "FigS1_LINCS_reversal", 8.5, 5.0)

message(
  "\nFinal plates written to:\n", fig_dir,
  "\n\nGenerated:\n",
  "Fig1, Fig2, Fig3, Fig4, Fig5, Fig6, FigS1_LINCS_reversal\n"
)