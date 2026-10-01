rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(patchwork)
library(stringr)

ROOT <- "C:/endometriosis_immunoediting"
RES  <- file.path(ROOT, "results")
OUT  <- file.path(ROOT, "figures", "final_plates")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

find1 <- function(pattern) {
  x <- list.files(
    RES,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(x) == 0) stop("Missing file: ", pattern)
  x[1]
}

save_fig <- function(p, name, w, h) {
  ggsave(
    file.path(OUT, paste0(name, ".tiff")),
    p, width = w, height = h, units = "in",
    dpi = 600, compression = "lzw", bg = "white"
  )

  ggsave(
    file.path(OUT, paste0(name, ".pdf")),
    p, width = w, height = h, units = "in",
    device = "pdf", bg = "white"
  )
}

theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey35"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "grey92"),
      panel.grid.minor = element_blank()
    )
}

blue <- "#235784"
red <- "#C00000"
purple <- "#5B2C92"
grey <- "#9B9B9B"
green <- "#2E7D32"

meta <- read_csv(
  find1("^S4_meta_all_genes.*csv$"),
  show_col_types = FALSE
)

tier <- read_csv(
  find1("^S4_meta_Tier1_core.*csv$"),
  show_col_types = FALSE
)

sens <- read_csv(
  find1("^S4_sensitivity_summary.*csv$"),
  show_col_types = FALSE
)

path <- read_csv(
  find1("^S5_meta_pathways.*csv$"),
  show_col_types = FALSE
)

pair <- read_csv(
  find1("^S5_concordance_GSE11691.*csv$"),
  show_col_types = FALSE
)

s6 <- read_csv(
  find1("^S6_01_tests.*csv$"),
  show_col_types = FALSE
)

hyp <- read_csv(
  find1("^S8_01_hypotheses.*csv$"),
  show_col_types = FALSE
)

names(meta) <- tolower(names(meta))
names(tier) <- tolower(names(tier))
names(sens) <- tolower(names(sens))
names(path) <- tolower(names(path))
names(pair) <- tolower(names(pair))
names(s6) <- tolower(names(s6))
names(hyp) <- tolower(names(hyp))

pick <- function(dat, choices) {
  hit <- choices[choices %in% names(dat)]
  if (length(hit) == 0) {
    stop("Column missing: ", paste(choices, collapse = ", "))
  }
  hit[1]
}

gene <- pick(meta, c("gene", "symbol"))
M <- pick(meta, c("m", "pooled_m", "effect"))
FDR <- pick(meta, c("fdr", "padj"))
I2 <- pick(meta, c("i2"))

tier_gene <- pick(tier, c("gene", "symbol"))
tier_M <- pick(tier, c("m", "pooled_m", "effect"))

robust <- "robust_all_7"

if (!robust %in% names(tier)) {
  stop("Column robust_all_7 is missing from S4_meta_Tier1_core file.")
}

meta <- meta %>%
  mutate(
    tier_plot = case_when(
      .data[[FDR]] < 0.05 & abs(.data[[M]]) >= 0.5 ~ "Tier1",
      .data[[FDR]] < 0.05 ~ "Tier2",
      TRUE ~ "Tier3"
    ),
    neglog = -log10(pmax(.data[[FDR]], 1e-300))
  )

Fig2A <- ggplot(
  meta,
  aes(x = .data[[M]], y = neglog, colour = tier_plot)
) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_point(size = 0.65, alpha = 0.55) +
  scale_colour_manual(
    values = c(
      Tier1 = red,
      Tier2 = "#88A8D8",
      Tier3 = "grey75"
    )
  ) +
  labs(
    title = "Cross-cohort random-effects meta-analysis",
    subtitle = "Dashed line: BH-FDR = 0.05",
    x = "Pooled meta effect M (log2 fold change)",
    y = "-log10 BH-FDR",
    colour = NULL
  ) +
  theme_pub()

Fig2B <- ggplot(meta, aes(x = .data[[I2]])) +
  geom_histogram(binwidth = 5, fill = "grey75", colour = "white") +
  geom_vline(
    xintercept = median(meta[[I2]], na.rm = TRUE),
    colour = purple,
    linewidth = 1
  ) +
  labs(
    title = "Between-cohort heterogeneity",
    subtitle = paste0(
      "Median I2 = ",
      round(median(meta[[I2]], na.rm = TRUE), 1),
      "%"
    ),
    x = "I2 (%)",
    y = "Genes"
  ) +
  theme_pub()

Fig2 <- Fig2A | Fig2B
save_fig(Fig2, "Fig2", 13, 5.5)

rob <- tier %>%
  filter(.data[[robust]] == TRUE) %>%
  mutate(
    direction = ifelse(
      .data[[tier_M]] > 0,
      "Up-regulated",
      "Down-regulated"
    )
  ) %>%
  arrange(.data[[tier_M]]) %>%
  mutate(rank = row_number())

Fig3A <- ggplot(
  sens,
  aes(x = reorder(scenario, tier1_n), y = tier1_n)
) +
  geom_col(fill = blue) +
  geom_hline(yintercept = 485, linetype = "dashed", colour = red) +
  coord_flip() +
  labs(
    title = "Tier-1 stability across sensitivity scenarios",
    subtitle = "Dashed line: base Tier-1 set (n = 485)",
    x = NULL,
    y = "Tier-1 genes"
  ) +
  theme_pub()

Fig3B <- ggplot(
  rob,
  aes(x = rank, y = .data[[tier_M]], colour = direction)
) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_point(size = 2.2) +
  scale_colour_manual(
    values = c(
      "Up-regulated" = red,
      "Down-regulated" = blue
    )
  ) +
  labs(
    title = paste0("Robust core: ", nrow(rob), " genes"),
    x = "Robust-core genes ranked by pooled M",
    y = "Pooled M (log2 fold change)",
    colour = NULL
  ) +
  theme_pub()

Fig3 <- Fig3A | Fig3B
save_fig(Fig3, "Fig3", 13, 5.5)

path <- path %>%
  mutate(
    label = str_replace_all(pathway, "_", " "),
    lo = M - 1.96 * SE,
    hi = M + 1.96 * SE,
    significance = ifelse(
      FDR < 0.05,
      "BH-FDR < 0.05",
      "Not significant"
    )
  ) %>%
  arrange(M) %>%
  mutate(label = factor(label, levels = label))

Fig4A <- ggplot(
  path,
  aes(x = M, y = label, colour = significance)
) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      "BH-FDR < 0.05" = red,
      "Not significant" = grey
    )
  ) +
  labs(
    title = "Pathway-level random-effects meta-analysis",
    subtitle = "Only proliferation meets BH-FDR < 0.05",
    x = "Pooled pathway score (95% CI)",
    y = NULL,
    colour = NULL
  ) +
  theme_pub()

pair <- path %>%
  select(pathway, M) %>%
  left_join(pair, by = "pathway") %>%
  mutate(
    label = str_replace_all(pathway, "_", " "),
    direction = ifelse(
      same_sign,
      "Same direction",
      "Opposite direction"
    )
  )

Fig4B <- ggplot(
  pair,
  aes(
    x = M,
    y = logfc_gse11691,
    colour = direction,
    label = label
  )
) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_text(nudge_y = 0.06, size = 3, check_overlap = TRUE) +
  scale_colour_manual(
    values = c(
      "Same direction" = blue,
      "Opposite direction" = red
    )
  ) +
  labs(
    title = "Reserved paired cohort concordance",
    subtitle = "GSE11691 was not pooled into the meta-analysis",
    x = "Six-cohort pooled M",
    y = "Paired GSE11691 log2 fold change",
    colour = NULL
  ) +
  theme_pub()

Fig4 <- Fig4A | Fig4B
save_fig(Fig4, "Fig4", 13, 5.5)

tier1 <- s6 %>%
  filter(
    status == "TESTED",
    score == "tier1_signature",
    contrast == "C1"
  ) %>%
  mutate(lineage = reorder(lineage, delta_median))

Fig6A <- ggplot(
  tier1,
  aes(
    x = delta_median,
    y = lineage,
    colour = FDR < 0.05
  )
) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c(
      `TRUE` = red,
      `FALSE` = grey
    )
  ) +
  labs(
    title = "Independent external validation",
    subtitle = "Tier-1 directional signature across lineages",
    x = "Median delta: lesion minus eutopic",
    y = NULL,
    colour = "BH-FDR < 0.05"
  ) +
  theme_pub()

serp <- s6 %>%
  filter(
    status == "TESTED",
    score == "SERPINE1_logCPM",
    lineage == "Endothelial cells"
  ) %>%
  slice(1)

mural <- s6 %>%
  filter(
    status == "TESTED",
    score == "tier1_signature",
    lineage == "Smooth muscle cells"
  ) %>%
  slice(1)

Fig6B <- data.frame(
  endpoint = c(
    "Tier-1 signature\nsmooth muscle",
    "SERPINE1\nendothelium"
  ),
  delta = c(mural$delta_median, serp$delta_median),
  FDR = c(mural$FDR, serp$FDR),
  verdict = c("Validated", "Not replicated")
) %>%
  ggplot(aes(x = endpoint, y = delta, fill = verdict)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0("FDR=", sprintf("%.3f", FDR))),
    vjust = -0.35,
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      Validated = green,
      `Not replicated` = grey
    )
  ) +
  labs(
    title = "Positive and negative validation endpoints",
    x = NULL,
    y = "Median delta"
  ) +
  theme_pub() +
  theme(legend.position = "none")

Fig6 <- Fig6A | Fig6B
save_fig(Fig6, "Fig6", 13, 5.5)

FigS1 <- ggplot(hyp, aes(x = hypothesis, y = bh_fdr)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = red) +
  geom_col(fill = grey, width = 0.65) +
  geom_text(
    aes(label = sprintf("FDR=%.3f", bh_fdr)),
    vjust = -0.35,
    size = 3.3
  ) +
  labs(
    title = "Pre-registered reversal hypotheses: none confirmed",
    subtitle = "Dashed line: BH-FDR = 0.05",
    x = NULL,
    y = "BH-FDR"
  ) +
  ylim(0, 1.1) +
  theme_pub()

save_fig(FigS1, "FigS1_LINCS", 8, 5)

message("DONE: figures saved in ", OUT)