# S4 knock summary -- cross-cohort random-effects meta-analysis

Generated: 2026-09-22 09:50:27

## Design
- Random-effects REML meta across the 6 locked primary disease-vs-control
  contrasts; GSE11691 (paired ectopic-vs-eutopic) analysed as a separate
  directional concordance arm (estimand separation).
- Presence rule k >= 4; Tier-1 = BH-FDR < 0.05 AND |M| >= 0.5 AND sign
  concordance >= 4 cohorts. Heterogeneity (tau2, I2, Q) reported for all.

## Primary results (measured)
- Genes entering meta (k >= 4): 21352
- Tier-1 core genes: 485
- Tier-2 (FDR-only) genes: 1292
- Tier-1 up / down: 121 / 364
- Median tau2 across Tier-1: 0.1832

## GSE11691 concordance arm
- Common genes: 13036 | spearman rho = 0.324 | Tier-1 sign-match = 0.803 (391 genes)

## Sensitivity battery
- sans-GSE7307: Tier-1 retention vs base = 0.2 | criterion >= 0.70 | FAIL
- LOCO:B-G7305: Tier-1 retention vs base = 0.305154639175258 | criterion >= 0.50 | FAIL
- LOCO:B-G6364: Tier-1 retention vs base = 0.942268041237113 | criterion >= 0.50 | PASS
- LOCO:B-G25628: Tier-1 retention vs base = 0.494845360824742 | criterion >= 0.50 | FAIL
- LOCO:B-G7307: Tier-1 retention vs base = 0.2 | criterion >= 0.50 | FAIL
- LOCO:B-G51981: Tier-1 retention vs base = 0.274226804123711 | criterion >= 0.50 | FAIL
- LOCO:B-G120103: Tier-1 retention vs base = 0.82680412371134 | criterion >= 0.50 | PASS
- FE-vs-RE: spearman(M_RE, M_FE) | Tier-1 Jaccard = 0.912 | 0.305 | criterion rho > 0.95 | FAIL
- QC-swap:B-G51981-S1: Tier-1 retention vs base = 1 | criterion >= 0.50 (reported) | PASS
- QC-swap:B-G51981-S2: Tier-1 retention vs base = 1 | criterion >= 0.50 (reported) | PASS

## scRNA projection (bulk core -> scRNA pseudobulk legs)
- Legs analysed: 27 | median sign-match rate: 0.507
- Legs with FDR < 0.05 (Wilcoxon): 20 of 27

## Immunoediting synthesis
- cytotoxicity: 0/6 in Tier-1 (up NA / down NA)
- antigen_presentation: 0/6 in Tier-1 (up NA / down NA)
- senescence_dormancy: 1/4 in Tier-1 (up  1 / down  0)
- proliferation: 0/3 in Tier-1 (up NA / down NA)
- stromal_ecm: 0/5 in Tier-1 (up NA / down NA)
- Escape Score: NA

## Artifact MD5s
- S4_meta_all_genes.csv: 703b7a15c96c797ecba435920c197052
- S4_meta_Tier1_core.csv: 3afbc9d85bec0b4409137cb73cc854bb
- S4_meta_low_presence.csv: 0a9a06b3ec26c0971815aa382491e2a0
- S4_gse11691_concordance.csv: b2d7bcb252c9935f335842ef4c94f8f6
- S4_sensitivity_summary.csv: da259446b8dacc2229b4e40be7c76e93
- S4_scrna_projection.csv: f99445f0625880efddd98ccda95ebd8c
- S4_scrna_lineage_readout.csv: 741550f4a21fadcff817d9e1c9c88bc0
- S4_immunoediting_synthesis.csv: fdc9187e4483b30bbb6603ccdb96f80d
