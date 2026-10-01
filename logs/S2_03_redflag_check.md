# S2_03 Red-Flag Check — focused measurement

**Date:** 2026-09-18 · **Protocol:** one focused check per stage, graded 🔴 Blocker / 🔵 Note
**Focused question:** Do the 5 arrays flagged by the registered ≥2-of-3 rule coincide with
independent geometric evidence (PCA of the same normalized matrices), and is every visually
extreme array explained by its per-metric statistics? (Two independent instruments: AQM
outlierDetection statistics vs. PCA embedding of the same data.)

## Measurement (from uploaded result files: 6 per-cohort CSVs + 7 PCA PNGs)

All 5 flagged arrays are cohort-level geometric extremes in their PCA:

| flagged array | PCA position | stat vs threshold (boxplot KS) | stat vs threshold (heatmap dist) |
|---|---|---|---|
| GSE6364 GSM150191 | extreme right (PC1 ≈ +310) | 0.0339 > 0.0189 | 15.27 > 13.74 |
| GSE6364 GSM150214 | extreme right (PC1 ≈ +270) | 0.0336 > 0.0189 | 15.31 > 13.74 |
| GSE7307 GSM176240 | extreme left (PC1 ≈ −235) | 0.0421 > 0.0393 | 22.16 > 21.91 |
| GSE51981 GSM1256659_134 | extreme top (PC2 ≈ +245) | 0.0579 > 0.0464 | 107.62 > 107.22 |
| GSE11691 GSM296885 | extreme top-left | 0.0496 > 0.0224 (cohort max KS) | 9.69 > 6.43 (cohort max dist) |

Every visually extreme NON-flagged array is numerically explained (fails ≤1 metric):

| array | PCA position | fails | margin |
|---|---|---|---|
| GSE7305 GSM175775 | extreme top (PC2 ≈ +230) | heatmap only | KS 0.0083 < 0.0124 |
| GSE7305 GSM175783 | extreme right | boxplot only | dist 7.353 < 7.490 |
| GSE7307 GSM176127 | left periphery | heatmap only | KS 0.0329 < 0.0393 |
| GSE51981 GSM1256655_103 | extreme top (with _134) | boxplot only | dist 106.71 < 107.22 |
| GSE51981 GSM1256685_72 | far right (PC1 ≈ +220) | heatmap only (120.9 > 107.2, cohort max dist) | KS below threshold |
| GSE6364 GSM150196 | left periphery | 0 | — |

## Grading

- 🔵 **Borderline pair GSE51981 _103/_134:** the heatmap threshold (107.22) splits the two
  top-of-PCA arrays by <0.5 units; _134 crosses on both metrics, _103 on one. The registered
  rule is mechanical and pre-registered — it stands as applied; the pair is recorded here so
  S4 can treat them as a sensitivity pair if either is ever removed (removal needs its own
  documented reason per design §6).
- 🔵 **GSE25628 platform mislabel** (series tagged GPL570; CELs measured HG-U133A_2, 22277
  probe sets) — see verification report §5-F1; manifest correction delivered.
- 🔵 **Documentation correction** A.1.3 (neqc drops, not retains, control probes) — report §5-F3.
- 🔵 Cosmetic items (warnings counts, renv namespace notice, geom_text legend glyph) — report §5-F4.

**Zero 🔴 blockers. Concordance between the two independent instruments is complete.**
