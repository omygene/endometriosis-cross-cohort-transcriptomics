# S2 Pre-Registered Design — QC + Integration

**Registered:** 2026-09-11, BEFORE any S2 execution (charter rule R5: criteria fixed before seeing results).
**Scope:** QC, filtering, integration, and annotation for all single-cell and bulk datasets acquired in S1.
**Binding:** deviations from this document require a logged amendment with justification; post-hoc threshold tuning after inspecting outcomes is forbidden (R5/R6).

---

## 1. Entry checks (from S1 report, pre-registered)

| # | Check | Pass criterion |
|---|---|---|
| E1 | Re-count tar members at first local load (GSE51981, GSE179640, GSE214411) | GSE51981 = 148; GSE179640 = 63 files (see Amendment S2-A below); GSE214411 = 39 files (13 samples × 3) |
| E2 | GSE7305 ↔ GSE7307 subset overlap | compare GSM titles/annotations; if overlapping donors suspected → flag affected samples, exclude one copy in S4 sensitivity analysis |
| E3 | GSE25628 metadata parse | 22 samples assigned case/control labels from GEO characteristics; ambiguities documented in S2 report |

### Amendment S2-A (2026-09-12): E1 correction for GSE179640 — documented BEFORE QC execution

**Trigger:** first local load returned 63 tar members vs the pre-registered 59 → pre-registered stop condition fired; execution halted and escalated.
**Root cause (documentation error, not data error):** the value 59 was the GSM *sample* count, mistakenly registered as the *file* count. The official GEO filelist (`logs/S2_GSE179640_GEO_filelist.txt`, fetched 2026-09-12 from `ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/suppl/filelist.txt`) lists **63 files for 59 samples**: 57 samples × 1 file (33 single-cell H5 + 24 bulk `featurecounts.txt.gz`) + 2 single-cell samples (GSM6102564_EOR02, GSM6102566_EOR04) × 3 files (10x barcodes/features/matrix trios) = 63. The local archive matches the official record exactly (63/63). **E1 for GSE179640 therefore PASSES against the corrected criterion 63.**
**Consequences registered here (before any QC output is seen):**
1. GSE179640 single-cell arm = 35 samples: 33 H5 + 2 mtx-trio samples (EOR02, EOR04) handled by the same trio-staging loader as GSE214411/GSE183837.
2. The 24 bulk `featurecounts` files are an author-provided bulk subseries **not part of the pre-registered analysis plan**; their presence is recorded in the manifest and they are not analysed in S2–S4 (any future use requires a new amendment).
3. **Load-time exclusion rule (fixed now):** GEO's `features.tsv.gz` for EOR02 (124 B) and EOR04 (140 B) are implausibly small (normal: ~10^5–10^6 B). Both samples are loaded with a guard: if the features file is <1000 bytes or `Read10X` fails, the sample is **excluded at load time** and the exclusion is written to `results/S2_GSE179640_load_exclusions.csv`. This is a data-integrity exclusion, not an outcome-based one.
4. Windows `bsdtar` failed to list `GSE7305_RAW.tar` (status 1) during the first E2 run, so that overlap listing was inconclusive; E2 is re-run with R's internal tar reader (`scripts/S2_00b_entry_checks_fix.R`) before S2_01.

## 2. Single-cell QC thresholds (per dataset, fixed now)

Global gene filter: genes detected in ≥3 cells. Cells outside all limits are removed; counts before/after reported per sample.

| Dataset | min genes/cell | max genes/cell | max mito % | Hb filter | Rationale |
|---|---|---|---|---|---|
| GSE213216 (processed rds) | report-only (authors filtered) | report-only | report-only | report-only | Author-processed objects; we audit and report distributions, no re-filtering. Stated limitation. |
| GSE179640 | 200 | 6000 | 20% | <5% | 10x tissue standard; paper's own QC used similar bounds |
| GSE214411 | 200 | 6000 | 10% | <5% | Paper used mito <10%; adopted as-is for comparability |
| GSE183837 (Lai RIF) | 200 | 6000 | 20% | <5% | 10x endometrium standard |

Doublets: scDblFinder per sample (expected doublet rate from 10x loading, homotypic adjustment on). Doublets removed; rate per sample reported. If scDblFinder fails on a sample, fall back to per-sample nCount-based cutoff and FLAG in report (pre-registered fallback).

## 3. Integration design (fixed)

- HVGs: 2000 (vst) per dataset before merge.
- Method: **Harmony** (theta=2, default) on top 30 PCs; batch variable = sample; secondary covariate = dataset.
- Sensitivity arm (pre-registered): repeat with scVI if Harmony mixing fails (see gate §5). scVI: 30 latent dims, 400 epochs max, early stopping; seed=42.
- GSE179640 hormonal treatment: included as covariate in downstream comparisons; sensitivity analysis with treated samples excluded — pre-registered, results of both reported.
- Seeds: set.seed(42) globally; all stochastic steps record their seed in the script header.

## 4. Annotation (fixed marker panels)

Major lineages validated against canonical markers AND the original publications' labels:
- Stromal: PDGFRB, DCN, COL1A1
- Epithelial: EPCAM, KRT8, KRT19 (ciliated: FOXJ1, PIFO; secretory: PAEP, CXCL14)
- Endothelial: PECAM1, VWF
- Perivascular: RGS5, MCAM
- NK/T: PTPRC, CD3D, NCAM1, KLRD1, NKG7, GZMB
- Myeloid: CD68, CD14, LYZ, ITGAX
- B/plasma: MS4A1, IGHG1
- Mast: TPSB2, KIT

Validation rule: annotation accepted per dataset only if ≥90% of clusters carry marker support AND are consistent with the source publication's labels at major-lineage level; conflicts documented cluster-by-cluster in the validation table (not silently resolved).

## 5. S2 QC gate (charter Section 4) — pass criteria fixed now

1. Per-dataset QC metrics (cells retained, genes/cell, mito %) reported for every sample.
2. Batch mixing: cLISI (batch) and iLISI (cell type) reported; integration accepted if major lineages are iLISI-separable and no sample forms an exclusive cluster (visual + LISI). If failed → scVI arm; if both fail → analyze datasets separately and report honestly.
3. Annotation validation table complete with ≥90% rule (§4).
4. All figures (UMAPs, QC violins) and the QC report saved to `figures/` and `results/S2_qc_report.pdf`.

## 6. Bulk cohorts QC (fixed)

- GPL570 cohorts (GSE51981, GSE7305, GSE7307-subset, GSE25628, GSE6364): RMA (oligo), arrayQualityMetrics; sample flagged outlier if ≥2 of 3 AQ metrics fail (pre-registered); outliers reported, removed from S4 only with documented reason.
- GPL96 (GSE11691): RMA separately; never merged with GPL570 at raw level (cross-platform integration happens at score level in S4, pre-registered decision).
- Agilent GSE120103: limma neqc normalization; kept as advanced-stage replication only.
- Probe annotation: fixed platform brainarray/standard GPL mappings recorded with versions in the S2 report.

## 7. Honest-reporting clause

Any dataset failing its QC is not silently dropped: the failure, metrics, and decision (exclude/repair/downgrade) are written into `logs/S2_verification_report.md` before proceeding.
