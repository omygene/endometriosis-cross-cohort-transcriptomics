# Benign Immunoediting in Endometriosis — Analysis Record & Methods Evidence Dossier

**Project:** An Elimination–Equilibrium–Escape (3E) Immunoediting Framework for Endometriosis Resolved by Single-Cell Analysis and Computational Reversal
**Compiled:** 2026-09-29 · **Sources:** every value below is traceable to a saved file (Charter, dated amendments, verification reports, versioned scripts with MD5 manifests, and locked result tables) per rules R1/R2. No value is reconstructed from memory.
**Companion documents:** `Endometriosis Project Charter.docx` (approved 2026-09-11); `analysis_log.md`; `logs/charter_amendment_01_dataset_table.md`; `logs/charter_amendment_02_data_residency.md`; `logs/S0_pre_S1_verification_report.md`; `logs/S1_verification_report.md`; `logs/S2_03_amendment_A.md` (+Addendum A.1); `logs/S2_03_redflag_check.md`; `logs/S2_03_verification_report.md`; `logs/S3_design.md`; `logs/S3_to_submission_master_plan.md`.

---

## 1. Governance, Preregistration, and Audit Trail

**The analysis charter was approved by the lead researcher on 2026-09-11, before any data analysis began.** All analyses were conducted under binding rules R1–R9 (no fabricated numbers; full traceability; confirmation gates; reproducibility with locked seeds; pre-registered QC gates; honest reporting; supporting files; division of labor; periodic accuracy review).

Binding design statements:
- **Sole statistical claim gate: BH-FDR < 0.05.** Effect-size tier |log2FC| ≥ 0.5 is interpretive only and never used for exclusion.
- **Single-cell integration:** Harmony (θ = 2 per covariate; covariates = sample and dataset), npcs = 50, dims = 1:30, clustering resolution = 0.5, **random seed = 42** (also set immediately before every scDblFinder call per sample).
- **External validation:** strictly patient-level (unit = patient); ≥20 cells per patient-by-lineage pseudobulk; arms with <4 patients declared NOT_TESTABLE and dropped symmetrically.
- **Every stage change after the Charter is a dated, numbered amendment approved in writing** (Charter Amendment 01/02; S2_03-A + Addendum A.1; S2-A; S2-B/B3; dated in-script amendments). Nothing was silently changed.

**Recommended wording (Methods):**
> "Analyses were conducted according to a time-stamped, version-controlled analysis charter finalized on 11 September 2026 before data analysis, with all subsequent changes recorded as dated, numbered charter amendments approved in writing by the lead researcher."

**Recommended wording (Data & Code Availability):**
> "The analysis charter, its dated amendments, and the complete analysis record (scripts, manifests, checksum registries, verification reports, and locked result files) are deposited at [repository/DOI]. Depositing the charter with an external timestamp (e.g., OSF/Zenodo) is recommended at submission; until then the claim 'preregistered' in the strict sense is intentionally not made."

---

## 2. Execution Timeline (code- and report-anchored; every date verifiable in a saved file)

| Date | Event | Source |
|---|---|---|
| 2026-09-11 | Charter + rules approved; accession verification (10/11 items, CLUE deferred to S8 by design); Amendment 01 (dataset table) + 02 (data residency); **S1 acquisition COMPLETE — QC gate PASS, all 13 items** | Charter; S0/S1 verification reports |
| 2026-09-12 | S0 environment lock (R 4.6.1 / Bioconductor 3.23) | analysis_log.md |
| 2026-09-13 | S2_01 v4: three QC defects fixed (scDblFinder seed policy); S2_05 rescue audit (negative) | script headers |
| 2026-09-14 | Propagation diagnostic lifts GSE179640 blocker; GSE183837 scoped separately (S2_01b) | script headers |
| 2026-09-16 | S2_02 RAM-isolated staged integration; smoke test measured | script headers |
| 2026-09-17 | GSE120103 Agilent fix; journal figure request (PI) | verification reports |
| 2026-09-18 | **S2_03 bulk QC complete** (286 AQM arrays; Amendment S2_03-A registered before analysis); S3 design frozen | S2_03 reports; S3_design.md |
| 2026-09-19 → 21 | S3 hardening: hormonal census; gene/object/join censuses; presto gate | script headers |
| 2026-09-22 | S3_01 rewritten to SYMBOL-level DE (**failure ledger #18**); S4 audit gate; **S6_00 structural audit (PI written approval)** | S6_00 manifest |
| 2026-09-23 | **S8_00 LINCS audit 17:15:52 → S8_01 confirmation 17:17:49; S6_01/S6_02 locked (manifests)** | S8 manifests |
| 2026-09-25 | **S8b_01 dienogest (v1.0) — last in-silico reversal attempt** | S8b manifest |

---

## 3. Datasets (all verified against primary sources before analysis; S1 gate PASS)

All files downloaded **2026-09-11**; MD5 recorded per dataset (per-file MD5s for GSE213216 objects and the GSE7307 subset); GSE51981 MD5 identical on two independent machines (strongest verification).

| Cohort | Accession | Role | n analyzed | Key verification notes |
|---|---|---|---|---|
| Fonseca atlas | GSE213216 | Primary single-cell atlas | 6 author objects; 373,851 cells (structural audit OK) | per-file MD5s; authors' processed objects; no re-QC (audit-only by design) |
| Lai RIF cohort | GSE183837 | scRNA reference (RIF vs Ctrl at WOI) | 9 samples (3 Ctrl + 6 RIF) | **Charter mislabel corrected same day: Lai et al. 2022, not Marečková (manifest v2)** |
| Tan/Flynn | GSE179640 | Independent scRNA validation | 33 samples (GEO 59 records incl. bulk/organoid arms) | hormonally treated cohort — hormonal covariate arm mandatory |
| Huang | GSE214411 | Early-stage (equilibrium candidate) | 13 samples (6 EMS + 7 controls) | eutopic endometrium across cycle phases |
| Zou | PRJNA713993 | Qualitative supporting evidence ONLY | 2 SRA runs (n=1+1 scRNA) | **no GSE exists**; never used for score derivation (Amendment 01-D1) |
| GSE51981 | GSE51981 | Primary bulk (severity-graded) | 148 samples | dual-machine MD5 match; severity-NA retained |
| GSE11691 | GSE11691 | Paired bulk | 9 women (18 samples) | GPL96; never merged raw with GPL570 |
| GSE7305 / GSE7307 | — | Bulk validation | 20 (10/10); 41 (18/23 subset) | GSE7307 subset GSM list frozen; overlap with GSE7305 flagged as sensitivity item (E2) |
| GSE25628 | GSE25628 | Bulk validation | 22 CEL | **CELs self-identify as HG-U133A_2 (22,277 probe sets) though series tagged GPL570 — analyzed with GPL571 geometry (measured correction F1)** |
| GSE6364 | GSE6364 | Bulk validation | 37 CEL | GEO official 37 |
| GSE120103 | GSE120103 | Severe-stage replication only | 36 downloaded; **31 analyzed** | 5 GE1-v5 chip-design exclusions (fail-loud, documented); all stage IV ovarian |
| SenMayo | MSigDB M45803 | Senescence panel | 124 symbols in GMT (page lists 125 — discrepancy recorded, R6) | — |

---

## 4. Locked Computational Environment

R 4.6.1; Bioconductor 3.23; renv.lock (frozen; covers 416 packages). Key analysis packages: Seurat 5.5.1 / SeuratObject 5.4.0, harmony 2.0.5, lisi 1.0, edgeR 4.10.5, limma (incl. neqc), GSVA 2.6.6 (ssgseaParam), metafor (REML), GEOquery 2.80.0, oligo + pd.hg.u133.plus.2 / pd.hg.u133a / pd.hg.u133a.2, affy/affyPLM/gcrma, arrayQualityMetrics 3.68.0, scDblFinder, bluster, SingleCellExperiment, HDF5Array. **Post-S2_03 packages added under registered finding F2 (pd.hg.u133a.2 3.12.0) — renv snapshot updated; recorded.**

---

## 5. Single-Cell QC (S2_01 / S2_01b; thresholds frozen in `S2_preregistered_design.md` §2)

**Per-cohort thresholds (THRESH block, verbatim):** min.genes = 200, max.genes = 6,000 for all three cohorts; max.mt = 20% (GSE179640, GSE183837) and 10% (GSE214411); **max.hb = 5%** (all). `CreateSeuratObject(min.cells = 3, min.features = 0)`. **No nCount threshold; no ribosomal threshold; no ambient-RNA correction (no SoupX/CellBender/decontX anywhere in the pipeline — verified by full code search).** Doublets: scDblFinder (default expected rate; seed 42 immediately before each sample). Load-time exclusion: files < 1,000 bytes (S2-A #3).

**QC results (locked CSVs; sums independently reproduce the integration cell count):**

| Cohort | samples | before QC | after QC | retention | doublets (weighted) |
|---|---|---|---|---|---|
| GSE179640 | 33 | 159,834 | 109,398 | 68.4% | 7.33% |
| GSE183837 (Lai RIF) | 9 | 67,949 | 40,801 | 60.0% | 7.80% |
| GSE214411 (Huang) | 13 | 155,361 | 83,260 | 53.6% | 10.93% |
| **Total** | **55** | **383,144** | **233,459** | **60.9%** | — |

All 55 samples retained (no sample-level exclusions). Stratified retention by condition is archived (`S2_qc_stratified.csv`). GSE213216: structural audit of the 6 author objects (373,851 cells; authors' doublets 28,060 = 7.5%; metadata classes: Endometriosis 118,625 / Eutopic 118,144 / Endometrioma 59,709 / Unaffected ovary 51,927 / No endometriosis detected 25,446).

---

## 6. Integration & Annotation (S2_02, v2.4)

Harmony integration (θ = 2 per covariate; sample + dataset), npcs = 50, dims 1:30, resolution 0.5, seed 42, union-HVG; RAM-isolated staged parts (p2a–p2d); multi-layer Seurat 5 objects.

**Final object:** **29 clusters (0–28), 233,459 cells, 8 lineages** (locked 8-lineage panel):

| Lineage | Cells | | Lineage | Cells |
|---|---|---|---|---|
| Stromal | 116,352 | | Myeloid | 14,062 |
| Epithelial | 38,951 | | Perivascular | 9,448 |
| NKT | 33,780 | | B_plasma | 2,290 |
| Endothelial | 17,252 | | Mast | 1,324 |

**Batch-mixing metrics:** iLISI_lineage = 1.000 (perfect lineage mixing); cLISI_sample = 4.636; cLISI_dataset = 1.393 (documented as partially tautological; score-level integration used accordingly). **Annotation support: 28/29 clusters supported**; cluster 8 (Perivascular, n = 9,448) support = FALSE (ratio 1.017; ambiguous with Stromal — carried as a documented limitation). Cluster 28 (n = 87) is **exclusive** (96.6% from one patient sample, GSM6102556_E09_EcO; stromal panel supported) — retained as a documented patient-specific state.

---

## 7. Bulk QC & Normalization (S2_03; Amendment S2_03-A registered before analysis)

- Affymetrix GPL570/GPL96: `oligo::rma` → arrayQualityMetrics **on the RMA-normalized ExpressionSet** with `do.logtransform = FALSE` (input-object amendment; all thresholds unchanged). Agilent GSE120103: `limma::neqc` only (no AQM registered); genes$Status marshalled from ControlType (measured 0/−1 coding; fail-loud re-check); **neqc drops control probes (45,015 → 43,376)** — documentation correction registered (F3).
- **286 arrays across 6 cohorts** (GSE7305 20, GSE6364 37, GSE25628 22, GSE7307 41, GSE51981 148, GSE11691 18). Outlier rule: **≥2 of 3** AQM methods (between-array distribution = Kolmogorov–Smirnov; MA-plot = Hoeffding's D, fixed 0.15; between-array distances = sum of L1). Outliers extracted programmatically from module objects (not HTML).
- **Exactly 5 flagged arrays:** GSE6364 GSM150191, GSM150214; GSE7307 GSM176240; GSE51981 GSM1256659_134; GSE11691 GSM296885. **Reported, not removed** (removal requires a separate documented reason); independent PCA of the same matrices confirmed all five as geometric extremes, and every visually extreme non-flagged array failed ≤1 metric (complete concordance of two independent instruments — `S2_03_redflag_check.md`). Borderline pair _103/_134 documented (margin < 0.5 distance units); sensitivity arms S1/S2 later showed **Jaccard = 1** (no effect on conclusions).
- GSE120103: 5 GE1-v5 chip-design arrays excluded by fail-loud column-binding stop (36 → 31); frozen primary contrast 13 vs 9.
- GSE25628: measured HG-U133A_2 geometry (GPL571 annotation used in S3; F1).

---

## 8. Differential Analysis & Contrasts (S3; frozen in `S3_contrasts_frozen.tsv`, MD5 c23feac0…)

18 contrasts specified in the frozen design (19 rows executed including the three GSE25628 exploratory arms). Primary analyses: **limma eBayes on RMA/neqc log2 values per cohort; unpaired except GSE11691 (paired, duplicateCorrelation blocking per woman); GSE6364 with cycle-phase covariate; GSE51981 severity-NA retained (no documented reason to exclude).**

| Contrast | n (test vs ref) | Design |
|---|---|---|
| B-G7305 | 10 vs 10 | unpaired |
| B-G6364 | **21 vs 16** | unpaired + phase covariate |
| B-G25628 | 16 vs 6 | unpaired |
| B-G7307 | 18 vs 23 | unpaired (+ sans-GSE7307 sensitivity E2) |
| B-G51981 | 77 vs 34 | unpaired; cleanest control group |
| B-G11691 | 9 vs 9 | **paired (the only paired contrast)** |
| B-G120103 | 13 vs 9 | unpaired (post chip-design exclusions) |
| B-G51981-MILD / -SEV | 27 vs 34 / 48 vs 34 | secondary |
| B-G51981-S1 / -S2 / -S3 | 77 vs 34 / 77 vs 34 / 77 vs 71 | sensitivity |
| B-G120103-2A / -2B / -IC | 9 vs 9 / 4 vs 9 / 13 vs 9 | secondary / secondary / sensitivity |
| B-G25628-EXPL / -ECT / -EUT | 8 vs 8 / 8 vs 6 / 8 vs 6 | exploratory |
| B-G120103-FCIC | 9 vs 9 | exploratory (fertility effect itself) |

Gene mapping: SYMBOL-level (S3_01 rewritten after failure ledger #18); Affy via hgu133plus2.db / hgu133a.db / hgu133a2.db; Agilent via genes symbol column; probe collapse = probe with highest mean expression per symbol (max-mean). GSE11691 pairing file self-audited 9×2 (any mismatch = FATAL).

**scRNA arms:** per-dataset pseudobulk (unit = sample/patient) ≥10 cells per (lineage × sample); CPM>1 in ≥ half samples before voom; **dual arms per Amendment S2-B/B3 — arm A includes GSE214411 sample N1 (27 legs), arm B excludes N1 (9 legs); median arm-A/B Spearman = 0.969.** GSE179640 hormonal covariate arm mandatory + no-treated sensitivity arm (results of both arms reported). Wilcoxon per-lineage (24 legs) = sensitivity only (pseudoreplication caveat flagged on every row). SenMayo z-mean defined per cell on log1p atlas-wide z-scores.

**Locked knock summary (S3_knock_summary.md):** 19/19 contrasts completed; **11,496 genes FDR < 0.05 in ≥2 primary cohorts, 6,720 with unanimous direction (58.5%)**; sans-GSE7307 sensitivity (E2): 10,411 genes, 55.4% unanimous. Per-contrast gene counts archived (e.g., B-G51981: 17,353; B-G6364: 608; B-G11691: 615).

---

## 9. Cross-Cohort Meta-Analysis (S4)

Score-level integration only (no raw cross-platform merge). **REML random-effects meta-analysis (metafor::rma) per gene** with fixed-effects sensitivity, τ²/I²/Q reported (never used as filters). BH-FDR. Inputs registry-locked; **Tier-1 core = 485 genes** (`S4_meta_Tier1_core.csv`); full per-gene table 21,355 genes (`S4_meta_all_genes.csv`). Tier labels: interpretive only.

---

## 10. Pathway Meta-Analysis (S5)

Six pathways, k = 6 cohorts each, REML with full heterogeneity statistics; BH-FDR across the 6-test family; k_min = 4 presence rule; effect tiers interpretive. Locked table `S5_meta_pathways.csv`. (Tier per pathway archived; proliferation pathway = Tier1, others Tier3 with exact statistics in the locked CSV.)

---

## 11. External Validation (S6; patient-level by design, S6_00 structural audit PI-approved 2026-09-22)

- Unit = **patient**; ≥20 cells per patient-by-lineage pseudobulk; **arms with <4 patients declared NOT_TESTABLE and dropped symmetrically**; BH-FDR < 0.05 sole gate; six scores (5 locked axes + SenMayo) × patient × lineage; verdicts P1–P4/E3.
- GSE213216 (authors' object; 373,851 cells) used **score-level only** — no reprocessing, no reclustering, no filtering (declared event; structural audit only).
- S6_02 endothelial subcluster analysis (18 subclusters; 23,226 cells) explicitly labeled **EXPLORATORY** (dilution hypothesis).
- Locked manifests: S6_00 v1.0 (MD5 712e13e…), S6_01 v1.0 (7f42c9e…), S6_02 v1.0 (a1c07…).

---

## 12. Connectivity Mapping & Reversal (S8; scripts v1.0, 2026-09-23)

**Source:** iLINCS live API (audit 2026-09-23 17:15:52; confirmation 17:17:49); L1000 space 12,328 genes (978 landmarks); gene info from GSE92742. **Stage S8 was prespecified in the Charter (deferred by design, verified 2026-09-11); the specific knock-down confirmation queries were defined after primary results — stated explicitly.**

- **SERPINE1-KD:** trt_sh.cgs, library LIB_6, **96 h** — **7 cell lines: HA1E, HEPG2, HT29, A375, MCF7, PC3, A549** (7 signatures; n_overlap = 46 Tier-1 genes each).
- **Robust core = 66 genes (14 up / 52 down)**; 48/66 in L1000; **8 landmarks (EGR1, FOS, CCND1, LYPLA1, PLA2G4A, GADD45B, NNT, UGDH)**; Tier-1 ∩ landmarks = 46. **MKI67 absent from the 978 landmarks → proliferation arm = TOP2A + PCNA only** (documented).
- **Estradiol:** trt_cp, library LIB_5 (modern L1000); 136 signatures across **14 source cell lines**; aggregation = **median per cell line**; statistical unit = cell line; acceptance rule ≥3 signatures per line → **10 aggregated lines plotted: MCF7 (56), VCAP (11), A549 (9), HA1E (8), HCC515 (6), HEPG2 (5), HT29 (6), PC3 (5), A375 (3), NKDBA (3)**; lines dropped: FIBRNPC (2), NEU (2), NPC (2), ASC (1) — **this is the complete 14→10 accounting**. No dose/duration filter (documented limitation; KD arm fixed at 96 h). Legacy CMAP (LIB_2) retained as an audit arm only.
- **Statistics:** Pearson r_tier1 (locked Tier-1 M vs signature); Wilcoxon one-sided vs 0; **BH-FDR across the 4-hypothesis family; claim gate = FDR < 0.05 + direction matching the preregistration.**
- **Compounds with zero signatures documented:** tiplaxtinin, TM5441, PAI-039, forskolin.
- **Results (negative; reported per R6):** P1a p = 0.234, FDR = 0.313 · P1b p = 1.000, FDR = 1.000 · P2a p = 0.216, FDR = 0.313 · P2b p = 0.0322, FDR = 0.129 — **all NOT CONFIRMED.** Median r_tier1: KD = −0.057; estradiol = −0.049. Median proliferation arm: KD = −0.238; estradiol = +0.291. Locked outputs: `S8_01_kd_scores.csv`, `S8_01_estradiol_scores.csv`, `S8_01_hypotheses.csv`.

---

## 13. Dienogest Pharmacological Test (S8b_01 v1.0, 2026-09-25; fully exploratory, no claim gate)

**Background census (measured 2026-09-23 across GEO/eutils):** danazol, GnRH agonists, leuprolide, norethindrone, gestrinone — **zero usable datasets**. Dienogest was the last available in-silico option; **registered STOP RULE: this is the final in-silico reversal attempt — next step is in vitro validation.**

- **Cohorts:** GSE75423 = endometriotic stromal cells (ECSC); GSE75425 = normal endometrial stromal cells (NESC); treatment = dienogest + dibutyryl-cAMP; platform Agilent GPL13497.
- **Design:** paired — 4 patients × (vehicle, treated) = 8 samples per cohort (columns 1–4 vehicle, 5–8 treated, same patient order); paired logFC = mean(treated) − mean(vehicle) per patient.
- **Test:** Pearson r between locked Tier-1 signature and paired logFC; **measured Tier-1 = 290 genes**; bootstrap **gene-level resampling, B = 2000, seed = 42**, 95% CI by quantiles.
- **Results (exploratory):** ECSC r = 0.199 [0.071, 0.328]; NESC r = 0.224 [0.120, 0.316] — same direction as the lesion signature → **NO reversal**. ECSC-vs-NESC response correlation r = 0.627 across 11,439 genes. Key-gene table archived (incl. TOP2A, PCNA, IL6, SOCS3, FOS, FOSB, EGR1, C1QB).
- **Dose/duration are read from GEO series-matrix metadata at runtime (not hard-coded); exact values must be quoted from the GSE75423/GSE75425 GEO records in the manuscript.**

---

## 14. Supporting Files Inventory (deposit with the paper)

1. `dataset_manifest.tsv` (v2, all MD5s + download dates + corrections log)
2. `S3_contrasts_frozen.tsv` (MD5 c23feac0880fa6f832ee21a1b0c647f2) + `S3_GSE11691_pairs.tsv` + `S1_GSE7307_subset.txt` (frozen 41-GSM list)
3. QC: `S2_qc_GSE179640.csv`, `S2_qc_GSE183837.csv`, `S2_qc_GSE214411.csv`, `S2_qc_stratified.csv`, `S2_qc_GSE213216_audit.csv`, `S2_02_annotation_validation.csv`, `S2_02_lisi.csv`, `S2_02_exclusive_clusters.csv`, `S2_GSE179640_load_exclusions.csv`, `S2_rescue_audit_summary.csv`, `S2_bulk_outliers*.csv` (6 + merged, 286 arrays)
4. Locked results: `S3_knock_summary.md`, `S4_meta_Tier1_core.csv` (485 genes), `S4_meta_all_genes.csv` (21,355 genes), `S5_meta_pathways.csv` (6 pathways), `S6_01_scores_long.csv`, S8/S8b verdicts and scores CSVs
5. Registries + MD5 manifests: S4/S5/S6 input registries; S6_00/S6_01/S6_02 manifests; S8_00/S8_01/S8b_01 manifests
6. Governance: Charter, Amendments 01/02, S2-A, S2-B/B3, S2_03-A + Addendum A.1, verification reports (S0, S1, S2_03), `S3_design.md`, `S3_to_submission_master_plan.md`
7. Code: 74 R scripts + frozen design docs (MD5 manifests per stage)
8. `renv.lock` + sessionInfo (to be added at submission)

---

## 15. Honest Limitations Carried into the Manuscript

1. Novelty claim based on a web-level screen (2026-09-11); formal documented screen (PubMed/Scopus with exported query logs) scheduled at S9; manuscript wording "to our knowledge, first application" with search date.
2. Zou et al. scRNA n=1+1 → qualitative evidence only (never used for statistics).
3. GSE120103 covers only stage IV (severe-end replication by design; 2B n=4 flagged).
4. GSE179640 hormonally treated (covariate + sensitivity arms reported).
5. GSE25628 unpaired geometry (U133A_2; correct measured platform used).
6. SenMayo GMT 124 vs MSigDB page 125 symbols (discrepancy recorded).
7. GSE6364 frozen contrast 21 vs 16 (authoritative) vs 22/15 label enumeration — one-sample label difference recorded for reconciliation.
8. cLISI_dataset partially tautological (condition confounded with dataset) → score-level integration.
9. S8 negative results reported as found (R6); estradiol arm without dose/duration filter; MKI67 absent from L1000 landmarks.
10. `failure_ledger` file not included in this dossier (reference #18 verified from scripts); dienogest dose/duration to be quoted from GEO records.

---

*Every number in this dossier traces to a saved file. Nothing herein is simulated or reconstructed from memory (R1/R2).*

---

## Addendum (2026-09-29) — Dienogest Treatment Parameters (from GEO primary records)

**Verified directly from the GSE75423 GEO series page (primary source; the value the S8b script reads at runtime):**

> "Decidualization of cultured ECSCs was induced by **12 days-culture with a combination of 0.5 mM dibutyryl-cAMP and 100 nM dienogest** in 10% charcoal-stripped heat-inactivated FBS."

- **GSE75423 (ECSC mRNA):** 0.5 mM dibutyryl-cAMP + 100 nM dienogest, **12 days**; platform GPL13497 (Agilent-026652 Whole Human Genome 4x44K v2); 8 samples (4 non-decidualized + 4 decidualized); submission 2016; citation PMID 27412773 (Aoyagi Y, Nasu K).
- **GSE75425 (NESC mRNA):** sister subseries of the same SuperSeries **GSE75427** ("Decidualized and non-decidualized ECSCs and NESCs") under the same citation (PMID 27412773) and the same experimental design — same regimen (0.5 mM dibutyryl-cAMP + 100 nM dienogest, 12 days) applies; its individual series page should be cited as the confirming record at submission.
- Terminology note for the manuscript: GEO frames the arms as "decidualized vs non-decidualized"; the S8b script labels them treated vs vehicle (paired per patient). Both describe the same paired contrast and should be worded consistently.

This closes the last open factual item of the Methods (dose and duration).
