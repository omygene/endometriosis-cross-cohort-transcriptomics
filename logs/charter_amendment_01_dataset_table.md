# Charter Amendment 01 — Dataset Table Update

**Date:** 2026-09-11
**Approved by:** Lead researcher (instruction: "اتبع توصياتك التي تحافظ على قوة البحث وقابلة للدفاع أمام المراجعين")
**Amends:** Charter Section 3 (Data Sources)
**Basis:** Verification report `logs/S0_pre_S1_verification_report.md` + accession verification of the two candidate additions (this session).

---

## 1. Decisions adopted

### D1 — Zou et al. peritoneal fluid scRNA-seq (PRJNA713993): downgraded to qualitative supporting evidence
- **Reason:** scRNA-seq was performed on n=1 endometriosis + n=1 control only (10,280 + 7,250 cells). Cannot support statistical claims or score derivation.
- **Role going forward:** narrative/mechanistic support for peritoneal NK-cell dysfunction (decreased GZMB/GNLY/GZMH, elevated XCL1/CCL3) — cited in Discussion, never used as a discovery or validation cohort.
- **Defensibility value:** prevents a reviewer from killing the paper on a single-cell n=1 comparison; shows we audited our sources.

### D2 — GSE120103: retained as advanced-stage replication only
- **Reason:** all endometriosis cases are stage IV ovarian; cannot inform the severity gradient.
- **Role going forward:** replication of the Immunoediting Score at the severe end of the spectrum; explicitly labeled as such in methods.

### D3 — GSE179640 added to the charter dataset table (single-cell arm)
- **Verified identity:** Tan Y, Flynn WF, Robson P, Luciano DE, Courtois ET — "Single cell analysis of endometriosis reveals a coordinated transcriptional program driving immunotolerance and angiogenesis across eutopic and ectopic tissues." GEO public 2022-05-05; BioProject PRJNA744463; SRA SRP327335; processed data available (GSE179640_RAW.tar, 758.5 MB, H5/MTX/TSV).
- **Content:** >122,000 cells across 14 individuals; 31 single-cell biopsies = 3 control endometrium + 9 eutopic endometrium + 8 ectopic peritoneal + 6 ectopic peritoneal-adjacent + 4 ectopic ovary; 10x Genomics, Illumina NovaSeq 6000 (GPL24676); plus bulk RNA-seq arm (24 biopsies) and organoid cell-hashing arm.
- **Role:** independent single-cell validation arm; peritoneal-lesion and peritoneal-adjacent samples directly serve the equilibrium-to-escape mapping; compensates for the Zou et al. limitation.
- **Caveat recorded:** cohort is hormonally treated (stated in the original study) — must be handled as a covariate/sensitivity analysis at S2, and stated in limitations.
- **Strategic note (novelty protection):** this is the closest conceptual neighbor in the literature — its title claims an "immunotolerant peritoneal niche." It does NOT use the elimination–equilibrium–escape formalism, compute any phase score, or map a transition. The manuscript's Positioning section must cite it explicitly and state the difference: descriptive immunotolerance vs. a formal, quantifiable, three-phase process with a validated score and a reversal strategy. Confronting it head-on raises defensibility; ignoring it is the main reviewer risk.

### D4 — GSE7307 added to the charter dataset table (bulk arm)
- **Verified identity:** Roth R. (Neurocrine Biosciences), Affymetrix U133 Plus 2.0 (GPL570) multi-tissue "human body index" series (677 samples, >90 tissue types), public 2007-04-09. The endometriosis-relevant subset (18 endometriosis + 23 normal endometrium, per Chen et al. 2022, BMC Genomics) is extracted at S1.
- **Role:** additional bulk validation; same platform (GPL570) as GSE51981 and GSE7305 → clean cross-cohort processing.
- **Caveat recorded:** it is a multi-tissue series, not a dedicated endometriosis study; the subset extraction (exact GSM list) must be documented in the dataset manifest with the criteria used.

---

## 2. Updated dataset table (supersedes Charter Section 3)

| Role in study | Dataset | Accession | Status |
|---|---|---|---|
| Primary single-cell atlas | Fonseca et al. — endometriomas, peritoneal lesions, eutopic, unaffected ovary, disease-free peritoneum | GSE213216 | Verified 2026-09-11 |
| Integrated endometrial atlas (normal reference) | Marečková et al. / HECA | GSE183837 | Verified 2026-09-11 |
| Independent single-cell validation arm | Tan/Flynn et al. — eutopic, ectopic peritoneal, peritoneal-adjacent, ectopic ovary (+bulk arm) | GSE179640 | **Added by Amendment 01** (hormonally treated cohort — covariate at S2) |
| Early-stage (equilibrium-candidate) gradient | Huang et al. minimal/mild eutopic endometrium | GSE214411 | Resolved 2026-09-11 |
| Peritoneal immune compartment | Zou et al. peritoneal fluid scRNA-seq | PRJNA713993 | Resolved 2026-09-11 — **qualitative supporting evidence only (n=1+1)** |
| Bulk validation (primary, severity-graded) | Tamaresis/Burney et al., n=148, minimal/mild vs moderate/severe | GSE51981 | Verified 2026-09-11 |
| Bulk validation (paired) | Paired eutopic/ectopic, 9 women | GSE11691 | Verified 2026-09-11 |
| Bulk validation (additional) | Ovarian endometriosis vs normal | GSE7305 | Verified 2026-09-11 |
| Bulk validation (additional) | Body-index subset: 18 endometriosis + 23 normal endometrium | GSE7307 | **Added by Amendment 01** (subset extraction documented at S1) |
| Bulk validation (additional) | Eutopic endometrium cohort | GSE25628 | Verified 2026-09-11 |
| Bulk validation (additional) | Eutopic endometrium cohort | GSE6364 | Verified 2026-09-11 |
| Bulk replication (advanced stage only) | Stage IV ovarian, fertile/infertile | GSE120103 | Verified 2026-09-11 — severe-end replication only |
| Senescence gene panel | SenMayo (Saul et al. 2022) | MSigDB SAUL_SEN_MAYO (M45803) | Verified 2026-09-11 |
| Connectivity mapping | LINCS L1000 / CLUE | clue.io; GSE92742 | Deferred to S8 |
| Structural follow-up | PDB / PubChem | To be selected at S8 | Unchanged |

---

## 3. Effect on the workflow

- S1 manifest now covers 10 downloadable datasets (9 GEO/BioProject + SenMayo gene-set file) instead of 8.
- No other charter section is modified. Rules R1–R9 and QC gates (Section 4) are unchanged.
- The revised-positioning argument against Tan/Flynn (GSE179640) "immunotolerant niche" is logged for the S9 manuscript outline.
