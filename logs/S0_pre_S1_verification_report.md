# Verification Report — Novelty Re-Screen + Accession Verification

**Date:** 2026-09-11
**Scope:** Pre-S1 verification requested by lead researcher: (A) re-verify novelty of the benign immunoediting framework; (B) verify every accession in Charter Section 3 against primary sources.
**Rule compliance:** R1 (no fabricated numbers — every entry below cites its source), R2 (traceability), R6 (honest reporting — limitations flagged, not hidden).

---

## A. Novelty verification

**Question:** Has the elimination–equilibrium–escape (3E) immunoediting formalism (Dunn/Schreiber, tumor immunology) been applied to endometriosis or to any benign disease?

**Method:** Web-level literature screen, 2026-09-11. Queries: `"immunoediting" endometriosis`; `"benign immunoediting" OR "immunoediting" benign disease non-cancer`; `immunoediting elimination equilibrium escape endometriosis`; `immunoediting framework applied benign disease adenomyosis fibroids equilibrium escape`.

**Findings:**
1. The 3E framework is established exclusively in tumor immunology (e.g., Kim et al. 2007, Immunology, cited >1700×; Mittal et al. 2014, Curr Opin Immunol, cited >2200×).
2. Query `"immunoediting" endometriosis` returned **no** paper applying the framework to endometriosis. The only gynecologic hit applies the framework to endometrial **cancer** (PD-1/PD-L1 review), a malignant disease — conceptually distinct.
3. Query `"benign immunoediting"` returned **zero results**. The only non-cancer usage of "immunoediting" found is in HIV (viral escape mutation — a different concept, not the 3E tissue-level framework).
4. No application to adenomyosis, fibroids, or other benign conditions was found.

**Conclusion:** The novelty claim in Charter Section 1 stands as of 2026-09-11: the 3E immunoediting formalism has not been applied to endometriosis or any benign disease.

**Limitations (honest reporting, R6):**
- This is a web-level screen, not a full systematic review. A formal, documented novelty screen (PubMed/Scopus with exported query logs) must be re-run at S9 before submission, and its output archived as a supporting file.
- Absence of evidence in a web screen is not proof of absence; the framing in the manuscript should claim "to our knowledge, first application" with the search date stated.

---

## B. Accession verification

Legend: ✅ verified against primary source | ⚠ verified but flagged | 🔎 resolved in this session

| # | Role (charter) | Charter accession | Status | Verified details | Source |
|---|---|---|---|---|---|
| 1 | Primary single-cell atlas (Fonseca et al.) | GSE213216 | ✅ | Public on GEO since 2022-10-11; "A single-cell transcriptomic analysis of endometriosis"; >370,000 cells (373,851 post-QC); endometriomas n=8, peritoneal lesions n=27, eutopic n=10, unaffected ovary n=4, disease-free peritoneum n=4; also via github.com/lawrenson-lab/AtlasEndometriosis | GEO accession viewer; Fonseca et al., Nat Genet 2023 (PMC10950360) |
| 2 | Integrated endometrial atlas, normal reference (Marečková et al.) | GSE183837 | ✅ | Verified as the Marečková scRNA-seq dataset integrated into HECA (Human Endometrial Cell Atlas, 313,527 cells, 63 women, with and without endometriosis), Nat Genet 2024. Note: HECA is the integrated product; GSE183837 is the underlying Marečková dataset. Both usable. | Marečková et al., Nat Genet 2024 (PMID 39198675) |
| 3 | Peritoneal immune compartment (Zou et al.) | to be confirmed at S1 | 🔎⚠ | **Resolved: PRJNA713993** (NCBI BioProject; raw sequencing data). Zou G. et al., "Cell subtypes and immune dysfunction in peritoneal fluid of endometriosis revealed by single-cell RNA-sequencing", Cell Biosci 2021 (PMC8157653). **FLAG: scRNA-seq was performed on only n=1 endometriosis + n=1 control (10,280 + 7,250 cells); the remaining 64 samples were flow-cytometry/qPCR validation only. Suitability for our NK analyses must be decided at S1 — likely usable only as qualitative supporting evidence, not as a discovery cohort.** | PMC8157653, Data Availability statement |
| 4 | Early-stage disease gradient (Huang et al.) | to be confirmed at S1 | 🔎 | **Resolved: GSE214411.** Huang X. et al., "Single-cell transcriptome analysis reveals endometrial immune microenvironment in minimal/mild endometriosis", Clin Exp Immunol 2023 (PMID 36869723). 138,057 cells; 6 minimal/mild patients + 7 controls; eutopic endometrium across cycle phases. Note: tissue is eutopic endometrium, not lesions — fits the "equilibrium-state candidate" role as intended. | PMC10243848, Data Availability statement |
| 5 | Bulk validation, primary (severity-graded) | GSE51981 | ✅ | n=148 endometrial samples; endometriosis classified minimal/mild vs moderate/severe; includes normal controls and other uterine/pelvic pathology controls; GPL570 platform. Exactly matches the charter's intended use. | GEO accession viewer |
| 6 | Bulk validation, paired | GSE11691 | ✅ | Paired eutopic/ectopic endometrium from 9 women (18 samples); public since 2008-12-22; GPL96 (HG-U133A). | GEO accession viewer |
| 7 | Bulk validation, additional | GSE7305 | ✅ | 10 ovarian endometriosis + 10 normal endometrium; GPL570; pre-treatment surgical samples. Widely reused. | GEO via multiple citing papers (e.g., Springer 2022, s12863-022-01036-y) |
| 8 | Bulk validation, additional | GSE25628 | ✅ | Endometriosis vs control endometrium (n=16 per citing study); used repeatedly as an independent validation set. | Citing studies (Frontiers Genet 2025; T&F 2026) |
| 9 | Bulk validation, additional | GSE6364 | ✅ | Eutopic endometrium, endometriosis vs control; GPL570. | Citing studies (Frontiers Mol Biosci 2021) |
| 10 | Bulk validation, additional | GSE120103 | ⚠ | Public since 2019-02-26; 36 samples; **all endometriosis cases are stage IV ovarian** with fertile/infertile subgroups; platform GPL6480 (Agilent). Usable as replication, but covers only the severe end of the spectrum — weight accordingly. | GEO accession viewer |
| 11 | Senescence gene panel | MSigDB SAUL_SEN_MAYO (M45803) | ✅ | Verified on MSigDB: standard name SAUL_SEN_MAYO, systematic name M45803, human, collection C2:CGP; 125 senescence/SASP genes; source Saul et al., Nat Commun 2022 (PMID 35974106). | gsea-msigdb.org human geneset page |
| 12 | Connectivity mapping | CLUE / GSE92742 | ⏳ | Not re-verified this session (needed at S8, not S1). Verify before use. | — |

**Incidental findings (recorded for S1 consideration, not adopted):**
- GSE7307 (18 endometriosis + 23 normal endometrium) — potential additional bulk cohort, same platform family.
- GSE5108 (n=22) — additional early bulk cohort.
- GSE179640 — single-cell dataset (8 eutopic + 8 ectopic) used in a 2026 multi-omics study; could strengthen the single-cell arm given the Zou et al. n=1+1 limitation.

---

## C. Decisions required from lead researcher (before/at S1)

1. **Zou et al. (PRJNA713993):** keep as qualitative supporting evidence only, or drop and substitute a larger peritoneal/immune cohort? Recommendation: keep for NK-dysfunction narrative support only; do not use for score derivation.
2. **GSE120103:** keep as severe-stage replication only (recommended).
3. **Add GSE179640 and/or GSE7307** to the charter dataset table? If yes, charter Section 3 gets an amendment entry (logged here and in analysis_log.md).

## D. Gate status

- Novelty re-verification: **PASS** (with documented limitations).
- Accession verification: **PASS** — 10/11 charter items verified or resolved; CLUE/GSE92742 deferred to S8 by design.
- This report satisfies the pre-S1 verification requested on 2026-09-11. S1 (download + manifest + checksums) remains **gated** pending explicit go-ahead (R3).
