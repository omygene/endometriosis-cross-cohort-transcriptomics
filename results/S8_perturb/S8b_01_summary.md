# S8b_01 -- exploratory drug positive-control (machine-printed)

Date: 2026-09-25 13:16:26.147652

## Census verdict (pre-measured, see S8b design doc)
- Usable standard-of-care drug expression datasets on human endometriotic tissue: ONE
  (GSE75423, dienogest+cAMP on ECSCs, paired 4v4); GSE75425 normal-ESC counterpart.
- danazol / GnRH-agonist / leuphin: no usable human endometriosis expression dataset.

## Exploratory results (no claim gate; single small dataset)

| dataset | Tier-1 measured | r | boot 95% CI | direction |
|---|---|---|---|---|
| GSE75423_ECSC_dienogest | 290 | 0.199 | [0.071, 0.328] | NO reversal |
| GSE75425_NESC_dienogest | 290 | 0.224 | [0.12, 0.316] | NO reversal |

- correlation of dienogest responses ECSC vs NESC (11439 genes): r = 0.627
- key genes (see S8b_01_key_genes.csv)

## Machine verdict (exploratory, machine-printed)
**EXPLORATORY (no claim gate): dienogest+cAMP on endometriotic stromal cells (GSE75423, paired 4v4, single dataset) shifts the Tier-1 signature TOWARD the lesion direction (r = 0.199, bootstrap 95% CI 0.071 to 0.328), i.e. drug-reversal is NOT observed; normal stromal cells show the same pattern (r = 0.224). Documented as exploratory; the positive-control question for standard-of-care drugs is CLOSED in silico (census: no usable danazol/GnRH datasets exist).**

*STOP RULE (locked): this was the LAST in-silico attempt regardless of outcome. Next stage = experimental in-vitro perturbation (independent study).*

*Every number measured at run time from GEO downloads and the locked Tier-1 file.
 Exploratory: single dataset, 4 patient pairs, platform covers only part of Tier-1.*
