# S2_03 Verification Report — independent re-check (R2)

**Date:** 2026-09-18
**Verifier:** session agent (R2: verified from uploaded files only, never from chat text)
**Stage under review:** S2_03 bulk QC + normalization (v3.1 pipeline + lead-authored GSE120103 fix)
**Inputs verified:** 7 cohort logs, collect log, agilent diagnostic log, 6 per-cohort outlier CSVs,
merged outlier CSV, 7 PCA PNGs, lead-authored fix script + its log, amendment copy,
S2_02 red-flag measurement log.

## 1. Governance chain — VERIFIED

- Amendment copy on the lead's machine is **byte-identical** to the registered version:
  `logs/S2_03_amendment_A.md` MD5 `bfd31bfb80303a63df6af810a89ea9c0` (includes Addendum A.1).
- Lead's written approval of Amendment S2_03-A + Addendum A.1 and of the 5-array
  GSE120103 exclusion was given explicitly in the lead's stage-completion message (item ③),
  same session. Recorded here per R3/R5.
- Lead-authored fix script `scripts/S2_03_GSE120103_fix.R` MD5 `8cd35692c14b5b31fdabc570fd464a8b`
  matches the entry in `analysis_log.md` (2026-09-17). Code review: mirrors v3.1 Agilent leg,
  fail-loud guards (stops if exclusion count ≠ expected or ControlType coding ≠ 0/-1),
  writes only the registered outputs. Compliant with R4/R9.

## 2. Execution integrity — VERIFIED

Timeline is a coherent single-machine chain of short-lived processes (B.6/B.7):
GSE7305 21:05→21:07 · GSE6364 21:08→21:11 · GSE25628 21:11→21:13 · GSE7307 21:14→21:17 ·
GSE51981 21:18→21:30 · GSE11691 21:31→21:32 · GSE120103 failed-loud 21:33 ·
fix 23:20→23:20 · collect 23:22.

| cohort | files read | matrix | platform pkg loaded | AQM outliers (box/ma/heat) | flagged |
|---|---|---|---|---|---|
| GSE7305 | 20/20 ✓ | 54675 × 20 | pd.hg.u133.plus.2 | 1 / 0 / 1 | 0 |
| GSE6364 | 37/37 ✓ | 54675 × 37 | pd.hg.u133.plus.2 | 4 / 0 / 3 | 2 |
| GSE25628 | 22/22 ✓ | **22277 × 22** | **pd.hg.u133a.2** (see §4-F1) | 1 / 0 / 0 | 0 |
| GSE7307 | 41/41 ✓ | 54675 × 41 | pd.hg.u133.plus.2 | 1 / 0 / 3 | 1 |
| GSE51981 | 148/148 ✓ | 54675 × 148 | pd.hg.u133.plus.2 | 2 / 0 / 3 | 1 |
| GSE11691 | 18/18 ✓ | 22283 × 18 | pd.hg.u133a | 2 / 0 / 1 | 1 |
| GSE120103 | 31 (5 excluded, §3) | 43376 × 31 | limma/neqc | not registered (design §6) | — |

- File counts match the pre-registered CFG at every cohort (fail-loud guard never fired falsely).
- GSE7307 groups attached from the frozen S1 list: endometriosis 18 / normal_endometrium 23 ✓.
- R10 check: GSE51981 peak Vcells 3.26 GB vs 4–5 GB estimate — within envelope. No OOM anywhere.
- Every cohort log ends with its DONE marker; collect log reports 286 arrays / 6 cohorts ✓.

## 3. GSE120103 design-mismatch exclusion — VERIFIED (measured, documented)

- First run failed loudly at `GSM3393522_..._GE1-v5_95_Feb07_...`: *"fewer rows than files
  previously read"* — direct measurement that 5 arrays are a different chip design
  (GE1-v5_95_Feb07, fewer probes) that cannot be column-bound with the 31 GE1_107_Sep09 arrays.
- Fix log: exactly 5 GE1-v5 files excluded (names listed), 31 read, `RGList dim 45015 × 31`.
- ControlType coding measured on the full set: `-1=153 / 0=43376 / 1=1486 / NA=0` — identical
  to the 2-file diagnostic (21:04) that preceded it. Status marshalling applied per Addendum A.1.
- neqc output: 43376 × 31, E range 4.62 .. 18.85. RAM trivial (117 MB max).
- Decision recorded by the lead in `analysis_log.md` with the script MD5. **Governance complete.**

## 4. Outlier rule — VERIFIED numerically (from the CSVs themselves)

- Merged `S2_bulk_outliers.csv`: 286 rows = 20+37+22+41+148+18 ✓; column contract as registered.
- Merged table is value-identical to the concatenation of the 6 per-cohort CSVs ✓.
- `flagged ⟺ n_fail ≥ 2` holds on all 286 rows; `n_fail` equals the sum of the three `out_*`
  columns on all rows ✓.
- Exactly 5 flagged: GSE6364 GSM150191 + GSM150214 · GSE7307 GSM176240 ·
  GSE51981 GSM1256659_134 · GSE11691 GSM296885 — matches the lead's report ✓.
- maplot fired 0 times in all cohorts; threshold constant 0.15 (fixed, as locked-binary-verified).
  Consistent with the conservative fixed threshold on normalized data; KS-metric's reduced
  post-normalization power was already registered in Amendment §3.
- Per design §6: flagged arrays are **reported** here; removal from S4 requires a separate
  documented reason. No removal happened in S2_03.

## 5. Findings carried forward

### 🔵 F1 — GSE25628 chip geometry differs from the series platform tag (MEASURED)
The 22 CELs self-identify as **HG-U133A_2** (oligo loaded `pd.hg.u133a.2`; matrix 22277 probe
sets) — not HG-U133_Plus_2 (GPL570, 54675) as the GEO series tag / manifest line implies.
Impact on S2_03: none (oligo auto-detects from CEL headers; RMA used the correct design package).
Impact on S4: gene-space mapping for GSE25628 must use the U133A 2.0 annotation; cohort stays
score-level-only in cross-platform integration (rule already in force).
**Action:** manifest line for GSE25628 platform corrected to
`Affymetrix HG-U133A_2 (GPL571 geometry; series tagged GPL570)` with this measurement cited.

### 🔵 F2 — Environment changed on the fly (to be locked)
oligo fetched and installed `pd.hg.u133a.2` 3.12.0 from Bioconductor during the GSE25628 leg
(download completed, install DONE, no hang). The renv library now differs from the lockfile.
**Action:** `renv::snapshot()` after this stage (joins the already-pending harmony 2.0.5 lock).

### 🔵 F3 — Documentation correction to Addendum A.1.3
A.1.3 carried the note that neqc *retains* control probes. Measurement refutes it:
45015 rows in → 43376 out (= regular probes only; the 153 negative + 1486 positive controls
are dropped by neqc). No downstream filtering decision needed. Correct the note at next
amendment-file touch; no design impact.

### 🔵 F4 — Cosmetic/log-level notes (R9: never blockers)
- Per-cohort "There were N warnings" lines: contents not captured in logs; consistent with
  known ks.test tie warnings post-quantile. To be read at next session if any doubt arises.
- GSE120103: one `normexp.signal` numerical-accuracy notice (documented limma behavior for
  very-low-intensity probes; adjusted intensities floored).
- renv namespace notice (BiocManager loaded pre-activation from user cache): cosmetic; all
  analysis packages resolved from the renv library.
- GSE7307 quick-PCA legend shows a stray "a" key glyph (geom_text legend artifact) — cosmetic;
  superseded by the journal figures script.

## 6. Archive completeness

- `S2_02_redflag_check.txt` received and reviewed: 233459 cells ✓, 29 clusters ✓,
  exclusive-cluster candidates 26/27/28 all n < 500 ✓ (0 exclusive clusters under the
  registered rule), cluster 28 stromal-panel supported (median 2.67 vs 0 in non-stromal) ✓.
  Archive gap closed.

## Verdict

All S2_03 outputs verified against files; governance chain complete; zero blockers
(see `logs/S2_03_redflag_check.md`). **S2_03 is recommended for closure.**
