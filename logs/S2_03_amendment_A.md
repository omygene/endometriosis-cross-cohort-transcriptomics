# Amendment S2_03-A — Bulk QC: AQM input object and outlier extraction

**Registered:** 2026-09-18, BEFORE any S2_03 analysis run (charter R5).
**Parent design:** `logs/S2_preregistered_design.md` Section 6 (bulk cohorts QC) — all
thresholds and the outlier rule remain unchanged; this Amendment changes only the
*input object* of the QC instrument and registers the extraction method.

## 1. Trigger (measured, per R9/B.4 — not theoretical)

1. **Locked-source verification of arrayQualityMetrics 3.68.0** (the installed binary,
   `reference/locked_sources/arrayQualityMetrics`): the package has NO method for
   oligo's `ExpressionFeatureSet`. Its data-preparation methods cover exactly
   `AffyBatch`, `ExpressionSet`, `MAList`, `NChannelSet`, `RGList`. The v2 script's
   call `arrayQualityMetrics(<oligo raw object>, ...)` would therefore error inside its
   `tryCatch` and the stage would continue **without any AQM output** — the registered
   outlier rule would silently never be evaluated.
2. **Preflight measurement on the lead's machine** (`logs/S2_03_preflight.txt`,
   2026-09-18): `pd.hg.u133.plus.2` 3.12.0 and `pd.hg.u133a` 3.12.0 ARE installed;
   `hgu133plus2cdf` and `hgu133acdf` are NOT installed. Locked-source verification of
   affy 1.90.0 shows that `affy::ReadAffy` with a missing CDF package attempts a
   network install — the documented non-interactive hang risk. oligo's own behavior
   for a missing `pd.*` package is a clean immediate `stop()` (verified in the locked
   binary), so the oligo normalization route is safe as-is.

## 2. Registered change

**AQM input = the RMA-normalized `ExpressionSet`** (oligo route, unchanged), with
`do.logtransform = FALSE` (RMA output is already log2). Normalization is unchanged:
`oligo::rma` for GPL570 and GPL96; `limma::neqc` for Agilent GSE120103.

The outlier rule is **unchanged**: an array is flagged if **≥2 of 3** arrayQualityMetrics
outlier methods call it an outlier. The three methods remain:
between-array distribution (Kolmogorov–Smirnov), MA-plot (Hoeffding's D), and
between-array distances (sum of L1 distances). Outlier calls are extracted
**programmatically from the module objects** (`outlierDetection` slots `statistic`,
`threshold`, `which`, `description` — slot structure verified in the locked binary),
not by parsing the HTML report; the HTML report is still produced as the
human-readable QC deliverable. Flagged arrays are reported and are removed from S4
**only with a documented reason** (parent design §6, unchanged).

## 3. Registered limitations (to be carried into the Methods/Limitations text)

- Post-normalization, the distribution-based (KS) metric has reduced discriminative
  power because quantile normalization matches array distributions by construction.
  The MA-plot and distance metrics retain full power on normalized data.
- The raw-data-only AQM modules (RLE, NUSE, RNA degradation, PM/MM) are
  `AffyBatch`-only and are therefore not produced. Recovering them would require
  installing `hgu133plus2cdf`/`hgu133acdf` (an environment amendment) — **not done**;
  may be revisited later by a new amendment.

## 4. Registered additions and clarifications

- **GSE7307-subset** receives the same AQM treatment as the other GPL570 cohorts
  (the parent design §6 lists it among them; v2 had omitted it).
- **GSE11691 (GPL96)** also receives AQM on its RMA ExpressionSet for symmetric QC
  evidence. It is still **never merged with GPL570 at raw level**; cross-platform
  integration remains score-level in S4 (unchanged).
- **GSE120103 (Agilent)**: `neqc` only, per design §6 (no AQM registered for it).
- **GSE7307 group labels** are attached by GSM-accession regex extraction with a
  hard stop on any unmapped sample. (v2's strip-extension matching was measured on
  the real filenames in the preflight: 41/41 exact — kept as documented fallback
  behavior; the regex form is the delivered one.)
- **Execution model**: one cohort per short-lived R process chained by
  `scripts\run_S2_03.bat` (Addendum B.6/B.7 pattern — hangs on this machine occurred
  only in long-lived processes). Missing inputs are **fail-loud**, never silently
  skipped (GSE183837 lesson).

## 5. Scope guard

No threshold, gate, cohort list, or biological comparison is changed by this
Amendment. S2_03 remains a QC/normalization stage; no differential analysis is
performed here.

---

# Addendum S2_03-A.1 — Agilent `genes$Status` marshalling + classifier fixes

**Registered:** 2026-09-18, AFTER the smoke test (`logs/S2_03_smoke.txt`) and BEFORE
any S2_03 analysis run (charter R5). This Addendum changes **no design element** —
no threshold, no cohort list, no outlier rule, no comparison. It registers one
data-preparation detail for GSE120103 and two implementation-level bug fixes in
`S2_03_bulk_qc.R`, each forced by a measurement (R9).

## A.1.1 Trigger (measured)

Smoke leg 5 (GSE120103, 2 files): `read.maimages` succeeded, `genes$Status` was
absent, and `limma::neqc` stopped with **"Detection p values not found in the
data."** — the run verdict was "SMOKE DONE WITH FAILURES".

## A.1.2 Root cause (verified in the locked binary, B.4)

In the installed limma 3.68.5 (`reference/locked_sources/limma`):

- `neqc` calls `nec` for background correction; `nec` routes on `genes$Status`:
  - **Status present** → `normexp.fit.control`, which matches `tolower(Status)`
    against the defaults `negctrl = "negative"`, `regular = "regular"` and stops
    cleanly with "No regular probes found" / "Fewer than two negative control
    probes found" if the labelling is wrong — i.e. it fails loudly, never
    silently.
  - **Status absent** → `normexp.fit.detection.p`, which requires a detection
    p-value column (`x$other`, default name "Detection"). These Agilent FE files
    carry none, producing the exact measured stop.
- `read.maimages(source = "agilent", green.only = TRUE)` on these files returns
  `genes` columns `Row, Col, ControlType, ProbeName, SystematicName` (measured in
  the smoke log) — i.e. no automatic Status marshalling happens for this FE
  format/version combination.

## A.1.3 Registered data-preparation detail (GSE120103 only)

`rg$genes$Status` is marshalled from `genes$ControlType`:
`ControlType == 0 → "regular"`, `ControlType == -1 → "negative"`, anything else
`→ "other"`. The **0/-1 coding is not assumed**: it is measured first by
`scripts/S2_03_agilent_diag.R` (table of ControlType on 2 real files) and the main
script re-checks it fail-loud (hard stop if either code is absent or if
`ControlType` itself is missing). If the diagnostic refutes the coding, the stage
stops and escalates — no fallback guessing.

Note carried forward: `neqc` retains control probes in its output object; whether
to drop non-regular probes is an S2_04 (downstream analysis) decision, recorded
here so it is not mistaken for an oversight.

## A.1.4 Implementation-level fixes in the outlier extractor (no design change)

Both were exposed by the smoke run's module dump:

- **F-A** — `outlierDetection@description` is a **length-2 vector** (metric text +
  `"data-driven"`/`"fixed"`), so `grepl` on it returns length 2 and `if (grepl(...))`
  errors under R ≥ 4.2 ("condition has length > 1"). Fix: collapse to a single
  string before pattern matching.
- **F-B** — the MA-plot module's description is literally `"<i>D<sub>a</sub></i>"`
  (Hoeffding's D, threshold 0.15 fixed); the string "Hoeffding" never occurs in
  it, so the v3 classifier could not identify the maplot module and would have
  hit its own FATAL stop on every cohort. Fix: identification pattern widened to
  `"Hoeffding" | "D<sub>a</sub>" | "D_a"`. The measured descriptions (smoke log)
  remain the anchor: "Kolmogorov-Smirnov" → boxplot, `"<i>D<sub>a</sub></i>"` →
  maplot, "sum of distances" → heatmap.

Both fixes touch only how the **already-produced** AQM module objects are read;
the ≥2-of-3 rule, the three metrics, and all thresholds are unchanged.

## A.1.5 Scope guard

Unchanged from Amendment S2_03-A §5: no threshold, gate, cohort list, or
biological comparison is modified. S2_03 remains a QC/normalization stage.
