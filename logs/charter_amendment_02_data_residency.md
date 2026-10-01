# Charter Amendment 02 — Data Residency & Storage Architecture

**Date:** 2026-09-11
**Status:** Adopted 2026-09-11 (lead researcher delegated the decision: "اتبع توصياتك" / no preference — assistant recommendation adopted per standing instruction).
**Amends:** Charter Section 6 ("All files live under /mnt/agents/output/endometriosis_immunoediting/")

## Reason
Two platform constraints discovered during S1 acquisition:

1. **100 MiB per-file cap** on the shared project storage (measured: writes stop/fail at exactly 104,857,600 bytes). Every archive >100 MB must be stored as 90 MB split parts. For GSE213216 (16.8 GB raw) this would mean ~187 parts — impractical and fragile.
2. **Outbound bandwidth throttling** of the sandbox after heavy downloads (~30–60 KB/s observed after ~6 GB), making multi-GB acquisition from the sandbox unreliable. The lead researcher's own network/machine is unaffected.

## Adopted architecture (if confirmed)

- **Heavy raw data** (single-cell archives, microarray archives >100 MB, future large intermediate objects such as integrated Seurat .rds): primary copy lives on the **lead researcher's local machine** under `data/raw/` and `data/processed/`.
- **Shared project storage** (`/mnt/agents/output/endometriosis_immunoediting/`): keeps (a) all small raw inputs (≤100 MiB: GSE11691, GSE7305, GSE25628, SenMayo GMT), (b) split parts for mid-size archives already acquired, (c) all scripts, tables, figures, logs, manifests, and reports, (d) the analysis log.
- **Traceability is preserved regardless of residency (R2):** `tables/dataset_manifest.tsv` records accession, source URL, download date, exact byte size, and MD5 of every input file wherever it lives. For locally held files, the lead researcher runs `scripts/S1_download_datasets.sh`, which produces a local MD5 log that is uploaded and merged into the manifest.
- **Verification rule:** no dataset enters S2 unless its size and MD5 in the manifest match a downloaded file — on either storage tier.

## No other change
Rules R1–R9, QC gates, stage workflow S0–S9, and Amendment 01 dataset decisions are unchanged.
