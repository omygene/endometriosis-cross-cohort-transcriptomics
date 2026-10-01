# Analysis Log — Benign Immunoediting in Endometriosis

Project root: `/mnt/agents/output/endometriosis_immunoediting/`
Charter: `/mnt/agents/output/Endometriosis Project Charter.docx` (approved by lead researcher on 2026-09-11)
Companion framing document reviewed: `Revised Proposals (Post-Novelty Screen) (1).docx` (Proposal B = this study; consistent with charter, no conflicts).
Outstanding document: `Data Availability and Feasibility Audit.docx` (not yet uploaded; accession verification performed independently on 2026-09-11 and documented in `logs/S0_pre_S1_verification_report.md`).

## Execution model

- Lead researcher runs R locally (R8). Assistant writes/reviews code, verifies numbers, maintains documentation and this log.
- No stage proceeds without explicit confirmation (R3). Verification report closes each stage.

---

## 2026-09-11 — Session 1

### Charter confirmation

- Charter rules R1–R9 and workflow S0–S9 approved explicitly by lead researcher.

### S0 (partial — infrastructure)

- Folder tree created: `data/raw`, `data/processed`, `scripts`, `results`, `figures`, `tables`, `logs`.
- `analysis_log.md` opened (this file).
- PENDING: environment lock. renv.lock + sessionInfo.txt must be generated on the lead researcher's R installation (script to be provided by assistant in S0 completion step).

### Novelty re-verification (requested by lead researcher, 2026-09-11)

- Web-level screen for applications of the elimination–equilibrium–escape immunoediting framework outside cancer.
- Result: no application of the immunoediting framework to endometriosis or any benign disease found. "immunoediting" + endometriosis returns only endometrial *cancer* literature. "benign immunoediting" returns zero hits.
- Conclusion: novelty claim stands as of 2026-09-11. Full details + sources in `logs/S0_pre_S1_verification_report.md`.
- Caveat: this is a web-level screen; repeat a formal screen before manuscript submission (S9).

### Accession verification (pre-S1 component, requested by lead researcher)

- All charter accessions checked against GEO/MSigDB/publications. Two previously unconfirmed accessions resolved:
- Zou et al. peritoneal fluid scRNA-seq → PRJNA713993 (BioProject). **Caveat flagged: scRNA-seq performed on n=1 endometriosis + n=1 control only.**
- Huang et al. minimal/mild endometriosis scRNA-seq → GSE214411.
- One flag: GSE120103 contains only stage IV ovarian endometriosis (usable, but limited severity range).
- Details + sources in `logs/S0_pre_S1_verification_report.md`.

### Charter Amendment 01 (approved by lead researcher, 2026-09-11)

- D1: Zou et al. (PRJNA713993) downgraded to qualitative supporting evidence only (scRNA n=1+1).
- D2: GSE120103 retained as advanced-stage replication only.
- D3: GSE179640 (Tan/Flynn et al.) added as independent single-cell validation arm; flagged as closest conceptual neighbor ("immunotolerant niche") — must be confronted in manuscript positioning.
- D4: GSE7307 added as bulk validation (GPL570; endometriosis subset extraction documented at S1).
- Full amendment: `logs/charter_amendment_01_dataset_table.md`.

### S1 progress — acquisition (2026-09-11)

- Platform constraint discovered: the project storage mount caps single files at 100 MiB (104,857,600 bytes). All archives >100 MB are stored as 90 MB split parts (`*.part_00`, `part_01`, ...); reassembly = `cat part_* > file`. Original-file MD5 is recorded in the manifest (computed on the intact archive before splitting).
- NCBI throttles single long connections; chunked HTTP range downloads used. The `ftp.ncbi.nlm.nih.gov` and `www.ncbi.nlm.nih.gov/geo/download/?acc=...&format=file` endpoints serve byte-different archives (same size, same 27 members, different tar container bytes) — each dataset's manifest MD5 names its endpoint.
- Completed + verified on storage:
- GSE11691 RAW.tar (62,914,560 B) — direct file.
- GSE25628 RAW.tar (43,048,960 B) — direct file.
- GSE7305 RAW.tar (91,064,320 B) — direct file.
- GSE120103 RAW.tar (109,291,520 B, MD5 9482b104a86f051ce86a1bca597f69b6) — 2 parts, reassembly MD5 verified; 36 tar members.
- GSE6364 RAW.tar (190,904,320 B, MD5 959690329da456c1f42a645c97dd3d09) — 3 parts, reassembly MD5 verified; 37 tar members.
- GSE183837 RAW.tar (517,140,480 B, MD5 ced92ae3a4c43b8c030a143faa4a280b, www endpoint) — 6 parts, each part individually MD5-verified, full reassembly MD5 verified; 27 tar members (3 Ctrl + 6 RIF samples × 3 files).
- SenMayo GMT (794 B, 124 gene symbols in downloaded row; MSigDB page lists 125 — discrepancy recorded).
- GSE51981 (in progress): full archive fetched via www endpoint (746,424,320 B, MD5 c3bc305129ce43ca97d93acad05effa4, 148 tar members = exactly the published n=148). Split into 8× 90 MB parts. Parts 00–04 copied and individually MD5-verified on storage. Parts 05–07 pending: outbound network throttled to ~30 KB/s after ~6 GB cumulative downloads (affects ftp.ncbi, www.ncbi and ftp.ebi.ac.uk alike — EBI mirror E-GEOD-51981 confirmed to exist as fallback, individual CEL layout). Resume when throughput recovers.
- Pending: GSE179640, GSE214411, GSE7307 subset (41 CELs), GSE213216 (processed objects via Cedars Box — blocked from sandbox; needs lead researcher browser download), PRJNA713993 (reference-only, no download by design).
- Incident note: storage mount crashed 3× under heavy write bursts (each crash wipes the local /tmp staging area). Protocol adopted: sequential single-part copies with per-part MD5 verification; range-request chunked downloads; full reassembly MD5 check per dataset.

### S1 — manifest + residency amendment (2026-09-11)

- `tables/dataset_manifest.tsv` created (v1): 14 rows; MD5s recorded for all completed downloads (GSE11691 0c0a95ce..., GSE7305 212fd3ae..., GSE25628 4b161cf1..., GSE6364 95969032..., GSE120103 9482b104..., GSE183837 ced92ae3..., SenMayo 322420eb...).
- Tar member counts verified against publications: GSE11691 = 18 CEL (9 paired women), GSE7305 = 20 CEL, GSE120103 = 36 Agilent txt (36 samples), GSE51981 = 148 members (= published n=148), GSE183837 = 27 members, GSE25628 = 22 CEL (vs n=16 in citing papers — flagged for S2 metadata reconciliation).
- `scripts/S1_GSE7307_subset.txt`: frozen 41-GSM subset list (18 disease + 23 normal) extracted from GSE7307 series matrix.
- `scripts/S1_download_datasets.sh`: lead-researcher-side downloader with resume + checksums (covers GSE51981, GSE179640, GSE214411, GSE7307 subset, GSE213216 options).
- `logs/charter_amendment_02_data_residency.md` PROPOSED: heavy raw data lives on lead's machine; shared storage keeps ≤100 MiB files + parts + all deliverables; manifest tracks MD5 regardless of residency. Awaiting confirmation.

### Amendment 02 adopted (2026-09-11)

- Lead researcher delegated the decision; assistant recommendation adopted per standing instruction ("اتبع توصياتك").
- Consequence: remaining heavy acquisition (GSE51981 parts 05-07 equivalent, GSE179640, GSE214411, GSE7307 subset, GSE213216) moves to the lead researcher's machine via `scripts/S1_download_datasets.sh`. Sandbox throttle persisted across probes (30-60 KB/s on ftp.ncbi, www.ncbi, ftp.ebi alike).
- On return of the local MD5 log, manifest will be completed and S1 QC gate evaluated: per-file size+MD5 match, sample counts vs publications, GSM list frozen.

### Windows adaptation (2026-09-11)

- Lead researcher works on Windows. Provided `scripts/S1_download_datasets.ps1` (native PowerShell + built-in curl.exe + Get-FileHash MD5). The bash version remains for reference. R/renv workflow (S0_environment_lock.R) is OS-independent and unchanged.

### S1 CLOSED (2026-09-11)

- Lead's local logs received and audited. All 13 acquisition items verified; QC gate S1 PASS.
- Manifest v2 written with corrections: GSE183837 = Lai et al. 2022 RIF cohort (not Mareckova); HECA = composite reference-source-only; Zou = PRJNA713993 with SRR13962189/SRR13962190 (no GSE exists); GSE214411 = 6 EMS + 7 controls (GEO official 13); GSE25628 = 22 per GEO official; GSE6364 = 37 per GEO official.
- GSE213216: 6/6 processed Seurat objects on lead's machine, per-file MD5s recorded.
- Superseded artifacts deleted (incomplete portal GSE51981 parts; over-appended sandbox GSE214411 copy).
- Process note: authoritative MD5s taken from uploaded log files, not from retyped chat text (two transposed values caught).
- Final report: `logs/S1_verification_report.md`. S2 remains gated pending go-ahead.

### S2 preparation package delivered (2026-09-11)
- `logs/S2_preregistered_design.md`: QC thresholds, integration design, annotation panels, and S2 gate criteria — all fixed BEFORE execution (R5).
- Scripts for local execution: `S2_00_entry_checks.R` (E1-E3), `S2_01_qc_singlecell.R` (Seurat + scDblFinder QC per dataset), `S2_02_integration_annotation.R` (Harmony + LISI + marker validation), `S2_03_bulk_qc.R` (RMA/AQM for GPL570 & GPL96, neqc for Agilent).
- `S0_environment_lock.R` updated with S2 dependencies (scDblFinder, affy, oligo, arrayQualityMetrics, lisi).
- Execution order on lead's machine: S0_environment_lock.R → S2_00 → S2_01 → S2_02 → S2_03 → return outputs for verification report.

### Next gated actions (updated 2026-09-11, post-S1 closure)

1. S0 completion (still open): lead researcher runs `scripts/S0_environment_lock.R` locally → returns `renv.lock` + `logs/sessionInfo_S0.txt`.
2. S2 (QC + integration): gated — awaiting written go-ahead. On approval, assistant delivers S2 script package (per-dataset QC with pre-registered thresholds + integration design) for local execution.
3. S2 pre-registered entry checks (from S1 report): re-count tar members at first local load (GSE51981/GSE179640/GSE214411); GSE7305↔GSE7307 overlap check; GSE25628 metadata parse (22 samples).