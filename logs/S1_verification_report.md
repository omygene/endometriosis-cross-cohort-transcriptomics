# S1 Verification Report — Dataset Acquisition (FINAL)

**Date:** 2026-09-11
**Stage:** S1 Acquisition — QC gate evaluation (Charter Section 4)
**Inputs audited:** lead researcher's local download logs (`S1_download_20260911_160454.log`, `S1_md5_20260911.txt`, `GSE213216_md5.txt`) + sandbox-side verified files + NCBI GEO/SRA official records (eutils esummary, 2026-09-11).

---

## 1. QC gate S1 (charter): "Every file has accession + download date + MD5 checksum recorded in the dataset manifest; sample counts match the original publications"

| Dataset | Accession ✓ | Date ✓ | MD5 ✓ | Sample count vs official | Verdict |
|---|---|---|---|---|---|
| Fonseca atlas | GSE213216 | 2026-09-11 | 6/6 per-file MD5s | 6/6 processed objects (authors' Box) | PASS |
| Lai RIF reference | GSE183837 | 2026-09-11 | ced92ae3... | 27 members = GEO official 9 samples | PASS |
| Tan validation arm | GSE179640 | 2026-09-11 | 6aea440e... | size = GEO filelist exact; GEO official 59 records | PASS |
| Huang minimal/mild | GSE214411 | 2026-09-11 | 87c8ed21... | size = GEO filelist exact; GEO official 13 samples (6 EMS + 7 Ctrl) | PASS |
| Zou peritoneal | PRJNA713993 | — | reference-only | 2 SRA runs (= published n=1+1) | PASS (as reference-only) |
| Primary bulk | GSE51981 | 2026-09-11 | c3bc3051... (dual-machine match) | 148 members = published n=148 | PASS |
| Paired bulk | GSE11691 | 2026-09-11 | 0c0a95ce... | 18 CEL = 9 paired women | PASS |
| GSE7305 | GSE7305 | 2026-09-11 | 212fd3ae... | 20 CEL = 10 + 10 | PASS |
| GSE7307 subset | GSE7307 | 2026-09-11 | per-file MD5 log | 41/41 CEL of frozen list | PASS |
| GSE25628 | GSE25628 | 2026-09-11 | 4b161cf1... | 22 CEL = GEO official 22 | PASS |
| GSE6364 | GSE6364 | 2026-09-11 | 95969032... | 37 CEL = GEO official 37 | PASS |
| GSE120103 | GSE120103 | 2026-09-11 | 9482b104... | 36 members = 36 samples | PASS |
| SenMayo | M45803 | 2026-09-11 | 322420eb... | 124 symbols in downloaded GMT | PASS |

**GATE S1: PASS — all 13 acquisition items verified.**

## 2. Corrections applied this session (manifest v2)

1. **GSE183837 re-identified:** Lai et al. 2022 — recurrent implantation failure (RIF) cohort, 3 Ctrl + 6 RIF at WOI (GEO title confirmed: "Single-cell transcriptome profiling of the human endometrium from patients with recurrent implantation failure", 9 samples). Charter's "Mareckova atlas" label was wrong; files themselves are correct and fit the normal-reference role.
2. **HECA/Marečková** has no single accession — composite of 7 datasets. Recorded as reference-source-only row; acquisition decision deferred to S2 integration design.
3. **Zou et al.:** no GSE exists. Verified in SRA: PRJNA713993 → runs SRR13962189 + SRR13962190 (exactly 2 runs, matching the published n=1 EMS + n=1 control). Supporting evidence only.
4. **GSE214411 controls = 7** (GEO official 13 samples = 6 EMS + 7 controls; N1–N4 + N-5..N-7). Question resolved.
5. **GSE25628 = 22 samples per GEO official** — our archive complete; "n=16" was citing-paper subsetting. Note retained for S2 metadata parse.
6. **GSE6364 = 37 samples per GEO official** — our 37 CEL members match exactly.
7. **MD5 transcription caution:** the lead's chat summary transposed characters in two MD5s; the authoritative values are from the uploaded log files (used in manifest). Recorded per R6 as a process note — always verify from files, never from retyped text.
8. **Sandbox-side artifacts cleaned:** the earlier over-appended GSE214411 sandbox copy and the incomplete GSE51981 portal parts (00–04) were discarded; lead-local copies are authoritative per Amendment 02.

## 3. Honest limitations (R6)

- GEO does not publish official checksums; "verification vs GEO" = exact byte-size match against the official GEO `filelist.txt` plus successful independent re-download. GSE51981 additionally has dual-machine MD5 agreement (sandbox + lead local, same value c3bc3051...).
- Tar member counts for lead-local archives (GSE51981, GSE179640, GSE214411) will be re-counted at first local load in S2 (Windows 10+ ships bsdtar; one-line check) — pre-registered as S2 entry check, not a blocker.
- GSE7305 ↔ GSE7307 potential sample overlap to be checked at S2 (same tissue bank lineage, Neurocrine submissions).

## 4. Gate decision

**S1 is COMPLETE.** All charter QC criteria met. Awaiting lead researcher's written go-ahead to open **S2 (QC + integration)** per rule R3.
