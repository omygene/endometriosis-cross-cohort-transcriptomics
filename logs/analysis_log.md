2026-09-17: S2_03 GSE120103 completed via scripts/S2_03_GSE120103_fix.R (lead-authored, MD5 8cd35692c14b5b31fdabc570fd464a8b). Measurement found 5 arrays of chip design GE1-v5_95_Feb07 - all excluded, 31/36 processed with neqc. No AQM for this cohort by design (S2_03-A). 
2026-09-18: S3_common.R config filled from S3_00 census (lead-authored): n1_sample=GSM6605437_N1, group_col=condition, ems=EMS, ctrl=N, lineage_col=lineage. hormonal_col/treated_lbl left EMPTY -- atlas has no treatment column; escalated to session. MD5 = 32b19b664e3d20e1997b67d520c4be09 
2026-09-18: S3_common.R config filled from S3_00 census (lead-authored): n1_sample=GSM6605437_N1, group_col=condition, ems=EMS, ctrl=N, lineage_col=lineage. hormonal_col/treated_lbl left EMPTY -- atlas has no treatment column; escalated to session. MD5 = 32b19b664e3d20e1997b67d520c4be09 
2026-09-20: S3_common.R v2.6 config re-filled after NA-fix replacement (lead-authored): n1_sample=GSM6605437_N1, group_col=condition, lineage_col=lineage. hormonal_* left EMPTY pending final S3_03-05 (Amendment A-HCOV/A-NOTX approved by lead 2026-09-19). MD5 = 147b2745ced588d88ee510ab0fbc5466 
2026-09-20: S3_common.R v2.6 config re-filled after NA-fix replacement (lead-authored): n1_sample=GSM6605437_N1, group_col=condition, lineage_col=lineage. hormonal_* left EMPTY pending final S3_03-05 (Amendment A-HCOV/A-NOTX approved by lead 2026-09-19). MD5 = 147b2745ced588d88ee510ab0fbc5466 
2026-09-20: S3_common.R v2.7 config re-filled (lead-authored, 3rd fill -- same measured values): n1_sample=GSM6605437_N1, group_col=condition, lineage_col=lineage. Amendment A applied (hormonal arms locked). MD5 = e6b465112fbb6100f28dc82629d69014 
2026-09-20: S3_common.R v2.8 config re-filled (lead-authored, 4th fill, same measured values): n1_sample=GSM6605437_N1, group_col=condition, lineage_col=lineage. get_rna slot-to-layer fix applied. MD5 = 6fc033640abdc79ab478dd8610f27bb8 
Mon 09/21/2026  0:19:58.13 S3_common.R replaced with v3.0 (per-dataset expression fix) - config filled: n1_sample="GSM6605437_N1" group_col="condition" lineage_col="lineage" - MD5 8d9f8226fb4593ec2ee9b66d0032442 
Mon 09/21/2026  0:41:47.98 S3_common.R replaced with v3.0.1 (per-dataset fix, config pre-filled by agent from lead values) - MD5 344edd5f2ac4ffe6855ffa8e9554d6cf 
Mon 09/21/2026  0:41:54.92 S3_smoke.R replaced with v2.12.1 (smoke leg 11 fix - test gates not Seurat) - MD5 d6743e133fa06bbbbd7b13d43b021ee4 
Tue 09/22/2026  0:22:50.61 S4 package v1.0 installed (7 files verified) - metafor 5.2-1 locked in renv.lock (R5) 
Tue 09/22/2026  8:31:41.31 A3 approved and executed: S3_common v3.0.2 + S3_01 v2.0 (symbol-level regeneration, 10 files MD5-verified) - S4 registry deleted - MD5s: S3_common=21c4541b... S3_01=0e21c38c... run_S4.bat=b86f13a0... 
Tue 09/22/2026  8:59:25.30 S4_common.R replaced with v1.1.1 (audit namespace overlap fix, MD5 a7bc38ab512edcb19cc8f04528f69c12) 
Tue 09/22/2026 14:08:17.63 S4_05 + S5 package installed (A4, 5 files MD5-verified, pathway-level RE meta for pre-registered axes) 
Tue 09/22/2026 15:50:16.04 make_figures.R v1.0 run (A5, 5 figures from locked files, MD5 3997f6efd9a5711f0541faf86fd7de02) 
Tue 09/22/2026 18:19:56.40 make_figures.R v2.1 installed (fix #22: cairo_pdf replaced with pdf device, MD5 bf150e6769b2080d3c8ff7e4d00a6b34) 
Tue 09/22/2026 22:19:52.71 make_figures.R v2.2 installed (fix #24: fig4b title clipping, MD5 28492734551d6655f2b65728c1a0cae4) 
Tue 09/22/2026 22:24:21.12 S6_00_inspect.R run (S6 preflight diagnostic, MD5 8850f994845dbe77adc61c443be5a478) 
Tue 09/22/2026 23:22:28.42 S6_01_validate.R run (external validation GSE213216, pre-registered verdicts locked in code, MD5 0f563fb8505661a0d23ec4fdfe8a1f2e) 
Wed 09/23/2026  0:07:45.96 make_figures.R final run (S6 Fig6 added, fault #25 fixed, MD5 99a68ffa071615c5085dac67faf94fe0) 
Wed 09/23/2026  0:29:22.59 S6_02 exploratory subcluster test run (SERPINE1 across 18 endothelial subclusters, hypothesis-generating only, MD5 eb2d4c4cd689d436698219a534f53409) 
Wed 09/23/2026  0:45:46.01 make_figures.R v3.1 installed (median line weight halved 0.6 to 0.3, S6b slots closed in manuscript, MD5 fed8d99e460cdb40e8fd8a2dd02dfada) 
Wed 09/23/2026 17:15:24.36 S8 reverse-confirmation run (LINCS, 4 pre-locked hypotheses, S8_00=88dcdc9c... S8_01=57df53d8...) 
Fri 09/25/2026 14:11:11.11 make_figures v3.3 + make_figS1 v1.0 installed (PI review layer: Fig7 colored, HA1E flagged, S8b supplementary figure with certified re-derivation gate) 
Fri 09/25/2026 15:08:52.58 make_figS1 v1.1 (legend restored: points/fit-CI; titles shortened vs combo clipping) 
