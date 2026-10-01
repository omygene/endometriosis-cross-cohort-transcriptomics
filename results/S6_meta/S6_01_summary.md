# S6_01 -- GSE213216 external validation (locked design v1.0)

Run: 2026-09-22 23:28:26.809338 | input md5 verified vs registry: 2218a6817475f69dfdeefa7451ada0ed
Unit = patient. Pseudobulks (>=20 cells): 271 | patients: 21
Tests: 195 | tested 109 | NOT_TESTABLE 86 | one BH family across all primary tests (FDR<0.05 sole gate)

## Patients per tissue (pseudobulk level)

- Disease-free peritoneum: 2 patients
- Endometrioma: 7 patients
- Eutopic Endometrium: 10 patients
- Peritoneal lesion: 12 patients
- Unaffected ovary: 4 patients

```
   Disease-free peritoneum Endometrioma Eutopic Endometrium Peritoneal lesion
1                        0            1                   0                 0
10                       1            0                   1                 1
11                       0            1                   1                 1
12                       0            0                   1                 0
13                       0            0                   1                 0
14                       0            0                   1                 1
15                       0            0                   0                 0
16                       0            0                   1                 1
17                       0            0                   1                 0
18                       0            0                   1                 0
19                       0            1                   0                 1
2                        0            1                   0                 0
20                       0            0                   0                 1
21                       0            0                   0                 1
3                        0            0                   1                 0
4                        0            0                   0                 1
5                        0            0                   0                 1
6                        0            0                   1                 0
7                        0            1                   0                 1
8                        1            1                   0                 1
9                        0            1                   0                 1
   Unaffected ovary
1                 0
10                0
11                0
12                0
13                0
14                1
15                1
16                0
17                0
18                1
19                0
2                 0
20                0
21                0
3                 0
4                 0
5                 0
6                 1
7                 0
8                 0
9                 0
```

## Pre-registered verdicts (locked rules, machine-printed)

- **P1 (SERPINE1 up, Endothelial, C1): NOT_CONFIRMED**
- **P2 (proliferation down, Mesenchymal, C1): NOT_CONFIRMED**
- **E3 (Tier-1 signature up, Mesenchymal, C1): NOT_CONFIRMED**
- P3 (directional only): senescence_dormancy delta=-0.654 (raw p=0.422) | SenMayo delta=-0.0229 (raw p=0.945) -- directional-only endpoint by design
- P4 (no direction claimed): cytotoxicity delta=-0.282 FDR=0.418 | antigen_presentation delta=-0.332 FDR=0.613 -- reported as-is, no direction claimed by design

## FDR < 0.05 hits (all, unfiltered)

- E3 | tier1_signature | Smooth muscle cells | C1 | delta=0.613 | FDR=0.0297 | n=17v10

Negative and NOT_TESTABLE results are reported in S6_01_tests.csv with the
same weight as positives (R1/R6). C1 sensitivity with patients shared across
arms removed: column p_indep_patient.

S6_01 DONE -- send S6_01_summary.md + S6_01_tests.csv back to the session.
