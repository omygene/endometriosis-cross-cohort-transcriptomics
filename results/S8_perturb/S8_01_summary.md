# S8_01 -- Reverse-confirmation summary (machine-printed)

Date: 2026-09-23 17:17:49.646169

## Pre-registered hypothesis family (one BH family, sole gate FDR < 0.05)

| H | test | n | p | BH-FDR | verdict |
|---|------|---|---|--------|---------|
| P1a | SERPINE1-KD r_tier1 < 0 (7 cell lines) | 7 | 0.234 | 0.313 | NOT CONFIRMED |
| P1b | SERPINE1-KD prolif(TOP2A,PCNA) > 0 | 7 | 1 | 1 | NOT CONFIRMED |
| P2a | estradiol per-line r_tier1 < 0 (cell-line medians) | 10 | 0.216 | 0.313 | NOT CONFIRMED |
| P2b | estradiol prolif > 0 | 10 | 0.0322 | 0.129 | NOT CONFIRMED |

## Median effects
- SERPINE1-KD: median r_tier1 = -0.057 | median prolif = -0.238
- estradiol:   median r_tier1 = -0.049 | median prolif = 0.291

## Story verdicts (machine-printed)
- **P1 NOT CONFIRMED: SERPINE1 knockdown does not significantly reverse the lesion signature -- reported as such**
- **P2 NOT CONFIRMED: estradiol reversal not significant -- reported as such**

*Negative r_tier1 = the perturbation moves landmark-measurable Tier-1 genes in the
 direction OPPOSITE to the lesion (reversal). Proliferation readout uses the two
 proliferation-axis genes measurable in L1000 landmark space (TOP2A, PCNA; MKI67
 is absent from the 978). The 66-gene robust core is represented by its 8 landmark
 genes in this space (coverage measured in S8_00). Nothing here is a claim beyond
 the locked FDR gate.*
