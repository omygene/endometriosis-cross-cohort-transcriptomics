# S5 sensitivity summary (pathway level)

Base significant set (BH-FDR < 0.05 across the 6 pathway tests): 1 of 6 -- proliferation

| scenario | n sig | overlap w/ base | sign match | max |dM| |
|---|---|---|---|---|
| sans-G7307 | 0 | 0/1 | 0.833 | 0.082 |
| LOCO-B-G7305 | 0 | 0/1 | 0.833 | 0.200 |
| LOCO-B-G6364 | 1 | 1/1 | 0.833 | 0.278 |
| LOCO-B-G25628 | 0 | 0/1 | 0.833 | 0.115 |
| LOCO-B-G7307 | 0 | 0/1 | 0.833 | 0.082 |
| LOCO-B-G51981 | 1 | 0/1 | 0.833 | 0.164 |
| LOCO-B-G120103 | 3 | 1/1 | 0.833 | 0.141 |
| swap-B-G51981-S1 | 1 | 1/1 | 1.000 | 0.000 |
| swap-B-G51981-S2 | 1 | 1/1 | 1.000 | 0.000 |
| FE-vs-RE | 5 | 1/1 | 0.833 | 0.180 |

Per-pathway FE FDRs and scenario FDRs: S5_sensitivity_scenarios.csv

Interpretation locks (design section 14): FDR < 0.05 is the only
significance claim; |M| >= 0.5 is interpretive; tau2/I2 reported, never
filtered; a scenario flipping a pathway across the FDR boundary is
reported as a boundary flip under heterogeneity, not hidden.
