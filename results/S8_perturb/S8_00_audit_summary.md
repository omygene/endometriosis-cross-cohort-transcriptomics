# S8_00 -- LINCS/iLINCS audit summary (machine-printed)

Date: 2026-09-23 17:15:52.825928

## Robust core re-derived from locked files
- robust core: **66** genes (14 up / 52 down) -- Tier-1 retained under all 6 LOCO + sans-G7307
- internal consistency gate vs locked value (66 = 14/52): PASS

## L1000 coverage (GSE92742 gene info)
- L1000 space: 12328 genes (978 landmarks)
- robust core in L1000: 48/66; landmark: 8/66 (EGR1, FOS, CCND1, LYPLA1, PLA2G4A, GADD45B, NNT, UGDH)
- Tier-1 in landmark space: 46 (13 up / 33 down) [criterion G2 >= 40: PASS]
- SERPINE1 in landmark space: TRUE

## iLINCS availability
- SERPINE1-KD (CGS, LIB_6): 7 signatures / 7 cell lines [criterion G1 >= 5: PASS]
- estradiol (LIB_5, L1000): 136 signatures [criterion G3 >= 3: PASS]
- pre-registered absent candidates (documented, not silently skipped):
  - tiplaxtinin: n=0
  - TM5441: n=0
  - PAI-039: n=0
  - forskolin: n=0

## VERDICT
**GO -- S8_01 reverse-confirmation may proceed (all locked criteria met)**

*Every number in this file was measured at run time from the locked inputs and
 the live iLINCS API. Nothing was copied from chat text (R1/R2).*
