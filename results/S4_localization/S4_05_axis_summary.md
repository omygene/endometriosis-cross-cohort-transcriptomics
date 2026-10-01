# S4_05 -- Axis localization summary (measurement only)

Locked criteria: signal flag = BH-FDR < 0.05 only (design section 13).
No exclusion, no weighting, no pooling. B-G11691 is its own paired arm.
'++/--' = FDR < 0.05; '+/-' = detected but not FDR-significant;
'.' = gene absent from that output (measured, never dropped silently).

## cytotoxicity (CD8A, CD8B, GZMB, PRF1, NKG7, GNLY)

### bulk: gene x primary contrast (direction of TEST vs REF)
| gene | B-G7305 | B-G6364 | B-G25628 | B-G7307 | B-G51981 | B-G120103 | B-G11691 |
|------|------|------|------|------|------|------|------|
| CD8A | ++ | - | + | ++ | ++ | - | + |
| CD8B | - | + | - | + | ++ | - | - |
| GZMB | - | + | + | - | ++ | + | - |
| PRF1 | - | + | + | - | ++ | - | - |
| NKG7 | -- | + | + | - | ++ | - | - |
| GNLY | -- | + | + | -- | ++ | + | - |

### scRNA: significant cells (FDR < 0.05) per lineage x arm
- GSE214411 | Mast | B: 1 of 6 axis genes

## antigen_presentation (HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2)

### bulk: gene x primary contrast (direction of TEST vs REF)
| gene | B-G7305 | B-G6364 | B-G25628 | B-G7307 | B-G51981 | B-G120103 | B-G11691 |
|------|------|------|------|------|------|------|------|
| HLA-A | - | + | ++ | - | ++ | - | + |
| HLA-B | + | + | ++ | - | ++ | - | + |
| HLA-C | + | + | + | - | ++ | - | + |
| B2M | + | + | - | + | -- | - | + |
| TAP1 | -- | - | - | -- | ++ | + | - |
| TAP2 | - | - | -- | -- | + | + | - |

### scRNA: significant cells (FDR < 0.05) per lineage x arm
- (no significant cells in any leg)

## senescence_dormancy (CDKN1A, CDKN2A, GADD45A, SERPINE1)

### bulk: gene x primary contrast (direction of TEST vs REF)
| gene | B-G7305 | B-G6364 | B-G25628 | B-G7307 | B-G51981 | B-G120103 | B-G11691 |
|------|------|------|------|------|------|------|------|
| CDKN1A | ++ | - | ++ | ++ | ++ | - | - |
| CDKN2A | -- | + | - | - | ++ | -- | - |
| GADD45A | ++ | - | + | ++ | + | - | - |
| SERPINE1 | + | + | + | ++ | ++ | + | - |

### scRNA: significant cells (FDR < 0.05) per lineage x arm
- GSE214411 | Endothelial | A: 1 of 4 axis genes
- GSE214411 | Endothelial | B: 1 of 4 axis genes
- GSE214411 | Mast | A: 1 of 4 axis genes
- GSE214411 | Mast | B: 1 of 4 axis genes
- GSE214411 | Perivascular | A: 1 of 4 axis genes
- GSE214411 | Perivascular | B: 1 of 4 axis genes

## proliferation (MKI67, TOP2A, PCNA)

### bulk: gene x primary contrast (direction of TEST vs REF)
| gene | B-G7305 | B-G6364 | B-G25628 | B-G7307 | B-G51981 | B-G120103 | B-G11691 |
|------|------|------|------|------|------|------|------|
| MKI67 | -- | + | - | -- | -- | - | - |
| TOP2A | -- | + | - | -- | -- | - | - |
| PCNA | -- | + | -- | -- | -- | - | - |

### scRNA: significant cells (FDR < 0.05) per lineage x arm
- GSE214411 | Endothelial | B: 2 of 3 axis genes
- GSE214411 | Mast | B: 1 of 3 axis genes

## stromal_ecm (VIM, COL1A1, COL3A1, FN1, ACTA2)

### bulk: gene x primary contrast (direction of TEST vs REF)
| gene | B-G7305 | B-G6364 | B-G25628 | B-G7307 | B-G51981 | B-G120103 | B-G11691 |
|------|------|------|------|------|------|------|------|
| VIM | + | + | + | ++ | -- | + | + |
| COL1A1 | - | + | ++ | + | ++ | - | + |
| COL3A1 | -- | + | + | -- | -- | - | + |
| FN1 | ++ | + | + | ++ | -- | - | + |
| ACTA2 | ++ | + | ++ | + | ++ | + | + |

### scRNA: significant cells (FDR < 0.05) per lineage x arm
- GSE214411 | Mast | A: 5 of 5 axis genes
- GSE214411 | Mast | B: 5 of 5 axis genes
- GSE214411 | Myeloid | B: 1 of 5 axis genes
