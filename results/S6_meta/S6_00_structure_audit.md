# S6_00 -- GSE213216 structural audit (pure measurement)

Run: 2026-09-22 22:25:41
SeuratObject available: TRUE

## auxiliary.seurat.shared.rds
- path: `data/raw/GSE213216_processed/auxiliary.seurat.shared.rds`
- class: Seurat | genes: 30354 | cells: 373851
- assays: SCT, RNA
- gene ID style: symbol-like
- axis genes present: 24/24 | Tier-1 core present: 463/485 | SenMayo present: 124/124
- metadata columns (24): nCount_RNA, nFeature_RNA, orig.ident, Patient.No., Stage, Fresh.Frozen, Major.Class, Index, percent.mito, scublet_doublet_score, scublet_predicted_doublet, DF.classifications, nCount_SCT, nFeature_SCT, SCT_snn_res.0.5, seurat_clusters, doublet_detection_score, doublet_detection_doublet, SCT_snn_res.3, S.Score, G2M.Score, Phase, active.cluster, Menstrual.Cycle

  - `Major.Class` (5 unique): Endometriosis=118625; Eutopic Endometrium=118144; Endometrioma=59709; Unaffected ovary=51927; No endometriosis detected=25446
  - `DF.classifications` (2 unique): Singlet=345791; Doublet=28060
  - `seurat_clusters` (114 unique): 0=16732; 1=15930; 2=14841; 3=13402; 4=11772; 5=9734; 6=9332; 7=8291; 8=7716; 9=7536; 10=7307; 11=7256; 12=7166; 13=7127; 14=7114; 15=7053; 16=7046; 17=6738; 18=6478; 19=6355; 20=6229; 21=5754; 22=5693; 23=5678; 24=5593; 25=5538; 26=4863; 27=4639; 28=4427; 29=4304; 30=4299; 31=4233; 32=3902; 33=3863; 34=3824; 35=3557; 36=3483; 37=3372; 38=3357; 39=3345
  - `active.cluster` (9 unique): Mesenchymal cells=149051; T/NK cells=101217; Epithelial cells=38456; Myeloid cells=27436; Endothelial cells=23226; Smooth muscle cells=18314; B/Plasma cells=8278; Erythrocytes=6186; Mast cells=1687

## endothelial_clusters.shared.rds
- path: `data/raw/GSE213216_processed/endothelial_clusters.shared.rds`
- class: Seurat | genes: 30354 | cells: 23226
- assays: SCT, RNA
- gene ID style: symbol-like
- axis genes present: 24/24 | Tier-1 core present: 463/485 | SenMayo present: 124/124
- metadata columns (24): nCount_RNA, nFeature_RNA, orig.ident, Patient.No., Stage, Fresh.Frozen, Major.Class, Index, percent.mito, scublet_doublet_score, scublet_predicted_doublet, DF.classifications, nCount_SCT, nFeature_SCT, SCT_snn_res.0.5, seurat_clusters, doublet_detection_score, doublet_detection_doublet, SCT_snn_res.3, S.Score, G2M.Score, Phase, active.cluster, Menstrual.Cycle

  - `Major.Class` (5 unique): Eutopic Endometrium=13816; Endometriosis=4733; Endometrioma=2430; Unaffected ovary=1818; No endometriosis detected=429
  - `DF.classifications` (2 unique): Singlet=21593; Doublet=1633
  - `seurat_clusters` (18 unique): 0=2845; 1=2690; 2=2365; 3=2352; 4=2111; 5=1837; 6=1680; 7=1679; 8=1163; 9=1054; 10=861; 11=722; 12=721; 13=322; 14=318; 15=285; 16=202; 17=19
  - `active.cluster` (1 unique): Endothelial cells=23226; Mesenchymal cells=0; Epithelial cells=0; Smooth muscle cells=0; Erythrocytes=0; Mast cells=0; Myeloid cells=0; T/NK cells=0; B/Plasma cells=0

## EnEpi_cells.shared.rds
- path: `data/raw/GSE213216_processed/EnEpi_cells.shared.rds`
- class: Seurat | genes: 30354 | cells: 4885
- assays: RNA, SCT
- gene ID style: symbol-like
- axis genes present: 24/24 | Tier-1 core present: 463/485 | SenMayo present: 124/124
- metadata columns (24): orig.ident, nCount_RNA, nFeature_RNA, Patient.No., Stage, Fresh.Frozen, Major.Class, Index, percent.mito, scublet_doublet_score, scublet_predicted_doublet, DF.classifications, nCount_SCT, nFeature_SCT, SCT_snn_res.0.5, seurat_clusters, doublet_detection_score, doublet_detection_doublet, SCT_snn_res.3, S.Score, G2M.Score, Phase, active.cluster, Menstrual.Cycle

  - `Major.Class` (3 unique): Eutopic Endometrium=3977; Endometriosis=769; Endometrioma=139
  - `DF.classifications` (2 unique): Singlet=4511; Doublet=374
  - `seurat_clusters` (5 unique): 7=1496; 13=993; 11=926; 14=891; 19=579
  - `active.cluster` (5 unique): MUC5B+=1496; Glandular secretory=993; SOX9+ LGR5+=926; IHH+ SPDEF+=891; Ciliated=579

## EnS_cells.shared.rds
- path: `data/raw/GSE213216_processed/EnS_cells.shared.rds`
- class: Seurat | genes: 30354 | cells: 4290
- assays: RNA, SCT
- gene ID style: symbol-like
- axis genes present: 24/24 | Tier-1 core present: 463/485 | SenMayo present: 124/124
- metadata columns (28): orig.ident, nCount_RNA, nFeature_RNA, Patient.No., Stage, Fresh.Frozen, Major.Class, Index, percent.mito, scublet_doublet_score, scublet_predicted_doublet, DF.classifications, nCount_SCT, nFeature_SCT, SCT_snn_res.0.5, seurat_clusters, doublet_detection_score, doublet_detection_doublet, SCT_snn_res.3, S.Score, G2M.Score, Phase, active.cluster, id.cells, select, selected.cells, Mesen.Class, Menstrual.Cycle

  - `Major.Class` (4 unique): Eutopic Endometrium=3158; Endometriosis=818; No endometriosis detected=238; Endometrioma=76
  - `DF.classifications` (2 unique): Singlet=3963; Doublet=327
  - `seurat_clusters` (2 unique): 18=2303; 19=1987
  - `active.cluster` (2 unique): EnS_Proliferative=2303; EnS_Secretory=1987
  - `id.cells` (4290 unique): AAACCCAAGACTCTTG-H12=1; AAACCCAAGGCTGAAC-FT_SA24031=1; AAACCCAAGTATGTAG-H12=1; AAACCCAAGTCGCGAA-FT_SA24031=1; AAACCCACACGAGAAC-FT_SA24031=1; AAACCCATCCATCCGT-S19286_CF=1; AAACCTGAGGGTATCG-FT_SA24024=1; AAACCTGGTAAGTGTA-FT_SA24024=1; AAACCTGTCAACGAAA-FT_SA24024=1; AAACCTGTCACAGGCC-FT_SA24024=1; AAACCTGTCCTCTAGC-FT_SA24024=1; AAACGAACACTGGATT-FT_SA24031=1; AAACGAAGTAACGATA-S19259_BC=1; AAACGAAGTGGCAACA-S19259_BC=1; AAACGAAGTTCCGCAG-H12=1; AAACGAATCAAGATAG-H12=1; AAACGAATCACCCTTG-BEME180=1; AAACGAATCTCCGTGT-FT_SA24031=1; AAACGCTAGCTGACAG-FT_SA24031=1; AAACGCTCAGCTGAAG-S19259_BC=1; AAACGCTGTGAAAGTT-BEME180=1; AAACGCTGTGATGGCA-H12=1; AAACGGGAGTTACGGG-FT_SA24024=1; AAACGGGCATATACGC-FT_SA24024=1; AAACGGGTCGTGACAT-FT_SA24024=1; AAAGAACAGGAAGTGA-BEME180=1; AAAGAACAGGATGTTA-FT_SA24031=1; AAAGAACAGGGAGTTC-BEME180=1; AAAGAACCACAACATC-H12=1; AAAGAACTCCCATAGA-G12=1; AAAGATGCAATAACGA-FT_SA24024=1; AAAGATGCATACTACG-FT_SA24024=1; AAAGCAAAGGCATTGG-FT_SA24024=1; AAAGCAAAGTCCGTAT-FT_SA24024=1; AAAGCAACACAGGTTT-FT_SA24024=1; AAAGCAACATCACGAT-FT_SA24024=1; AAAGCAATCTGAGTGT-FT_SA24024=1; AAAGGATAGCACTCTA-BEME180=1; AAAGGATAGCTACAAA-H12=1; AAAGGATAGTACCCTA-FT_SA24031=1
  - `selected.cells` (1 unique): Yes=4290
  - `Mesen.Class` (1 unique): Enometrial_type_stroma=4290

## epithelial.annotated.shared.rds
- path: `data/raw/GSE213216_processed/epithelial.annotated.shared.rds`
- class: Seurat | genes: 30354 | cells: 13771
- assays: RNA, SCT
- gene ID style: symbol-like
- axis genes present: 24/24 | Tier-1 core present: 463/485 | SenMayo present: 124/124
- metadata columns (24): orig.ident, nCount_RNA, nFeature_RNA, Patient.No., Stage, Fresh.Frozen, Major.Class, Index, percent.mito, scublet_doublet_score, scublet_predicted_doublet, DF.classifications, nCount_SCT, nFeature_SCT, SCT_snn_res.0.5, seurat_clusters, doublet_detection_score, doublet_detection_doublet, SCT_snn_res.3, S.Score, G2M.Score, Phase, active.cluster, Menstrual.Cycle

  - `Major.Class` (5 unique): Eutopic Endometrium=5163; Endometriosis=5034; No endometriosis detected=1846; Endometrioma=1106; Unaffected ovary=622
  - `DF.classifications` (2 unique): Singlet=12718; Doublet=1053
  - `seurat_clusters` (11 unique): 4=2684; 5=2049; 6=1633; 7=1509; 11=1097; 13=1000; 14=916; 16=899; 18=816; 19=732; 22=436
  - `active.cluster` (11 unique): KRT10/ACTA2 (1)=2684; KRT10/ACTA2 (2)=2049; Mesothelial (1)=1633; MUC5B+=1509; SOX9+ LGR5+=1097; Glandular secretory=1000; IHH+ SPDEF+=916; KRT10/ACTA2 (3)=899; Mesothelial (2)=816; Ciliated=732; Mesothelial (3)=436

## mesenchymal.annotated.shared.rds
- path: `data/raw/GSE213216_processed/mesenchymal.annotated.shared.rds`
- class: Seurat | genes: 30354 | cells: 82735
- assays: RNA, SCT
- gene ID style: symbol-like
- axis genes present: 24/24 | Tier-1 core present: 463/485 | SenMayo present: 124/124
- metadata columns (24): orig.ident, nCount_RNA, nFeature_RNA, Patient.No., Stage, Fresh.Frozen, Major.Class, Index, percent.mito, scublet_doublet_score, scublet_predicted_doublet, DF.classifications, nCount_SCT, nFeature_SCT, SCT_snn_res.0.5, seurat_clusters, doublet_detection_score, doublet_detection_doublet, SCT_snn_res.3, S.Score, G2M.Score, Phase, cluster.id, Menstrual.Cycle

  - `Major.Class` (5 unique): Extra-ovarian endometriosis=36365; Endometrioma=17555; No endometriosis detected=12838; Unaffected ovary=8848; Eutopic Endometrium=7129
  - `DF.classifications` (2 unique): Singlet=76360; Doublet=6375
  - `seurat_clusters` (22 unique): 0=9126; 1=7542; 2=6776; 3=5575; 4=5545; 5=5408; 6=5286; 7=5238; 8=4975; 9=4966; 10=4382; 11=3501; 12=2811; 13=2400; 14=2256; 15=1982; 16=1649; 17=1332; 18=991; 19=666; 20=243; 21=85
  - `cluster.id` (13 unique): Fibro (1)=13236; Fibro (2)=12126; C7_Activated_Fibro=10655; GAS5+=9999; C7_Fibro (1)=7713; C7_Fibro (2)=7068; Fibro (3)=6415; C7_Fibro (3)=5162; Fibro (4)=2822; EnS_Proliferative=2303; EnS_Secretory=1987; Smooth_Muscle=1735; Activated_Fibro=1514

## Next step (locked sequence)

Send this audit + S6_00_meta_values.csv back to the session. The S6_01
contrast/lineage/sample mapping is locked ONLY from these measured values
(R2/R9), then S6_01_validate.R is finalized for the PI to run (R8).
