# endometriosis-cross-cohort-transcriptomics
Reproducible code, frozen metadata and final results for a cross-cohort transcriptomic analysis of endometriosis-associated tissue contrasts.
Cross-cohort transcriptomic analysis of endometriosis-associated tissue contrasts
Overview
This repository contains the analysis code, frozen metadata, final result tables, figure-generation scripts and quality-control outputs for a cross-cohort transcriptomic study of endometriosis-associated tissue contrasts.
The study integrates six prespecified bulk transcriptomic cohorts using random-effects meta-analysis. A paired cohort is analysed separately for directional concordance. Three single-cell RNA-sequencing datasets are used for discovery-stage cellular localization, and an independent single-cell atlas is reserved for external validation using patient-level pseudobulk profiles. LINCS perturbational analyses and in-vitro decidualization-response analyses are reported separately.
Study design
•	Bulk discovery meta-analysis: GSE7305, GSE6364, GSE25628, GSE7307, GSE51981 and GSE120103.
•	Reserved paired concordance arm: GSE11691.
•	scRNA-seq discovery localization: GSE179640, GSE183837 and GSE214411.
•	Independent external validation atlas: GSE213216.
•	Exploratory decidualization-response datasets: GSE75423 and GSE75425.
The six bulk discovery contrasts are heterogeneous in tissue source, disease context and comparator definition. Therefore, all primary pooled analyses use a random-effects framework and report heterogeneity rather than filtering on it.
Key analytical rules
•	The primary statistical claim gate is Benjamini-Hochberg false-discovery rate (BH-FDR) < 0.05 within each prespecified testing family.
•	A Tier-1 gene satisfies BH-FDR < 0.05, absolute pooled log2 fold change >= 0.5, and directional concordance in at least four of six discovery cohorts.
•	The robust core is defined as genes retaining Tier-1 status in all six leave-one-cohort-out analyses and the prespecified analysis excluding GSE7307.
•	The paired GSE11691 cohort is never pooled with the six-cohort discovery meta-analysis.
•	Single-cell statistical analyses use sample- or patient-level pseudobulk profiles; individual cells are not treated as independent observations.
•	External validation uses GSE213216 patient-level pseudobulks and retains non-testable planned contrasts in the test matrix.
•	LINCS and decidualization analyses are interpreted according to their preregistered or explicitly exploratory status.
Repository structure
scripts/             Analysis and figure-generation R scripts
metadata/            Dataset manifest, frozen contrast registry and group mappings
results/final/       Final result tables used for manuscript claims
results/qc/          QC summaries, audits and diagnostics
results/sensitivity_full/  Full all-gene outputs for sensitivity scenarios
figures/main/        Main manuscript figures
figures/supplementary/ Supplementary manuscript figures
figures/qc_microarray/ Detailed microarray QC diagnostics
reproducibility/     Analysis logs, amendments, package checks and environment lockfile

Reproducing the environment
The analysis was performed in R 4.6.1 with Bioconductor 3.23. The exact package environment is captured in renv.lock.
install.packages("renv")
renv::restore()

Run the workflow from the repository root. Script stages follow the prefix convention S0 through S8b. The stage-specific scripts and their registered inputs/outputs provide the execution order and audit trail.
Primary outputs
The principal final outputs are provided as machine-readable supplementary data files:
•	Supplementary Data 1: all-gene random-effects meta-analysis.
•	Supplementary Data 2: Tier-1 genes and robust-core membership.
•	Supplementary Data 3: sensitivity-analysis summary.
•	Supplementary Data 4: pathway-level random-effects meta-analysis.
•	Supplementary Data 5: paired GSE11691 concordance.
•	Supplementary Data 6: single-cell pseudobulk localization results.
•	Supplementary Data 7: all planned external-validation tests.
•	Supplementary Data 8: exploratory endothelial SERPINE1 subcluster results.
•	Supplementary Data 9: LINCS coverage and perturbation results.
•	Supplementary Data 10: exploratory decidualization-response results.
Public data availability
Raw and processed public data are not redistributed in this repository. They remain available from the Gene Expression Omnibus (GEO), Sequence Read Archive (SRA), MSigDB and LINCS resources under the accessions and resource identifiers documented in metadata/dataset_manifest.tsv and Supplementary Table S1.
Quality control
Detailed bulk microarray QC diagnostics, including principal-component analyses, array-level outlier diagnostics, MA plots, RLE/boxplot diagnostics, array-distance heatmaps and intensity-variance plots, are retained under figures/qc_microarray/. Single-cell QC diagnostics are retained under results/qc/ and represented in Supplementary Figures S1-S5.
Reproducibility record
The repository includes frozen contrast definitions, input registries, audit logs, analysis amendments, package checks and an renv.lock environment lockfile. These files are intended to provide an auditable record of analytical inputs, decisions and software dependencies.
License and citation
Add a license before making the repository public. Recommended: MIT License for code. Please cite the associated manuscript and the Zenodo DOI for the archived release when available.
Contact
For questions regarding the analysis, contact the corresponding author listed in the manuscript.
