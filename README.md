# Integrative Single-Cell, Spatial and Multi-Cohort Transcriptomic Analysis of Plasma Membrane Tension Regulators in Clear Cell Renal Cell Carcinoma (`KIRC_MTRG`)

Complete reproducible R analysis pipeline accompanying the manuscript:

> **Integrative Single-Cell, Spatial and Multi-Cohort Transcriptomic Analysis of Plasma Membrane Tension Regulators Reveals Tumour–Stroma Mechanical Crosstalk and Predicts Therapy Response in Clear Cell Renal Cell Carcinoma**  
> Xiaoyang Guo, Mingjie Wu, Zhiyue Sun, Ke Zeng, Zhenghang Liao, Bo Wen*  
> *Department of Urology, Shenzhen Hospital of Integrated Traditional Chinese and Western Medicine, Shenzhen 518104, Guangdong, China*

- **Repository Identifier:** `https://github.com/xiaoyangguo2024/kirc-mtrg-mechanotransduction` (Zenodo metadata included in `.zenodo.json`)
- **License:** MIT License (`LICENSE`)

## Pipeline Architecture (`R/00_config.R` – `R/11_figures.R`)

All global parameters and random seeds (`seed = 20261108`) are centralized in `R/00_config.R`:
- `R/00_config.R`: Global parameters, directory structure, and statistical helper functions.
- `R/01_geneset.R`: Curation of the 97-gene core MTRG catalog (`MT_core`), 125-gene extended catalog (`MT_all`), 4 biophysical modules (`MT_cortex`, `MT_BAR`, `MT_channel`, `MT_ECM`), bidirectional Membrane Tension Index (`MTI`), and 10 published ccRCC prognostic signatures (`published_KIRC_signatures_sources.csv`).
- `R/02_download_bulk.R`: Retrieval and uniform standardization of 8 bulk transcriptomic cohorts (`TCGA-KIRC`, `GSE29609`, `GSE22541`, `GSE53757`, `GSE73731`, `GSE67501`, `IMvigor210`, `GSE176307`).
- `R/03_scRNA.R` & `R/03b_sc_cellcycle.R`: Single-cell RNA-seq processing (`GSE159115`, 14,208 cells across 7 ccRCC tumours and 6 normal kidneys; `200 < nFeature_RNA < 7000`, `percent.mt < 25%`), Harmony integration, `copyKAT` negative-control validation, sample-paired `edgeR` pseudobulk differential expression, G1-only Hallmark validation, and `CellChat` intercellular communication.
- `R/04_bulk_score.R`: Bulk `ssGSEA` pathway scoring, proliferation overlap audit, somatic mutation analysis (`maftools`), and module-level Cox regressions.
- `R/05_candidates.R`: Identification of 37 candidate genes and 16 prognostic genes by intersecting `MT_all`, `WGCNA` key-module hubs (`|MM| > 0.50`, `|GS| > 0.25`), `TCGA-KIRC` `DESeq2` DEGs, and single-cell pseudobulk DEGs (`P < 0.05`, `|log2FC| > 0.30`).
- `R/06_ml_signature.R`: Benchmark of 117 machine-learning survival configurations (`Mime1`) with fixed risk direction across cohorts (`StepCox[forward] + survival-SVM`, 16-gene `MTRS`) and head-to-head comparison against 10 published ccRCC signatures.
- `R/06b_specificity.R`: 2,000-iteration random gene-set permutation benchmark (`Null A`: 1,000 random expressed 16-gene sets; `Null B`: 1,000 random prognostic 16-gene sets) and proliferation-adjusted Cox models.
- `R/07_clinical_nomogram.R`: Multivariable clinical Cox regression, nomogram construction (`rms`), bootstrap calibration, and decision curve analysis (`dcurves`).
- `R/08_tme_therapy.R`: 7-algorithm immune/stromal deconvolution (`IOBR`), checkpoint expression, `TIDE`, immunotherapy clinical response, and `GDSC2` drug sensitivity (`oncoPredict`).
- `R/09_hub_gene_spatial.R`: Exact linear SHAP decomposition (`shapviz`, `NTN4`, `RNF152`, `GNG7`) and 10x Genomics Visium spatial transcriptomics (`GSE210041`, 7,079 spots across 2 ccRCC slices) with partial Spearman correlations controlling for proliferation and UMI depth.
- `R/10_tables.R`: Generation of `Table 1` and `Supplementary_Tables_S1-S10.xlsx` (`Tables_all.xlsx`).
- `R/11_figures.R`: Generation of publication-ready 600 dpi main (`Figure1`–`Figure9`) and supplementary (`FigureS1`–`FigureS3`) figures.

## Reproduction

Run individual steps from the project root directory via `./run_step.sh <STEP>` (for example, `./run_step.sh 01_geneset` through `./run_step.sh 11_figures`). Exact package versions are documented in `Table S10` (`results/10_tables/Table_S10_package_versions.csv`).
