# Analysis Workflow

This document summarizes the manuscript-associated analysis order. The repository preserves the final locked analysis logic and reference outputs; it is not a redesign of the study.

## 1. GSE18965

The manuscript-method main analysis uses the GEO deposited processed log2 Series Matrix. The workflow maps GPL96 probes to gene symbols, collapses multiple probes to the highest mean expression probe per gene, performs limma analysis with log2FC defined as atopic asthma minus non-atopic healthy controls, and applies the fixed DEG threshold FDR < 0.05 and absolute log2FC >= 0.5.

Raw CEL analyses using RMA and GCRMA are retained as preprocessing sensitivity analyses and do not replace the Series Matrix main analysis.

## 2. Gene-Set Scoring

Scores are computed for fixed T2, IFN, and repair-ECM gene sets. The repair-ECM list is fixed at 32 historical genes. External ECM gene sets are sensitivity definitions and are not used to replace the fixed historical manuscript score.

## 3. RNA-seq Cohorts

GSE152004 is used for T2/IFN/repair-ECM scoring, PCA, and k=3 epithelial-state clustering. GSE118761 nasal epithelium is the primary cross-cohort projection tissue, with tracheal epithelium treated as tissue sensitivity.

## 4. PCA and Clustering

PCA and k=3 clustering are run on locked score inputs. The k=3 state sizes in the locked manuscript are E1=338, E2=299, and E3=58. PCA factor labels should be interpreted from the current loadings/correlation tables rather than assumed to match unrecoverable historical numbering.

## 5. Cross-Cohort Projection

GSE152004 centroids are projected to GSE118761 epithelial samples without combining raw matrices across cohorts. The locked manuscript projection counts are nasal E1/E2/E3=30/18/7 and tracheal E1/E2/E3=15/29/5.

## 6. Hospital Clinical Models

Hospital clinical analyses require protected de-identified patient-level data and are not publicly reproducible from this repository alone. The public repository provides the expected input schema and a synthetic demonstration file. The synthetic file must not be used to reproduce manuscript clinical estimates.

## 7. Tables and Figures

Reference outputs and figure previews are stored under `outputs/`. These are manuscript-associated locked outputs or audit reference outputs, depending on the file. See `figure_table_code_map.csv` and `reproducibility_notes.md`.

