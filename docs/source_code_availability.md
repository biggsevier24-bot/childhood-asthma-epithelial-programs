# Source Code Availability

This repository was corrected to avoid placeholder analysis scripts. Recovered real R source code from the project workspace was copied into the public repository and sanitized for public paths.

## Recovered R Source Code

| public file | recovered source | status |
|---|---|---|
| `R/01_GSE18965_preprocessing.R` | `01_source_inputs/GSE18965_Codex_full_package/GSE18965_series_matrix_reproduction.R` | recovered real source |
| `R/02_GSE18965_DEG_and_ECM_scores.R` | `01_source_inputs/GSE18965_Codex_full_package/GSE18965_GO_reaudit.R` | recovered real source |
| `R/03_RNAseq_preprocessing.R` | `04_scripts/External_ECM_multicohort_analysis.R` | recovered real source, public paths sanitized |
| `R/04_gene_set_scoring.R` | `04_scripts/External_ECM_multicohort_analysis.R` | recovered real source, public paths sanitized |
| `R/05_sensitivity_analyses.R` | `04_scripts/ECM_third_axis_corrected_rerun.R` | recovered real source |
| `R/06_tables_and_figures.R` | `04_scripts/final_locked_analysis_20260714.R` | recovered real source, output/private-data paths sanitized |

## Python Analysis Source

Independent original Python analysis scripts for PCA, clustering, projection, statistical modeling, hospital modeling, or output generation were not found in the project workspace.

The manuscript Methods reported use of Python packages including scikit-learn and statsmodels for some analyses. However, independent original Python source files implementing those analyses were not recovered from the available project workspace.

The files in `python/` therefore explicitly raise `source script unavailable` and are not represented as executable analysis implementations. They are retained only to document the unavailable Python source paths and to avoid silently presenting placeholder scripts as real analysis code.

No placeholder or reconstructed Python code is represented as the original implementation.

## Duplicate R Mapping Note

`R/03_RNAseq_preprocessing.R` and `R/04_gene_set_scoring.R` intentionally preserve the same recovered comprehensive source script for workflow mapping. The original recovered implementation combined RNA-seq preprocessing, metadata harmonization, gene-set scoring, and external ECM multicohort analyses in one R script.

## Safety Modifications Applied

- Absolute local paths were replaced with relative public repository paths or environment-variable inputs.
- Rerun outputs were redirected to `outputs/user_rerun/`.
- Protected hospital input is controlled through `HOSPITAL_PRIVATE_DATA`.
- Manuscript-locked outputs are not overwritten by the public rerun launcher.
- No model formulas, thresholds, random seeds, phenotype definitions, or locked manuscript results were intentionally changed.
