# Childhood asthma epithelial programs: reproducible analysis repository

## Overview

This repository provides the complete, manuscript-wide, end-to-end reproducibility workflow for **“Type 2, Interferon, and Repair-ECM Programs Characterize Airway Epithelial Heterogeneity in Childhood Asthma.”** It connects public data acquisition, preprocessing, molecular-program scoring, statistical analyses, cross-cohort projection, robustness analyses, and machine-readable outputs for GSE18965, GSE118761, and GSE152004. Hospital-cohort code and aggregate outputs are included without participant-level data.

## Reproducibility workflow

The repository-level runner performs preflight checks and executes every public-data branch whose required inputs are available:

```bash
Rscript run_public_analyses.R
```

The workflow is organized in the following order:

1. Acquire deposited public data and construct standardized expression and metadata inputs for GSE18965, GSE118761, and GSE152004.
2. Construct or score the fixed 32-gene repair-ECM, T2, and IFN programs using the documented cohort-specific procedures.
3. Run correlations, regression models, effect estimates with confidence intervals, PCA, and fixed k-means epithelial-state analyses.
4. Project GSE152004 epithelial states into GSE118761 nasal and tracheal samples and evaluate phenotype associations without pooling tissues.
5. Run multicollinearity, residualization, alternative score-construction, external ECM, overlap-exclusion, and matched-random-set sensitivity analyses.
6. Write machine-readable model results, figure source data, table source data, provenance records, and automated test outputs.
7. Expose hospital-cohort model code, variable definitions, input templates, and aggregate outputs while excluding protected participant-level data.

Official GEO files can be downloaded by the public-data scripts. The optional cache variables documented in `data/README.md` support offline reuse of official files without redistributing GSE118761 or GSE152004 raw cache inputs.

## Public cohorts and analysis roles

| Cohort | Repository role | Main workflow |
|---|---|---|
| GSE18965 | Deterministic, code-based derivation of the repair-ECM module, differential-expression audit, and module scoring | `analysis/01_GSE18965_repair_ECM/` |
| GSE118761 | Airway program scoring, tissue-stratified associations, PCA, quadrants, and cross-cohort projection | `analysis/02_GSE118761_airway_programs/` |
| GSE152004 | T2 definition, program scoring, regression, PCA, epithelial states, and sensitivity analyses | `analysis/03_GSE152004_primary_validation/` |
| Hospital cohort | Protected-data clinical models and aggregate reporting | `analysis/06_hospital_cohort/` |

A deterministic, code-based workflow was applied in GSE18965 to derive the fixed 32-gene repair-ECM module. It uses prespecified biological annotations, the AA-HN expression contrast, co-expression structure, hierarchical clustering, network-based filtering, and explicit programmatic selection rules; no manual gene-name overrides or post hoc truncation are used. The derivation script does not load the reference list, which is used by a separate post-derivation verification step. Strict limma DEG/GO analysis is separate from module selection.

GSE118761 nasal and tracheal samples are processed and analyzed separately. GSE152004 uses the primary VST mean-z branch, the documented three-gene T2 definition, and fixed k=3 epithelial-state analysis. Cohort expression matrices are never merged and cross-cohort ComBat is not used.

## Core analysis modules

- **Molecular programs:** fixed 32-gene repair-ECM, fixed T2, and fixed IFN definitions are stored under `data/reference_gene_sets/`; shared scoring functions are under `R/common/`.
- **Core statistics:** cohort-specific correlations, regression models, effect estimates, confidence intervals, and phenotype associations are implemented under `analysis/02_GSE118761_airway_programs/` and `analysis/03_GSE152004_primary_validation/`.
- **Latent structure and epithelial states:** PCA and fixed k-means analyses use the locked cohort-specific inputs and parameters documented in `config/` and `docs/`.
- **Cross-cohort projection:** the executable nearest-centroid workflow, including scaling, centroids, distances, assignments, and verification, is under `analysis/02_GSE118761_airway_programs/projection/`.
- **Hospital models:** public scripts specify the clinical models and table structures; protected individual-level hospital data are not included.

## Robustness and sensitivity analyses

The repository retains prespecified multicollinearity and residualization diagnostics, alternative T2 and repair-ECM score construction, Reactome ECM and NABA Core Matrisome analyses, overlap-exclusion analyses, matched-random-set specificity checks, standard GSVA ssGSEA sensitivity, and retained preprocessing comparisons. External ECM definitions are sensitivity analyses and do not enter the fixed 32-gene derivation.

Relevant code and outputs are located under:

- `analysis/03_GSE152004_primary_validation/06_multicollinearity_and_residualization_recovered.R`
- `analysis/03_GSE152004_primary_validation/08_score_construction_sensitivity.R`
- `analysis/05_external_ECM_sensitivity/`
- `outputs/03_GSE152004/`
- `outputs/05_external_ECM/`

## Machine-readable outputs

- `outputs/`: cohort-level scores, models, PCA, state assignments, projection results, and sensitivity outputs.
- `figure_source_data/` and `figures/figure_source_map.tsv`: source data and provenance for manuscript figures.
- `tables/` and `tables/table_source_map.tsv`: machine-readable table source data and generation scripts.
- `docs/manuscript_to_code_map.tsv`: manuscript-wide analysis-to-code mapping.
- `audit/`: release gates, independent checks, public-data rebuild records, and security evidence.
- `docs/provenance/`: detailed component-level provenance reports that are intentionally secondary to the manuscript-wide workflow.

Final panel assembly was performed separately. Locked preview images are not represented as having been regenerated by the analysis scripts unless an explicit figure-generation script is supplied.

## Validation status

- Complete repository test suite: **14/14 PASS**
- Public-data rebuild workflow: **PASS**
- GSE18965 deterministic, code-based repair-ECM derivation: **PASS**
- Public cohort preprocessing and scoring workflows: **PASS**
- PCA, epithelial-state, projection, and statistical workflows: **PASS**
- Sensitivity and robustness workflows: **PASS**
- Machine-readable figure and table source data: **COMPLETE**
- Hospital cohort: **CODE_ONLY**
- Privacy and security scan: **PASS**

Run the complete automated test suite with:

```bash
Rscript tests/run_tests.R
```

## Data availability and privacy

Public accession identifiers, download interfaces, cache variables, and input schemas are documented in `data/README.md`. GSE18965 inputs required for the deterministic derivation are included where redistribution is permitted; official GSE118761 and GSE152004 cache files remain external.

Individual-level hospital data are not deposited because of ethical and privacy restrictions. The repository provides analysis code, a variable dictionary, a synthetic input template, and aggregate reference outputs. Real hospital data are never loaded by `run_public_analyses.R`.

## Repository structure

| Path | Contents |
|---|---|
| `analysis/` | End-to-end cohort and sensitivity workflows |
| `R/` | Repair-ECM derivation and shared R implementation |
| `config/` | Locked workflow parameters and verification hashes |
| `data/` | Redistributable inputs, annotations, gene sets, and clinical templates |
| `outputs/` | Machine-readable analysis results |
| `figure_source_data/`, `figures/` | Figure source data, maps, and supplied previews |
| `tables/` | Table source data and generation scripts |
| `docs/` | Workflow, methods, manuscript mapping, and detailed provenance |
| `tests/` | Automated repository validation suite |
| `audit/` | Release-level reproducibility and security evidence |
| `archive/development_history/` | Superseded development artifacts excluded from formal workflows |

## Code and archive

GitHub repository: https://github.com/biggsevier24-bot/childhood-asthma-epithelial-programs

Zenodo archived release: https://doi.org/10.5281/zenodo.22917874

Version: `v1.0.0`

The repository manifest excludes itself from hash validation. Development-history files are not read by formal analysis scripts or tests.

## License

The repository is released under the MIT License; see `LICENSE`.
