# Type 2, Interferon, and Repair-ECM Programs Characterize Airway Epithelial Heterogeneity in Childhood Asthma

This repository contains the manuscript-associated analysis code for a study of airway epithelial heterogeneity in childhood asthma. The analyses evaluate three epithelial programs: Type 2 inflammation, interferon response, and a fixed historical repair-ECM score.

## Overview

The project integrates public transcriptomic cohorts and a protected hospital clinical cohort to describe epithelial program heterogeneity. Public expression data are used for differential-expression analysis, module scoring, PCA, clustering, cross-cohort projection, and sensitivity analyses. The hospital cohort is used only for the locked clinical Type 2 analyses reported in the manuscript.

## Study Cohorts

- **GSE18965**: Affymetrix GPL96 airway epithelial microarray data; 9 atopic asthma and 7 non-atopic healthy controls in the locked comparison.
- **GSE118761**: nasal and tracheal epithelial RNA-seq data; 55 nasal and 49 tracheal epithelial samples in the locked projection analyses.
- **GSE152004**: airway epithelial RNA-seq data; 695 samples in the locked T2/IFN/repair-ECM scoring, PCA, and k=3 epithelial-state analyses.
- **Hospital clinical cohort**: 325 children, including 158 with physician-diagnosed asthma. The original patient-level data are protected and are not included in this repository.

## Repository Scope

This repository provides recovered R source code, fixed gene-set manifests, reference outputs, and analysis documentation. Independent original Python analysis source files could not be recovered and are explicitly documented as unavailable.

Public transcriptomic data can be obtained from GEO. Protected hospital patient-level data are not shared. Only a data dictionary and a small synthetic template are provided for the clinical input format.

## Analysis Workflow

1. Overall study design and cohort-level workflow summary.
2. GSE18965 preprocessing, differential expression, fixed 32-gene repair-ECM audit, external ECM definitions, and microarray preprocessing sensitivity.
3. GSE118761 nasal repair-ECM/IFN scoring, IFN-ECM quadrants, and atopic versus non-atopic wheeze discrimination.
4. GSE152004 and GSE118761 PCA-derived latent factors, with PCA run independently in each cohort.
5. GSE152004 k=3 epithelial-state clustering, clustering stability, external ECM sensitivity, and projection to GSE118761 nasal and tracheal samples.
6. Hospital clinical Type 2 models, using protected clinical data only when supplied by an authorized user.
7. GSE152004 molecular T2-high logistic models.
8. Locked tables, figure previews, and output checks.

## Figure and Table Provenance

Figure 1 summarizes the overall study design and analysis workflow. It is not a GSE18965 public-data statistical plot.

Figure 2 is the GSE18965 figure, covering the fixed 32-gene repair-ECM audit, differential expression, external ECM definitions, and microarray preprocessing sensitivity. The public code can reproduce the audit from public data, but the historical 32-gene original selection pathway cannot be fully reconstructed.

Figure 3 uses GSE118761 for nasal repair-ECM/IFN scores, IFN-ECM quadrants, and atopic versus non-atopic wheeze discrimination.

Figure 4 uses GSE152004 and GSE118761 for independently run PCA-derived latent factors and correlations with T2, IFN, repair-ECM, and clinical variables.

Figure 5 uses GSE152004 and GSE118761 for the locked three-state clustering and projection framework. The locked manuscript results are retained as GSE152004 E1/E2/E3=338/299/58, GSE118761 nasal E1/E2/E3=30/18/7, and GSE118761 tracheal E1/E2/E3=15/29/5.

Figure 6a depends on protected hospital clinical data and cannot be fully reproduced from public data. Figure 6b uses GSE152004 molecular T2-high logistic models and can be reproduced from public transcriptomic data and metadata.

Table 1 describes the characteristics and analytical roles of the study cohorts. Table 2 is the protected hospital-cohort baseline characteristics table. Table 3 combines public transcriptomic models and protected hospital clinical models. Detailed panel-level mapping is in `docs/figure_table_code_map.csv`.

## Reproducibility Categories

Analyses in this repository are classified into five categories:

1. **Workflow/reference documentation**: documentation or figure structure is provided, but no statistical rerun is implied.
2. **Partly rerunnable from public data**: recovered code and public data support part of the rerun, while specific historical branches or intermediate files remain incompletely reconstructable.
3. **Rerunnable only with required intermediate inputs**: recovered source code and audit outputs are provided, but complete end-to-end rerun additionally requires standardized expression matrices, processed metadata, fixed-score files, random-set memberships, or other intermediates listed in `data_manifest/required_intermediate_inputs.csv`.
4. **Reference output only / historical branch incompletely reconstructable**: reference outputs and analysis documentation are provided, but a historical source script or intermediate selection branch is unavailable.
5. **Not publicly reproducible**: analysis depends on protected individual-level hospital clinical data.

Figure 1 is workflow/reference documentation.

Figure 2 is partly rerunnable from public GSE18965 data; the historical 32-gene selection and some CEL intermediate branches are incompletely reconstructable.

Figures 3-5 and Figure 6b have recovered source code and audit outputs in this repository. Complete end-to-end public rerun additionally requires the intermediate standardized expression/metadata inputs listed in `data_manifest/required_intermediate_inputs.csv`.

Figure 6a and hospital analyses are not publicly reproducible because protected clinical data are required.

## Software Environment

The environment files in `environment/` were generated from the analysis workstation used to prepare the final locked branch:

- R: 4.2.3
- Python: 3.13.14
- R package versions: `environment/installed_R_packages.csv`
- R session details: `environment/R_sessionInfo.txt`
- Python package targets: `environment/requirements.txt`

Install R and Python dependencies in your preferred environment. Then place public GEO files and required intermediate inputs under the relative paths documented in `data_manifest/public_data_download_instructions.md` and `data_manifest/required_intermediate_inputs.csv`. The repository intentionally uses relative paths such as `data/public/GSE18965/`, `data/public/GSE118761/`, `data/public/GSE152004/`, and `outputs/user_rerun/`.

Run:

```bash
bash run_all.sh
```

`run_all.sh` calls recovered R analysis scripts only when preflight checks find the required inputs for each step. Missing public or intermediate inputs print `SKIP: required intermediate input unavailable` and do not trigger a failing analysis attempt. The script writes rerun outputs to `outputs/user_rerun/`, does not overwrite `outputs/manuscript_locked/`, and does not run protected hospital models unless `HOSPITAL_PRIVATE_DATA` points to an authorized local de-identified file.

## Locked Manuscript Outputs

This repository corresponds to the final manuscript-associated locked analysis version. The current figures, main tables, and supplementary tables are treated as manuscript-locked reference outputs. Manuscript-locked reference materials are provided in `outputs/manuscript_locked/`. Audit rerun outputs are provided separately in `outputs/audit_reruns/` and must not be described as manuscript source data.

If reruns with different software versions, random initialization behavior, or unrecoverable historical branches produce small numeric differences, use the locked outputs as the manuscript reference and see `docs/reproducibility_notes.md`.

The repository does not claim that all historical intermediate steps are fully reconstructable.

## Repair-ECM Historical Score Provenance

The 32-gene repair-ECM list was retained as a fixed historical discovery-derived score. The complete historical intermediate candidate list, correlation-clustering output, and prespecified connectivity threshold were unavailable. The exact historical selection pathway could not be fully reconstructed. The revision therefore audited the fixed score using differential-expression results, alternative preprocessing, external ECM definitions, and sensitivity analyses. This repository does not claim to reconstruct an unavailable historical selection procedure.

## Data Availability

Public data:

- GSE18965: GEO accession GSE18965
- GSE118761: GEO accession GSE118761
- GSE152004: GEO accession GSE152004

Protected clinical data:

The hospital patient-level data are not publicly shared because of privacy and ethics restrictions. De-identified data may be considered by the corresponding author upon reasonable request and institutional ethics approval.

## Code Availability

GitHub repository: pending

Archived release: pending

Code availability scope:

- recovered R analysis scripts
- fixed gene-set manifests
- public data and required intermediate input manifests
- manuscript-locked reference outputs and audit rerun outputs
- source-code availability notes documenting unavailable independent Python source files

## Contact

Corresponding author: Chuangli Hao  
Children's Hospital of Soochow University  
Email: hcl_md@sina.com

