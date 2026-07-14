# Public Repository Final Audit

This final audit documents the last minimal repository correction. No statistical analysis was rerun, and no locked Figure, Table, sample size, OR, IRR, P value, cluster size, or projection result was changed.

## 1. Safety Scan Status

The repository safety scan passed using `tests/test_no_private_information.py`.

No identifiable patient data, real hospital individual-level records, local absolute paths, passwords, tokens, API keys, GitHub credentials, temporary Office files, or private clinical data were detected.

## 2. Manuscript-Locked Results

Manuscript-locked results were not changed.

`outputs/manuscript_locked/` contains manuscript-locked reference materials and must not be overwritten by reruns. `outputs/audit_reruns/` contains audit rerun outputs that may differ from final formatted manuscript tables and must not be described as manuscript source data.

## 3. Analyses That Can Be Directly Rerun

No manuscript Figure 2-6 analysis is represented as a guaranteed one-command rerun from raw GEO downloads alone.

Repository checks and documentation/workflow validation can be run directly. Statistical reruns require the public and intermediate inputs listed in `data_manifest/required_intermediate_inputs.csv`.

## 4. Analyses Requiring Missing Intermediate Inputs

Figures 3-5 and Figure 6b have recovered source code and audit outputs, but complete end-to-end public rerun requires standardized expression matrices, processed metadata, fixed score files, external gene-set resources, random-set membership files, and corrected external FDR outputs listed in `data_manifest/required_intermediate_inputs.csv`.

When those files are absent, `run_all.sh` prints:

`SKIP: required intermediate input unavailable`

and does not attempt to run a failing analysis.

## 5. Reference Output Only

The following are reference-output/documentation only or incompletely reconstructable:

- Figure 1: workflow/reference documentation.
- Historical fixed 32-gene repair-ECM selection pathway: original candidate list, correlation clustering output, and connectivity threshold were unavailable.
- Supplementary Table S4A: standalone IFN-ECM quadrant profile source CSV unavailable.
- Supplementary Table S4B: standalone cross-validation fold/AUC source CSV unavailable.
- Supplementary Table S12: standalone mean-z versus official ssGSEA model output unavailable in the public archive.

## 6. Analyses Requiring Protected Clinical Data

The following are not publicly reproducible because protected individual-level clinical data are required:

- Figure 6a.
- Table 2.
- Hospital components of Table 3.
- Supplementary Table S10.
- Any hospital clinical model rerun, including asthma, previous wheeze, frequent exacerbation, exacerbation count, and ICS-adjusted sensitivity models.

The repository provides code, a data dictionary, and a synthetic input structure only.

## 7. Original Python Source Code

Independent original Python source files were not recovered.

The manuscript Methods reported use of Python packages including scikit-learn and statsmodels, but the original Python scripts implementing those analyses were unavailable in the project workspace. The `python/` files therefore explicitly report `source script unavailable`; no placeholder or reconstructed Python code is represented as the original implementation.

## 8. S4A, S4B, S12, And Table 3 Mapping Corrections

The mapping file `docs/figure_table_code_map.csv` was corrected as follows:

- Supplementary Table S4A now maps to IFN-ECM quadrant profile reference output and no longer points to projection assignments.
- Supplementary Table S4B now maps to cross-validation fold/AUC reference output and no longer points to PCA correlations.
- Supplementary Table S12 now maps to mean-z versus official ssGSEA ECM scoring comparison and no longer points to `ECM_random_set_specificity.csv`.
- Table 3 now explicitly includes GSE118761 nasal atopy/asthma regression models, GSE152004 molecular T2-high models, and protected hospital models.

## 9. Source Code Availability

Recovered real R source code was copied into `R/` and sanitized for public paths. `R/03_RNAseq_preprocessing.R` and `R/04_gene_set_scoring.R` intentionally derive from the same recovered comprehensive source script and are duplicated only for workflow mapping.

See `docs/source_code_availability.md` for the source-code availability table.

## 10. Required Intermediate Inputs

`data_manifest/required_intermediate_inputs.csv` lists required inputs for `R/01` through `R/06`, including GEO files, standardized expression matrices, processed metadata, fixed score files, external gene-set resources, random-set membership files, corrected external FDR files, and protected hospital input.

## 11. Citation And Zenodo Metadata

`CITATION.cff` and `.zenodo.json` include all five manuscript authors in order:

1. Yu-xiang Zhang
2. Yu Miao
3. Jiating Xue
4. Nana Wang
5. Chuangli Hao

README records GitHub repository and archived release as pending. No invalid angle-bracket GitHub URL or Zenodo DOI placeholders remain.

## Final Status

The public repository is ready for user-side manual review and upload. The user still needs to add the real GitHub URL and Zenodo DOI after creating those records.
