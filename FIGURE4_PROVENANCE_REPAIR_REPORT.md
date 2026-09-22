> SUPERSEDED BY `FIGURE4_GSE118761_PCA_REPRODUCIBILITY_CLOSURE.md`. Retained for audit history.

# Figure 4 provenance repair report

## Scope

This repair addresses only Figure 4 source-data provenance. No statistical analysis, PCA, correlation analysis, module scoring, figure generation, or manuscript editing was performed.

## Repository changes

The repository now contains `figure_source_data/Figure4/`, which separates four provenance classes:

1. **Authoritative current GSE152004 source**: copied byte-for-byte from `outputs/03_GSE152004/`.
2. **Historical-like GSE118761 machine-readable source**: a retained 30-row nasal PCA correlation table that contains the historical repair-ECM, atopy, and asthma/wheeze signals.
3. **Locked GSE118761 historical values**: five key Figure 4 correlations, including two values whose exact standalone historical machine-readable source remains pending.
4. **Current GSE118761 repository comparison**: copied byte-for-byte from the current reproducible outputs and clearly labeled as comparison data rather than historical Figure 4 source.

The canonical current file `outputs/02_GSE118761/GSE118761_nasal_PCA_correlations.csv` was not modified.

## Locked GSE118761 values

| PC | Variable | Locked rho | Provenance status |
|---|---|---:|---|
| PC1 | IFN | 0.715 | Historical locked value; standalone historical machine-readable source not yet isolated |
| PC3 | T2 | -0.609 | Historical locked value; standalone historical machine-readable source not yet isolated |
| PC6 | repair-ECM | 0.299 | Historical-like table contains 0.299350649350649 |
| PC7 | atopy | 0.451 | Historical-like table contains 0.451341281866493 |
| PC9 | wheeze/asthma | 0.254 | Historical-like table contains 0.254309047143049 |

## Files added under Figure 4 source data

- `figure_source_data/Figure4/FIGURE4_PROVENANCE_NOTE.md`
- `figure_source_data/Figure4/Figure4_core_summary.tsv`
- `figure_source_data/Figure4/Figure4_locked_key_values.tsv`
- `figure_source_data/Figure4/GSE152004_PCA_correlations.csv`
- `figure_source_data/Figure4/GSE152004_PCA_scores.csv`
- `figure_source_data/Figure4/GSE152004_processed_metadata.csv`
- `figure_source_data/Figure4/GSE152004_module_scores_recomputed.csv`
- `figure_source_data/Figure4/GSE118761_nasal_PCA_correlations_historical_like.csv`
- `figure_source_data/Figure4/GSE118761_nasal_PCA_correlations_current_repo.csv`
- `figure_source_data/Figure4/GSE118761_nasal_PCA_scores_current_repo.csv`
- `figure_source_data/Figure4/GSE118761_processed_metadata.csv`
- `figure_source_data/Figure4/GSE118761_all_sample_module_scores.csv`

## Other repository files

- Added: `FIGURE4_PROVENANCE_REPAIR_REPORT.md`
- Modified: `README.md` (added only the `Figure 4 provenance` section)
- Mechanically refreshed: `repository_file_manifest.csv`

## Validation

- All 12 required Figure 4 source-data/provenance files are present.
- The four GSE152004 copies are byte-identical to their authoritative repository outputs.
- The four current GSE118761 comparison copies are byte-identical to their current repository outputs.
- The historical-like GSE118761 table contains 30 correlation rows and exactly reproduces the PC6-repair-ECM, PC7-atopy, and PC9-asthma locked values after rounding.
- `Figure4_core_summary.tsv` uses only the four allowed `source_type` values.
- No file under `outputs/` was modified.
- No file under `figures/` was modified.
- No manuscript file was modified.
- No figure was generated or changed.

## Interpretation boundary

The current GSE118761 PCA output is retained for reproducibility, comparison, and debugging. It is not represented as the historical Figure 4 source and must not be silently substituted for the locked historical values.
