# Figure 4 provenance

## Status

GSE118761 Figure 4 PCA provenance is computationally closed from the preserved repository input. The canonical implementation is `analysis/02_GSE118761_airway_programs/03_PCA.R`.

## Canonical GSE118761 implementation

For each tissue separately (nasal and tracheal):

1. rank genes by within-tissue variance;
2. retain the top 2,000 variable genes;
3. run PCA with `prcomp(..., center=TRUE, scale.=TRUE)` and retain PC1-PC10;
4. calculate Spearman correlations with T2, IFN, fixed 32-gene repair-ECM, atopy, asthma, and wheeze;
5. write sample-level PCA scores and the complete correlation table as machine-readable outputs.

The nasal analysis reproduces the manuscript Figure 4 values:

- PC1-IFN: rho = 0.715079365079365 (display 0.715)
- PC3-T2: rho = -0.608585858585859 (display -0.609)
- PC6-repair-ECM: rho = 0.299350649350649 (display 0.299)
- PC7-atopy: rho = 0.451341281866493 (display 0.451)
- PC9-asthma: rho = 0.254309047143049 (display 0.254)
- PC9-wheeze: rho = 0.254309047143049; asthma and wheeze are identically coded in the analyzed nasal samples.

## Canonical source files

- `GSE152004_PCA_correlations.csv`
- `GSE152004_PCA_scores.csv`
- `GSE118761_nasal_PCA_correlations.csv`
- `GSE118761_nasal_PCA_scores.csv`
- `GSE118761_tracheal_PCA_correlations.csv`
- `GSE118761_tracheal_PCA_scores.csv`
- `Figure4_core_summary.tsv`
- `Figure4_locked_key_values.tsv`

Pre-closure current-repository comparison files and the pre-closure `03_PCA.R` are retained only under `archive/development_history/Figure4_GSE118761_PCA_preclosure/` for audit history. They are not formal analysis inputs.

## GSE152004

The GSE152004 PCA source remains authoritative and unchanged. The copies in this directory are byte-identical to the corresponding canonical outputs under `outputs/03_GSE152004/`.
