# Public release checklist

Status: READY

## Manuscript-wide workflow

- [x] GSE18965, GSE118761, and GSE152004 public-data acquisition and preprocessing interfaces are documented and executable.
- [x] Fixed repair-ECM, T2, and IFN program definitions and scoring code are included.
- [x] Correlation, regression, effect-estimate, confidence-interval, PCA, epithelial-state, and phenotype-association workflows are included.
- [x] Cross-cohort projection is executable from sample-level public-derived inputs with scaling, centroid, distance, assignment, and verification records.
- [x] Multicollinearity, residualization, alternative score construction, external ECM, overlap-exclusion, matched-random-set, and ssGSEA sensitivity resources are retained.
- [x] Machine-readable cohort results, figure source data, table source data, and manuscript-to-code maps are included.

## Validation and provenance

- [x] Complete repository test suite passed 14/14.
- [x] Public-data rebuild workflow completed using official GEO inputs supplied through documented cache interfaces.
- [x] Deterministic repair-ECM derivation and post-derivation verification passed.
- [x] Public cohort preprocessing, scoring, PCA, epithelial-state, projection, and statistical workflows passed.
- [x] Detailed component provenance is retained below the repository root under `docs/provenance/`, `figure_source_data/`, and `audit/`.
- [x] Superseded development artifacts are isolated under `archive/development_history/` and are not formal pipeline inputs.

## Data protection and release integrity

- [x] Hospital-cohort analysis remains code-only with a variable dictionary, synthetic template, and aggregate outputs.
- [x] Protected participant-level hospital data are excluded.
- [x] GSE118761 and GSE152004 raw GEO cache files are excluded.
- [x] Privacy, credential, secret, and absolute-path scans passed.
- [x] Repository manifest covers all release files except itself.
