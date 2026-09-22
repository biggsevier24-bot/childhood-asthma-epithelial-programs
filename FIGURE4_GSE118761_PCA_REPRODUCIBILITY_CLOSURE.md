# Figure 4 / GSE118761 PCA reproducibility closure

## Decision

`FIGURE4_GSE118761_PCA_REPRODUCIBILITY = PASS`

The historical Figure 4 GSE118761 PCA recipe was tested against the preserved repository input and reproduced all five locked manuscript relationships without parameter tuning or hard-coded analysis targets.

## Recovered implementation

- nasal and tracheal analyzed separately
- genes ranked by within-tissue variance
- top 2,000 variable genes retained
- PCA: `prcomp(..., center=TRUE, scale.=TRUE)`, PC1-PC10
- Spearman association with T2, IFN, fixed repair-ECM, atopy, asthma, and wheeze

## Validation results

| Relationship | Recomputed rho | Manuscript display | Status |
|---|---:|---:|---|
| PC1-IFN | 0.715079365079365 | 0.715 | PASS |
| PC3-T2 | -0.608585858585859 | -0.609 | PASS |
| PC6-repair-ECM | 0.299350649350649 | 0.299 | PASS |
| PC7-atopy | 0.451341281866493 | 0.451 | PASS |
| PC9-asthma/wheeze | 0.254309047143049 | 0.254 | PASS |

In the analyzed nasal samples, asthma and wheeze are identically coded, so PC9 has the same correlation with both variables.

## Repository changes

- canonical `analysis/02_GSE118761_airway_programs/03_PCA.R` updated to the recovered historical implementation
- full PC1-PC10 correlation tables now include T2, IFN, repair-ECM, atopy, asthma, wheeze, P value, and n
- sample-level PCA score outputs are generated
- Figure 4 source-data files are generated directly from canonical code
- `tests/test_GSE118761_Figure4_PCA.R` added and registered in `tests/run_tests.R`
- pre-closure PCA script/comparison outputs archived under `archive/development_history/Figure4_GSE118761_PCA_preclosure/`
- README and Figure 4 provenance documentation updated

## Release boundary

This closure resolves the Figure 4 PCA provenance issue. The recovered PCA script and dedicated test passed in the final R-enabled runtime validation, the complete repository suite passed 14/14, and the public-data rebuild completed successfully using official GEO files supplied through the documented external cache mechanism.

## Final runtime validation

An independent numerical preflight against retained machine-readable outputs passed all 21 checked items, including exact32 status, GSE152004 primary sample/event counts and core ORs, state counts, alternative repair-score OR, GSE118761 sample/projection counts, all Figure 4 key correlations, source-data dimensions, asthma/wheeze identity in the nasal subset, and the absence of obvious protected hospital raw data. The detailed table is `audit/Figure4_GSE118761_PCA_independent_validation.tsv`.

Final validation used R 4.2.3. `analysis/02_GSE118761_airway_programs/03_PCA.R` passed, `tests/test_GSE118761_Figure4_PCA.R` passed, the complete suite passed 14/14, and `run_public_analyses.R` completed with exit code 0. No PCA parameter, test threshold, Figure 4 value, or manuscript value was changed during runtime validation.
