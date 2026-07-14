# Reproducibility Notes

This public archive is designed for transparent code availability and provenance. It does not claim that all historical intermediate steps can be exactly reconstructed.

## Known Limitations and Analysis Provenance

1. The repair-ECM 32-gene score is a fixed historical discovery-derived score. The unavailable historical candidate list, clustering output, and connectivity threshold are not reconstructed.
2. GSE18965 Series Matrix analysis is the manuscript-method main analysis. CEL-RMA and CEL-GCRMA are sensitivity analyses and can differ because of preprocessing, probe summarization, and annotation choices.
3. Current software versions may produce small numerical differences in PCA, clustering, bootstrap, or model estimates.
4. The locked manuscript Figure/Table values must not be overwritten by reruns. In particular, the final revision instructions preserve hospital primary asthma OR=2.71 and GSE118761 projection counts nasal E1/E2/E3=30/18/7 and tracheal E1/E2/E3=15/29/5.
5. Some machine-readable audit outputs included under `outputs/audit_reruns/machine_readable_tables/` may reflect final locked audit reruns rather than the final formatted manuscript tables. They are retained for transparency and should be interpreted together with `docs/analysis_decisions.md`.
6. Hospital clinical results cannot be publicly reproduced without protected de-identified patient-level data. The synthetic hospital file is a format demonstration only.
7. Figure 1 is an overall study design and analysis workflow figure, not a standalone GSE18965 statistical analysis. Figure 2 is the GSE18965 fixed repair-ECM/differential-expression/preprocessing-sensitivity figure.
8. Table 1 is a cohort characteristics and analytical-role table, not the hospital baseline characteristics table. Table 2 contains hospital baseline characteristics and requires protected individual-level clinical data. Table 3 combines public transcriptomic model components with protected hospital clinical model components.
9. Before the first Zenodo archive is created, the repository intentionally records the archived release as pending and does not include an invalid placeholder DOI in `CITATION.cff` or `.zenodo.json`.

## Software Versions

Actual package records are provided in `environment/`. Older historical software versions that were not recoverable are not fabricated.

