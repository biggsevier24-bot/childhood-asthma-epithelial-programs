# Locked method parameters

- Global random seed: `20260626`.
- GSE18965 effect: AA minus HN; primary limma eBayes defaults; DEG threshold FDR < 0.05 and absolute log2FC >= 0.5.
- GSE118761: gene-wise standardization is performed across the 104-sample cohort; nasal is primary and tracheal is a tissue sensitivity analysis; tissues are not pooled in ordinary regression.
- GSE152004 T2 score: mean within-cohort gene-wise z score of CST1, CLCA1, and SERPINB2; T2-high is score greater than or equal to the median (348/347).
- Primary program score: mean gene-wise z score within cohort or within the prespecified tissue stratum; regression effects are reported per one standard deviation of the resulting program score.
- GSE152004 clustering: k=3, seed 20260626, nstart 500.
- ssGSEA: official GSVA implementation only.
- External ECM definitions never enter the repair-ECM 32-gene derivation.

The complete frozen GSE18965 deterministic derivation parameters are in `config/repair_ECM_final_derivation.yaml`. The explicit selection rule was calibrated against the locked manuscript module membership; the derivation script does not load the reference list.
