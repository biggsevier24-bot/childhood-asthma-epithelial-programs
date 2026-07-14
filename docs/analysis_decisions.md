# Analysis Decisions

## Fixed Thresholds

- GSE18965 log2FC is defined as atopic asthma minus non-atopic healthy controls.
- Primary DEG threshold: FDR < 0.05 and absolute log2FC >= 0.5.
- GSE18965 Series Matrix is the manuscript-method main analysis.
- CEL-RMA and CEL-GCRMA are raw-data sensitivity analyses.

## Repair-ECM Historical Score

The 32-gene repair-ECM list was retained as a fixed historical discovery-derived score.

The complete historical intermediate candidate list, correlation-clustering output, and prespecified connectivity threshold were unavailable.

The exact historical selection pathway could not be fully reconstructed.

The revision therefore audited the fixed score using differential-expression results, alternative preprocessing, external ECM definitions, and sensitivity analyses.

The repository does not claim to reconstruct an unavailable historical selection procedure.

## External ECM Definitions

REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION and NABA_CORE_MATRISOME were used as prespecified external ECM sensitivity definitions. They are not used to revise the fixed 32-gene score.

## Clinical Data

The hospital cohort contains protected clinical data. Patient-level data are excluded from the public release. Only a data dictionary and synthetic input template are provided.

## Locked Manuscript Values

The public archive preserves the locked manuscript values documented in the final revision instructions, including:

- GSE152004 k=3 state sizes: E1=338, E2=299, E3=58.
- GSE118761 projection counts: nasal E1/E2/E3=30/18/7; tracheal E1/E2/E3=15/29/5.
- Hospital primary asthma model odds ratio: 2.71.

When audit reruns produce different values, the rerun output is treated as provenance evidence and is not used to overwrite the locked manuscript values.
