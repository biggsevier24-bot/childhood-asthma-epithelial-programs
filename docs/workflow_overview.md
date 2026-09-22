# Workflow overview

1. GSE18965 rebuilds the frozen revision-era repair-ECM reconstruction and verifies the result against the locked 32 genes only after derivation.
2. GSE118761 computes cohort/tissue-specific mean gene-wise z scores, then performs nasal/tracheal analyses separately.
3. GSE152004 uses VST expression, a three-gene T2 score, the `>= median` T2-high rule, fixed k=3 clustering, and prespecified logistic diagnostics.
4. External ECM definitions are sensitivity analyses only and are isolated from GSE18965 derivation.
5. Hospital scripts expose model specifications but require protected inputs outside this repository.

Final panel assembly was performed separately. Source maps identify the machine-readable statistics supporting each panel and table.
