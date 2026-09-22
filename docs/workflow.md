# Workflow

1. Read the deposited GSE18965 log2 Series Matrix and assign 9 AA and 7 HN samples from deposited titles.
2. Parse GPL96 and retain one probe per gene: the highest mean expression across all 16 samples, with lexical probe ID as the exact-tie rule.
3. Build the 140-gene candidate pool using fixed GO Biological Process and curated gene-family annotations.
4. Calculate AA-HN effects and signed/absolute Pearson correlations.
5. Construct the fixed feature cluster and the fixed 115-gene network-prefilter set.
6. Calculate the 11 rule features.
7. Apply the frozen depth-5 explicit decision rule in `R/05_repair_ECM_final_selection_rule.R`.
8. Write 32 genes before loading the manuscript reference list.
9. Verify exact membership separately with `R/07_verify_final_module.R`.

The workflow is a revision-era calibrated reconstruction. It is not the historical discovery algorithm.
