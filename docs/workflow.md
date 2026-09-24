# Workflow

1. Read the deposited GSE18965 log2 Series Matrix and assign 9 AA and 7 HN samples from deposited titles.
2. Parse GPL96 and retain one probe per gene: the highest mean expression across all 16 samples, with lexical probe ID as the exact-tie rule.
3. Build the 140-gene candidate pool using prespecified GO Biological Process and curated gene-family annotations.
4. Calculate AA-HN effects and signed/absolute Pearson correlations.
5. Construct the fixed feature cluster and the fixed 115-gene network-prefilter set.
6. Calculate the 11 rule features.
7. Apply the frozen depth-5 explicit decision rule in `R/05_repair_ECM_final_selection_rule.R`.
8. Write 32 genes before loading the manuscript reference list.
9. Verify exact membership separately with `R/07_verify_final_module.R`.

A deterministic, code-based workflow was applied in GSE18965 to derive the fixed 32-gene repair-ECM module. It uses the AA-HN expression contrast, co-expression structure, hierarchical clustering, network-based filtering, and explicit programmatic selection rules. No manual gene-name overrides or post hoc truncation are used. The explicit selection rule was calibrated against the locked manuscript module membership; the derivation script does not load the reference list, which is used only by the separate post-derivation verification step.
