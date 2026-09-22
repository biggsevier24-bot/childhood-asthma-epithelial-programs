# External sensitivity inputs

External ECM definitions are deliberately independent of the fixed 32-gene derivation. To run `R/09_external_ECM_sensitivity.R`, provide:

- `standardized_gene_expression.tsv`: gene symbols in rows and samples in columns.
- `external_ECM_gene_sets.gmt`: fixed external definitions in GMT format.

These inputs are not required by `R/06_run_full_derivation.R` and are not used to force 32-gene membership.
