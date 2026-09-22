# Reconstruction provenance

Simple biologically constrained B1 rules did not recover exact membership. B2 therefore used the fixed 32-gene manuscript membership as development labels while excluding gene symbols and all reference-derived frequency/near-solution features from the predictors.

The prescribed decision-tree search evaluated depths 3-10, minimum leaf sizes 1-3, gini/entropy criteria, and unweighted/balanced classes, with and without P/FDR. The selected exact solution excluded P and FDR and had depth 5, 15 leaves, 14 internal conditions, 8 positive leaf clauses, and 11 unique features.

Thresholds were simplified only when exact 32/32 membership was retained. The final rule was frozen, converted to ordinary R if/else logic, and rerun in an isolated reference-free directory. Only after the output file was written was the reference reopened for exact-set verification.
