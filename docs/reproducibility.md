# Reproducibility

The GSE18965 derivation is base-R and directly rerunnable from included files. GSE118761 and GSE152004 are rebuilt from official GEO count matrices and Series Matrix metadata by the download and preparation scripts listed in `data/README.md`.

Reference outputs are retained so reviewers can inspect reported values. The final GSE118761 projection is executable from the formal standardized inputs and preserves sample-level scores, tissue-specific scaling parameters, GSE152004 reference centroids, all distances, final assignments, and dependent source data. Historical unreproduced aggregate values are retained only under `archive/development_history/GSE118761_projection_historical_unreproduced/` and are never read by the formal pipeline.

R 4.2.3 was used for the frozen exact32 workflow. Package and operating-system details from the release validation are stored in `sessionInfo.txt`.
