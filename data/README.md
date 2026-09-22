# Data inputs

## Included public inputs

`data/GSE18965/` contains the deposited GSE18965 log2 Series Matrix and GPL96 annotation snapshot needed by the exact32 workflow. Fixed annotation resources and reference gene-set definitions are also included. The fixed manuscript 32-gene membership is not read by derivation scripts and is used only for post-derivation verification.

## Public RNA-seq inputs

The repository rebuilds GSE118761 and GSE152004 standardized inputs from their official GEO count and Series Matrix files. Run:

```bash
Rscript analysis/02_GSE118761_airway_programs/00_download_GSE118761.R
Rscript analysis/02_GSE118761_airway_programs/01_prepare_GSE118761.R
Rscript analysis/03_GSE152004_primary_validation/00_download_GSE152004.R
Rscript analysis/03_GSE152004_primary_validation/00_download_raw_counts.R
Rscript analysis/03_GSE152004_primary_validation/01_prepare_GSE152004.R
```

For offline or cached runs, set `GSE118761_COUNTS_CACHE`, `GSE118761_SERIES_CACHE`, `GSE152004_COUNTS_CACHE`, and `GSE152004_SERIES_CACHE` to the corresponding downloaded files. Generated expression files use `gene_symbol` in the first column and GEO sample IDs in subsequent columns. GSE118761 phenotype parsing explicitly protects `non-wheezer`, `non-atopic`, and `non-asthmatic` labels from substring miscoding.

The downloaded raw files and generated matrices are intentionally excluded from the release ZIP to keep the repository compact; the scripts and exact official URLs are retained. The GSE152004 raw-count file is deposited at:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152004/suppl/GSE152004_695_raw_counts.txt.gz`

Place it at `data/public/GSE152004/raw/GSE152004_695_raw_counts.txt.gz`, or set `GSE152004_COUNTS_CACHE` to the downloaded file, to rerun the secondary mean-log2(CPM+1) sensitivity branch.

## Hospital cohort

No patient-level hospital data are included. `data/clinical_template/` contains a data dictionary and synthetic structural example only.
