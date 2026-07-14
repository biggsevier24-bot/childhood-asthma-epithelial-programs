# Public Data Download Instructions

This repository does not include full GEO expression files or raw sequencing/microarray files.

Place downloaded public data under:

- `data/public/GSE18965/`
- `data/public/GSE118761/`
- `data/public/GSE152004/`

Download pages:

- GSE18965: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18965
- GSE118761: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118761
- GSE152004: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152004

Expected public inputs:

- GSE18965: GEO Series Matrix, GPL96 annotation, and optional CEL files for CEL-RMA/CEL-GCRMA sensitivity analysis.
- GSE118761: expression matrix and metadata sufficient to identify nasal and tracheal epithelial samples.
- GSE152004: expression matrix and metadata sufficient to compute T2, IFN, and repair-ECM scores and molecular T2-high status.

Protected hospital data are not public data and must not be placed in these folders.
