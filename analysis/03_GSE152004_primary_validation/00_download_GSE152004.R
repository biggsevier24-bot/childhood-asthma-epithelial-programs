options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/public_geo_inputs.R"))
raw_dir <- file.path(repo, "data/public/GSE152004/raw")
source(file.path(repo, "analysis/03_GSE152004_primary_validation/00_download_raw_counts.R"))
download_or_copy_geo(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152004/matrix/GSE152004_series_matrix.txt.gz",
  file.path(raw_dir, "GSE152004_series_matrix.txt.gz"), "GSE152004_SERIES_CACHE"
)
cat("GSE152004 official GEO inputs ready\n")
