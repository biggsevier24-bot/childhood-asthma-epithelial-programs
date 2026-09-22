options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/public_geo_inputs.R"))
raw_dir <- file.path(repo, "data/public/GSE118761/raw")
download_or_copy_geo(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118761/suppl/GSE118761_genecounts.csv.gz",
  file.path(raw_dir, "GSE118761_genecounts.csv.gz"), "GSE118761_COUNTS_CACHE"
)
download_or_copy_geo(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118761/matrix/GSE118761_series_matrix.txt.gz",
  file.path(raw_dir, "GSE118761_series_matrix.txt.gz"), "GSE118761_SERIES_CACHE"
)
cat("GSE118761 official GEO inputs ready\n")
