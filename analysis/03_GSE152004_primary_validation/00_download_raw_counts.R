options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(repo, "R/common/public_geo_inputs.R"))

destination <- file.path(
  repo, "data/public/GSE152004/raw/GSE152004_695_raw_counts.txt.gz"
)
download_or_copy_geo(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE152nnn/GSE152004/suppl/GSE152004_695_raw_counts.txt.gz",
  destination,
  "GSE152004_COUNTS_CACHE"
)
cat("GSE152004 official raw-count matrix ready: ", destination, "\n", sep = "")
