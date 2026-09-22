options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
required_scripts <- c(
  "R/common/public_geo_inputs.R",
  "analysis/02_GSE118761_airway_programs/00_download_GSE118761.R",
  "analysis/02_GSE118761_airway_programs/01_prepare_GSE118761.R",
  "analysis/03_GSE152004_primary_validation/00_download_GSE152004.R",
  "analysis/03_GSE152004_primary_validation/01_prepare_GSE152004.R"
)
stopifnot(all(file.exists(file.path(repo, required_scripts))))
txt <- paste(vapply(file.path(repo, required_scripts), function(x) paste(readLines(x, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
stopifnot(grepl("ftp.ncbi.nlm.nih.gov", txt, fixed = TRUE), grepl("varianceStabilizingTransformation", txt, fixed = TRUE))
meta118 <- file.path(repo, "data/public/GSE118761/GSE118761_metadata.csv")
meta152 <- file.path(repo, "data/public/GSE152004/GSE152004_metadata.csv")
if (file.exists(meta118)) {
  m <- read.csv(meta118, check.names = FALSE)
  stopifnot(nrow(m) == 104L, sum(m$tissue == "nasal") == 55L, sum(m$tissue == "tracheal") == 49L)
  stopifnot(all(m$asthma[grepl("non-asthmatic", m$diagnosis, ignore.case = TRUE)] == 0L))
  stopifnot(all(m$atopy[grepl("^Non-atopic", m$diagnosis, ignore.case = TRUE)] == 0L))
  stopifnot(all(m$wheeze[grepl("non-wheezer", m$sample_title, ignore.case = TRUE)] == 0L))
}
if (file.exists(meta152)) stopifnot(nrow(read.csv(meta152, check.names = FALSE)) == 695L)
cat("PASS: official GEO public-input rebuild interface\n")
