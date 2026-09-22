options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
pca_script <- file.path(repo, "analysis/02_GSE118761_airway_programs/03_PCA.R")
run <- system2(rscript, shQuote(pca_script), stdout = TRUE, stderr = TRUE)
status <- attr(run, "status")
if (!is.null(status) && status != 0L) stop(paste(run, collapse = "\n"))

x <- read.csv(file.path(repo, "outputs/02_GSE118761/GSE118761_nasal_PCA_correlations.csv"), check.names = FALSE)
get_rho <- function(pc, variable) {
  z <- x[x$PC == pc & x$variable == variable, "rho"]
  stopifnot(length(z) == 1L)
  as.numeric(z)
}

# Manuscript display targets (3 d.p.) with tolerance chosen to require the
# recomputed values to round to the locked Figure 4 values.
targets <- c(
  PC1_IFN = 0.715,
  PC3_T2 = -0.609,
  PC6_repair_ECM = 0.299,
  PC7_atopy = 0.451,
  PC9_asthma = 0.254,
  PC9_wheeze = 0.254
)
observed <- c(
  PC1_IFN = get_rho("PC1", "IFN"),
  PC3_T2 = get_rho("PC3", "T2"),
  PC6_repair_ECM = get_rho("PC6", "repair_ECM"),
  PC7_atopy = get_rho("PC7", "atopy"),
  PC9_asthma = get_rho("PC9", "asthma"),
  PC9_wheeze = get_rho("PC9", "wheeze")
)
stopifnot(all(abs(observed - targets) < 5e-4))
stopifnot(nrow(x) == 60L, all(c("PC", "variable", "rho", "P", "n") %in% names(x)))

m <- read.csv(file.path(repo, "outputs/02_GSE118761/GSE118761_recomputed_module_scores.csv"), check.names = FALSE)
nasal <- m[tolower(m$tissue) == "nasal", , drop = FALSE]
stopifnot(nrow(nasal) == 55L, identical(as.integer(nasal$asthma), as.integer(nasal$wheeze)))
stopifnot(file.exists(file.path(repo, "figure_source_data/Figure4/GSE118761_nasal_PCA_correlations.csv")))
stopifnot(file.exists(file.path(repo, "figure_source_data/Figure4/GSE118761_nasal_PCA_scores.csv")))
cat("PASS: GSE118761 historical Figure 4 PCA recipe reproduces all five locked correlations\n")
