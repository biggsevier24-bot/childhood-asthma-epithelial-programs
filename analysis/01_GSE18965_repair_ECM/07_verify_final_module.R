options(stringsAsFactors = FALSE)
suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "R/07_verify_final_module.R"
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
derived_file <- file.path(repo_root, "outputs/final_derived_genes.tsv")
reference_file <- file.path(repo_root, "data/reference/manuscript_reference32.csv")
stopifnot(file.exists(derived_file), file.exists(reference_file))
derived <- sort(unique(toupper(trimws(read.delim(derived_file, check.names = FALSE)[[1]]))))
reference <- sort(unique(toupper(trimws(read.csv(reference_file, check.names = FALSE)[[1]]))))
missing <- setdiff(reference, derived)
extra <- setdiff(derived, reference)
overlap <- intersect(reference, derived)
result <- data.frame(
  REFERENCE_N = length(reference), DERIVED_N = length(derived), OVERLAP_N = length(overlap),
  MISSING_N = length(missing), EXTRA_N = length(extra),
  MISSING = paste(missing, collapse = ";"), EXTRA = paste(extra, collapse = ";"),
  JACCARD = length(overlap) / length(union(reference, derived)),
  EXACT_MATCH_32 = setequal(reference, derived), stringsAsFactors = FALSE
)
write.table(result, file.path(repo_root, "outputs/verification.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (name in names(result)) cat(name, " = ", result[[name]][1], "\n", sep = "")
stopifnot(result$REFERENCE_N == 32L, result$DERIVED_N == 32L, result$EXACT_MATCH_32)
