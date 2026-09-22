options(stringsAsFactors = FALSE)
suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "tests/test_full_derivation.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")

run <- system2(rscript, c("--vanilla", shQuote(file.path(repo_root, "R/06_run_full_derivation.R"))),
               stdout = TRUE, stderr = TRUE)
if (!is.null(attr(run, "status")) && attr(run, "status") != 0L) stop(paste(run, collapse = "\n"))
verify <- system2(rscript, c("--vanilla", shQuote(file.path(repo_root, "R/07_verify_final_module.R"))),
                  stdout = TRUE, stderr = TRUE)
if (!is.null(attr(verify, "status")) && attr(verify, "status") != 0L) stop(paste(verify, collapse = "\n"))
result <- read.delim(file.path(repo_root, "outputs/verification.tsv"), check.names = FALSE)
stopifnot(result$REFERENCE_N == 32L, result$DERIVED_N == 32L, result$OVERLAP_N == 32L,
          result$MISSING_N == 0L, result$EXTRA_N == 0L, result$EXACT_MATCH_32)
cat("REPAIR_ECM_FINAL_PASS = TRUE\nFINAL_N = 32\nEXACT_MATCH_32 = TRUE\n")
