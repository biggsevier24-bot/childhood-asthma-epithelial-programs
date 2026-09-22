options(stringsAsFactors = FALSE)
suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "tests/test_no_hardcoding.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
formal_scripts <- file.path(repo_root, "R", sprintf("%02d_%s.R", 0:6, c(
  "prepare_GSE18965", "map_GPL96_probes", "build_candidate_pool",
  "generate_expression_features", "generate_network_features",
  "repair_ECM_final_selection_rule", "run_full_derivation"
)))
stopifnot(all(file.exists(formal_scripts)))
text <- paste(unlist(lapply(formal_scripts, readLines, warn = FALSE)), collapse = "\n")
reference <- toupper(trimws(read.csv(file.path(repo_root, "data/reference/manuscript_reference32.csv"), check.names = FALSE)[[1]]))
symbol_hits <- reference[vapply(reference, function(symbol) {
  grepl(paste0("(?<![A-Z0-9])", symbol, "(?![A-Z0-9])"), toupper(text), perl = TRUE)
}, logical(1))]
reference_file_used <- grepl("manuscript_reference32|repair_ECM_32_genes|data[/\\\\]reference", text, ignore.case = TRUE)
top32 <- grepl("top[_ -]?32|head\\s*\\([^)]*,\\s*32\\s*\\)|rank\\s*<=?\\s*32", text, ignore.case = TRUE, perl = TRUE)
manual <- grepl("manual[_ -]?(include|exclude|override|add|remove)", text, ignore.case = TRUE, perl = TRUE)
gene_specific <- length(symbol_hits) > 0L
result <- c(
  paste0("GENE_NAME_SPECIFIC_RULE = ", toupper(as.character(gene_specific))),
  paste0("REFERENCE_FILE_USED_IN_DERIVATION = ", toupper(as.character(reference_file_used))),
  paste0("TOP32_TRUNCATION = ", toupper(as.character(top32))),
  paste0("MANUAL_OVERRIDE = ", toupper(as.character(manual))),
  paste0("MATCHED_REFERENCE_SYMBOLS = ", paste(symbol_hits, collapse = ";"))
)
dir.create(file.path(repo_root, "audit"), showWarnings = FALSE)
writeLines(result, file.path(repo_root, "audit/no_hardcoding_audit.txt"))
cat(paste(result, collapse = "\n"), "\n")
stopifnot(!gene_specific, !reference_file_used, !top32, !manual)
