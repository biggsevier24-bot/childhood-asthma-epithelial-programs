options(stringsAsFactors = FALSE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/", mustWork = TRUE)
roots <- file.path(repo, c("analysis", "R"))
files <- unlist(lapply(roots, list.files, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE))
text <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
dependency_pattern <- "(read|source|load|readRDS|readLines|read\\.csv|read\\.delim).{0,120}archive[/\\\\]"
stopifnot(!any(grepl(dependency_pattern, text, ignore.case = TRUE, perl = TRUE)))
cat("PASS: formal analysis code has no archive input dependency\n")
