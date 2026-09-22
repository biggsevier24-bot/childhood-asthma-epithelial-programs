source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 05_summarize_projection_counts.R")

d <- projection_read_tsv(file.path(projection_out, "GSE118761_projection_assignments_FINAL.tsv"))
grid <- expand.grid(tissue = c("nasal", "tracheal"), state = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
obs <- as.data.frame(table(tissue = d$tissue, state = d$assigned_state), stringsAsFactors = FALSE)
names(obs)[3] <- "n"
counts <- merge(grid, obs, by = c("tissue", "state"), all.x = TRUE, sort = FALSE)
counts$n[is.na(counts$n)] <- 0L
counts <- counts[order(match(counts$tissue, c("nasal", "tracheal")), match(counts$state, c("E1", "E2", "E3"))), ]
counts$denominator <- ave(counts$n, counts$tissue, FUN = sum)
counts$percentage <- 100 * counts$n / counts$denominator
projection_write_tsv(counts, file.path(projection_out, "GSE118761_projection_counts_FINAL.tsv"))

stopifnot(sum(counts$n) == as.integer(projection_config$expected_sample_counts$total))
stopifnot(unique(counts$denominator[counts$tissue == "nasal"]) == as.integer(projection_config$expected_sample_counts$nasal))
stopifnot(unique(counts$denominator[counts$tissue == "tracheal"]) == as.integer(projection_config$expected_sample_counts$tracheal))
for (i in seq_len(nrow(counts))) {
  expected <- projection_expected_count(counts$tissue[i], counts$state[i])
  if (counts$n[i] != expected) stop("Projection count mismatch for ", counts$tissue[i], " ", counts$state[i])
}
projection_log("FINAL_COUNTS", paste(counts$tissue, counts$state, counts$n, sep = ":", collapse = ";"))
projection_log("END 05_summarize_projection_counts.R")
