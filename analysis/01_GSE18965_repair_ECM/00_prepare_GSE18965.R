read_gse18965_series <- function(series_file) {
  lines <- readLines(gzfile(series_file, "rt"), warn = FALSE)
  extract_tag <- function(tag) {
    hit <- grep(paste0("^", tag, "\\t"), lines, value = TRUE)
    if (length(hit) != 1L) stop("Could not uniquely extract ", tag)
    gsub('^"|"$', "", strsplit(hit, "\t", fixed = TRUE)[[1]][-1])
  }
  table_begin <- grep("^!series_matrix_table_begin", lines)
  table_end <- grep("^!series_matrix_table_end", lines)
  stopifnot(length(table_begin) == 1L, length(table_end) == 1L, table_end > table_begin)
  metadata <- data.frame(
    gsm = extract_tag("!Sample_geo_accession"),
    title = extract_tag("!Sample_title"),
    platform_id = extract_tag("!Sample_platform_id"),
    stringsAsFactors = FALSE
  )
  metadata$group <- ifelse(grepl("^AA_", metadata$title), "AA",
                           ifelse(grepl("^HN_", metadata$title), "HN", NA_character_))
  stopifnot(nrow(metadata) == 16L, sum(metadata$group == "AA") == 9L,
            sum(metadata$group == "HN") == 7L, !anyNA(metadata$group),
            all(metadata$platform_id == "GPL96"), !anyDuplicated(metadata$gsm))
  series_table <- read.delim(
    text = paste(lines[(table_begin + 1L):(table_end - 1L)], collapse = "\n"),
    check.names = FALSE, quote = "\""
  )
  probe_expression <- as.matrix(series_table[, metadata$gsm, drop = FALSE])
  storage.mode(probe_expression) <- "double"
  rownames(probe_expression) <- series_table$ID_REF
  list(metadata = metadata, probe_expression = probe_expression)
}
