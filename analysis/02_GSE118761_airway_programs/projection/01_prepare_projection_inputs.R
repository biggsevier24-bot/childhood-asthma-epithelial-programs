source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "00_projection_common.R"))
projection_log("START 01_prepare_projection_inputs.R")

paths <- list(
  gse118_expression = projection_input_path("gse118_expression", "GSE118761_EXPR"),
  gse118_metadata = projection_input_path("gse118_metadata", "GSE118761_META"),
  gse152_scores = projection_input_path("gse152_scores", "GSE152004_SCORES"),
  gse152_assignments = projection_input_path("gse152_assignments", "GSE152004_ASSIGNMENTS"),
  t2_gene_set = projection_input_path("t2_gene_set", ""),
  ifn_gene_set = projection_input_path("ifn_gene_set", ""),
  repair_gene_set = projection_input_path("repair_gene_set", "")
)

expr118 <- projection_read_expression(paths$gse118_expression)
meta118 <- read.csv(paths$gse118_metadata, check.names = FALSE, stringsAsFactors = FALSE)
scores152 <- read.csv(paths$gse152_scores, check.names = FALSE, stringsAsFactors = FALSE)
assign152 <- read.csv(paths$gse152_assignments, check.names = FALSE, stringsAsFactors = FALSE)

required_meta <- c("sample_id", "tissue", "atopy", "asthma", "wheeze")
if (!all(required_meta %in% names(meta118))) stop("GSE118761 metadata is missing required fields")
if (anyDuplicated(meta118$sample_id)) stop("Duplicate GSE118761 sample_id")
if (!setequal(colnames(expr118), meta118$sample_id)) stop("GSE118761 expression and metadata IDs differ")
meta118 <- meta118[match(colnames(expr118), meta118$sample_id), , drop = FALSE]
meta118$tissue <- tolower(meta118$tissue)

expected <- projection_config$expected_sample_counts
stopifnot(nrow(meta118) == as.integer(expected$total))
stopifnot(sum(meta118$tissue == "nasal") == as.integer(expected$nasal))
stopifnot(sum(meta118$tissue == "tracheal") == as.integer(expected$tracheal))
stopifnot(!anyDuplicated(scores152$sample_id), !anyDuplicated(assign152$sample_id))
stopifnot(setequal(scores152$sample_id, assign152$sample_id), nrow(scores152) == 695L)

saveRDS(list(expr118 = expr118, meta118 = meta118, scores152 = scores152,
             assign152 = assign152, paths = paths), file.path(projection_out, "projection_inputs.rds"))

sample_audit <- data.frame(
  sample_id = meta118$sample_id,
  dataset = "GSE118761",
  tissue = meta118$tissue,
  included = TRUE,
  exclusion_reason = NA_character_,
  atopy = meta118$atopy,
  asthma = meta118$asthma,
  wheeze = meta118$wheeze,
  age = if ("age" %in% names(meta118)) meta118$age else NA,
  sex = if ("sex" %in% names(meta118)) meta118$sex else NA_character_
)
projection_write_tsv(sample_audit, file.path(projection_out, "projection_sample_audit.tsv"))

dimensions <- list(
  gse118_expression = c(nrow(expr118), ncol(expr118)),
  gse118_metadata = c(nrow(meta118), ncol(meta118)),
  gse152_scores = c(nrow(scores152), ncol(scores152)),
  gse152_assignments = c(nrow(assign152), ncol(assign152)),
  t2_gene_set = c(length(projection_read_gene_set(paths$t2_gene_set)), 1),
  ifn_gene_set = c(length(projection_read_gene_set(paths$ifn_gene_set)), 1),
  repair_gene_set = c(length(projection_read_gene_set(paths$repair_gene_set)), 1)
)
created_from <- c(
  gse118_expression = "GSE118761 GEO gene counts; nonzero >=20%; DESeq2 blind VST; gene symbols collapsed by rowsum",
  gse118_metadata = "GSE118761 deposited Series Matrix metadata",
  gse152_scores = "Formal GSE152004 VST program-score pipeline",
  gse152_assignments = "Formal GSE152004 k=3 pipeline; seed 20260626; nstart 500",
  t2_gene_set = "Locked manuscript T2 three-gene definition",
  ifn_gene_set = "Locked manuscript IFN definition",
  repair_gene_set = "Locked fixed 32-gene repair-ECM definition"
)
used_by <- c(
  gse118_expression = "02_calculate_projection_features.R",
  gse118_metadata = "01_prepare_projection_inputs.R;04_assign_nearest_centroid.R;06_summarize_projection_phenotypes.R",
  gse152_scores = "03_build_reference_centroids.R",
  gse152_assignments = "03_build_reference_centroids.R",
  t2_gene_set = "02_calculate_projection_features.R",
  ifn_gene_set = "02_calculate_projection_features.R",
  repair_gene_set = "02_calculate_projection_features.R"
)
manifest <- do.call(rbind, lapply(names(paths), function(nm) data.frame(
  file = basename(paths[[nm]]),
  relative_path = projection_config$inputs[[nm]],
  description = created_from[[nm]],
  nrow = dimensions[[nm]][1],
  ncol = dimensions[[nm]][2],
  sha256 = projection_sha256(paths[[nm]]),
  created_from = created_from[[nm]],
  used_by_script = used_by[[nm]]
)))
projection_write_tsv(manifest, file.path(projection_out, "provenance/input_manifest.tsv"))
for (i in seq_len(nrow(manifest))) projection_log("INPUT", manifest$relative_path[i], manifest$sha256[i])
projection_log("SAMPLES total=104 nasal=55 tracheal=49 duplicates=0")
projection_log("END 01_prepare_projection_inputs.R")
