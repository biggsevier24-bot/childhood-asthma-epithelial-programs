options(stringsAsFactors = FALSE)
if (file.exists("config/use_local_r_libs.R")) source("config/use_local_r_libs.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(clusterProfiler)
  library(ggplot2)
  library(openxlsx)
})

out_dir <- "GSE18965_GO_reaudit"
series_dir <- "GSE18965_series_matrix_reproduction"
gpl_file <- "GPL96-57554.txt"

subdirs <- c(
  "00_logs", "01_current_annotation_all_DEGs",
  "02_current_annotation_upregulated", "03_current_annotation_downregulated",
  "04_legacy_GPL96_sensitivity", "05_targeted_term_audit",
  "06_comparison", "07_report"
)

safe_reset_dir <- function(path) {
  target <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(path) && !startsWith(target, paste0(root, "/"))) {
    stop("Refusing to delete output outside working directory: ", target)
  }
  if (file.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
safe_reset_dir(out_dir)
for (d in subdirs) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "00_logs", "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}
log_msg("Started GO reaudit.")

required <- c(
  file.path(series_dir, "04_limma_primary", "series_complete_limma_primary_default_eBayes_results.csv"),
  file.path(series_dir, "04_limma_primary", "series_DEG_flow_counts.csv"),
  file.path(series_dir, "04_limma_primary", "series_upregulated_FDR0.05_log2FC0.5.csv"),
  file.path(series_dir, "04_limma_primary", "series_downregulated_FDR0.05_log2FC_minus0.5.csv"),
  file.path(series_dir, "03_expression", "series_gene_level_expression.csv"),
  file.path(series_dir, "GSE18965_series_matrix_reproduction_report.md"),
  file.path(series_dir, "GSE18965_series_vs_CEL_comparison_report.md"),
  gpl_file
)
input_manifest <- data.frame(
  file = required,
  exists = file.exists(required),
  size_bytes = ifelse(file.exists(required), file.info(required)$size, NA_real_),
  md5 = ifelse(file.exists(required), tools::md5sum(required), NA_character_),
  stringsAsFactors = FALSE
)
write.csv(input_manifest, file.path(out_dir, "00_logs", "input_manifest.csv"), row.names = FALSE)
if (any(!input_manifest$exists)) stop("Missing required inputs: ", paste(input_manifest$file[!input_manifest$exists], collapse = ", "))

limma_all <- read.csv(required[1], check.names = FALSE)
gene_expr <- read.csv(required[5], check.names = FALSE)
background <- unique(limma_all$gene_symbol)
all_deg <- unique(limma_all$gene_symbol[limma_all$FDR < 0.05 & abs(limma_all$log2FC) >= 0.5])
up_deg <- unique(limma_all$gene_symbol[limma_all$FDR < 0.05 & limma_all$log2FC >= 0.5])
down_deg <- unique(limma_all$gene_symbol[limma_all$FDR < 0.05 & limma_all$log2FC <= -0.5])

stopifnot(length(background) == 12548, length(all_deg) == 772, length(up_deg) == 496, length(down_deg) == 276)

clusterProfiler_available <- requireNamespace("clusterProfiler", quietly = TRUE)
method_note <- if (clusterProfiler_available) {
  paste0(
    "Current-annotation main analysis used clusterProfiler::enrichGO ",
    "with OrgDb=org.Hs.eg.db, keyType=SYMBOL, ont=BP, all 12,548 tested genes ",
    "as universe, and BH multiple-testing correction."
  )
} else {
  "clusterProfiler was attempted but unavailable; fallback used org.Hs.eg.db, GO.db, AnnotationDbi, hypergeometric ORA, and BH correction."
}

get_go_desc <- function(go_ids) {
  con <- AnnotationDbi::dbconn(GO.db)
  if (!length(go_ids)) return(character())
  quoted <- paste(DBI::dbQuoteString(con, unique(go_ids)), collapse = ",")
  tab <- DBI::dbGetQuery(con, paste0("select go_id, term from go_term where go_id in (", quoted, ") and ontology = 'BP'"))
  out <- tab$term[match(go_ids, tab$go_id)]
  unname(out)
}

build_current_mapping <- function(symbols) {
  con <- AnnotationDbi::dbconn(org.Hs.eg.db)
  quoted <- paste(DBI::dbQuoteString(con, unique(symbols)), collapse = ",")
  sql <- paste0(
    "select distinct gene_info.symbol as gene_symbol, go_bp_all.go_id as GO_ID ",
    "from gene_info join go_bp_all on gene_info._id = go_bp_all._id ",
    "where gene_info.symbol in (", quoted, ")"
  )
  rows <- DBI::dbGetQuery(con, sql)
  rows <- unique(rows)
  rows$Description <- get_go_desc(rows$GO_ID)
  rows <- rows[!is.na(rows$Description), ]
  rows[, c("gene_symbol", "GO_ID", "Description")]
}

ora_from_mapping <- function(input_genes, background_genes, mapping, label) {
  input_genes <- intersect(unique(input_genes), background_genes)
  mapping_bg <- unique(mapping[mapping$gene_symbol %in% background_genes, c("gene_symbol", "GO_ID", "Description")])
  bg_genes <- intersect(unique(background_genes), unique(mapping_bg$gene_symbol))
  mapping_input <- mapping_bg[mapping_bg$gene_symbol %in% input_genes, , drop = FALSE]
  input_annotated <- intersect(input_genes, unique(mapping_input$gene_symbol))
  N <- length(bg_genes)
  n <- length(input_annotated)
  if (n == 0) {
    return(data.frame(
      analysis = label, GO_ID = character(), Description = character(),
      input_gene_count = integer(), background_gene_count = integer(),
      GeneRatio = character(), BgRatio = character(), enrichment_fold = numeric(),
      raw_P = numeric(), BH_FDR = numeric(), gene_symbols = character(),
      stringsAsFactors = FALSE
    ))
  }
  bg_counts <- as.data.frame(table(mapping_bg$GO_ID), stringsAsFactors = FALSE)
  names(bg_counts) <- c("GO_ID", "background_gene_count")
  input_counts <- as.data.frame(table(mapping_input$GO_ID), stringsAsFactors = FALSE)
  names(input_counts) <- c("GO_ID", "input_gene_count")
  res <- merge(input_counts, bg_counts, by = "GO_ID", all.x = TRUE, sort = FALSE)
  if (!nrow(res)) res <- data.frame(
    analysis = label, GO_ID = character(), Description = character(),
    input_gene_count = integer(), background_gene_count = integer(),
    GeneRatio = character(), BgRatio = character(), enrichment_fold = numeric(),
    raw_P = numeric(), gene_symbols = character(), stringsAsFactors = FALSE
  )
  term_desc <- unique(mapping_bg[, c("GO_ID", "Description")])
  term_genes <- split(mapping_input$gene_symbol, mapping_input$GO_ID)
  if (nrow(res)) {
    res$analysis <- label
    res$Description <- term_desc$Description[match(res$GO_ID, term_desc$GO_ID)]
    res$GeneRatio <- paste0(res$input_gene_count, "/", n)
    res$BgRatio <- paste0(res$background_gene_count, "/", N)
    res$enrichment_fold <- (res$input_gene_count / n) / (res$background_gene_count / N)
    res$raw_P <- phyper(res$input_gene_count - 1, res$background_gene_count, N - res$background_gene_count, n, lower.tail = FALSE)
    res$gene_symbols <- vapply(res$GO_ID, function(id) paste(sort(unique(term_genes[[id]])), collapse = "; "), character(1))
  }
  res$BH_FDR <- p.adjust(res$raw_P, method = "BH")
  res <- res[, c("analysis", "GO_ID", "Description", "input_gene_count", "background_gene_count", "GeneRatio", "BgRatio", "enrichment_fold", "raw_P", "BH_FDR", "gene_symbols")]
  res[order(res$BH_FDR, res$raw_P, -res$input_gene_count), ]
}

ratio_to_numeric <- function(x) {
  vapply(strsplit(x, "/", fixed = TRUE), function(parts) {
    if (length(parts) != 2) return(NA_real_)
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
}

run_clusterprofiler_go <- function(input_genes, background_genes, label) {
  input_genes <- intersect(unique(input_genes), background_genes)
  if (!length(input_genes)) {
    return(data.frame(
      analysis = label, GO_ID = character(), Description = character(),
      input_gene_count = integer(), background_gene_count = integer(),
      GeneRatio = character(), BgRatio = character(), enrichment_fold = numeric(),
      raw_P = numeric(), BH_FDR = numeric(), gene_symbols = character(),
      stringsAsFactors = FALSE
    ))
  }
  eg <- clusterProfiler::enrichGO(
    gene = input_genes,
    universe = unique(background_genes),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = FALSE
  )
  tab <- as.data.frame(eg)
  if (!nrow(tab)) {
    return(data.frame(
      analysis = label, GO_ID = character(), Description = character(),
      input_gene_count = integer(), background_gene_count = integer(),
      GeneRatio = character(), BgRatio = character(), enrichment_fold = numeric(),
      raw_P = numeric(), BH_FDR = numeric(), gene_symbols = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(
    analysis = label,
    GO_ID = tab$ID,
    Description = tab$Description,
    input_gene_count = tab$Count,
    background_gene_count = as.integer(vapply(strsplit(tab$BgRatio, "/", fixed = TRUE), `[`, character(1), 1)),
    GeneRatio = tab$GeneRatio,
    BgRatio = tab$BgRatio,
    enrichment_fold = ratio_to_numeric(tab$GeneRatio) / ratio_to_numeric(tab$BgRatio),
    raw_P = tab$pvalue,
    BH_FDR = tab$p.adjust,
    gene_symbols = gsub("/", "; ", tab$geneID, fixed = TRUE),
    stringsAsFactors = FALSE
  )
  out[order(out$BH_FDR, out$raw_P, -out$input_gene_count), ]
}

parse_gpl_go <- function(x) {
  if (is.na(x) || !nzchar(x)) return(data.frame(GO_ID = character(), Description = character(), stringsAsFactors = FALSE))
  parts <- unlist(strsplit(x, " /// ", fixed = TRUE))
  out <- do.call(rbind, lapply(parts, function(part) {
    fields <- strsplit(part, " // ", fixed = TRUE)[[1]]
    if (length(fields) < 2) return(NULL)
    data.frame(GO_ID = paste0("GO:", fields[1]), Description = fields[2], stringsAsFactors = FALSE)
  }))
  if (is.null(out)) return(data.frame(GO_ID = character(), Description = character(), stringsAsFactors = FALSE))
  unique(out)
}

build_legacy_mapping <- function() {
  gpl <- read.delim(gpl_file, sep = "\t", comment.char = "#", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
  ann <- data.frame(
    probe_id = gpl$ID,
    gene_symbol = trimws(gpl$`Gene Symbol`),
    go_bp_raw = gpl$`Gene Ontology Biological Process`,
    stringsAsFactors = FALSE
  )
  ann <- ann[!is.na(ann$gene_symbol) & ann$gene_symbol != "" & !grepl("///", ann$gene_symbol, fixed = TRUE), ]
  rows <- do.call(rbind, lapply(seq_len(nrow(ann)), function(i) {
    pgo <- parse_gpl_go(ann$go_bp_raw[i])
    if (!nrow(pgo)) return(NULL)
    data.frame(gene_symbol = ann$gene_symbol[i], pgo, stringsAsFactors = FALSE)
  }))
  rows <- unique(rows)
  rows[rows$gene_symbol %in% background, ]
}

log_msg("Building current org.Hs.eg.db/GO.db BP mapping.")
current_mapping <- build_current_mapping(background)
log_msg("Current mapping built: ", nrow(current_mapping), " gene-term rows.")
saveRDS(current_mapping, file.path(out_dir, "00_logs", "current_orgHs_GO_BP_gene_to_term_mapping.rds"), compress = FALSE)
write.csv(data.frame(
  mapping = "current_orgHs_GOdb",
  rows = nrow(current_mapping),
  genes_mapped = length(unique(current_mapping$gene_symbol)),
  GO_terms = length(unique(current_mapping$GO_ID))
), file.path(out_dir, "00_logs", "current_orgHs_GO_BP_gene_to_term_mapping_summary.csv"), row.names = FALSE)
log_msg("Building legacy GPL96 BP mapping.")
legacy_mapping <- build_legacy_mapping()
log_msg("Legacy mapping built: ", nrow(legacy_mapping), " gene-term rows.")
saveRDS(legacy_mapping, file.path(out_dir, "04_legacy_GPL96_sensitivity", "legacy_GPL96_GO_BP_gene_to_term_mapping.rds"), compress = FALSE)
write.csv(data.frame(
  mapping = "legacy_GPL96",
  rows = nrow(legacy_mapping),
  genes_mapped = length(unique(legacy_mapping$gene_symbol)),
  GO_terms = length(unique(legacy_mapping$GO_ID))
), file.path(out_dir, "04_legacy_GPL96_sensitivity", "legacy_GPL96_GO_BP_gene_to_term_mapping_summary.csv"), row.names = FALSE)
log_msg("Running current annotation ORA.")
if (clusterProfiler_available) {
  log_msg("Using clusterProfiler::enrichGO for current-annotation main analysis.")
  cur_all <- run_clusterprofiler_go(all_deg, background, "current_clusterProfiler_all_772_DEGs")
  cur_up <- run_clusterprofiler_go(up_deg, background, "current_clusterProfiler_upregulated_496")
  cur_down <- run_clusterprofiler_go(down_deg, background, "current_clusterProfiler_downregulated_276")
} else {
  cur_all <- ora_from_mapping(all_deg, background, current_mapping, "current_all_772_DEGs")
  cur_up <- ora_from_mapping(up_deg, background, current_mapping, "current_upregulated_496")
  cur_down <- ora_from_mapping(down_deg, background, current_mapping, "current_downregulated_276")
}
log_msg("Writing current annotation ORA outputs.")
write.csv(cur_all, file.path(out_dir, "01_current_annotation_all_DEGs", "GO_all_DEGs_complete_current_annotation.csv"), row.names = FALSE)
write.csv(cur_up, file.path(out_dir, "02_current_annotation_upregulated", "GO_upregulated_complete_current_annotation.csv"), row.names = FALSE)
write.csv(cur_down, file.path(out_dir, "03_current_annotation_downregulated", "GO_downregulated_complete_current_annotation.csv"), row.names = FALSE)
write.csv(head(cur_all, 20), file.path(out_dir, "01_current_annotation_all_DEGs", "GO_all_DEGs_top20_current_annotation.csv"), row.names = FALSE)
write.csv(head(cur_up, 20), file.path(out_dir, "02_current_annotation_upregulated", "GO_upregulated_top20_current_annotation.csv"), row.names = FALSE)
write.csv(head(cur_down, 20), file.path(out_dir, "03_current_annotation_downregulated", "GO_downregulated_top20_current_annotation.csv"), row.names = FALSE)

log_msg("Running legacy GPL96 ORA.")
leg_all <- ora_from_mapping(all_deg, background, legacy_mapping, "legacy_GPL96_all_772_DEGs")
leg_up <- ora_from_mapping(up_deg, background, legacy_mapping, "legacy_GPL96_upregulated_496")
leg_down <- ora_from_mapping(down_deg, background, legacy_mapping, "legacy_GPL96_downregulated_276")
log_msg("Writing legacy GPL96 ORA outputs.")
write.csv(leg_all, file.path(out_dir, "04_legacy_GPL96_sensitivity", "GO_legacy_GPL96_all_DEGs.csv"), row.names = FALSE)
write.csv(leg_up, file.path(out_dir, "04_legacy_GPL96_sensitivity", "GO_legacy_GPL96_upregulated.csv"), row.names = FALSE)
write.csv(leg_down, file.path(out_dir, "04_legacy_GPL96_sensitivity", "GO_legacy_GPL96_downregulated.csv"), row.names = FALSE)

target_terms <- data.frame(
  target_label = c(
    "extracellular matrix organization", "collagen fibril organization",
    "wound healing", "cell adhesion", "focal adhesion assembly",
    "extracellular matrix disassembly", "collagen catabolic process",
    "tissue remodeling"
  ),
  target_GO_ID = c(
    "GO:0030198", "GO:0030199", "GO:0042060", "GO:0007155",
    "GO:0048041", "GO:0022617", "GO:0030574", "GO:0048771"
  ),
  stringsAsFactors = FALSE
)

extract_targets <- function(res, method, direction) {
  if (!nrow(res)) {
    return(data.frame(
      method = method, direction = direction, target_label = target_terms$target_label,
      target_GO_ID = target_terms$target_GO_ID,
      matched = FALSE, GO_ID = NA_character_, Description = NA_character_,
      input_gene_count = NA_integer_, background_gene_count = NA_integer_,
      GeneRatio = NA_character_, BgRatio = NA_character_, enrichment_fold = NA_real_,
      raw_P = NA_real_, BH_FDR = NA_real_, gene_symbols = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(seq_len(nrow(target_terms)), function(i) {
    term <- target_terms$target_label[i]
    target_id <- target_terms$target_GO_ID[i]
    idx <- which(res$GO_ID == target_id)
    if (!length(idx)) idx <- which(tolower(res$Description) == tolower(term))
    if (!length(idx)) {
      data.frame(method = method, direction = direction, target_label = term, matched = FALSE,
                 target_GO_ID = target_id,
                 GO_ID = NA_character_, Description = NA_character_,
                 input_gene_count = NA_integer_, background_gene_count = NA_integer_,
                 GeneRatio = NA_character_, BgRatio = NA_character_, enrichment_fold = NA_real_,
                 raw_P = NA_real_, BH_FDR = NA_real_, gene_symbols = NA_character_,
                 stringsAsFactors = FALSE)
    } else {
      rows <- res[idx, ]
      data.frame(method = method, direction = direction, target_label = term, matched = TRUE,
                 target_GO_ID = target_id,
                 rows[, c("GO_ID", "Description", "input_gene_count", "background_gene_count", "GeneRatio", "BgRatio", "enrichment_fold", "raw_P", "BH_FDR", "gene_symbols")],
                 stringsAsFactors = FALSE)
    }
  }))
}

target_audit <- rbind(
  extract_targets(cur_all, if (clusterProfiler_available) "current_clusterProfiler_orgHs_GOdb" else "current_orgHs_GOdb", "all_DEGs"),
  extract_targets(cur_up, if (clusterProfiler_available) "current_clusterProfiler_orgHs_GOdb" else "current_orgHs_GOdb", "upregulated"),
  extract_targets(cur_down, if (clusterProfiler_available) "current_clusterProfiler_orgHs_GOdb" else "current_orgHs_GOdb", "downregulated"),
  extract_targets(leg_all, "legacy_GPL96", "all_DEGs"),
  extract_targets(leg_up, "legacy_GPL96", "upregulated"),
  extract_targets(leg_down, "legacy_GPL96", "downregulated")
)
write.csv(target_audit, file.path(out_dir, "05_targeted_term_audit", "GO_targeted_ECM_terms_audit.csv"), row.names = FALSE)

sig_count <- function(x) sum(x$BH_FDR < 0.05, na.rm = TRUE)
top_overlap <- function(a, b) length(intersect(head(a$GO_ID, 20), head(b$GO_ID, 20)))
comparison <- data.frame(
  analysis = c("all_DEGs", "upregulated", "downregulated"),
  current_GO_terms = c(nrow(cur_all), nrow(cur_up), nrow(cur_down)),
  legacy_GO_terms = c(nrow(leg_all), nrow(leg_up), nrow(leg_down)),
  current_FDR_lt_0_05 = c(sig_count(cur_all), sig_count(cur_up), sig_count(cur_down)),
  legacy_FDR_lt_0_05 = c(sig_count(leg_all), sig_count(leg_up), sig_count(leg_down)),
  top20_GO_ID_overlap = c(top_overlap(cur_all, leg_all), top_overlap(cur_up, leg_up), top_overlap(cur_down, leg_down)),
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(out_dir, "06_comparison", "GO_annotation_method_comparison.csv"), row.names = FALSE)

target_summary <- subset(target_audit, grepl("^current_", method))
fdr_for <- function(term, direction) {
  x <- target_summary[tolower(target_summary$target_label) == tolower(term) & target_summary$direction == direction, ]
  if (!nrow(x) || all(is.na(x$BH_FDR))) NA_real_ else min(x$BH_FDR, na.rm = TRUE)
}
fmt_fdr <- function(term) paste(signif(c(fdr_for(term, "all_DEGs"), fdr_for(term, "upregulated"), fdr_for(term, "downregulated")), 4), collapse = " / ")
fig_repro <- "No"
fig_action <- "Redraw or remove Figure 1B unless original analysis provenance can justify it; current ORA should replace it if GO is retained."
fig_implication <- "Figure 1B cannot be retained as showing high FDR-significant ECM/wound-healing enrichment unless the exact original input and multiple-testing output are supplied."

summary_lines <- c(
  paste0("- Tested background genes: ", length(background)),
  paste0("- All DEGs: ", length(all_deg)),
  paste0("- Upregulated DEGs: ", length(up_deg)),
  paste0("- Downregulated DEGs: ", length(down_deg)),
  paste0("- Significant GO BP terms for all DEGs at FDR<0.05: ", sig_count(cur_all)),
  paste0("- Significant GO BP terms for upregulated genes: ", sig_count(cur_up)),
  paste0("- Significant GO BP terms for downregulated genes: ", sig_count(cur_down)),
  paste0("- ECM organization FDR (all/up/down): ", fmt_fdr("extracellular matrix organization")),
  paste0("- Collagen fibril organization FDR (all/up/down): ", fmt_fdr("collagen fibril organization")),
  paste0("- Wound healing FDR (all/up/down): ", fmt_fdr("wound healing")),
  paste0("- Cell adhesion FDR (all/up/down): ", fmt_fdr("cell adhesion")),
  paste0("- Original Figure 1B reproducible: ", fig_repro),
  paste0("- Recommended action for Figure 1B: ", fig_action),
  paste0("- Critical implication for repair-ECM derivation: ", fig_implication)
)
writeLines(summary_lines, file.path(out_dir, "00_logs", "terminal_summary.txt"))

report <- c(
  "# GSE18965 GO Reaudit Report",
  "",
  "## Fixed Summary",
  summary_lines,
  "",
  "## Scope",
  "This third-stage audit only evaluates GO Biological Process overrepresentation and whether the original Figure 1B can be reproduced. No manuscript, figure, response-letter, threshold, group, or gene-list edits were made.",
  "",
  "## Main Method",
  paste0("- ", method_note),
  "- Main input for Methods reproduction: all 772 DEGs at FDR<0.05 and |log2FC|>=0.5.",
  "- Directional supplementary inputs: 496 upregulated and 276 downregulated DEGs.",
  "- Background: all 12,548 genes entering the Series Matrix gene-level limma analysis.",
  "- Current annotation source: org.Hs.eg.db + GO.db through clusterProfiler for the main analysis; explicit org.Hs.eg.db/GO.db mappings are also archived for auditability.",
  "- Test: hypergeometric overrepresentation test through clusterProfiler::enrichGO, followed by BH FDR across all tested GO BP terms for each analysis.",
  "- Legacy GPL96 annotation is reported only as sensitivity analysis.",
  "",
  "## Figure 1B Reproducibility",
  paste0("- Original Figure 1B reproducible: ", fig_repro),
  paste0("- Recommended action: ", fig_action),
  "Nominal P values were not used as FDR. Target ECM/wound-healing/cell-adhesion terms must be judged by BH-FDR.",
  "",
  "## Targeted Terms",
  "All requested target terms are reported in `GO_targeted_ECM_terms_audit.xlsx` with raw P and BH-FDR for all/up/down inputs under current and legacy annotations.",
  "Targeted-term matching uses fixed GO IDs; `matched=FALSE` and NA statistics mean that the exact requested GO term was not returned in that analysis, and similar child/parent labels were not substituted.",
  "",
  "## Legacy Sensitivity",
  "The 2014 GPL96 built-in GO annotations were parsed as a legacy sensitivity only. These results are not the primary GO analysis and should not be used to preserve Figure 1B if current annotation does not support it.",
  "",
  "## Files",
  "- Complete current annotation results: `01_current_annotation_all_DEGs`, `02_current_annotation_upregulated`, `03_current_annotation_downregulated`.",
  "- Legacy GPL96 sensitivity: `04_legacy_GPL96_sensitivity`.",
  "- Method comparison: `06_comparison/GO_annotation_method_comparison.csv`.",
  "- Consolidated workbook: `GSE18965_GO_reaudit_all_results.xlsx`."
)
writeLines(report, file.path(out_dir, "07_report", "GSE18965_GO_reaudit_report.md"))

changelog <- c(
  "# CHANGELOG_GO_reaudit",
  "",
  "- Reaudited GO Biological Process ORA using all 772 Series Matrix DEGs as the primary Methods-reproduction input.",
  "- Added directional supplementary ORA for 496 upregulated and 276 downregulated DEGs.",
  "- Installed and used Bioconductor annotation packages org.Hs.eg.db 3.16.0 and GO.db 3.16.0.",
  "- Installed clusterProfiler 4.6.2 after updating htmltools and installing ggiraph/shadowtext dependencies; current-annotation main analysis now uses clusterProfiler::enrichGO.",
  "- Used all 12,548 limma-tested genes as background for every current-annotation analysis.",
  "- Kept GPL96 built-in GO annotations only as legacy sensitivity analysis.",
  "- Reported raw P and BH-FDR for complete GO results; did not substitute nominal P for FDR.",
  "- Did not modify manuscript, Figure 1, response letter, thresholds, groups, gene list, or probe rules."
)
writeLines(changelog, file.path(out_dir, "CHANGELOG_GO_reaudit.md"))

write.csv(data.frame(package = c("AnnotationDbi", "org.Hs.eg.db", "GO.db", "clusterProfiler"),
                     available = c(TRUE, TRUE, TRUE, clusterProfiler_available),
                     version = c(as.character(packageVersion("AnnotationDbi")),
                                 as.character(packageVersion("org.Hs.eg.db")),
                                 as.character(packageVersion("GO.db")),
                                 if (clusterProfiler_available) as.character(packageVersion("clusterProfiler")) else NA_character_)),
          file.path(out_dir, "00_logs", "annotation_package_versions.csv"), row.names = FALSE)

add_sheet <- function(wb, name, df) {
  addWorksheet(wb, substr(name, 1, 31))
  writeData(wb, substr(name, 1, 31), df)
  freezePane(wb, substr(name, 1, 31), firstRow = TRUE)
}

target_wb <- createWorkbook()
add_sheet(target_wb, "targeted_ECM_terms", target_audit)
saveWorkbook(target_wb, file.path(out_dir, "05_targeted_term_audit", "GO_targeted_ECM_terms_audit.xlsx"), overwrite = TRUE)

comparison_wb <- createWorkbook()
add_sheet(comparison_wb, "annotation_method_summary", comparison)
add_sheet(comparison_wb, "current_all_top20", head(cur_all, 20))
add_sheet(comparison_wb, "legacy_all_top20", head(leg_all, 20))
add_sheet(comparison_wb, "targeted_terms", target_audit)
saveWorkbook(comparison_wb, file.path(out_dir, "06_comparison", "GO_annotation_method_comparison.xlsx"), overwrite = TRUE)

all_wb <- createWorkbook()
add_sheet(all_wb, "README_summary", data.frame(item = sub("^- ", "", summary_lines), stringsAsFactors = FALSE))
add_sheet(all_wb, "GO_all_current", cur_all)
add_sheet(all_wb, "GO_up_current", cur_up)
add_sheet(all_wb, "GO_down_current", cur_down)
add_sheet(all_wb, "GO_legacy_all", leg_all)
add_sheet(all_wb, "GO_legacy_up", leg_up)
add_sheet(all_wb, "GO_legacy_down", leg_down)
add_sheet(all_wb, "targeted_ECM_terms", target_audit)
add_sheet(all_wb, "method_comparison", comparison)
add_sheet(all_wb, "package_versions", read.csv(file.path(out_dir, "00_logs", "annotation_package_versions.csv")))
saveWorkbook(all_wb, file.path(out_dir, "07_report", "GSE18965_GO_reaudit_all_results.xlsx"), overwrite = TRUE)

file.copy(file.path(out_dir, "07_report", "GSE18965_GO_reaudit_report.md"),
          file.path(out_dir, "GSE18965_GO_reaudit_report.md"), overwrite = TRUE)
file.copy(file.path(out_dir, "07_report", "GSE18965_GO_reaudit_all_results.xlsx"),
          file.path(out_dir, "GSE18965_GO_reaudit_all_results.xlsx"), overwrite = TRUE)
file.copy(file.path(out_dir, "05_targeted_term_audit", "GO_targeted_ECM_terms_audit.xlsx"),
          file.path(out_dir, "GO_targeted_ECM_terms_audit.xlsx"), overwrite = TRUE)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
log_msg("Completed GO reaudit.")

