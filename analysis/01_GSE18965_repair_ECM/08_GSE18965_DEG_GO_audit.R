# Independent audit; this script is not called by the 32-gene derivation.
options(stringsAsFactors = FALSE)
suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))
library_helper <- Sys.getenv("CODEX_R_LIB_HELPER", unset = "")
if (nzchar(library_helper) && file.exists(library_helper)) source(library_helper)
required <- c("limma", "clusterProfiler", "org.Hs.eg.db")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required audit packages: ", paste(missing, collapse = ", "))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else "R/08_GSE18965_DEG_GO_audit.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
expr <- as.matrix(read.delim(file.path(repo_root, "outputs/GSE18965_gene_expression.tsv"), row.names = 1, check.names = FALSE))
metadata <- read.delim(file.path(repo_root, "outputs/sample_metadata.tsv"), check.names = FALSE)
group <- factor(metadata$group, levels = c("HN", "AA"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- limma::lmFit(expr[, metadata$gsm, drop = FALSE], design)
fit <- limma::contrasts.fit(fit, limma::makeContrasts(AA_minus_HN = AA - HN, levels = design))
fit <- limma::eBayes(fit)
result <- limma::topTable(fit, number = Inf, sort.by = "none")
result$gene <- rownames(result)
result <- result[, c("gene", setdiff(names(result), "gene"))]
result$prespecified_DEG <- result$adj.P.Val < 0.05 & abs(result$logFC) >= 0.5
audit_dir <- file.path(repo_root, "outputs/independent_DEG_GO_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
write.table(result, file.path(audit_dir, "limma_all_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
deg <- result$gene[result$prespecified_DEG]
go <- clusterProfiler::enrichGO(
  gene = deg, universe = result$gene, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
  keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", readable = TRUE
)
write.table(as.data.frame(go), file.path(audit_dir, "GO_BP_all_DEGs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Independent DEG/GO audit completed; it was not used for 32-gene selection.\n")
