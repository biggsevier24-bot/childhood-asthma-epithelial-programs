codex_r_helper <- Sys.getenv("CODEX_R_LIB_HELPER", "")
if (nzchar(codex_r_helper) && file.exists(codex_r_helper)) source(codex_r_helper)
options(stringsAsFactors = FALSE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
src <- read.delim(file.path(repo, "figures/Figure5/Figure5f_source_data.tsv"), check.names = FALSE)
src$state <- factor(src$state, levels = c("E1", "E2", "E3"))
src$phenotype <- factor(src$phenotype, levels = c("atopy", "asthma"))

draw_panel <- function() {
  old <- par(mar = c(4.5, 4.5, 1.5, 0.8), las = 1)
  on.exit(par(old))
  m <- xtabs(percentage ~ phenotype + state, src)
  cols <- c("#D97841", "#477A9E")
  bp <- barplot(m, beside = TRUE, col = cols, border = NA, ylim = c(0, 100),
                ylab = "Participants (%)", xlab = "Projected state", legend.text = c("Atopy", "Asthma"),
                args.legend = list(x = "topright", bty = "n"))
  text(bp, m + 4, labels = sprintf("%.1f", m), cex = 0.8)
  box(bty = "l")
}

png(file.path(repo, "figures/Figure5/Figure5f_rebuilt.png"), width = 1800, height = 1300, res = 300)
draw_panel()
dev.off()
tiff(file.path(repo, "figures/Figure5/Figure5f_rebuilt.tiff"), width = 6, height = 4.33, units = "in", res = 300, compression = "lzw")
draw_panel()
dev.off()
