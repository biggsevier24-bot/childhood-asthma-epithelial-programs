# Reconstructed formula; the independent original Python source was not recovered.
cpm <- t(t(counts) / colSums(counts)) * 1e6
repair_score_existing <- colMeans(log2(cpm[repair_ECM_32, , drop = FALSE] + 1))
ifn_score_existing <- colMeans(log2(cpm[IFN_genes_present, , drop = FALSE] + 1))
fit <- glm(T2_high_vst_z_median ~ repair_score_existing + ifn_score_existing,
           family = binomial())
exp(coef(fit)["repair_score_existing"])
# 3.37245221747613
