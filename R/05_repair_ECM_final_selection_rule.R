# Frozen revision-era calibrated rule. It uses only biological, expression,
# and co-expression features; it does not load a fitted model or gene list.
select_repair_ECM_gene <- function(feature_table) {
  required <- c(
    "AA_minus_HN", "mean_signed_correlation_to_selected_cluster",
    "aa_positive_integrin_family_absolute_median", "signed_degree_040",
    "aa_positive_go_collagen_signed_degree_060", "family_collagen",
    "aa_positive_matricellular_family_signed_mean", "signed_degree_060",
    "aa_positive_mmp_family_signed_max", "aa_positive_all_signed_median",
    "signed_degree_030"
  )
  missing <- setdiff(required, names(feature_table))
  if (length(missing)) stop("Missing frozen-rule features: ", paste(missing, collapse = ", "))
  with(feature_table, ifelse(
    AA_minus_HN < 0.125,
    ifelse(
      mean_signed_correlation_to_selected_cluster >= -0.50,
      ifelse(
        aa_positive_integrin_family_absolute_median >= 0.25,
        signed_degree_040 < 0.005,
        TRUE
      ),
      TRUE
    ),
    ifelse(
      aa_positive_go_collagen_signed_degree_060 < 0.766194332,
      ifelse(
        !family_collagen,
        ifelse(
          aa_positive_matricellular_family_signed_mean < 0.69,
          AA_minus_HN < 0.14,
          AA_minus_HN >= 0.19
        ),
        ifelse(
          signed_degree_060 >= 0.77,
          FALSE,
          aa_positive_mmp_family_signed_max < 0.939526755
        )
      ),
      ifelse(
        aa_positive_all_signed_median >= 0.838103855,
        mean_signed_correlation_to_selected_cluster >= 0.686,
        signed_degree_030 >= 0.90
      )
    )
  ))
}
