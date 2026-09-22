hospital_model_formulas <- list(
  asthma = asthma ~ clinical_T2_score + age + sex,
  frequent_exacerbations = frequent_exacerbations ~ clinical_T2_score + age + sex,
  wheeze = wheeze ~ clinical_T2_score + age + sex,
  exacerbation_count = exacerbation_count ~ clinical_T2_score + age + sex
)
cat("Model specifications loaded. Protected patient-level input is required for estimation.\n")
