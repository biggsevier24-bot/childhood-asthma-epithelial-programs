#!/usr/bin/env bash
set -euo pipefail

echo "childhood-asthma-epithelial-programs public rerun launcher"
echo "This repository does not claim to regenerate all manuscript results directly from raw GEO downloads."
echo "Default rerun output root: outputs/user_rerun/"
echo "Manuscript-locked outputs will not be overwritten."

mkdir -p outputs/user_rerun

if ! command -v Rscript >/dev/null 2>&1; then
  echo "FAIL: Rscript not found. Install R and required packages before rerunning analyses."
  exit 0
fi

check_required() {
  local label="$1"
  shift
  local missing=0
  echo "PREFLIGHT: ${label}"
  for required_file in "$@"; do
    if [ ! -f "$required_file" ]; then
      echo "MISSING: ${required_file}"
      missing=1
    fi
  done
  if [ "$missing" -ne 0 ]; then
    echo "SKIP: required intermediate input unavailable"
    return 1
  fi
  echo "PREFLIGHT OK: ${label}"
  return 0
}

run_r_step() {
  local label="$1"
  local script="$2"
  echo "START: ${label}"
  if Rscript "$script"; then
    echo "SUCCESS: ${label}"
  else
    echo "FAIL: ${label}"
    return 1
  fi
}

if check_required "Figure 2 / GSE18965 Series Matrix branch" \
  "data/public/GSE18965/GSE18965_series_matrix.txt.gz" \
  "data/public/GSE18965/GSE18965_family.soft.gz" \
  "data/public/GSE18965/GPL96-57554.txt" \
  "gene_sets/repair_ECM_fixed32.txt"; then
  run_r_step "GSE18965 preprocessing / Series Matrix branch" "R/01_GSE18965_preprocessing.R"
else
  echo "Download public GSE18965 files or provide the required intermediate files listed in data_manifest/required_intermediate_inputs.csv."
fi

if check_required "GSE18965 GO and repair audit" \
  "outputs/user_rerun/GSE18965_series_matrix_reproduction/04_limma_primary/series_complete_limma_primary_default_eBayes_results.csv" \
  "outputs/user_rerun/GSE18965_series_matrix_reproduction/03_expression/series_gene_level_expression.csv" \
  "data/public/GSE18965/GPL96-57554.txt"; then
  run_r_step "GSE18965 DEG, enrichment, and repair-ECM audit" "R/02_GSE18965_DEG_and_ECM_scores.R"
fi

if check_required "Figures 3-5 / RNA-seq standardized input generation" \
  "data/public/GSE118761/GSE118761_genecounts.csv.gz" \
  "data/public/GSE118761/GSE118761_series_matrix.txt.gz" \
  "data/public/GSE118761/GSE118761_IFN_T2_original32_PCA_cluster_scores.csv" \
  "data/public/GSE152004/GSE152004_695_raw_counts.txt.gz" \
  "data/public/GSE152004/GSE152004_series_matrix.txt.gz" \
  "data/public/GSE152004/GSE152004_metadata_with_T2_IFN_cluster.csv" \
  "gene_sets/repair_ECM_fixed32.txt" \
  "gene_sets/IFN_fixed50.txt" \
  "gene_sets/T2_fixed3.txt" \
  "gene_sets/Reactome_ECM_manifest.csv" \
  "gene_sets/NABA_core_matrisome_manifest.csv"; then
  run_r_step "RNA-seq preprocessing and metadata matching" "R/03_RNAseq_preprocessing.R"
  run_r_step "Gene-set scoring and external ECM multicohort analysis" "R/04_gene_set_scoring.R"
else
  echo "Provide the standardized score/metadata inputs listed in data_manifest/required_intermediate_inputs.csv before rerunning Figures 3-5 end-to-end."
fi

if check_required "External ECM corrected sensitivity" \
  "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_VST_gene_expression.csv" \
  "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_metadata.csv" \
  "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_VST_gene_expression.csv" \
  "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_metadata.csv" \
  "outputs/audit_reruns/random_gene_set_membership.csv.gz"; then
  run_r_step "External ECM sensitivity analyses" "R/05_sensitivity_analyses.R"
fi

if [ -n "${HOSPITAL_PRIVATE_DATA:-}" ]; then
  if check_required "Protected hospital models and final figure/table assembly" \
    "$HOSPITAL_PRIVATE_DATA" \
    "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_VST_gene_expression.csv" \
    "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE152004_metadata.csv" \
    "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_VST_gene_expression.csv" \
    "outputs/user_rerun/External_ECM_multicohort_analysis/02_gene_sets/standardized_inputs/GSE118761_metadata.csv"; then
    run_r_step "Tables, figures, and protected hospital models" "R/06_tables_and_figures.R"
  fi
else
  echo "SKIP: required intermediate input unavailable"
  echo "Protected hospital models require HOSPITAL_PRIVATE_DATA pointing to an authorized local de-identified file."
fi

if command -v python >/dev/null 2>&1; then
  echo "START: repository safety and layout checks"
  python tests/test_required_files.py
  python tests/test_output_schema.py
  python tests/test_no_private_information.py
  echo "SUCCESS: repository checks"
else
  echo "FAIL: python not found; repository checks were not run."
fi

echo "Run launcher finished."
