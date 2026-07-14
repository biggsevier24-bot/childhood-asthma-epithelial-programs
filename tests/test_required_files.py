#!/usr/bin/env python3
"""Required-file check for the public repository layout."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
required = [
    "README.md",
    "LICENSE",
    "CITATION.cff",
    ".zenodo.json",
    "CHANGELOG.md",
    "config/analysis_config.yml",
    "config/paths_example.yml",
    "R/01_GSE18965_preprocessing.R",
    "R/02_GSE18965_DEG_and_ECM_scores.R",
    "R/03_RNAseq_preprocessing.R",
    "R/04_gene_set_scoring.R",
    "R/05_sensitivity_analyses.R",
    "R/06_tables_and_figures.R",
    "python/01_PCA.py",
    "python/02_clustering.py",
    "python/03_cross_cohort_projection.py",
    "python/04_statistical_models.py",
    "python/05_hospital_models.py",
    "python/06_output_checks.py",
    "gene_sets/repair_ECM_fixed32.txt",
    "gene_sets/IFN_fixed50.txt",
    "gene_sets/T2_fixed3.txt",
    "gene_sets/Reactome_ECM_manifest.csv",
    "gene_sets/NABA_core_matrisome_manifest.csv",
    "data_manifest/public_data_manifest.csv",
    "data_manifest/required_intermediate_inputs.csv",
    "data_manifest/public_data_download_instructions.md",
    "data_manifest/hospital_data_dictionary.csv",
    "data_manifest/hospital_data_synthetic_example.csv",
    "data_manifest/private_data_not_shared.md",
    "environment/R_sessionInfo.txt",
    "environment/installed_R_packages.csv",
    "environment/requirements.txt",
    "environment/environment.yml",
    "docs/analysis_workflow.md",
    "docs/analysis_decisions.md",
    "docs/figure_table_code_map.csv",
    "docs/phenotype_definitions.md",
    "docs/reproducibility_notes.md",
    "docs/source_code_availability.md",
    "docs/public_repository_final_audit.md",
    "outputs/README.md",
    "outputs/manuscript_locked/README.md",
    "outputs/audit_reruns/README.md",
    "outputs/manuscript_locked/machine_readable_tables/manuscript_locked_reference_summary.csv",
    "run_all.sh",
]
missing = [p for p in required if not (ROOT / p).exists()]
if missing:
    raise SystemExit("Missing required files: " + ", ".join(missing))
print("Required-file check passed.")
