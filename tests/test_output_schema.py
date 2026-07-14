#!/usr/bin/env python3
"""Reference-output schema checks."""

from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]

checks = {
    "outputs/audit_reruns/machine_readable_tables/GSE18965_DEG_threshold_audit.csv": ["filter_definition", "total", "up", "down"],
    "outputs/audit_reruns/machine_readable_tables/GSE152004_final_state_profiles.csv": ["state_label", "n"],
    "outputs/audit_reruns/machine_readable_tables/hospital_analysis_population_audit.csv": ["full_cohort", "asthma_subgroup"],
    "outputs/audit_reruns/machine_readable_tables/final_consistency_tests.csv": ["test", "passed", "observed", "expected"],
    "outputs/manuscript_locked/machine_readable_tables/manuscript_locked_reference_summary.csv": ["item", "value", "source_note"],
}

for rel, cols in checks.items():
    path = ROOT / rel
    if not path.exists():
        raise SystemExit(f"Missing schema target: {rel}")
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in cols if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit(f"{rel} missing columns {missing}")

for fig in range(1, 7):
    if not (ROOT / f"outputs/manuscript_locked/figure_previews/Figure{fig}_final_locked.png").exists():
        raise SystemExit(f"Missing Figure{fig} preview")

print("Output-schema check passed.")
