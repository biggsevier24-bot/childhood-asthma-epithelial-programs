#!/usr/bin/env python3
"""Security scan for public release readiness."""

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {
    ".md", ".txt", ".csv", ".tsv", ".yml", ".yaml", ".json", ".cff",
    ".R", ".r", ".py", ".sh", ".gitignore"
}
ALLOWED_FILENAMES = {
    "hospital_data_dictionary.csv",
    "hospital_data_synthetic_example.csv",
    "private_data_not_shared.md",
    "test_no_private_information.py",
}

local_path = re.compile(r"((?<![A-Za-z])[A-Za-z]:[\\/]|/home/[^/\s]+|/Users/[^/\s]+)")
id_card = re.compile(r"\b\d{17}[\dXx]\b")
phone = re.compile(r"(?<!\d)(?:\+?86[- ]?)?1[3-9]\d{9}(?!\d)")
secret_assignment = re.compile(r"(?i)(password|api[_-]?key|secret|token)\s*[:=]\s*['\"]?[A-Za-z0-9_\-]{8,}")
raw_private_filename = re.compile(r"(?i)(hospital_raw|clinical_identifiable|medical_record|mrn|id_number)")
excel_temp = re.compile(r"~\$")

findings = []

for path in ROOT.rglob("*"):
    rel = path.relative_to(ROOT).as_posix()
    if path.is_dir():
        continue
    if rel == "tests/test_no_private_information.py":
        continue
    if any(part in {".git"} for part in path.parts):
        continue
    if path.name not in ALLOWED_FILENAMES and raw_private_filename.search(path.name):
        findings.append((rel, "suspicious private-data filename"))
    if excel_temp.search(path.name):
        findings.append((rel, "temporary office file"))
    if path.suffix not in TEXT_SUFFIXES and path.name != ".gitignore":
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="utf-8", errors="ignore")
    scan_personal_numeric_patterns = not rel.startswith("outputs/audit_reruns/machine_readable_tables/")
    for lineno, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if local_path.search(stripped):
            findings.append((f"{rel}:{lineno}", "absolute local path"))
        if scan_personal_numeric_patterns and id_card.search(stripped):
            findings.append((f"{rel}:{lineno}", "possible national ID number"))
        if scan_personal_numeric_patterns and phone.search(stripped):
            findings.append((f"{rel}:{lineno}", "possible phone number"))
        if secret_assignment.search(stripped):
            findings.append((f"{rel}:{lineno}", "possible secret assignment"))

report = ["# Security Scan Report", ""]
if findings:
    report.append("Status: FAIL")
    report.append("")
    for item, reason in findings:
        report.append(f"- {item}: {reason}")
else:
    report.append("Status: PASS")
    report.append("")
    report.append("No identifiable patient data, local absolute paths, secret assignments, raw protected clinical-data files, or temporary Office files were detected by the repository scan.")

(ROOT / "security_scan_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
if findings:
    sys.exit(1)
print("Security scan passed.")
