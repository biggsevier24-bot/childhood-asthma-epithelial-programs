# Security scan report

Status: PASS

Scope: complete, end-to-end, manuscript-wide v1.0.0 reproducibility workflow, including public cohort code and outputs, sensitivity analyses, provenance records, tests, and privacy-preserving hospital-cohort materials.

- No user-specific `C:/Users` or project-specific drive path was detected in public text files.
- Local paths in retained historical audit evidence were replaced with neutral provenance placeholders.
- No API key, password, access token, bearer token, or client-secret pattern was detected.
- No protected hospital participant-level data were detected. Hospital analysis remains code-only with a synthetic input template and aggregate outputs.
- No GSE118761 or GSE152004 raw GEO cache file is included; those files remain external inputs referenced through documented cache variables.
- Public sample-level records are limited to deposited GEO identifiers and public-derived expression, metadata, score, state, and projection outputs.
- Formal derivation scripts passed the separate gene-name, reference, top-N, and manual-override audit.
- The frozen release scan found zero local absolute-path hits, zero credential-pattern hits, zero Office-internal path hits, zero protected hospital-file hits, and zero redistributed GSE118761/GSE152004 cache files.
- Analysis code, locked numerical results, cohort definitions, thresholds, and gene sets were not changed by the repository-presentation update.
- Machine-readable security details are retained in `audit/security_scan_report.md`.
