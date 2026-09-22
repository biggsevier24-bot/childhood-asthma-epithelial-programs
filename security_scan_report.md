# Security scan report

Status: PASS

Scope: frozen final v1.0.0 local release after Figure 4 closure, 14/14 tests, and the successful cached public-data rebuild.

- No user-specific `C:/Users` or project-specific drive path was detected in public text files.
- Local paths in retained historical audit evidence were replaced with neutral provenance placeholders.
- No API key, password, access token, bearer token, or client-secret pattern was detected.
- No protected hospital participant-level data were detected. Hospital analysis remains code-only with a synthetic input template.
- No GSE118761 or GSE152004 raw GEO cache file is included; those files remain external inputs referenced only through documented environment variables.
- Public sample-level records are limited to deposited GEO identifiers and public-derived expression/metadata outputs.
- The projection assignment file contains 104 unique public samples: 55 nasal and 49 tracheal.
- Formal derivation scripts passed the separate gene-name/reference/top-N/manual-override audit.
- The frozen release scan found zero local absolute-path hits, zero credential-pattern hits, zero Office-internal path hits, zero protected hospital-file hits, and zero redistributed GSE118761/GSE152004 cache files.
- The locked manuscript results were not changed.
- Machine-readable details are retained in `audit/security_scan_report.md`.
