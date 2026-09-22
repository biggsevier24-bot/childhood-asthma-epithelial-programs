# Frozen final release security scan

SECURITY_SCAN_PASS = TRUE

- No local absolute Windows, macOS, sandbox, or Linux user path was detected in non-archive public text files.
- No API key, GitHub token, password assignment, bearer token, or client secret assignment was detected in non-archive public text files.
- No protected hospital participant-level input or real patient metadata is included.
- No GSE118761 or GSE152004 raw GEO cache file is included in the release.
- Hospital materials are limited to code, aggregate audit outputs, a data dictionary, and explicitly synthetic example rows.
- Public sample-level records use deposited GEO identifiers and public-derived metadata only.
- Historical/pre-closure audit artifacts are isolated under `archive/development_history/` and are not formal pipeline inputs.
