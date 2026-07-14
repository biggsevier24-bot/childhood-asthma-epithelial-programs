# Outputs

This folder separates manuscript-locked reference materials from audit rerun outputs.

- `manuscript_locked/`: materials that correspond to confirmed locked manuscript figures/tables or user-confirmed locked reference values. These files must not be overwritten by reruns.
- `audit_reruns/`: machine-readable audit rerun outputs retained for provenance. These files may differ from formatted manuscript tables and must not be described as manuscript source data.
- `user_rerun/`: default destination for user-initiated reruns.

Public code should write new results to `outputs/user_rerun/` unless the user intentionally prepares a new manuscript release.
