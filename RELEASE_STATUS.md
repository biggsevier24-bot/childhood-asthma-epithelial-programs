# Release status

Status: `FINAL_MANUSCRIPT_WIDE_RELEASE_VALIDATED`

The v1.0.0 release provides the complete, end-to-end, manuscript-wide reproducibility workflow.

- Complete repository test suite: 14/14 PASS
- Public-data rebuild workflow: PASS
- Deterministic repair-ECM workflow: PASS
- Public cohort preprocessing and molecular-program scoring: PASS
- Correlation, regression, PCA, epithelial-state, and projection workflows: PASS
- Robustness and sensitivity workflows: PASS
- Machine-readable figure and table source data: COMPLETE
- Hospital cohort: expected `CODE_ONLY`; protected participant-level data are not included
- Privacy and security scan: PASS

Official GSE118761 and GSE152004 cache files remain external and are accessed through the documented cache interfaces. Detailed component-level validation records are retained under `docs/provenance/` and `audit/` without displacing the manuscript-wide workflow from the repository root.
