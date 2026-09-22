# Release status

Status: `FINAL_LOCAL_RELEASE_FROZEN`

Final local v1.0.0 validation completed on 2026-09-22:

- Figure 4 / GSE118761 PCA runtime reproduction: PASS
- Figure 4 locked correlations: PASS
- complete repository test suite: 14/14 PASS
- `run_public_analyses.R`: PASS, exit code 0
- GSE118761 and GSE152004 public-input rebuilds: PASS using official GEO files supplied through the documented external cache variables
- GSE118761 projection preparation and deliverables: PASS
- GSE152004 score-construction sensitivity: PASS
- external ECM branch: PASS
- hospital cohort: expected `CODE_ONLY`; protected patient-level data are not included
- privacy and security scan: PASS

The external GEO cache files are not redistributed. GitHub upload is intentionally outside this local freeze task.
