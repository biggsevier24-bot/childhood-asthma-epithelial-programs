# GSE118761 projection provenance audit

## Status

`GSE118761_PROJECTION_METHOD_RECOVERED = NO`

The manuscript-locked projection counts are nasal E1/E2/E3 = 30/18/7 and tracheal E1/E2/E3 = 15/29/5. No retained executable code, sample-level assignment file, workbook, R object, archive, or documented finite method-grid configuration reproduced both count vectors. These values are therefore not represented as regenerated outputs.

## Evidence reviewed

- The retained `final_locked_analysis_20260714.R` branch uses tissue-specific standardization, GSE152004 centroids, and Euclidean nearest-centroid assignment. Its retained assignment files are archived under `archive/development_history/GSE118761_projection_conflict/` because they conflict with the manuscript counts.
- A separate conservative audit branch used GSE152004 reference scaling and did not recover the manuscript counts.
- The local project search covered R/Python scripts, CSV/TSV/XLSX/RDS files, ZIP archives, and the originating Codex session history.
- The session record states that the manuscript counts were to be locked only if reproduced; a later record explicitly notes that the executable rerun differed and did not overwrite the manuscript values.
- The finite, prespecified forensic grid compared retained score branches, cohort/reference/no scaling, Euclidean/Manhattan/correlation distance, and retained centroid sources. Zero configurations recovered both manuscript count vectors. The full grid and assignments are in this directory.

## Release consequence

The repository does not hardcode sample assignments or counts. Formal projection outputs are withheld pending recovery of the original executable method or a manuscript correction. This unresolved item prevents public-release readiness but does not alter the independently reproduced GSE18965 exact32 and GSE152004 locked results.
