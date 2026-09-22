# Public release checklist

Status: READY

- [x] GSE118761 projection rerun from sample-level public-derived inputs.
- [x] Nasal E1/E2/E3 counts verified as 21/27/7.
- [x] Tracheal E1/E2/E3 counts verified as 21/23/5.
- [x] Input hashes, score definitions, scaling parameters, centroids, distances, margins, assignments, and output hashes recorded.
- [x] Figure 5f source data rebuilt from the final assignment file.
- [x] Supplementary Table S9 rebuilt from the final assignment file.
- [x] Supplementary Table S8 audited and unchanged.
- [x] Manuscript change map generated; manuscript Word file not edited.
- [x] Clean projection rerun passed.
- [x] Figure 4 PCA runtime reproduction passed and all locked GSE118761 correlations were recovered.
- [x] Full public analysis runner completed with exit code 0 using official GEO files supplied through external cache variables.
- [x] GSE118761 and GSE152004 public-input rebuilds passed without redistributing the cache files.
- [x] All 14 repository tests passed, including Figure 4 PCA, sample-count, archive-independence, and security tests.
- [x] GSE152004 score-construction sensitivity reran successfully from the official cached raw-count input.
- [x] Historical unreproduced projection artifacts isolated under `archive/development_history/`.
- [x] Privacy, secret, and absolute-path scan passed.
- [x] Repository manifest regenerated from the frozen contents.
- [x] Final release ZIP structure, privacy exclusions, and manifest integrity independently verified.
