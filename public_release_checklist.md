# Public Release Checklist

- [x] No identifiable patient data detected by `tests/test_no_private_information.py`.
- [x] No real hospital-level records included.
- [x] No passwords or tokens detected.
- [x] No absolute local paths detected.
- [x] README complete.
- [x] Public data accessions documented.
- [x] Required intermediate inputs documented.
- [x] Gene sets included.
- [x] Actual software versions recorded.
- [x] CITATION.cff included.
- [x] .zenodo.json included.
- [x] LICENSE included.
- [x] GitHub repository status recorded as pending; no invalid GitHub URL placeholder remains.
- [x] Zenodo DOI placeholder removed before first archive; README records archived release as pending.
- [x] Current manuscript results were not recalculated or altered.
- [x] Protected hospital patient-level data excluded; only a synthetic format template is included.
- [x] Security scan passed and wrote `security_scan_report.md`.
- [x] Final public repository audit included.

Manual steps before public release:

1. Add the real GitHub URL to README, CITATION.cff, .zenodo.json if needed, and the manuscript after creating the repository.
2. Add the real Zenodo DOI to README, CITATION.cff, .zenodo.json if needed, and the manuscript after Zenodo creates it.
3. Review all figure/table captions in the manuscript against locked values before upload.
4. Confirm institutional data-sharing language for protected clinical data.
