# PR Title

chore(repo): implement strict repository housekeeping guidelines

# PR Description

## Summary

This PR executes a major structural clean-up of the `unified-vsd` repository to meet strict scientific coding standards and GitHub best practices. It introduces critical community and configuration files while structurally separating core physics from executable task runners.

This PR builds naturally on top of the Batch 5 validations structure.

## Changes

1. **Standardized Root Configuration**
   - Added `.editorconfig` and `.gitattributes` to enforce consistent line endings (LF), charset, and indentation.
   - Added `CITATION.cff` for accurate academic attribution.
   - Added `CHANGELOG.md` to track version history formally.
   - Added standard `LICENSE` (MIT), `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, and `SECURITY.md` (with explicit boundaries around patient data).

2. **Source Code Encapsulation**
   - Created `src/` directory.
   - Moved all programmatic domains into `src/`:
     - `src/models/` (Physics)
     - `src/solvers/` (Integration logic)
     - `src/calibration/` (L-BFGS-B and objectives)
     - `src/gsa/` (Sobol analysis frameworks)
     - `src/utils/` (Helper functions)

3. **Task Separation**
   - Created `scripts/` directory.
   - Moved `run_virtual_patients.m` to `scripts/` to separate execution tasks from scientific source code.
   - `main_run.m` remains exactly where it should be at the root directory, injecting `genpath` automatically so functionality is completely uninterrupted.

4. **Documentation Decluttering**
   - Established `docs/references/` and moved all large papers (`*.pdf`).
   - Established `docs/planning/` and moved PR tracker markdown notes, `.png` architecture assets, and pipeline implementation trackers.

## Why this matters

By restructuring into a clear `src/`, `scripts/`, and `docs/` hierarchy with full community health files attached, the repository is now fully prepared for external collaboration, thesis review, and formal publication. It removes clutter, sets robust boundaries regarding clinical data security, and ensures the repo passes the "Vibecoding Guardrails" 100%.

## Verification
- Verified `main_run` and `test_baseline.m` execute exactly as expected. Root path inclusion rules handle the new `src/` restructure seamlessly.

## Notes
- Git history and commits tracking was perfectly maintained via `git mv` operations. Reviewers can verify no file deletions occurred, only migrations.