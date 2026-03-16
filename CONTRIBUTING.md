# Contributing to the Unified VSD Project

First off, thank you for considering contributing to the Unified VSD cardiovascular simulation repository.

### Philosophy and Guardrails
This codebase adheres to strict scientific standards as described in `AGENTS.md` and `docs/VIBECODING_GUARDRAILS.md`. 
Every contribution must be:
- **Physically traceable**: Mapped to real physical properties or clinical metrics.
- **Numerically defensible**: Documenting solver choices, conditions, and assumptions.
- **Reproducible**: Seeded, stateless, and fully versioned.
- **Peer-reviewable**: Ready to be included as a scientific appendix.

### Branch Workflow
- `main` is stable and validated.
- Develop features or bug fixes in `feature/name` or `fix/name` branches.
- Use explicit commit messages explaining *what changed physically or numerically*.

### Testing
- Run `tests/test_baseline.m` to ensure physiological outputs fall within baseline healthy target ranges before submitting any changes.
- Ensure the GSA and optimization pipeline runs successfully end-to-end.

### Pull Requests
Your PR description must include:
1. What was changed.
2. Why it matters medically/numerically.
3. Verification steps taken (e.g., metric limit checks < 5%).
