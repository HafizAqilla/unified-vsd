# Unified VSD Lumped-Parameter Model

Patient-specific 0-D cardiovascular simulation for pediatric Ventricular
Septal Defect (VSD) analysis before and after surgical closure.

This repository is a scientific MATLAB codebase, not a throwaway script. The
main model combines a Valenti-style 14-state lumped-parameter circulation,
time-varying chamber elastance, VSD shunt physiology, evidence-aware clinical
target handling, and reproducible calibration outputs.

Last updated: 2026-05-19.

## Current Status

The current `main` branch includes the latest systemic-flow calibration work.
The model now produces auditable pre-surgery and post-surgery run folders with
clinical consistency checks, parameter plausibility tables, baseline provenance,
candidate snapshots, and validation reports.

Scientific interpretation remains important:

- The model is physiologically meaningful and reproducible.
- The main remaining bottleneck is simultaneous patient-specific matching of
  systemic flow, systemic pressure, vascular resistance, and chamber volumes.
- For unrepaired VSD, clinical cardiac output is treated as systemic flow
  (`Qs_Lmin`) rather than raw LV outflow (`LVCO_Lmin`).
- Derived quantities such as EF, SVR, Qp/Qs, and stroke volumes are audited so
  they are not silently double-counted as independent measurements.

## Latest Changes On Main

- Added a centralized parameter registry in `config/build_parameter_registry.m`
  with scenario-aware bounds, units, source notes, and plausibility anchors.
- Added baseline provenance export through
  `config/build_baseline_provenance.m` and `scripts/audit_baseline_provenance.m`.
- Added evidence-aware case profiles in
  `src/calibration/build_case_calibration_profile.m`, including `full_data`,
  `sparse_cath`, and `synthetic_benchmark` governance paths.
- Added target-tier handling in `src/calibration/build_target_tiers.m` so direct
  measurements, derived consistency checks, and excluded targets are explicit.
- Added clinical consistency audits in `src/validation/`, including stroke
  volume, flow, pressure, and target availability checks.
- Added vascular R-C coupling helpers so coupled resistance/compliance movement
  is tracked and penalized more defensibly during calibration.
- Added parameter plausibility scoring and calibration classification through
  `src/calibration/evaluate_parameter_plausibility.m` and
  `src/calibration/classify_calibration_run.m`.
- Added best, scientific, and accepted candidate snapshots to `main_run.m`.
  Conservative rollback behavior is now visible instead of hiding strong
  but rejected candidates.
- Added age-validity and scaling annotations, including `lundquist_bsa` as the
  preferred pediatric scaling mode and `zhang` as a comparator.
- Added post-surgery warm-start support through
  `src/utils/apply_post_surgery_warm_start.m` and `run_post_surgery.m`.
- Added focused Reyna analysis scripts for uncertainty, targeted polish,
  pressure-flow grid polish, local sensitivity, and systemic-output audits.
- Added regression and governance tests for parameter registry, plausibility,
  target tiers, scaling modes, vascular R-C coupling, post-surgery warm start,
  and Reyna systemic-flow profile behavior.
- Route generated outputs into per-run archives under
  `results/runs/<timestamp>_<patient>_<scenario>/`.
- Ignore generated `results/` artifacts, MATLAB logs, and LaTeX helper files so
  GitHub keeps code and documentation instead of large run dumps.

## Scenarios

| Scenario | Physiology | VSD handling | Primary use |
|---|---|---|---|
| `pre_surgery` | Unrepaired VSD circulation | Active shunt resistance and Qp/Qs accounting | Fit pre-operative catheter/echo targets and evaluate systemic-flow consistency |
| `post_surgery` | Closed VSD / repair state | `R.vsd` closed numerically and optionally warm-started from pre-op fit | Evaluate post-operative hemodynamics from a pre-op calibrated seed |

## Requirements

Core runs require MATLAB and the Optimization Toolbox.

Optional GSA/PCE workflows use external uncertainty/sensitivity tooling such as
UQLab or SoBioS. Those packages are not committed to this repository. Keep
licensed toolboxes outside Git, for example under a local `Toolbox/` directory.

## Quick Start

From MATLAB:

```matlab
cd('D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd')
addpath(genpath(pwd))

setenv('UNIFIED_VSD_DO_GSA', '0')       % faster first run
setenv('UNIFIED_VSD_DO_PLOTS', '1')
setenv('UNIFIED_VSD_DO_OVERLAY', '1')

clinical = patient_template();          % fill fields before scientific use
main_run('pre_surgery', clinical)
```

Run an existing profile:

```matlab
clinical = patient_reyna();
main_run('pre_surgery', clinical)
```

Run the post-surgery handoff after a pre-surgery seed exists:

```matlab
run run_post_surgery
```

## Useful Environment Flags

| Variable | Effect |
|---|---|
| `UNIFIED_VSD_DO_GSA` | `0` disables GSA for faster runs; `1` enables GSA/PCE passes |
| `UNIFIED_VSD_DO_PLOTS` | `0` disables figure generation; `1` enables figures |
| `UNIFIED_VSD_DO_OVERLAY` | `0` disables overlay figures; `1` enables overlays |
| `UNIFIED_VSD_SCALING_MODE` | Selects scaling mode, commonly `lundquist_bsa` or `zhang` |
| `UNIFIED_VSD_FAST_CALIBRATION` | Enables shorter calibration settings for triage |
| `UNIFIED_VSD_FMINCON_PARALLEL` | Enables parallel fmincon behavior when appropriate |
| `UNIFIED_VSD_USE_PARPOOL` | Allows `main_run.m` to start a MATLAB parallel pool |
| `UNIFIED_VSD_MAX_FUN_EVALS` | Overrides calibration function-evaluation budget |
| `UNIFIED_VSD_MAX_ITERATIONS` | Overrides calibration iteration budget |

## Common Workflows

Fast no-plot regression sweep:

```matlab
run_quick_regression_suite
```

Audit baseline parameter provenance:

```matlab
audit_baseline_provenance
```

Compare pediatric scaling modes:

```matlab
run scripts/compare_scaling_modes.m
```

Run Reyna-focused polish or sensitivity helpers after a candidate package
exists:

```matlab
run_reyna_best_fit_search
run_reyna_15pct_gate_polish
run_reyna_targeted_grid_polish
run_reyna_pressure_flow_grid_polish
scan_reyna_gate_local_sensitivity
scan_reyna_targeted_sensitivity
run_reyna_uncertainty_ensemble
```

Run the cohort triage helper only after adding an authorized local
`patient_cohort_cases.m` manifest:

```matlab
run_patient_cohort_fast_pipeline
```

## Output Layout

Each `main_run` creates a self-contained archive:

```text
results/runs/<timestamp>_<patient>_<scenario>/
  figures/       publication-oriented plots
  gsa/           PCE/GSA checkpoints and summaries
  logs/          console diary
  mat/           parameter packages and run packages
  tables/        validation, plausibility, target, and audit tables
  run_manifest.txt
```

Important generated files include:

- `mat/params_best_candidate_<scenario>.mat`
- `mat/params_scientific_candidate_<scenario>.mat`
- `mat/params_accepted_candidate_<scenario>.mat`
- `mat/pre_to_post_seed_latest.mat` for pre-to-post handoff
- `tables/validation_best_candidate_<scenario>.csv`
- `tables/validation_accepted_candidate_<scenario>.csv`
- `tables/parameter_plausibility_<scenario>.csv`
- `tables/baseline_provenance_<scenario>.csv`
- `tables/clinical_consistency_*.csv`
- `tables/co_definition_audit_<scenario>.csv`
- `run_manifest.txt`

`results/` is intentionally ignored by Git. Commit code, documentation, and
figure-generation logic; do not commit regenerated run dumps or patient data
without documented authorization.

## Calibration And Validation Design

The current calibration pipeline is structured around reproducibility and
physiological review:

1. Load reference parameters from `config/default_parameters.m`.
2. Build baseline provenance from `config/build_baseline_provenance.m`.
3. Apply pediatric scaling with `src/utils/apply_scaling.m`.
4. Build a case profile from the available clinical evidence.
5. Map clinical values into seeded model parameters with
   `src/utils/params_from_clinical.m`.
6. Build the active parameter vector and bounds from the registry.
7. Run calibration with target-tier weights, plausibility penalties, and
   optional GSA-informed masks.
8. Simulate best and accepted candidates.
9. Export validation, clinical consistency, CO-definition, plausibility,
   baseline provenance, and manifest records.

The accepted candidate may roll back to baseline or another conservative state
when a lower-RMSE solution violates physiological or plausibility gates. The
best candidate remains saved for scientific inspection.

## Project Layout

```text
unified_vsd/
  main_run.m                         main pre/post scenario runner
  run_post_surgery.m                 post-op handoff runner
  run_patient_case.m                 patient batch entry point
  config/                            reference params, profiles, registries
  src/models/                        ODE RHS, elastance, valves, VSD shunt
  src/solvers/                       integration wrapper
  src/utils/                         scaling, ICs, indices, plotting, reports
  src/calibration/                   calibration profiles, bounds, objectives
  src/gsa/                           Sobol/PCE sensitivity helpers
  src/validation/                    clinical data and consistency audits
  scripts/                           operational audits and experiment runners
  tests/                             regression and governance checks
  docs/                              theory, data dictionary, decision memos
  results/                           generated output; ignored by Git
```

## Recommended Checks Before Sharing Results

Run the focused governance tests after calibration logic changes:

```matlab
run('tests/test_parameter_registry.m')
run('tests/test_parameter_plausibility.m')
run('tests/test_clinical_consistency_target_tiers.m')
run('tests/test_reyna_systemic_flow_profile.m')
run('tests/test_vascular_rc_coupling.m')
run('tests/test_scaling_modes.m')
run('tests/test_post_surgery_warm_start.m')
```

Run the classic physiology checks when touching model dynamics:

```matlab
run('tests/test_baseline.m')
run('tests/test_valve_logic.m')
run('tests/test_ic_perturbation.m')
run('tests/test_vsd_modes.m')
```

For broad smoke testing:

```matlab
run_quick_regression_suite
```

## Patient Data Policy

Patient-specific profiles and clinical records must be handled deliberately.
The template file is safe to version:

```text
config/patient_template.m
```

Real clinical files, raw patient records, large `.mat` outputs, and generated
run archives must not be committed without documented authorization. See
`AGENTS.md` and `docs/clinical_data_dictionary.md` for the project rules.

## Key Documentation

- `AGENTS.md` - MATLAB/physiology coding guardrails.
- `docs/theory_notes.md` - governing equations and assumptions.
- `docs/clinical_data_dictionary.md` - clinical field mapping.
- `docs/calibration_data_governance_notes.md` - target-tier policy.
- `docs/pak_dipo_co_objective_memo.md` - CO/Qs interpretation.
- `docs/reyna_pak_dipo_systemic_flow_revision.md` - Reyna systemic-flow notes.
- `docs/model_progress_bounds_and_results_20260512.md` - latest progress memo.
- `docs/scaling_method_comparison.md` - Zhang vs Lundquist workflow notes.
- `docs/zhang_vs_lundquist_audit_20260516.md` - scaling audit.
- `docs/zhang_vs_lundquist_target_values_20260516.md` - target comparison.

## References

1. Valenti (2023). Full-order 0-D cardiovascular model thesis, Table 3.3 and
   Eqs. 2.1-2.7.
2. Lundquist et al. (2025). Patient-specific pediatric cardiovascular lumped
   parameter modeling. ASAIO Journal.
3. Saltelli et al. (2010). Variance based sensitivity analysis of model output.
   Computer Physics Communications, 181, 259-270.
4. Jansen (1999). Analysis of variance designs for model output. Computer
   Physics Communications, 117, 35-43.
5. Byrd, Lu, and Nocedal (1995). A limited memory algorithm for bound
   constrained optimization. SIAM Journal on Scientific Computing, 16(5).
6. Gorlin and Gorlin (1951). Hydraulic formula for valve area. American Heart
   Journal, 41(1), 1-29.
