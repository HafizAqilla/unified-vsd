# Batch 4 Validation Checklist

## Scope

Validate fixed sequential execution in main_run:

baseline -> initial Sobol -> mask -> calibration -> final Sobol.

## Static checks

- Problems panel: no errors in
  - main_run.m
  - calibration/run_calibration.m
  - calibration/calibration_param_sets.m
  - utils/create_optimization_mask.m

## Functional checks

1. End-to-end run (pre_surgery)
- MATLAB call:
  - main_run('pre_surgery', clinical)

Expected:
- Baseline simulation completes.
- Initial Sobol GSA executes.
- Active mask is printed with active count and parameter names.
- Masked calibration executes.
- Final Sobol GSA executes.
- Results are saved in results/gsa and results/tables.

2. End-to-end run (post_surgery)
- MATLAB call:
  - main_run('post_surgery', clinical)

Expected:
- Same stage ordering and outputs as pre_surgery.
- Parameter mapping succeeds for post_surgery free-parameter list.

3. Mask consistency check
- Inspect saved MAT artefact fields:
  - sobol_threshold
  - optMask
  - calib_names_all

Expected:
- numel(optMask) == numel(calib_names_all)
- nnz(optMask) >= 1

4. GSA output consistency check
- Verify both gsa_init_out and gsa_final_out are present in saved artefact.

Expected:
- Both structs include per-metric ST vectors and cfg metadata.

## Failure conditions to monitor

- main_run:missingGsaMapping (calibration names not found in Sobol names)
- main_run:missingPrimaryMetricST (missing primary metric ST in Sobol output)
- run_calibration:missingFmincon (Optimization Toolbox unavailable)

## Reproducibility notes

- Maintain fixed Sobol seed in gsa_sobol_setup.
- Record threshold used for each run (saved as sobol_threshold).
