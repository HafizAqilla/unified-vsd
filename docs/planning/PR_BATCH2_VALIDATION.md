# Batch 2 Validation Checklist

## Scope

Validate Sobol-to-optimization masking behavior and compatibility with calibration flow.

## Static checks

- Pylance/Problems panel: no errors in
  - calibration/calibration_param_sets.m
  - calibration/run_calibration.m
  - utils/create_optimization_mask.m
  - tests/test_optimization_mask.m

## Functional checks

1. Run unit test script
- Command in MATLAB:
  - test_optimization_mask

Expected:
- All tests pass.

2. Legacy calibration path (no mask)
- Command path:
  - main_run with DO_CALIBRATION = true and no mask argument

Expected:
- run_calibration executes with full active set.
- No API break in existing workflow.

3. Masked calibration path
- Manual smoke test:
  - Construct a logical mask matching scenario free-parameter count.
  - Call run_calibration(params0, clinical, scenario, optMask).

Expected:
- Only masked parameters are optimized.
- Non-active parameters remain at baseline values.

## Failure conditions to monitor

- Mask length mismatch should throw calibration_param_sets:maskSizeMismatch.
- All-false mask should throw calibration_param_sets:emptyActiveSet.
- Threshold outside [0, 1] should error in create_optimization_mask.

## Reproducibility notes

- Keep Sobol seed fixed in GSA setup for deterministic mask generation.
- Store threshold value used to derive optMask in run metadata.
