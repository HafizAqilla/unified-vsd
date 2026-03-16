# Batch 3 Validation Checklist

## Scope

Validate migration from Nelder-Mead (`fminsearch`) to bounded `fmincon` interior-point with L-BFGS approximation.

## Static checks

- Problems panel: no errors in
  - calibration/run_calibration.m
  - calibration/calibration_param_sets.m
  - utils/create_optimization_mask.m

## Functional checks

1. Smoke test (legacy call, no mask)
- MATLAB call path:
  - run_calibration(params0, clinical, scenario)

Expected:
- fmincon runs successfully with active set = all scenario parameters.
- Bounds are respected without penalty wrapping.

2. Smoke test (masked call)
- MATLAB call path:
  - run_calibration(params0, clinical, scenario, optMask)

Expected:
- Only masked parameters are optimized.
- Frozen parameters remain at baseline values in `params_best`.
- `calib_out.xbest_all` reflects merged active + frozen values.

3. Objective improvement check
- Compare `calib_out.improvement`.

Expected:
- `calib_out.improvement >= 0` in normal successful runs.

4. Solver metadata check
- Verify `calib_out.exitflag` and `calib_out.output` are populated.

Expected:
- Exit reason is interpretable for diagnostics.

## Failure conditions to monitor

- Missing Optimization Toolbox should raise `run_calibration:missingFmincon`.
- Mask size mismatch should raise `calibration_param_sets:maskSizeMismatch`.
- All-false mask should raise `calibration_param_sets:emptyActiveSet`.

## Reproducibility notes

- Keep Batch 1 solver tolerances (rtol/atol) active in objective calls.
- Keep Batch 2 mask threshold documented with each calibration run.
