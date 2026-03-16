# PR Title

feat(calibration): add Sobol-based masking bridge for Batch 2 optimization

# PR Description

## Summary

This PR implements Batch 2 of the GSA-informed optimization roadmap by introducing a Sobol-driven parameter screening bridge between global sensitivity analysis and calibration.

The core objective is to freeze low-influence parameters and optimize only the active subset, improving calibration robustness and efficiency while preserving physiological interpretability.

## Changes

1. Added Sobol mask utility
- New file: utils/create_optimization_mask.m
- Converts Sobol total-order indices (ST) into a logical active-parameter mask.
- Supports both vector ST and multi-metric ST matrices.
- Uses conservative aggregation across metrics: max(ST_i,j).
- Sanitizes non-finite and negative ST values.
- Includes safety fallback to prevent empty optimization sets.

2. Updated calibration parameter configuration for masking
- Modified: calibration/calibration_param_sets.m
- Added optional input optMask.
- Added full parameter bookkeeping:
  - names_all
  - x0_all
  - lb_all
  - ub_all
- Produces active-only optimization vectors:
  - names
  - x0
  - lb
  - ub
- Added strict validation:
  - mask size must match parameter count
  - all-false masks are rejected

3. Added calibration runner compatibility
- Modified: calibration/run_calibration.m
- Added optional optMask argument.
- Preserves backward compatibility when no mask is provided.
- Passes optMask into calibration_param_sets.

4. Added unit tests for Batch 2 logic
- New file: tests/test_optimization_mask.m
- Verifies:
  - vector threshold behavior
  - matrix aggregation behavior
  - empty-mask safety fallback
  - invalid threshold error handling
  - calibration mask-size mismatch error handling

## Why this matters

This bridge operationalizes research best practice for high-dimensional calibration workflows:
- Perform global sensitivity analysis first.
- Screen out low-influence variables.
- Optimize only influential parameters.

This reduces optimizer burden and strengthens parameter identifiability while keeping the workflow reproducible and traceable.

## Backward compatibility

- Existing calls to run_calibration(params0, clinical, scenario) remain valid.
- If no mask is supplied, all scenario parameters remain active (legacy behavior).

## Next step (Batch 3)

Wire fmincon L-BFGS-B style optimization to the active subset produced by this mask pipeline.
