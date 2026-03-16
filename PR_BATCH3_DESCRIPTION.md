# PR Title

feat(calibration): migrate Batch 3 optimizer to fmincon interior-point with L-BFGS

# PR Description

## Summary

This PR implements Batch 3 of the GSA-informed optimization roadmap by replacing the previous Nelder-Mead (`fminsearch`) calibration engine with bounded `fmincon` using interior-point and L-BFGS Hessian approximation.

The implementation is fully compatible with Batch 2 masking, so only active parameters are optimized and frozen parameters remain fixed at baseline.

## Changes

1. Replaced Nelder-Mead optimizer with fmincon
- Modified: calibration/run_calibration.m
- Removed `fminsearch` path and bound-penalty wrapper objective.
- Added native bound-constrained optimization via `fmincon`.

2. Configured L-BFGS-oriented options
- `Algorithm = 'interior-point'`
- `HessianApproximation = 'lbfgs'`
- `FiniteDifferenceStepSize = 1e-5`
- Added practical solver controls:
  - `MaxFunctionEvaluations = 4000`
  - `MaxIterations = 300`
  - `OptimalityTolerance = 1e-6`
  - `StepTolerance = 1e-8`

3. Preserved and formalized active-subset optimization
- The objective remains defined on active parameters only (`calib.names`).
- Bounds are applied on active vectors (`calib.lb`, `calib.ub`).
- Full-parameter traceability is retained by reconstructing `xbest_all` from frozen + optimized subsets.

4. Extended calibration outputs for auditability
- Added fields in `calib_out`:
  - `names_all`
  - `mask`
  - `x0_all`
  - `xbest_all`
  - `exitflag`
  - `output`

5. Added explicit toolbox guard
- Throws a clear error if `fmincon` is unavailable (Optimization Toolbox missing).

## Why this matters

Compared with simplex-based Nelder-Mead, this step improves optimization quality for bounded, multi-parameter calibration by:
- enforcing physiological bounds natively,
- using a gradient-aware search strategy suitable for Batch 1 smoothing,
- scaling better as active parameter dimension grows.

## Backward compatibility

- Existing call style remains valid:
  - `run_calibration(params0, clinical, scenario)`
- Optional mask call remains supported:
  - `run_calibration(params0, clinical, scenario, optMask)`

## Next step (Batch 4)

Re-plumb `main_run.m` into a fixed sequential pipeline:
Phase 0 -> Initial GSA -> Mask creation -> Calibration -> Final GSA.
