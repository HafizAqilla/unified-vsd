# PR Title

feat(pipeline): implement Batch 4 sequential Sobol-to-calibration workflow

# PR Description

## Summary

This PR implements Batch 4 by rewiring the main execution flow into a fixed sequential pipeline:

Phase 0 baseline -> Initial Sobol GSA -> Sobol mask generation -> Masked calibration -> Final Sobol GSA.

The previous independent toggles for calibration and GSA are removed from control flow so the optimization is always informed by sensitivity screening.

## Changes

1. Sequential pipeline enforced in main entrypoint
- Modified: main_run.m
- Removed calibration/GSA branch gating logic from runtime flow.
- Introduced mandatory execution order for all core stages.

2. Initial Sobol GSA integrated before calibration
- Uses gsa_sobol_setup + gsa_run_sobol on pre-calibration parameters.
- Extracts ST from primary metrics for parameter screening.

3. Dynamic optimization mask generation wired
- Uses create_optimization_mask with threshold = 0.10.
- Maps Sobol parameter ordering to calibration parameter ordering safely.
- Throws explicit errors for missing parameter mappings or missing ST fields.

4. Masked calibration call enabled
- Calls run_calibration(params0, clinical, scenario, optMask).
- Retains active-subset reporting and calibrated metrics computation.

5. Final Sobol GSA runs on calibrated parameters
- Runs gsa_sobol_setup + gsa_run_sobol after calibration.
- Displays initial and final Sobol summary tables.
- Exports final Sobol matrix plot via make_gsa_matrix_table.

6. Persistent artefact saving upgraded
- Saves one pipeline MAT package with:
  - initial/final Sobol config + outputs
  - threshold
  - optimization mask
  - baseline/calibrated params and metrics
  - calibration outputs
- Keeps calibrated parameter table save in results/tables.

## Why this matters

This change operationalizes the intended research workflow:
- use global sensitivity to decide what is identifiable,
- optimize only influential variables,
- then quantify final post-calibration sensitivities.

It removes accidental mode drift from toggle-based execution and improves reproducibility.

## Notes

- DO_PLOTS is retained as a visualization-only toggle.
- Sobol runs are computationally heavy; this is expected by design in Batch 4.

## Next step (Batch 5)

Add explicit <5% target pass/fail reporting and final-vs-initial GSA export overlays in validation outputs.
