# PR Title

feat(validation): implement Batch 5 5-percent gate and Sobol overlay exports

# PR Description

## Summary

This PR implements Batch 5 by adding explicit <5% acceptance checks for primary calibration targets and exporting final-vs-initial Sobol sensitivity overlays for reporting.

The implementation is aligned with the sequential Batch 4 pipeline and AGENTS guardrails for traceability and reproducibility.

## Changes

1. Added strict primary-metric gate checks in validation
- Modified: utils/validation_report.m
- Introduced explicit 5% absolute error gate for:
  - QpQs
  - SAP_mean (MAP)
  - LVEF
- Emits clear warning lines in console when any primary metric exceeds 5%.
- Stores gate results in report.primary_gate.

2. Added optional validation export inputs
- Modified: utils/validation_report.m
- Added optional name/value args:
  - ResultsDir
  - GsaInitOut
  - GsaFinalOut
- Kept backward compatibility with prior call signature.

3. Added Sobol overlay table export (initial vs final)
- Modified: utils/validation_report.m
- Exports long-format overlay table with:
  - Metric
  - Parameter
  - ST_Initial
  - ST_Final
  - Delta_ST
- Output files:
  - results/tables/gsa_overlay_<scenario>.csv
  - results/tables/gsa_overlay_<scenario>.mat

4. Wired Batch 5 exports into main pipeline
- Modified: main_run.m
- Passes initial and final Sobol outputs plus results directory to validation_report.
- Keeps calibrated parameter package export in results/tables.

5. Updated implementation tracker
- Modified: LBFGSB_IMPLEMENTATION_BATCHES.md
- Marked Batch 5 as implemented.

## Why this matters

Batch 5 closes the loop on acceptance and reporting:
- calibration quality is now explicitly checked against the 5% target on the primary clinical endpoints,
- sensitivity changes from Phase 0 to final calibrated state are exported in a reproducible tabular form.

## Outputs introduced

- report.primary_gate
- results/tables/gsa_overlay_<scenario>.csv
- results/tables/gsa_overlay_<scenario>.mat

## Notes

- Console warning coloring depends on MATLAB terminal capabilities; warnings are explicit regardless of color support.
- Existing callers of validation_report remain compatible.
