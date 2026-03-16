# Batch 5 Validation Checklist

## Scope

Validate strict primary-metric gate checks and initial-vs-final Sobol overlay exports.

## Static checks

- Problems panel: no errors in
  - utils/validation_report.m
  - main_run.m

## Functional checks

1. Run full pipeline
- MATLAB call:
  - main_run('pre_surgery', clinical)

Expected:
- Validation report prints primary gate table with QpQs, SAP_mean, LVEF.
- Any >5% absolute error is flagged as warning.

2. Verify report struct
- Inspect `report.primary_gate`.

Expected:
- Columns: Metric, Clinical, Model, AbsError_pct, Pass_5pct
- Pass_5pct true only for <=5% absolute error.

3. Verify overlay exports
- Check files in results/tables:
  - gsa_overlay_pre_surgery.csv
  - gsa_overlay_pre_surgery.mat

Expected:
- Both files exist after run.
- CSV contains ST_Initial, ST_Final, Delta_ST.

4. Repeat on post_surgery
- MATLAB call:
  - main_run('post_surgery', clinical)

Expected:
- Same gate behavior and scenario-specific overlay files created.

## Failure conditions to monitor

- Missing Sobol metric tables should skip overlay rows and warn clearly.
- NaN clinical targets should not trigger false pass states.

## Reproducibility notes

- Keep fixed Sobol seed in gsa_sobol_setup.
- Preserve saved MAT and CSV overlays with run artifacts for auditability.
