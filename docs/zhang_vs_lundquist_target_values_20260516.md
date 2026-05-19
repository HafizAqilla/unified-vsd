# Zhang vs Lundquist Target Values - 2026-05-16

This table compares the best same-pipeline Reyna pre-surgery candidates from
the 2026-05-14 scaling comparison.

Source files:

- Zhang: `results/runs/20260514_160222_reyna_pre_surgery/tables/validation_best_candidate_pre_surgery.csv`
- Lundquist: `results/runs/20260514_155253_reyna_pre_surgery/tables/validation_best_candidate_pre_surgery.csv`

The gate is a simple absolute error check against 10%. The raw record labels
the function measurement as generic `EF`; it is treated as model-side `LVEF`
because it is documented with LV dimensions and LV volume derivation, not as
a separately measured right-ventricular EF. `RVEDV` is retained as visible
evidence but remains consistency-check-only because the echo RV stroke volume
conflicts with catheter-derived pulmonary stroke volume.

| Metric | Clinical | Zhang model | Zhang error | Zhang 10% gate | Lundquist model | Lundquist error | Lundquist 10% gate | Better |
|---|---:|---:|---:|---|---:|---:|---|---|
| RAP_mean | 5.000 | 5.328 | +6.56% | pass | 5.009 | +0.17% | pass | Lundquist |
| LAP_mean | 8.000 | 3.948 | -50.64% | fail | 6.336 | -20.80% | fail | Lundquist |
| PAP_mean | 15.000 | 14.247 | -5.02% | pass | 14.419 | -3.87% | pass | Lundquist |
| SAP_min | 57.000 | 55.457 | -2.71% | pass | 56.181 | -1.44% | pass | Lundquist |
| SAP_max | 100.000 | 86.864 | -13.14% | fail | 87.530 | -12.47% | fail | Lundquist |
| SAP_mean | 71.300 | 70.985 | -0.44% | pass | 71.209 | -0.13% | pass | Lundquist |
| Qp/Qs | 1.194 | 1.056 | -11.53% | fail | 1.123 | -5.99% | pass | Lundquist |
| SVR | 19.370 | 20.678 | +6.75% | pass | 21.205 | +9.47% | pass, barely | Zhang |
| CO_Lmin | 3.423 | 3.175 | -7.24% | pass | 3.122 | -8.79% | pass | Zhang |
| LVEDV | 41.000 | 46.196 | +12.67% | fail | 45.920 | +12.00% | fail | Lundquist |
| LVESV | 19.300 | 17.475 | -9.46% | pass | 17.342 | -10.15% | fail, barely | Zhang |
| RVEDV | 30.500 | 40.607 | +33.14% | fail, consistency-check only | 38.393 | +25.88% | fail, consistency-check only | Lundquist |
| RVESV | 12.000 | 14.078 | +17.32% | fail | 10.792 | -10.07% | fail, barely | Lundquist |
| EF reported, treated as LVEF | 0.528 | 0.622 | +17.75% | fail | 0.622 | +17.87% | fail | Zhang, tiny |

Summary:

| Method | Pass count, 10% gate | Fail count | Primary RMSE | Full RMSE |
|---|---:|---:|---:|---:|
| Zhang | 7 / 14 | 7 / 14 | 0.1734 | 0.1891 |
| Lundquist | 7 / 14 | 7 / 14 | 0.1065 | 0.1238 |

Interpretation:

- The pass/fail count is tied, but Lundquist has smaller errors on the clinically
  important pressure-shunt-right-heart cluster: `RAP_mean`, `PAP_mean`,
  `SAP_mean`, `Qp/Qs`, `RVEDV`, and `RVESV`.
- Zhang is slightly better for `CO_Lmin`, `SVR`, and `LVESV`.
- `RVEF` is not present as a direct clinical target. It should not be counted
  as a target/gate in this comparison. If someone asks, the RV volume pair
  mathematically implies `(30.5 - 12.0) / 30.5 = 0.607`, but that is a derived
  consistency note, not a measured RVEF field.
- The largest shared problem is not the scaling formula alone; it is the
  mixed-modality volume conflict, especially `RVEDV`, `LAP_mean`, and the
  reported EF treated as `LVEF`.
- The later accepted Lundquist polish improves beyond this table, but there is
  no equivalent later Zhang polish yet, so this table is the fair scaling-method
  comparison.
