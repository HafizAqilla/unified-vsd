# Reyna and Razka Best Pre-Surgery Error/RMSE Summary

Generated: 2026-05-21

This note compares the best pre-surgery results we have run so far for Reyna
and Razka, separated by clinical evidence class:

- **Protocol direct**: clinical value was directly reported in the protocol or
  catheter record and used as a comparator.
- **Protocol derived**: clinical comparator was calculated from available
  protocol values, for example mean systemic pressure or shunt flow.
- **No comparator**: model output exists, but the protocol did not provide a
  trusted pre-surgery comparator in the selected result table. These rows do
  not have error percentage and are not included in RMSE.

## Source Runs

| Patient | Selected run | Why selected | Source table |
|---|---|---|---|
| Reyna | `results/runs/20260520_100315_reyna_shunt_flow_finetune` | Best shunt-flow finetune run: lowest good primary/full RMSE with physiologic validity and near-zero shunt-flow error | `tables/reyna_shunt_flow_finetune_metrics.csv` |
| Razka | `results/runs/20260520_110545_razka_pre_surgery` | Better calibrated RMSE than the later Razka run and no rollback | `tables/validation_best_candidate_pre_surgery.csv` |

The selected runs are the best result folders from the runs reviewed here,
not necessarily the newest timestamp.

## Error and RMSE Definitions

Error percentage:

```text
error_pct = 100 * (model_value - clinical_value) / clinical_value
```

Group RMSE percentage:

```text
RMSE_pct = sqrt(mean(error_pct^2))
```

The project RMSE values in the run manifests are the same quantity as a
fraction, so `0.0624` means approximately `6.24%`.

## Best-Run Summary

| Patient | Baseline/start RMSE | Best calibrated RMSE | Best hard-target RMSE | Best soft-target RMSE | Status |
|---|---:|---:|---:|---:|---|
| Reyna | 0.0776 | **0.0624** | 0.0459 | 0.0846 | physiologic, shunt-flow tuned |
| Razka | 0.2290 | **0.0687** | 0.0367 | 0.0811 | `PROMISING_NEAR_MISS` |

Notes:

- Reyna `0.0459` is the best direct objective from the finetune manifest.
- Reyna primary/full RMSE is `0.0624`; the all-comparator table below also
  gives `6.24%`, which is the same value in percent units.
- Razka hard-target RMSE is much better than the full RMSE because the larger
  residuals are mostly pressure extrema rather than mean pressure/QpQs.

## RMSE by Evidence Class

| Patient | Evidence class | n with comparator | Start/baseline RMSE % | Best RMSE % | Best RMSE fraction | Comment |
|---|---:|---:|---:|---:|---:|---|
| Reyna | Protocol direct | 7 | 5.70 | **6.47** | 0.0647 | Direct pressure/QpQs rows; worsened slightly after shunt-flow polish because the polish prioritized derived shunt consistency |
| Reyna | Protocol derived | 3 | 11.18 | **5.66** | 0.0566 | Big improvement, driven by `Q_shunt_Lmin` |
| Reyna | All with comparator | 10 | 7.76 | **6.24** | 0.0624 | Main transparent clinical RMSE |
| Reyna | No comparator | 14 | N/A | N/A | N/A | Model-only outputs; no error percentage |
| Razka | Protocol direct | 7 | 20.57 | **7.30** | 0.0730 | Large improvement after calibration |
| Razka | Protocol derived | 1 | 35.14 | **2.14** | 0.0214 | Only `SAP_mean` has a derived clinical comparator in the selected table |
| Razka | All with comparator | 8 | 22.90 | **6.87** | 0.0687 | Main transparent clinical RMSE |
| Razka | No comparator | 16 | N/A | N/A | N/A | CO/shunt/PVR/SVR/volumes/EF are model-only for this protocol |

## RMSE by Calibration Tier

| Patient | Tier | n | Start/baseline RMSE % | Best RMSE % | Best RMSE fraction |
|---|---|---:|---:|---:|---:|
| Reyna | hard | 5 | 4.71 | **4.59** | 0.0459 |
| Reyna | soft | 3 | 12.41 | **8.46** | 0.0846 |
| Reyna | validation_only | 2 | 3.88 | **5.88** | 0.0588 |
| Razka | hard | 4 | 19.28 | **3.67** | 0.0367 |
| Razka | soft | 2 | 34.43 | **8.11** | 0.0811 |
| Razka | validation_only | 2 | 12.98 | **9.80** | 0.0980 |

## Reyna Detailed Metric Table

| Metric | Category | Unit | Clinical | Start model | Start err % | Best model | Best err % | Tier | Calib | Primary RMSE |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| RAP_min | No comparator | mmHg | N/A | 2.549 | N/A | 2.605 | N/A | unavailable | 0 | 0 |
| RAP_mean | Protocol direct | mmHg | 5.000 | 5.003 | 0.0622 | 5.113 | 2.267 | hard | 1 | 1 |
| RAP_max | No comparator | mmHg | N/A | 10.161 | N/A | 10.482 | N/A | unavailable | 0 | 0 |
| LAP_min | No comparator | mmHg | N/A | 4.494 | N/A | 4.575 | N/A | unavailable | 0 | 0 |
| LAP_mean | No comparator | mmHg | N/A | 6.403 | N/A | 6.523 | N/A | unavailable | 0 | 0 |
| LAP_max | No comparator | mmHg | N/A | 8.394 | N/A | 8.584 | N/A | unavailable | 0 | 0 |
| PAP_min | Protocol direct | mmHg | 10.000 | 9.860 | -1.398 | 10.114 | 1.138 | validation_only | 0 | 1 |
| PAP_max | Protocol direct | mmHg | 20.000 | 21.062 | 5.308 | 21.647 | 8.235 | validation_only | 0 | 1 |
| PAP_mean | Protocol direct | mmHg | 15.000 | 14.771 | -1.529 | 15.198 | 1.320 | hard | 1 | 1 |
| SAP_min | Protocol direct | mmHg | 57.000 | 55.339 | -2.913 | 54.630 | -4.158 | soft | 1 | 1 |
| SAP_max | Protocol direct | mmHg | 100.000 | 86.388 | -13.612 | 85.957 | -14.043 | soft | 1 | 1 |
| SAP_mean | Protocol derived | mmHg | 71.300 | 70.161 | -1.597 | 69.460 | -2.580 | hard | 1 | 1 |
| QpQs | Protocol direct | - | 1.194 | 1.181 | -1.114 | 1.214 | 1.634 | hard | 1 | 1 |
| Q_shunt_Lmin | Protocol derived | L/min | 0.6640 | 0.5553 | -16.374 | 0.6618 | -0.3368 | soft | 1 | 1 |
| PVR | No comparator | WU | N/A | 2.306 | N/A | 2.306 | N/A | unavailable | 0 | 0 |
| SVR | No comparator | WU | N/A | 21.205 | N/A | 20.761 | N/A | unavailable | 0 | 0 |
| CO_Lmin | Protocol derived | L/min | 3.423 | 3.073 | -10.230 | 3.099 | -9.453 | hard | 1 | 1 |
| VSD_frac_pct | No comparator | % | N/A | 15.303 | N/A | 17.591 | N/A | unavailable | 0 | 0 |
| LVEDV | No comparator | mL | N/A | 46.161 | N/A | 46.806 | N/A | unavailable | 0 | 0 |
| LVESV | No comparator | mL | N/A | 16.948 | N/A | 16.730 | N/A | unavailable | 0 | 0 |
| RVEDV | No comparator | mL | N/A | 38.685 | N/A | 39.386 | N/A | unavailable | 0 | 0 |
| RVESV | No comparator | mL | N/A | 10.900 | N/A | 11.014 | N/A | unavailable | 0 | 0 |
| LVEF | No comparator | - | N/A | 0.6328 | N/A | 0.6426 | N/A | unavailable | 0 | 0 |
| RVEF | No comparator | - | N/A | 0.7182 | N/A | 0.7204 | N/A | unavailable | 0 | 0 |

## Razka Detailed Metric Table

| Metric | Category | Unit | Clinical | Baseline model | Baseline err % | Calibrated model | Calibrated err % | Tier | Calib | Primary RMSE |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| RAP_min | No comparator | mmHg | N/A | 3.674 | N/A | 3.514 | N/A | unavailable | 0 | 0 |
| RAP_mean | Protocol direct | mmHg | 5.000 | 5.577 | 11.544 | 5.277 | 5.548 | hard | 1 | 1 |
| RAP_max | No comparator | mmHg | N/A | 9.001 | N/A | 8.019 | N/A | unavailable | 0 | 0 |
| LAP_min | No comparator | mmHg | N/A | 3.257 | N/A | 3.101 | N/A | unavailable | 0 | 0 |
| LAP_mean | No comparator | mmHg | N/A | 4.632 | N/A | 4.351 | N/A | unavailable | 0 | 0 |
| LAP_max | No comparator | mmHg | N/A | 6.860 | N/A | 6.272 | N/A | unavailable | 0 | 0 |
| PAP_min | Protocol direct | mmHg | 10.000 | 11.827 | 18.273 | 10.799 | 7.989 | validation_only | 0 | 1 |
| PAP_max | Protocol direct | mmHg | 22.000 | 21.617 | -1.741 | 19.510 | -11.319 | validation_only | 0 | 1 |
| PAP_mean | Protocol direct | mmHg | 15.000 | 16.597 | 10.645 | 15.031 | 0.2066 | hard | 1 | 1 |
| SAP_min | Protocol direct | mmHg | 59.000 | 75.970 | 28.763 | 52.773 | -10.555 | soft | 1 | 1 |
| SAP_max | Protocol direct | mmHg | 82.000 | 114.219 | 39.292 | 85.671 | 4.477 | soft | 1 | 1 |
| SAP_mean | Protocol derived | mmHg | 70.000 | 94.599 | 35.141 | 68.503 | -2.139 | hard | 1 | 1 |
| QpQs | Protocol direct | - | 1.210 | 1.239 | 2.430 | 1.158 | -4.287 | hard | 1 | 1 |
| Q_shunt_Lmin | No comparator | L/min | N/A | 0.5602 | N/A | 0.3376 | N/A | unavailable | 0 | 0 |
| PVR | No comparator | WU | N/A | 4.126 | N/A | 4.319 | N/A | unavailable | 0 | 0 |
| SVR | No comparator | WU | N/A | 38.046 | N/A | 29.613 | N/A | unavailable | 0 | 0 |
| CO_Lmin | No comparator | L/min | N/A | 2.340 | N/A | 2.135 | N/A | unavailable | 0 | 0 |
| VSD_frac_pct | No comparator | % | N/A | 19.317 | N/A | 13.653 | N/A | unavailable | 0 | 0 |
| LVEDV | No comparator | mL | N/A | 39.082 | N/A | 37.518 | N/A | unavailable | 0 | 0 |
| LVESV | No comparator | mL | N/A | 12.851 | N/A | 15.151 | N/A | unavailable | 0 | 0 |
| RVEDV | No comparator | mL | N/A | 45.543 | N/A | 38.870 | N/A | unavailable | 0 | 0 |
| RVESV | No comparator | mL | N/A | 23.724 | N/A | 19.167 | N/A | unavailable | 0 | 0 |
| LVEF | No comparator | - | N/A | 0.6712 | N/A | 0.5962 | N/A | unavailable | 0 | 0 |
| RVEF | No comparator | - | N/A | 0.4791 | N/A | 0.5069 | N/A | unavailable | 0 | 0 |

## Clinical Interpretation

### Reyna

The strongest result is the shunt-flow fit. `Q_shunt_Lmin` improved from
`-16.37%` error to `-0.34%`, while all-comparator RMSE improved from `7.76%`
to `6.24%`. The main remaining residual is systemic systolic pressure:
`SAP_max = -14.04%`. This is why Reyna's direct-protocol RMSE is slightly
worse after the shunt-flow polish, even though the overall full RMSE improved.

Rows marked `No comparator` are not validation failures. They are model outputs
without trusted clinical pre-surgery comparator values in the selected result
table, so they must be reported as simulation-derived outputs only.

### Razka

Razka improves strongly from baseline to calibrated: all-comparator RMSE drops
from `22.90%` to `6.87%`. Mean pulmonary pressure and derived mean systemic
pressure fit very well:

- `PAP_mean`: `+0.21%`
- `SAP_mean`: `-2.14%`
- `QpQs`: `-4.29%`

The largest residuals are pressure extrema:

- `PAP_max`: `-11.32%`
- `SAP_min`: `-10.55%`
- `PAP_min`: `+7.99%`

CO, shunt flow, PVR, SVR, ventricular volumes, and EF are model-only for Razka
because the protocol did not provide CO/Fick or echo volume comparators.

## Bottom Line

| Patient | Best all-comparator RMSE | Best hard-target RMSE | Strongest fit | Main limitation |
|---|---:|---:|---|---|
| Reyna | **0.0624** | **0.0459** | Shunt flow and mean pressures | Systemic systolic pressure remains low vs protocol |
| Razka | **0.0687** | **0.0367** | Mean pressures and QpQs | No CO/volume comparators; pressure extrema still residual |

Both patients support the same main conclusion: the model can reproduce the
available pre-surgery hemodynamic targets for two VSD patients with roughly
`6-7%` transparent all-comparator RMSE, but missing CO/volume data limits how
strongly we can validate derived flow and chamber-volume outputs.
