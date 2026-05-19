# Zhang vs Lundquist Scaling Audit - 2026-05-16

## Scope

This audit compares the Reyna pre-surgery results for:

- `zhang` scaling: weight-based pediatric allometry.
- `lundquist_bsa` scaling: BSA-based pediatric allometry, accepted by aliases
  `lundquist` and `lundqvist`.

The audit separates:

1. same-pipeline bounded comparison from 2026-05-14, and
2. current best accepted Lundquist result after targeted local polish.

## Sources

- Same-pipeline comparison:
  `results/scaling_method_comparison/summary_20260514.csv`
- Lundquist run:
  `results/runs/20260514_155253_reyna_pre_surgery`
- Zhang run:
  `results/runs/20260514_160222_reyna_pre_surgery`
- Current accepted Lundquist polish:
  `results/runs/20260515_173742_reyna_targeted_grid_polish`
- Scaling implementation:
  `src/utils/apply_physiological_scaling.m`
  `src/utils/apply_scaling.m`

## Same-Pipeline Result

Both methods were run through the same bounded fast comparison. Both were
classified `REJECT` by strict gates, so the comparison uses the saved
best/scientific candidate rather than the rolled-back artifact.

| Method | Baseline primary RMSE | Best primary RMSE | Best full RMSE | Hard RMSE | Soft RMSE | Param warnings |
|---|---:|---:|---:|---:|---:|---:|
| Lundquist | 0.2363 | 0.1065 | 0.1238 | 0.0932 | 0.0934 | 7 |
| Zhang | 0.2317 | 0.1734 | 0.1891 | 0.1012 | 0.1146 | 5 |

Interpretation:

- Zhang starts very slightly better at baseline primary RMSE.
- Lundquist calibrates much better: `0.1065` vs `0.1734`.
- Lundquist full RMSE is also better: `0.1238` vs `0.1891`.
- Zhang has fewer plausibility warnings, so it is not "invalid"; it is just
  less successful for Reyna's pressure-flow-volume target set.

## Key Residuals After Same-Pipeline Calibration

| Metric | Lundquist error | Zhang error | Better |
|---|---:|---:|---|
| `CO_Lmin` | -8.79% | -7.24% | Zhang slightly |
| `QpQs` | -5.99% | -11.53% | Lundquist |
| `PAP_mean` | -3.87% | -5.02% | Lundquist |
| `SAP_mean` | -0.13% | -0.44% | Lundquist |
| `RAP_mean` | +0.17% | +6.56% | Lundquist |
| `RVEDV` consistency check | +25.88% | +33.14% | Lundquist |
| `RVESV` | -10.07% | +17.32% | Lundquist |
| `LVESV` | -10.15% | -9.46% | Zhang slightly |
| `LVEF` | +17.87% | +17.75% | essentially tied |

Zhang's main failure is not CO. It is the coupled shunt/right-heart block:
`QpQs`, `RAP_mean`, `RVEDV`, and `RVESV` all move worse than Lundquist.

## Current Best Accepted Result

After targeted local polishing of the Lundquist candidate:

| Quantity | Same-pipeline Lundquist | Current accepted Lundquist |
|---|---:|---:|
| Primary RMSE | 0.1065 | 0.0884 |
| Full RMSE | 0.1238 | 0.1093 |
| Hard RMSE | 0.0932 | 0.0826 |
| Soft RMSE | 0.0934 | 0.0891 |
| Primary targets above 10% | 4-ish in strict comparison | 3 |
| Physiology valid | yes | yes |

Current accepted Lundquist target errors:

| Metric | Error |
|---|---:|
| `RAP_mean` | -0.44% |
| `PAP_mean` | +0.36% |
| `SAP_min` | -7.43% |
| `SAP_max` | -9.44% |
| `SAP_mean` | -1.34% |
| `QpQs` | -5.34% |
| `SVR` | +9.47% |
| `CO_Lmin` | -9.95% |
| `LVEDV` | +14.65% |
| `LVESV` | -2.93% |
| `RVESV` | -9.15% |
| `LVEF` | +13.91% |
| `RVEDV` consistency check | +25.65% |

## Why Lundquist Wins Here

For Reyna:

- Weight = 13.4 kg
- BSA = 0.588 m2
- Adult BSA reference = 1.73 m2
- Adult weight reference = 70 kg

This gives two very different size anchors:

- Zhang weight ratio: `13.4 / 70 = 0.191`
- Lundquist BSA ratio: `0.588 / 1.73 = 0.340`

Because Zhang starts from a smaller weight anchor, it tends to produce smaller
volume/compliance anchors and stronger pediatric resistance/elastance shifts.
That can be useful in some pediatric settings, but in Reyna it leaves less room
to match the measured pressure-flow pattern while also keeping chamber volumes
visible.

The code difference is also asymmetric:

- Lundquist scales vascular inertance.
- Zhang currently does not scale inertance.
- Lundquist is now the default in `apply_scaling.m`.

## Data Governance Caveat

Neither scaling method solves the clinical inconsistency:

- Fick/QpQs implies pulmonary stroke volume about `34.35 mL/beat`.
- Echo-derived LV stroke volume is `21.7 mL/beat`.
- Echo-derived RV stroke volume is `18.5 mL/beat`.

Therefore the persistent `RVEDV`, `LVEDV`, and `LVEF` errors are not purely a
scaling-method problem. They are a mixed-modality target conflict.

## Verdict

Use Lundquist/BSA scaling as the default for Reyna and the current cohort.

Reason:

- It gives the best calibrated hemodynamic fit.
- It keeps all direct pressure-flow targets under 10% after targeted polish.
- It outperforms Zhang in Qp/Qs, PAP, RAP, RVESV, and full RMSE.

Keep Zhang as a sensitivity analysis, not as the primary Reyna result.

Zhang is useful to show that the conclusion is scaling-prior-sensitive, but the
current evidence does not support Zhang as the main scaling method for Reyna.
