# Reyna 15% Gate Polish Experiment

Date: 2026-05-13

## Purpose

This experiment tested whether the latest Reyna scientific candidate could be polished so the remaining high-error metrics fall below a 15% absolute-error gate without damaging already acceptable targets.

Starting candidate:

`results/runs/20260511_223816_reyna_pre_surgery/mat/params_best_candidate_pre_surgery.mat`

Starting remaining offenders:

| Metric | Error |
|---|---:|
| CO_Lmin | -18.17% |
| SAP_max | -18.14% |
| RVEDV | +24.85% |
| LVEF | +15.45% |
| LVESV | -14.69% |

## Implemented Method

Two reproducible MATLAB scripts were added:

| Script | Purpose |
|---|---|
| `scripts/run_reyna_15pct_gate_polish.m` | Direct-ODE 15% hinge polish starting from the saved best candidate. |
| `scripts/scan_reyna_gate_local_sensitivity.m` | One-at-a-time local sensitivity scan around the best candidate. |

The polish objective:

- targets `CO_Lmin`, `SAP_max`, `SAP_min`, `RVEDV`, `LVEF`, and `LVESV`
- protects `PAP_mean`, `QpQs`, `RAP_mean`, `SAP_mean`, `SVR`, `LVEDV`, `RVESV`, and `LAP_mean`
- starts from the saved scientific candidate, not the rolled-back baseline
- uses direct ODE simulation rather than a surrogate
- tries the requested small active set first
- optionally tests a minimal `vsd.Cd` expansion
- includes an explicit RV-bound relaxation probe for `E.RV.EA` and `V0.RV`

## Best Saved Result

Best saved run:

`results/runs/20260513_001445_reyna_15pct_gate_polish`

This reproduces the same selected best candidate first observed in:

`results/runs/20260512_233056_reyna_15pct_gate_polish`

| Metric | Start Error | Polished Error | Status |
|---|---:|---:|---|
| CO_Lmin | -18.17% | -15.88% | near miss |
| SAP_max | -18.14% | -16.92% | near miss |
| SAP_min | -14.88% | -7.60% | pass |
| RVEDV | +24.85% | +24.31% | persistent failure |
| LVEF | +15.45% | +14.23% | pass |
| LVESV | -14.69% | -11.48% | pass |

RMSE improved:

| Run | RMSE |
|---|---:|
| Starting scientific candidate | 0.1287 |
| Best saved 15% polish | 0.1153 |

The best saved polish reduced overall RMSE and fixed `SAP_min`, `LVEF`, and `LVESV`, but it did not make every measured metric pass the 15% gate.

## Local Sensitivity Finding

The local scan showed that the RV volume error is constrained by RV parameters:

- increasing `E.RV.EA` reduces `RVEDV`, but `E.RV.EA` reaches its upper bound
- decreasing `V0.RV` reduces `RVEDV`, but `V0.RV` reaches its lower bound
- changes that improve systemic pressure can worsen CO or RVEDV
- changes that improve RVEDV alone do not solve CO and SAP_max simultaneously

This means the remaining error is not just an optimizer failure. It is a pressure-flow-volume consistency issue.

## Clinical Consistency Check

The clinical flow and volume targets are internally difficult to satisfy together.

Using Reyna's HR = 119 bpm:

| Quantity | Value |
|---|---:|
| Cath Qs | 3.423 L/min |
| Qp from Qs x QpQs | 4.087 L/min |
| Stroke volume implied by Qs | 28.76 mL |
| Stroke volume implied by Qp | 34.35 mL |
| LV volume SV from LVEDV-LVESV | 21.70 mL |
| RV volume SV from RVEDV-RVESV | 18.50 mL |
| RVEDV required for Qp if RVESV = 12 mL | 46.35 mL |
| 15% upper gate for clinical RVEDV | 35.07 mL |

Therefore, the cath pulmonary flow target implies an RV stroke volume far larger than the RV volume measurements allow. A model matching Qp/Qs and CO cannot also keep RVEDV below the 15% gate unless it violates flow-volume consistency.

## Conclusion

The current best engineering result is the 15% gate polish with RMSE 0.1153. It is scientifically better than the previous candidate for several metrics, but forcing every measured target below 15% is not currently defensible because `RVEDV` conflicts with the measured flow block.

Recommended next step:

Treat `RVEDV` as a lower-confidence or consistency-check target, not as a hard 15% calibration gate, unless the clinical RV volume measurement is re-verified from the source record.
