# Reyna Pak Dipo Systemic-Flow Revision

Date: 2026-05-11
Case: Reyna pre-surgery VSD calibration
Best pre-BSA-revision scientific run: `results/runs/20260510_225025_reyna_pre_surgery`
Corrected-BSA test run: `results/runs/20260511_134027_reyna_pre_surgery`
Current raw-PDF/PAP-mean-only run: `results/runs/20260511_203930_reyna_pre_surgery`

## Anthropometry Decision

The active Reyna patient file now follows the raw protocol anthropometry so the
main case remains strictly patient-specific:

```text
weight_kg = 13.4
height_cm = 95.0
BSA = 0.588 m^2
```

Keisya's 2026-05-11 anthropometry revision remains documented as an alternate
baseline-scaling experiment:

```text
weight_kg = 14.0
height_cm = 98.0
BSA = sqrt(14.0 * 98.0 / 3600) = 0.6173419726 m^2
```

Corrected-BSA run outcome:

| Quantity | Previous best | Corrected-BSA run |
|---|---:|---:|
| BSA | 0.5880 m^2 | 0.6173 m^2 |
| Baseline RMSE | 0.2513 | 0.2454 |
| Best scientific RMSE | 0.1195 | 0.1087 |
| Scientific improvement vs baseline | 52.4% | 55.7% |
| Operational status | REJECT / rollback | REJECT / rollback |

The corrected BSA moved the scientific candidate in the right direction but did
not by itself reach the target `RMSE < 0.10`. The candidate is still rejected by
the strict gate, mainly because `RAP_mean` and `SAP_mean` remain outside 5%.

The subsequent raw-PDF/PAP-mean-only run confirmed the remaining bottleneck:
systemic preload-flow-volume consistency, especially `CO_Lmin`, `SVR`,
`SAP_mean`, `RAP_mean`, `RVEDV`, `LVESV`, and `LVEF`.

Targeted systemic-flow polish run outcome:

| Quantity | Raw-PDF/PAP-mean-only run | Targeted systemic polish |
|---|---:|---:|
| Run folder | `20260511_203930_reyna_pre_surgery` | `20260511_223816_reyna_pre_surgery` |
| Baseline RMSE | 0.2699 | 0.2699 |
| Best scientific RMSE | 0.1702 | 0.1287 |
| Scientific improvement vs baseline | 36.9% | 52.3% |
| Operational status | REJECT / rollback | REJECT / rollback |

The targeted run improved the intended bottleneck family without reintroducing
PAP systolic/diastolic waveform targets. Best-candidate errors moved as follows:

| Metric | Previous error % | Targeted error % |
|---|---:|---:|
| `CO_Lmin` | -24.49 | -18.17 |
| `SVR` | +17.92 | +7.30 |
| `SAP_mean` | -11.39 | -10.86 |
| `RAP_mean` | -17.17 | +6.86 |
| `LVESV` | -26.93 | -14.69 |
| `LVEF` | +19.03 | +15.45 |
| `RVEDV` | +20.98 | +24.85 |

The remaining rejection is not primarily a pulmonary issue. It reflects an
unresolved systemic output and chamber-volume trade-off: `CO_Lmin` is still low
while `RVEDV` remains high and several active parameters sit near plausible
bounds. Further RMSE improvement should therefore be treated as a model
identifiability/preload-representation problem, not just another weight-tuning
problem.

## Purpose

This note documents the calibration cleanup made after Pak Dipo's advice on the Reyna case. The main change was not another round of blind tuning. The cleanup corrected how cardiac output is interpreted, treated systemic pressure-flow quantities as a coupled physiological target, and preserved the best scientific candidate even when the strict operational classifier rejected it.

The most important result is:

- Baseline RMSE: `0.2513`
- Best scientific candidate RMSE: `0.1195`
- Improvement from baseline: `52.4%`
- Operational accepted result in the run manifest: rolled back to baseline because the strict classifier rejected the candidate

This means the `validation_best_candidate_pre_surgery.csv` table is the scientifically informative result, while the manifest-level accepted calibration is conservative rollback behavior.

## Pak Dipo Advice and Implemented Changes

### 1. CO definition was locked to systemic flow

Pak Dipo's key point was that the clinical cardiac output value should be compared to effective systemic cardiac output, not raw LV outflow.

Implemented interpretation:

```text
clinical CO_Lmin = 3.423 L/min
model comparator = Qs_Lmin
not LVCO_Lmin
```

Rationale: in a VSD model, LV outflow can include flow that is not equivalent to effective systemic delivery. Therefore, systemic flow `Qs_Lmin` is the correct comparator for clinical CO in this case.

The latest CO audit supports this definition:

| Quantity | Best candidate value |
|---|---:|
| `Qs_Lmin` | 2.8520 L/min |
| `Qp_Lmin` | 3.2612 L/min |
| `Qao_Lmin` | 2.8522 L/min |
| `LVCO_Lmin` | 3.1266 L/min |
| `QpQs` | 1.1435 |
| Clinical CO target | 3.4230 L/min |
| CO comparator | `Qs_Lmin` |
| CO uncertainty | 0.5000 L/min |

### 2. The systemic quartet was treated as one coupled target

The systemic variables are physiologically linked:

```text
SVR = (SAP_mean - RAP_mean) / Qs
```

Therefore, `SAP_mean`, `RAP_mean`, `CO_Lmin`, and `SVR` should not be treated as four fully independent objectives. Optimizing them independently can double-count the same pressure-flow relationship and make the model improve pressure while damaging flow, or improve flow while producing unrealistic resistance.

Implemented direction:

- Keep `SAP_mean`, `RAP_mean`, `CO_Lmin`, and `SVR` coupled in the objective.
- Avoid rewarding pressure improvement if the model only gets there by suppressing systemic flow.
- Add guards against low-flow/high-SVR solutions.
- Preserve `SVR` and `CO_Lmin` as scientific validation quantities, not just tuning targets.

### 3. The best candidate was kept scientifically even when rejected operationally

The strict classifier rejected the best candidate because not all gates passed. That is useful operational caution, but it should not erase the scientific information from the run.

Interpretation:

- The accepted operational result remains conservative.
- The best scientific candidate shows the current model can substantially improve Reyna's hemodynamics.
- The remaining bottleneck is systemic pressure-flow consistency and ventricular volume plausibility, not a CO accounting bug.

## Best Scientific Candidate Result

Source artifact:

```text
results/runs/20260510_225025_reyna_pre_surgery/tables/validation_best_candidate_pre_surgery.csv
```

| Metric | Clinical | Model | Error % |
|---|---:|---:|---:|
| `RAP_mean` | 5.000 | 4.302 | -13.96 |
| `LAP_mean` | 8.000 | 7.330 | -8.38 |
| `PAP_min` | 10.000 | 11.045 | 10.45 |
| `PAP_max` | 20.000 | 20.730 | 3.65 |
| `PAP_mean` | 15.000 | 15.399 | 2.66 |
| `SAP_min` | 57.000 | 52.216 | -8.39 |
| `SAP_max` | 100.000 | 85.744 | -14.26 |
| `SAP_mean` | 71.300 | 68.038 | -4.57 |
| `QpQs` | 1.194 | 1.144 | -4.23 |
| `SVR` | 19.370 | 22.348 | 15.37 |
| `CO_Lmin` | 3.423 | 2.852 | -16.68 |
| `LVEDV` | 41.000 | 42.692 | 4.13 |
| `LVESV` | 19.300 | 16.418 | -14.93 |
| `RVEDV` | 30.500 | 37.654 | 23.46 |
| `RVESV` | 12.000 | 12.027 | 0.23 |
| `LVEF` | 0.528 | 0.615 | 16.56 |

Primary 5% target gate behavior:

| Metric | Clinical | Model | Absolute error % | Pass 5% |
|---|---:|---:|---:|---|
| `RAP_mean` | 5.000 | 4.302 | 13.96 | No |
| `PAP_max` | 20.000 | 20.730 | 3.65 | Yes |
| `PAP_mean` | 15.000 | 15.399 | 2.66 | Yes |
| `SAP_mean` | 71.300 | 68.038 | 4.57 | Yes |
| `QpQs` | 1.194 | 1.144 | 4.23 | Yes |

This is a strong improvement over the earlier behavior: pulmonary pressure, mean systemic pressure, and shunt ratio now land inside the 5% primary gate. The remaining primary miss is `RAP_mean`.

## Active Parameter Bounds and Plausibility

Sources:

```text
results/runs/20260510_225025_reyna_pre_surgery/tables/parameter_registry_active_pre_surgery.csv
results/runs/20260510_225025_reyna_pre_surgery/tables/parameter_plausibility_best_candidate_pre_surgery.csv
```

| Parameter | Unit | Fitted value | Lower bound | Upper bound | Ratio to baseline | Plausibility |
|---|---:|---:|---:|---:|---:|---|
| `group.R_sys_scale` | ratio | 0.5021 | 0.2500 | 2.8000 | 0.5021 | WARNING |
| `R.SVEN` | mmHg*s/mL | 0.1886 | 0.0588 | 0.4119 | 1.2817 | OK |
| `group.R_pul_scale` | ratio | 0.5526 | 0.4500 | 2.8000 | 0.5526 | WARNING |
| `C.SAR` | mL/mmHg | 0.5306 | 0.3785 | 0.6813 | 1.3015 | OK |
| `E.LV.EA` | mmHg/mL | 7.3484 | 7.0031 | 30.7620 | 0.5255 | WARNING |
| `E.LV.EB` | mmHg/mL | 0.2053 | 0.0840 | 0.3501 | 1.4661 | OK |
| `E.RV.EA` | mmHg/mL | 5.3484 | 1.3531 | 5.9045 | 2.1740 | OK |
| `E.RV.EB` | mmHg/mL | 0.1047 | 0.1041 | 0.4920 | 0.5532 | WARNING |
| `V0.LV` | mL | 5.7928 | 3.8238 | 7.1377 | 1.1362 | OK |
| `V0.RV` | mL | 8.4036 | 8.0524 | 18.3530 | 0.6181 | WARNING |

Summary:

- No active fitted parameter is outside its registry bounds.
- Several parameters carry plausibility warnings because they moved substantially from the baseline or sit close to a bound.
- The most important warnings are systemic resistance scale, pulmonary resistance scale, LV systolic elastance, RV diastolic elastance, and RV unstressed volume.

This means the candidate is physiologically possible inside the current registry, but it is not yet clean enough to call fully robust without explanation.

## Scientific Interpretation

The best candidate is a real step forward. It no longer looks like the model is merely fixing pressure by exploding systemic resistance. Compared with the previous problematic result, systemic flow recovered meaningfully and SVR became much less excessive:

| Metric | Earlier problematic pattern | Best scientific candidate |
|---|---:|---:|
| `CO_Lmin` | 2.174 L/min | 2.852 L/min |
| `SVR` | 27.575 WU | 22.348 WU |
| `SAP_mean` | 65.081 mmHg | 68.038 mmHg |
| `RAP_mean` | 5.129 mmHg | 4.302 mmHg |
| RMSE | 0.2001 | 0.1195 |

The model has moved in the desired physiological direction:

- Mean pulmonary pressures are close to clinical data.
- Mean systemic arterial pressure is within 5%.
- Qp/Qs is within 5%.
- Systemic flow is still low, but much less low than before.
- SVR is still high, but no longer exploding.

The remaining scientific problem is not "which CO variable should be used." That part is now resolved. The remaining problem is patient-specific pressure-flow-volume balance.

## Remaining Weak Points

The candidate should be presented carefully because several clinically relevant metrics remain outside ideal tolerance:

| Metric | Error % | Interpretation |
|---|---:|---|
| `RVEDV` | 23.46 | RV filling volume is too high |
| `CO_Lmin` | 16.68 | Effective systemic output is still low |
| `LVEF` | 16.56 | LV ejection fraction is too high |
| `SVR` | 15.37 | Systemic resistance is still high |
| `SAP_max` | 14.26 | Systolic systemic pressure remains low |
| `LVESV` | 14.93 | LV end-systolic volume is low |
| `RAP_mean` | 13.96 | Right atrial mean pressure is low |

For publication-quality calibration, the next refinement should be targeted, not broad. The active target set is now:

```text
CO_Lmin, SVR, SAP_mean, RAP_mean, RVEDV, LVESV, LVEF
```

Pulmonary waveform polishing is intentionally reduced. The preserved pulmonary/shunt targets are:

```text
PAP_mean, QpQs
```

## Publication Readiness

This result is not yet enough to claim a clinically validated patient-specific digital twin. It is closer to a publishable calibration workflow or model-development case study.

Publishable framing:

- A VSD hemodynamic model with explicit clinical-to-model CO mapping.
- A coupled systemic pressure-flow objective to avoid double-counting.
- A transparent distinction between best scientific candidate and strict accepted candidate.
- A documented preschool BSA extrapolation limitation.
- A reproducible run artifact with validation and parameter plausibility tables.

Not yet publishable as:

- A clinically validated predictive model.
- A final patient-specific surgical decision model.
- A parameter-identifiable physiology model without sensitivity/uncertainty follow-up.

## Recommended Next Step

Do not chase general RMSE blindly. The next step should be a targeted polish that tries to reduce the remaining left-tail errors while preserving the four already-good primary metrics.

Recommended guardrails for the next run:

- Keep `CO_Lmin` compared to `Qs_Lmin`.
- Keep the systemic quartet coupled.
- Penalize low systemic flow without allowing SVR inflation.
- Penalize RVEDV overshoot and LVESV/LVEF mismatch together.
- Preserve `PAP_mean` and `QpQs` while prioritizing systemic flow recovery.
- Reject candidates that improve RMSE only by pushing active parameters closer to implausible bounds.

## Reproducibility Artifacts

Use these files as the source of truth for this note:

```text
results/runs/20260510_225025_reyna_pre_surgery/run_manifest.txt
results/runs/20260510_225025_reyna_pre_surgery/tables/validation_best_candidate_pre_surgery.csv
results/runs/20260510_225025_reyna_pre_surgery/tables/parameter_registry_active_pre_surgery.csv
results/runs/20260510_225025_reyna_pre_surgery/tables/parameter_plausibility_best_candidate_pre_surgery.csv
results/runs/20260510_225025_reyna_pre_surgery/tables/co_definition_audit_pre_surgery.csv
```
