# Unified VSD Model Improvement Summary

Generated: 2026-05-21

This document summarizes the practical and scientific improvements made during
the recent Reyna/Razka/Ali work. The goal is to explain what changed, why it
matters physiologically, and what evidence we now have from the best available
runs.

## Executive Summary

The project moved from "can the optimizer make the numbers small?" toward
"can we trust why the numbers are small?" The main improvements are:

- calibration targets are now separated into direct protocol data, derived
  protocol data, unavailable data, and model-only outputs
- Reyna's pre-surgery interpretation is now treated more carefully because
  some volume/function values are suspected to come from post-operative or
  mixed-source records
- Razka is now a cleaner sparse-catheterization validation case with no forced
  CO, PVR, SVR, or volume targets when the record does not support them
- Ali Zhafran was converted into a data-governance/sensitivity case because
  the flow block is internally inconsistent
- best-result reporting now includes RMSE by evidence class, not only one
  aggregate number
- run outputs are easier to audit through manifests, CSV metric tables, and
  patient-specific run folders

The strongest current proof point is that the model can reproduce available
pre-surgery hemodynamic targets for two VSD patients:

| Patient | Best selected run | Best all-comparator RMSE | Hard-target RMSE | Main evidence |
|---|---|---:|---:|---|
| Reyna | `results/runs/20260520_100315_reyna_shunt_flow_finetune` | **0.0624** | **0.0459** | shunt flow, Qp/Qs, mean pressures |
| Razka | `results/runs/20260520_110545_razka_pre_surgery` | **0.0687** | **0.0367** | sparse cath pressures and Qp/Qs |

## 1. Data Governance Improvements

### 1.1 Direct, derived, unavailable, and model-only data are now separated

Before this cleanup, it was too easy to treat every number in a protocol as if
it had the same evidence strength. The current reporting separates:

| Evidence class | Meaning | Example |
|---|---|---|
| Direct protocol | measured/reported directly in the clinical protocol | `PAP_mean`, `SAP_sys`, `QpQs` |
| Derived protocol | calculated from available direct data | `SAP_mean`, `Q_shunt_Lmin`, `CO_Lmin` |
| Unavailable | not present or not reliable enough to compare | Razka `CO_Lmin`, `PVR`, volumes |
| Model-only output | produced by simulation but has no clinical comparator | simulated `LVEF`, `RVEF`, chamber volumes when not measured |

Why this matters:

- RMSE is no longer inflated or improved by unsupported targets.
- Missing clinical data are shown explicitly instead of silently ignored.
- A model-only number is not misrepresented as a validated number.

### 1.2 Clinical consistency checks are now part of the reasoning

The workflow now checks whether flow-derived stroke volumes and volume-derived
stroke volumes agree:

| Quantity | Formula | Interpretation |
|---|---|---|
| `SV_Qs` | `CO_Lmin * 1000 / HR` | systemic stroke volume implied by flow |
| `SV_Qp` | `CO_Lmin * QpQs * 1000 / HR` | pulmonary stroke volume implied by flow ratio |
| `SV_LV` | `LVEDV - LVESV` | LV stroke volume implied by echo volumes |
| `SV_RV` | `RVEDV - RVESV` | RV stroke volume implied by echo volumes |

This is important for Reyna because catheter-derived flow and echo-derived
volumes do not always imply the same stroke-volume scale. The project now
documents that mismatch as a clinical evidence issue instead of hiding it
inside calibration.

### 1.3 Ali Zhafran was downgraded from calibration target to governance case

For Ali, the catheter report is pre-release/pre-occluder, so the pressure data
are relevant to pre-surgery physiology. However, the flow block is inconsistent:

| Item | Issue |
|---|---|
| Protocol `Qp` | reported as `5.718 L/min`, but likely depends on indexed VO2 |
| Protocol `Qs` | reported as `3.15 L/min`, inconsistent with `Qp/Qs = 1.6` |
| Workbook `Qs` | gives about `3.50`, inconsistent with protocol `Qs` |
| PARI | reported, but unit/definition is not explicit |
| HR | cath/Fick line uses `111 bpm`, protocol/ECG line shows `150 bpm` |

Improvement:

- direct pressures and `QpQs` are retained
- ambiguous `CO_Lmin`, `Q_shunt_Lmin`, and `PVR_WU` are held out
- raw worksheet/protocol values are retained as audit fields
- Ali is treated as a sensitivity/governance case, not a primary validation
  proof case

## 2. Reyna Improvements

### 2.1 Shunt-flow finetuning improved the physiologic target that mattered most

Best selected Reyna run:

```text
results/runs/20260520_100315_reyna_shunt_flow_finetune
```

Key improvement:

| Quantity | Before finetune | After finetune |
|---|---:|---:|
| Direct objective RMSE | 0.0471 | **0.0459** |
| Primary/full RMSE | 0.0776 | **0.0624** |
| `Q_shunt_Lmin` error | -16.37% | **-0.34%** |
| Physiology validity | pass | pass |

Interpretation:

- The model now matches Reyna's derived shunt flow almost exactly.
- This is physiologically meaningful because VSD severity is strongly tied to
  the shunt volume and Qp/Qs behavior.
- The remaining main residual is systemic systolic pressure (`SAP_max`), not
  shunt flow.

### 2.2 Reyna is now handled as a mixed-evidence patient

Reyna has richer data than Razka, but richer does not automatically mean safer.
The improvement is that the pipeline now treats mixed-source evidence with more
care:

| Data group | Current handling |
|---|---|
| Mean pressures | high-priority calibration/validation targets |
| `QpQs` | high-priority shunt severity target |
| `Q_shunt_Lmin` | derived soft target |
| `CO_Lmin` | compared as systemic flow `Qs`, not LV outflow |
| Chamber volumes/function | reported carefully; not allowed to dominate when source timing or consistency is doubtful |

Why this matters:

- It prevents the model from overfitting a volume value that may not represent
  the same physiologic state as the catheter data.
- It makes the result easier to defend in a thesis or paper.

## 3. Razka Improvements

### 3.1 Razka now works as a sparse-catheterization validation case

Best selected Razka run:

```text
results/runs/20260520_110545_razka_pre_surgery
```

Main result:

| Quantity | Baseline | Calibrated |
|---|---:|---:|
| All-comparator RMSE | 0.2290 | **0.0687** |
| Hard-target RMSE | 0.1928 | **0.0367** |
| Soft-target RMSE | 0.3443 | **0.0811** |
| Calibration status | baseline | `PROMISING_NEAR_MISS` |
| Rollback applied | N/A | no |

Interpretation:

- Razka improved substantially after calibration.
- The model matches the high-value sparse targets well: mean pressures and
  Qp/Qs.
- The remaining errors are mostly pressure extrema (`PAP_max`, `SAP_min`),
  which are harder to fit without overfitting waveform shape.

### 3.2 Missing Razka data are now explicitly held out

The Razka protocol does not provide enough information to calculate some
derived hemodynamic quantities safely.

| Metric group | Current status |
|---|---|
| CO/Fick | unavailable, not used as target |
| Shunt flow in L/min | unavailable, not used as target |
| PVR/SVR | unavailable because CO is unavailable |
| LV/RV volumes | unavailable |
| EF | unavailable |

Why this matters:

- We can still validate Razka on the data that exist.
- We avoid inventing false precision from missing CO or echo volume data.
- Razka becomes a clean example of sparse patient calibration.

## 4. RMSE Reporting Improvements

A new result summary was added:

```text
docs/reyna_razka_best_error_rmse_summary_20260521.md
```

It reports:

- best run source folders
- error percentage formula
- RMSE formula
- RMSE by evidence class
- RMSE by calibration tier
- detailed metric-level error table for Reyna
- detailed metric-level error table for Razka
- model-only outputs with no clinical comparator

### Evidence-class RMSE snapshot

| Patient | Evidence class | Best RMSE % | Meaning |
|---|---|---:|---|
| Reyna | Protocol direct | 6.47% | direct pressure/QpQs data |
| Reyna | Protocol derived | 5.66% | derived `SAP_mean`, shunt flow, CO |
| Reyna | All with comparator | 6.24% | full transparent comparator set |
| Razka | Protocol direct | 7.30% | direct sparse cath data |
| Razka | Protocol derived | 2.14% | derived `SAP_mean` only |
| Razka | All with comparator | 6.87% | full transparent comparator set |

This is better than reporting only one RMSE because it shows whether good
performance comes from measured data or derived data.

## 5. Calibration Workflow Improvements

### 5.1 Case modes are more explicit

The pipeline now distinguishes different calibration contexts:

| Case mode | Purpose |
|---|---|
| `sparse_cath` | limited catheter record; restrict free parameters |
| `adaptive_patient` | richer patient record; use target governance and holdouts |
| synthetic/benchmark modes | test or reference use, not patient inference |

Why this matters:

- Razka should not be calibrated with the same freedom as a data-rich patient.
- Reyna needs more governance because more data also means more possible
  contradictions.

### 5.2 Target tiers are clearer

The validation/calibration target logic now uses categories such as:

| Tier | Meaning |
|---|---|
| `hard` | primary anchor; should be matched closely |
| `soft` | useful target, but lower certainty or waveform-sensitive |
| `validation_only` | reported for transparency, not optimized directly |
| `consistency_check_only` | used to expose data inconsistency, not fit |
| `unavailable` | no clinical comparator |

This makes the calibration defensible because the optimizer is not allowed to
treat every number as equally reliable.

### 5.3 Plausibility and parameter governance improved

The current workflow produces parameter plausibility tables and active
parameter registries. This helps answer:

- which parameters were allowed to change?
- how far did they move from baseline?
- were they inside bounds?
- did a low-RMSE solution require implausible parameters?

This is especially important for scientific simulation because a low RMSE is
not enough if the fitted physiology becomes unrealistic.

## 6. Reproducibility Improvements

Run outputs are now easier to audit because the project stores:

| Artifact | Purpose |
|---|---|
| `run_manifest.txt` | scenario, patient, scaling mode, RMSE, status, file paths |
| validation CSV tables | clinical vs model values and error percentages |
| parameter plausibility CSVs | fitted parameter values and bound checks |
| MAT parameter packages | reproducible calibrated parameter sets |
| patient-specific run folders | keeps results from different patients separate |

This matters because the result can be traced from:

```text
patient config -> calibration targets -> parameter set -> simulation metrics -> RMSE table
```

## 7. Documentation Improvements

Recent documentation now covers the pieces needed to explain the workflow:

| Document | Purpose |
|---|---|
| `docs/calibration_data_governance_notes.md` | explains target consistency and why some targets are held out |
| `docs/clinical_data_dictionary.md` | maps clinical record fields to MATLAB variables |
| `docs/reyna_razka_best_error_rmse_summary_20260521.md` | best Reyna/Razka RMSE and error tables |
| `docs/reyna_razka_status_update_2026-05-20.md` | broader status update and branch context |
| `docs/model_improvements_20260521.md` | this improvement summary |

The documentation now supports the core scientific claim:

> The model is not only fitted; the target selection, exclusions, and derived
> quantities are explicitly traceable.

## 8. What We Can Claim Now

Reasonable claim:

> The current workflow can reproduce available pre-surgery hemodynamic targets
> for Reyna and Razka with approximately 6-7% transparent all-comparator RMSE,
> while clearly separating measured targets from derived and unavailable data.

Stronger claim for Reyna:

> Reyna's shunt-flow-specific calibration reduced shunt-flow error from about
> `-16.37%` to `-0.34%`.

Stronger claim for Razka:

> Razka's sparse-catheter calibration reduced all-comparator RMSE from `0.2290`
> to `0.0687`, with hard-target RMSE `0.0367`.

What we should not overclaim:

- We should not claim Razka validates CO, PVR, SVR, or ventricular volumes,
  because those comparators are not available.
- We should not claim every Reyna volume/function value is a safe pre-surgery
  calibration target without resolving source timing and consistency.
- We should not treat Ali as a primary validation success until the ambiguous
  flow block is resolved.

## 9. Remaining Limitations

| Limitation | Why it matters | Current mitigation |
|---|---|---|
| Reyna mixed-source volume/flow consistency | may bias chamber mechanics if forced too hard | keep questionable targets out of primary RMSE or treat as consistency checks |
| Razka lacks CO/Fick and echo volumes | cannot validate flow magnitude or chamber volumes | validate only sparse cath pressures and Qp/Qs |
| Ali flow block inconsistent | CO/shunt/PVR cannot be trusted directly | retain raw audit fields, hold out ambiguous targets |
| Preschool scaling uncertainty | literature support is weaker in this age range | report age-validity regime and uncertainty |
| Pressure extrema residuals | systolic/diastolic waveform fit still imperfect | prioritize mean pressures and avoid overfitting extrema |

## 10. Recommended Next Steps

1. Keep Reyna and Razka as the two primary validation examples.
2. Use Ali as a data-governance example until the flow block is clarified.
3. Preserve the evidence-class RMSE tables in future reports.
4. Re-run Reyna only after deciding exactly which volume/function targets are
   valid pre-surgery comparators.
5. For Razka, avoid adding CO/PVR/SVR targets unless a source record provides
   the missing flow measurement.
6. When writing the thesis result section, report both:
   - all-comparator RMSE
   - evidence-class RMSE

## Bottom Line

The main improvement is not just better RMSE. The main improvement is that the
workflow is now more scientifically honest:

- missing data stay missing
- derived data are labeled as derived
- inconsistent data are held out or audited
- model-only outputs are not called validated outputs
- Reyna and Razka now show reproducible, traceable pre-surgery fits with
  roughly `6-7%` transparent RMSE

