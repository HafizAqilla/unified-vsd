# Reyna and Razka Model Status Update

Date: 2026-05-20
Branch for latest runs: `codex/baseline-integration`
Project root: `D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd`

## Executive Summary

This document captures the current state of the Unified VSD model for the two cases we are actively using as anchors:

- `Reyna` as the richer pre-op adaptive case
- `Razka` as the sparse-catheterization case

It also records what we improved in the model over the previous sessions, what modes and workflows now exist, what data are actually available for each case, what are derived versus directly reported, and how the latest runs behaved under the current baseline-integration branch.

The main headline from the latest branch is:

- `Razka` improved materially and now lands in `PROMISING_NEAR_MISS` with calibrated RMSE `0.0830` from baseline `0.1550`.
- `Reyna` does not currently survive the governance checks on this branch. A scientifically better candidate was found (`best candidate RMSE 0.0745`), but it was rejected and rolled back because the clinical consistency and target-governance checks failed. The accepted result therefore stayed at the baseline RMSE `0.2701`.

So, at the moment, the new integrated baseline looks favorable for `Razka`, but not yet safe for `Reyna` in the accepted scientific workflow.

## Scope of This Update

This report covers:

1. What we improved across recent sessions
2. What model modes and calibration modes currently exist
3. What baseline is currently active on this branch
4. Latest run results for `Reyna` and `Razka`
5. Data availability tables for both cases
6. Plausibility bounds and active calibration parameters
7. What is available, what is missing, and what still needs to be fixed

## Branch and Commit Context

### Current branch used for latest runs

- Branch: `codex/baseline-integration`
- Commit: `057f05a` - `Integrate Keisya baseline reference parameters`

### Relevant previous branch from the prior session

- Branch: `codex/reyna-shunt-flow-finetune`
- Commit: `3c515c3` - `Add Reyna hemodynamic shunt-flow finetune`

### Main branch head when this update was written

- `main` at `fe1f17b` - `Document systemic-flow calibration governance in README`

## What We Improved So Far

## 1. Core workflow and reproducibility improvements

These landed before the current baseline branch and are part of the current working model history:

| Commit | Improvement |
|---|---|
| `1818f4d` | Cleaned MATLAB startup path and aligned systemic output metrics |
| `70a89ba` | Routed simulation outputs into per-run patient archive folders |
| `e87aef5` | Added VSD shunt calibration diagnostics and regression helpers |
| `84f3951` | Improved preload initialization and added shunt-collapse guards |
| `1016511` | Updated Reyna anthropometry baseline and documented systemic CO interpretation |
| `1c54933` | Added systemic-flow calibration workflow and excluded generated outputs |
| `e58a1ee` | Merged `codex/reyna-systemic-flow-consistency` into `main` |
| `fe1f17b` | Documented systemic-flow calibration governance in README |

These changes are important because they made the pipeline much more reviewable and auditable:

- run folders are now timestamped and self-contained
- calibration behavior is easier to inspect through manifests and CSV tables
- systemic-flow consistency is treated explicitly rather than informally
- preload seeding and collapse-prone shunt states are less brittle

## 2. Reyna-specific finetune improvements from the previous session

On `codex/reyna-shunt-flow-finetune`, we introduced a more evidence-aware Reyna workflow:

- excluded H+1 post-op echo volume/function values from pre-op Reyna calibration
- added `Q_shunt_Lmin` as a derived soft target from `Qp - Qs`
- tightened sparse-cath/adaptive governance so chamber mechanics do not drift without chamber evidence
- added dedicated finetune scripts:
  - `scripts/run_reyna_shunt_flow_finetune.m`
  - `scripts/run_reyna_pulse_pressure_20pct_finetune.m`

### Reyna shunt-flow finetune result from previous session

Run folder:
`results/runs/20260520_210556_reyna_shunt_flow_finetune`

| Quantity | Before | After |
|---|---:|---:|
| Direct RMSE | `0.0471` | `0.0565` |
| Primary RMSE | `0.0776` | `0.0675` |
| Q_shunt error | `-16.37%` | `-0.15%` |
| Physiology validity | passed | passed |

Interpretation:

- this branch intentionally accepted a small direct-fit penalty
- in return, it nearly eliminated the shunt-flow mismatch
- for Reyna specifically, that was a meaningful physiological improvement

## 3. Baseline integration improvements on the current branch

On `codex/baseline-integration`, we integrated Keisya's new baseline package into the codebase:

- updated `config/default_parameters.m`
- added `config/baseline_reference_metrics.m`
- updated `config/build_baseline_provenance.m`
- updated baseline tests to respect the new relaxation timing behavior

### Main baseline changes introduced

Compared with the older main-branch baseline, the integrated Keisya baseline uses:

- smaller unstressed chamber volumes
- lower adult LV systolic elastance anchor (`E.LV.EA = 3.5`)
- higher `C.SAR = 1.33`
- higher `R.SC = 0.80`
- longer LV/RV relaxation fractions (`Tr_LV_frac = 0.40`, `Tr_RV_frac = 0.40`)
- lower RV initial volume seed (`120 mL`)

### Baseline integration validation that passed

- `test_baseline_reference_metrics.m`
- `test_baseline.m`
- `test_scaling_modes.m`
- `test_baseline_provenance_and_age_validity.m`
- `test_parameter_registry.m`
- `test_vascular_v0_blood_volume_consistency.m`
- `checkcode` on touched MATLAB files: `0` issues
- `git diff --check`: clean

## Model Modes and Workflows Currently Available

## 1. Clinical scenario modes

The code currently supports two top-level scenarios:

| Scenario | Meaning |
|---|---|
| `pre_surgery` | unrepaired VSD physiology |
| `post_surgery` | post-closure / post-repair physiology |

`main_run` supports both:

```matlab
main_run('pre_surgery', clinical)
main_run('post_surgery', clinical)
```

## 2. Scaling modes

Currently supported scaling modes:

| Scaling mode | Status |
|---|---|
| `lundquist_bsa` | current default |
| `zhang` | still available for comparison |

The latest `Reyna` and `Razka` runs on this branch used:

- `ScalingMode: lundquist_bsa`

### Scaling comparison snapshot from latest runs

For small preschool patients, the two scaling modes differ a lot in effective factors. Example:

| Quantity | Reyna Zhang | Reyna Lundquist | Razka Zhang | Razka Lundquist |
|---|---:|---:|---:|---:|
| `HR_bpm` | 121.55 | 105.38 | 120.78 | 104.73 |
| `R_SAR_factor` | 2.148 | 2.802 | 2.126 | 2.750 |
| `C_SAR_factor` | 0.200 | 0.357 | 0.204 | 0.364 |
| `E_LV_EA_factor` | 2.236 | 2.802 | 2.212 | 2.750 |
| `V0_LV_factor` | 0.276 | 0.357 | 0.281 | 0.364 |
| `V0_RV_factor` | 0.276 | 0.357 | 0.281 | 0.364 |

Interpretation:

- `lundquist_bsa` scales resistances and elastances upward more aggressively
- `zhang` shrinks compliance more aggressively
- this branch is therefore sensitive to which baseline is plugged into the scaling law

## 3. Calibration case modes

The case-governance layer currently distinguishes:

| Calibration mode | Meaning |
|---|---|
| `synthetic_benchmark` | synthetic benchmark only, not patient-specific inference |
| `sparse_cath` | sparse catheter case, restricted free-parameter set |
| `adaptive_patient` | richer real patient case with evidence-aware target and parameter selection |

Latest runs:

- `Reyna`: `adaptive_patient`
- `Razka`: `sparse_cath`

## 4. VSD modes

The code currently supports at least:

| VSD mode | Meaning |
|---|---|
| `linear_bidirectional` | symmetric linear resistance shunt |
| `orifice_bidirectional` | reduced-order orifice law with signed flow |

Current patient settings:

- `Reyna`: `orifice_bidirectional`
- `Razka`: `linear_bidirectional`

## 5. Runtime pipeline modes and toggles

`main_run` exposes these practical run controls:

| Toggle | Purpose |
|---|---|
| `DO_PLOTS` | generate hemodynamic figures |
| `DO_OVERLAY` | overlay baseline vs calibrated traces |
| `DO_GSA` | run pre/post calibration GSA |
| `DO_FAST_CALIBRATION` | faster calibration path |
| `DO_PARALLEL_FMINCON` | optional parallel optimization |
| `USE_PCE_IN_CALIBRATION` | optional PCE-assisted calibration path |

For the latest Reyna and Razka runs in this document:

- `DO_GSA = false`
- `DO_PLOTS = false`

So the latest results are direct calibration runs without the PCE GSA overhead.

## 6. Calibration stages currently implemented

The current codebase uses a staged calibration design:

- baseline scaling and clinical mapping
- baseline simulation
- optional pre-calibration GSA
- masked calibration
- systemic polish for supported cases
- plausibility polish
- validation and governance
- optional final GSA

For `adaptive_patient` cases like Reyna, the logic also includes:

- target tier governance
- validation holdouts
- plausibility penalties
- rollback when the best numerical fit is not scientifically acceptable

## Latest Run Summary

## 1. Reyna latest run

Run folder:
`results/runs/20260520_223546_reyna_pre_surgery`

Manifest highlights:

- `ScalingMode: lundquist_bsa`
- `CalibrationCaseMode: adaptive_patient`
- `TargetGovernance: direct_measurements_fit_derived_values_check`
- `ConsistencyOnlyTargets: RVEDV,LVEDV,LVESV`
- `ClinicalConsistencyAudit severity: critical`
- `CalibrationStatus: REJECT`
- `RollbackApplied: 1`

### Reyna RMSE summary

| RMSE type | Baseline | Accepted calibrated |
|---|---:|---:|
| Primary governed | `0.270094` | `0.270094` |
| Full transparent | `0.259541` | `0.259541` |
| Hard only | `0.228772` | `0.228772` |
| Soft only | `0.244781` | `0.244781` |

### Reyna best candidate versus accepted result

| Item | Value |
|---|---:|
| Best candidate RMSE | `0.074543` |
| Best candidate status | `REJECT` |
| Accepted result RMSE | `0.270094` |
| Rollback applied | `yes` |

Interpretation:

- the optimizer did find a much better numerical fit
- the scientific governance did not accept it
- because of rollback, the official accepted candidate remained the baseline-like solution

### Why Reyna was rejected on this branch

Manifest summary:

- `REJECT | primary_fail=4 | secondary_fail=2 | plausibility_warning=6 | plausibility_fail=0`

This means:

- there was no hard out-of-bounds implausibility failure
- the problem was not simple LB/UB violation
- the problem was failure of the target-governance and consistency logic

### Reyna primary gate metrics in accepted result

| Metric | Clinical | Model | Abs. error |
|---|---:|---:|---:|
| `RAP_mean` | `5` | `4.220` | `15.60%` |
| `PAP_mean` | `15` | `18.510` | `23.40%` |
| `SAP_mean` | `71.3` | `50.180` | `29.62%` |
| `QpQs` | `1.194` | `1.162` | `2.67%` |
| `CO_Lmin` | `3.423` | `2.373` | `30.68%` |

### Reyna clinical consistency audit

Stroke volume cross-checks from the run:

| Estimate | Value | Formula |
|---|---:|---|
| `SV_Qs` | `28.76 mL/beat` | `CO * 1000 / HR` |
| `SV_Qp` | `34.35 mL/beat` | `CO * QpQs * 1000 / HR` |
| `SV_LV` | `21.70 mL/beat` | `LVEDV - LVESV` |
| `SV_RV` | `18.50 mL/beat` | `RVEDV - RVESV` |

This mismatch is exactly why Reyna is delicate:

- catheter/systemic flow implies one stroke-volume scale
- echo-derived chamber volumes imply a meaningfully smaller one
- the pipeline now treats that as a scientific inconsistency, not just a fitting inconvenience

### Reyna accepted validation table

| Metric | Clinical | Accepted model | Error | Tier |
|---|---:|---:|---:|---|
| `RAP_mean` | `5` | `4.220` | `-15.60%` | `hard` |
| `LAP_mean` | `8` | `6.918` | `-13.52%` | `validation_only` |
| `PAP_min` | `10` | `15.946` | `59.46%` | `validation_only` |
| `PAP_max` | `20` | `20.948` | `4.74%` | `validation_only` |
| `PAP_mean` | `15` | `18.510` | `23.40%` | `hard` |
| `SAP_min` | `57` | `39.180` | `-31.26%` | `soft` |
| `SAP_max` | `100` | `62.419` | `-37.58%` | `soft` |
| `SAP_mean` | `71.3` | `50.180` | `-29.62%` | `hard` |
| `QpQs` | `1.194` | `1.162` | `-2.67%` | `hard` |
| `SVR` | `19.37` | `19.369` | `-0.006%` | `soft` |
| `CO_Lmin` | `3.423` | `2.373` | `-30.68%` | `hard` |
| `LVEDV` | `41` | `35.624` | `-13.11%` | `consistency_check_only` |
| `LVESV` | `19.3` | `13.078` | `-32.24%` | `consistency_check_only` |
| `RVEDV` | `30.5` | `34.096` | `11.79%` | `consistency_check_only` |
| `RVESV` | `12` | `12.318` | `2.65%` | `soft` |

### Reyna best candidate CO audit

| Candidate | `Qs_Lmin` | `Qp_Lmin` | `Qvsd_Lmin` | `QpQs` | `SAP_mean` | `RAP_mean` |
|---|---:|---:|---:|---:|---:|---:|
| baseline | `2.373` | `2.757` | `0.385` | `1.162` | `50.18` | `4.22` |
| scientific candidate | `3.469` | `3.974` | `0.505` | `1.145` | `71.48` | `5.25` |
| accepted candidate | `2.373` | `2.757` | `0.385` | `1.162` | `50.18` | `4.22` |

Interpretation:

- the rejected scientific candidate actually restored systemic flow and systemic pressure much better
- but the candidate still violated the governance logic enough to be rejected
- the accepted result is therefore not the best fit, but the best fit that satisfied the branch's acceptance rules

## 2. Razka latest run

Run folder:
`results/runs/20260520_225511_razka_pre_surgery`

Manifest highlights:

- `ScalingMode: lundquist_bsa`
- `CalibrationCaseMode: sparse_cath`
- `ClinicalConsistencyAudit: not_evaluable`
- `CalibrationStatus: PROMISING_NEAR_MISS`
- `RollbackApplied: 0`

### Razka RMSE summary

| RMSE type | Baseline | Calibrated |
|---|---:|---:|
| Primary governed | `0.155049` | `0.083017` |
| Full transparent | `0.155049` | `0.083017` |
| Hard only | `0.121241` | `0.042690` |
| Soft only | `0.161491` | `0.053866` |

### Razka primary gate metrics

| Metric | Clinical | Model | Abs. error |
|---|---:|---:|---:|
| `RAP_mean` | `5` | `4.675` | `6.49%` |
| `PAP_mean` | `15` | `15.624` | `4.16%` |
| `SAP_mean` | `70` | `72.552` | `3.65%` |
| `QpQs` | `1.21` | `1.215` | `0.39%` |

Interpretation:

- three of four primary anchors are already within 5%
- the only miss is `RAP_mean`, and even that is close
- this is why the status is `PROMISING_NEAR_MISS` rather than reject

### Razka accepted validation table

| Metric | Clinical | Accepted model | Error | Tier |
|---|---:|---:|---:|---|
| `RAP_mean` | `5` | `4.675` | `-6.49%` | `hard` |
| `PAP_min` | `10` | `11.785` | `17.85%` | `validation_only` |
| `PAP_max` | `22` | `19.779` | `-10.10%` | `validation_only` |
| `PAP_mean` | `15` | `15.624` | `4.16%` | `hard` |
| `SAP_min` | `59` | `58.172` | `-1.40%` | `soft` |
| `SAP_max` | `82` | `88.140` | `7.49%` | `soft` |
| `SAP_mean` | `70` | `72.552` | `3.65%` | `hard` |
| `QpQs` | `1.21` | `1.215` | `0.39%` | `hard` |

### Razka CO audit

| Candidate | `Qs_Lmin` | `Qp_Lmin` | `Qvsd_Lmin` | `QpQs` | `SAP_mean` | `RAP_mean` |
|---|---:|---:|---:|---:|---:|---:|
| baseline | `1.868` | `2.472` | `0.604` | `1.323` | `81.63` | `4.55` |
| scientific candidate | `2.018` | `2.451` | `0.433` | `1.215` | `72.55` | `4.68` |
| accepted candidate | `2.018` | `2.451` | `0.433` | `1.215` | `72.55` | `4.68` |

Interpretation:

- calibration pulled Razka toward the reported shunt ratio and systemic pressure cleanly
- there was no rollback
- the sparse-case governance seems to be doing what we want here

## Data Availability Tables

The goal here is to clearly separate:

- `Direct`: directly reported or directly measured
- `Derived`: back-calculated from reported values
- `Model setting`: a modeling choice, not clinical evidence
- `Unavailable`: currently absent

## 1. Reyna data availability

### Common data

| Field | Value | Status | Notes |
|---|---:|---|---|
| `patient_name` | `reyna` | model setting | reproducibility only |
| `age_years` | `3.17` | direct | protocol-derived age |
| `weight_kg` | `14.0` | direct | anthropometry revision |
| `height_cm` | `98.0` | direct | anthropometry revision |
| `sex` | `0` | direct | female code |
| `BSA` | `0.61734` | derived | Mosteller |
| `HR` | `119 bpm` | direct | protocol |

### Reyna pre-surgery data

| Field | Value | Status | Available | Notes |
|---|---:|---|---|---|
| `VSD_diameter_mm` | `3.025` | derived | yes | RV-side midpoint |
| `VSD_gradient_mmHg` | `69` | derived | yes | LV-RV systolic difference |
| `Q_shunt_Lmin` | `0.664` | derived | yes | `Qp - Qs` |
| `QpQs` | `1.194` | direct | yes | protocol |
| `VSD_mode` | `orifice_bidirectional` | model setting | yes | model choice |
| `PAP_sys_mmHg` | `20` | direct | yes | cath |
| `PAP_dia_mmHg` | `10` | direct | yes | cath |
| `PAP_mean_mmHg` | `15` | direct | yes | cath |
| `PVR_WU` | `NaN` | unavailable | no | not calculated |
| `SAP_sys_mmHg` | `100` | direct | yes | RFA cath |
| `SAP_dia_mmHg` | `57` | direct | yes | RFA cath |
| `SAP_mean_mmHg` | `71.3` | derived | yes | from sys/dia |
| `SVR_WU` | `19.37` | derived | yes | from MAP, RAP, Qs |
| `RAP_mean_mmHg` | `5` | direct | yes | cath |
| `LAP_mean_mmHg` | `8` | estimated / treated as direct input | yes | not directly measured |
| `LVEDP_mmHg` | `8` | estimated / treated as direct input | yes | not directly measured |
| `LVEDV_mL` | `41.0` | derived | yes | Teichholz from M-mode |
| `LVESV_mL` | `19.3` | derived | yes | Teichholz from M-mode |
| `RVEDV_mL` | `30.5` | direct | yes | protocol |
| `RVESV_mL` | `12.0` | direct | yes | protocol |
| `LVEF` | `0.528` | derived | yes | from LV volumes |
| `override_IC` | `true` | model setting | yes | calibration helper |
| `CO_comparator` | `Qs_Lmin` | model setting | yes | workflow control |
| `CO_uncertainty_Lmin` | `0.5` | model setting | yes | workflow control |
| `CO_Lmin` | `3.423` | derived | yes | Fick Qs |

### Reyna post-surgery data currently in mainline profile

| Field block | Status |
|---|---|
| post-op pressures | unavailable |
| post-op volumes | unavailable |
| post-op function | unavailable |
| post-op CO | unavailable |
| post-op Qp/Qs | unavailable |

Note: the current branch still has post-op fields largely empty in `config/patient_reyna.m`. The separate warm-start and post-op prediction workflow was discussed and partially developed, but the present `main`-compatible patient profile is still mostly `NaN` for post-op validation.

## 2. Razka data availability

### Common data

| Field | Value | Status | Notes |
|---|---:|---|---|
| `patient_name` | `razka` | model setting | reproducibility only |
| `patient_id` | `00948048` | model setting | identifier only |
| `age_years` | `4.64` | direct/derived from DOB | procedure date based |
| `weight_kg` | `14.3` | direct | measured |
| `height_cm` | `100.4` | direct | measured |
| `BSA` | `0.629` | direct report | stated in log |
| `sex` | `1` | placeholder | needs verification |
| `HR` | `111 bpm` | direct | monitor |
| `maturation_mode` | `normal` | model setting | workflow control |

### Razka pre-surgery data

| Field | Value | Status | Available | Notes |
|---|---:|---|---|---|
| `VSD_diameter_mm` | `2.5` | derived/communicated | yes | midpoint of 2-3 mm |
| `VSD_gradient_mmHg` | `NaN` | unavailable | no | not reported |
| `Q_shunt_Lmin` | `NaN` | unavailable | no | CO not measured |
| `QpQs` | `1.21` | direct | yes | cath report FR |
| `VSD_mode` | `linear_bidirectional` | model setting | yes | model choice |
| `PAP_sys_mmHg` | `22` | direct | yes | cath |
| `PAP_dia_mmHg` | `10` | direct | yes | cath |
| `PAP_mean_mmHg` | `15` | direct | yes | cath |
| `PVR_WU` | `NaN` | unavailable | no | CO not measured |
| `SAP_sys_mmHg` | `82` | direct | yes | DAO cath |
| `SAP_dia_mmHg` | `59` | direct | yes | DAO cath |
| `SAP_mean_mmHg` | `70` | direct report in record | yes | DAO cath |
| `SVR_WU` | `NaN` | unavailable | no | CO not measured |
| `RAP_mean_mmHg` | `5` | direct | yes | cath |
| `LAP_mean_mmHg` | `NaN` | unavailable | no | PCWP not measured |
| `LVEDV_mL` | `NaN` | unavailable | no | not reported |
| `LVESV_mL` | `NaN` | unavailable | no | not reported |
| `RVEDV_mL` | `NaN` | unavailable | no | not reported |
| `RVESV_mL` | `NaN` | unavailable | no | not reported |
| `LVEF` | `NaN` | unavailable | no | not reported |
| `CO_Lmin` | `NaN` | unavailable | no | explicitly not measured |
| `LVEDP_mmHg` | `9` | direct | yes | LV cath |
| `RVEDP_mmHg` | `8` | direct | yes | RV cath |
| `LVP_sys_mmHg` | `91` | direct | yes | LV cath |
| `RVP_sys_mmHg` | `26` | direct | yes | RV cath |
| `PARI_raw` | `0.71` | reported but unclassified | yes | definition not confirmed |
| `override_IC` | `false` | model setting | yes | default |

### Razka post-surgery data

| Field block | Status |
|---|---|
| post-op scenario | not applicable in current file |
| post-op pressures | unavailable |
| post-op volumes | unavailable |
| post-op function | unavailable |
| post-op CO | unavailable |

## Target Governance and Validation Tiers

## 1. Reyna target tier map

| Metric | Tier | Included in calibration | Included in primary RMSE | Note |
|---|---|---:|---:|---|
| `RAP_mean` | `hard` | 1 | 1 | anchor |
| `PAP_mean` | `hard` | 1 | 1 | anchor |
| `SAP_mean` | `hard` | 1 | 1 | anchor |
| `QpQs` | `hard` | 1 | 1 | anchor |
| `CO_Lmin` | `hard` | 1 | 1 | anchor |
| `SAP_min` | `soft` | 1 | 1 | waveform pressure |
| `SAP_max` | `soft` | 1 | 1 | waveform pressure |
| `SVR` | `soft` | 1 | 1 | derived check |
| `RVESV` | `soft` | 1 | 1 | chamber support |
| `LAP_mean` | `validation_only` | 0 | 1 | check only |
| `PAP_min` | `validation_only` | 0 | 1 | check only |
| `PAP_max` | `validation_only` | 0 | 1 | check only |
| `LVEDV` | `consistency_check_only` | 0 | 0 | excluded due inconsistency |
| `LVESV` | `consistency_check_only` | 0 | 0 | excluded due inconsistency |
| `RVEDV` | `consistency_check_only` | 0 | 0 | excluded due inconsistency |

## 2. Razka target tier map

| Metric | Tier | Included in calibration | Included in primary RMSE | Note |
|---|---|---:|---:|---|
| `RAP_mean` | `hard` | 1 | 1 | anchor |
| `PAP_mean` | `hard` | 1 | 1 | anchor |
| `SAP_mean` | `hard` | 1 | 1 | anchor |
| `QpQs` | `hard` | 1 | 1 | anchor |
| `SAP_min` | `soft` | 1 | 1 | waveform pressure |
| `SAP_max` | `soft` | 1 | 1 | waveform pressure |
| `PAP_min` | `validation_only` | 0 | 1 | check only |
| `PAP_max` | `validation_only` | 0 | 1 | check only |
| all chamber volumes | `unavailable` | 0 | 0 | not measured |
| `CO_Lmin` | `unavailable` | 0 | 0 | not measured |
| `SVR` | `unavailable` | 0 | 0 | not measurable from data |
| `PVR` | `unavailable` | 0 | 0 | not measurable from data |

## Plausibility and Bounds

Important point: in the latest results, the plausibility system is not mainly rejecting cases because parameters went outside hard bounds. The harder issue is whether a numerically good fit is scientifically acceptable given the data hierarchy and consistency checks.

## 1. Reyna active calibratable parameters and bounds

| Parameter | Baseline scaled | Seeded | LB | UB | Best candidate | Accepted | Within bounds |
|---|---:|---:|---:|---:|---:|---:|---|
| `group.R_sys_scale` | `1.0000` | `0.4608` | `0.25` | `2.80` | `0.3187` | `0.4608` | yes |
| `R.SVEN` | `0.1401` | `0.0646` | `0.0560` | `0.3923` | `0.3863` | `0.0646` | yes |
| `group.R_pul_scale` | `1.0000` | `1.0000` | `0.45` | `2.80` | `0.4596` | `1.0000` | yes |
| `C.SAR` | `0.4746` | `0.6689` | `0.5017` | `0.9031` | `0.5063` | `0.6689` | yes |
| `C.PAR` | `1.7842` | `3.4345` | `2.4042` | `4.9800` | `3.3354` | `3.4345` | yes |
| `E.LV.EA` | `9.8082` | `7.4501` | `5.8849` | `21.5780` | `6.9379` | `7.4501` | yes |
| `E.LV.EB` | `0.2242` | `0.2304` | `0.1345` | `0.5605` | `0.1346` | `0.2304` | yes |
| `E.RV.EA` | `2.3456` | `2.9349` | `1.2901` | `5.6294` | `2.8545` | `2.9349` | yes |
| `E.RV.EB` | `0.1970` | `0.1223` | `0.1084` | `0.5123` | `0.1090` | `0.1223` | yes |
| `E.LA.EA` | `0.9808` | `0.9808` | `0.4904` | `2.4520` | `1.8938` | `0.9808` | yes |
| `E.RA.EA` | `1.4073` | `1.4073` | `0.7037` | `3.5184` | `1.1638` | `1.4073` | yes |
| `V0.LV` | `1.2627` | `6.2800` | `0.9470` | `6.5940` | `5.2861` | `6.2800` | yes |
| `V0.RV` | `2.9999` | `5.4581` | `2.0999` | `5.7310` | `4.7205` | `5.4581` | yes |
| `vsd.Cd` | `0.7000` | `0.5227` | `0.4182` | `0.6272` | `0.5595` | `0.5227` | yes |

Interpretation:

- the rejected Reyna best candidate stayed inside every hard bound shown here
- the rejection was therefore governance-driven, not bound-driven
- the key conflict is between systemic-flow restoration and chamber/consistency evidence

## 2. Razka active calibratable parameters and bounds

| Parameter | Baseline scaled | Seeded | LB | UB | Accepted | Within bounds |
|---|---:|---:|---:|---:|---:|---|
| `group.R_sys_scale` | `1.0000` | `1.0000` | `0.50` | `2.50` | `0.7994` | yes |
| `R.SVEN` | `0.1375` | `0.1375` | `0.0550` | `0.3438` | `0.1496` | yes |
| `group.R_pul_scale` | `1.0000` | `1.0000` | `0.50` | `2.50` | `0.8961` | yes |
| `R.vsd` | `2.5950` | `2.5950` | `0.6487` | `10.3799` | `3.3827` | yes |

Interpretation:

- Razka is behaving like a well-governed sparse case
- only a few representative knobs are allowed to move
- the fit improves substantially without implausible drift

## What the Latest Results Seem to Mean

## 1. Reyna

The current integrated baseline made Reyna harder to accept scientifically, not easier.

What the latest run suggests:

- the optimizer can recover systemic flow and systemic pressure better than baseline
- but once the branch's consistency rules are applied, the candidate is rejected
- the current integrated baseline therefore does not yet give a trustworthy accepted Reyna solution

The likely reason is not a simple optimization bug. It is a model-evidence conflict:

- Reyna has mixed data sources
- some chamber volumes are derived from Teichholz/M-mode logic
- systemic Fick-style flow points imply a different stroke-volume scale
- the stronger governance now refuses to silently absorb that contradiction

## 2. Razka

Razka is much cleaner on this branch.

Why:

- data are sparse, but internally more coherent for the limited target set
- sparse-cath governance avoids overfitting unsupported chamber mechanics
- the integrated baseline gives a better starting point for the limited hemodynamic target set

So for Razka, the new baseline looks like a real improvement.

## 3. Why the two cases diverge

This is the important scientific point.

The new baseline is not "good" or "bad" in the abstract. It is interacting with:

- scaling law choice
- case governance mode
- whether chamber volumes are trusted as direct fit anchors
- whether CO is direct, derived, or inconsistent with volumes
- whether rollback is enabled when scientific criteria fail

That is why:

- `Razka` improves
- `Reyna` numerically improves but scientifically rejects

## Current Gaps and Open Problems

## 1. Reyna pre-op evidence conflict remains unresolved

Most important unresolved issue:

- Teichholz-derived LV volumes and catheter/systemic-flow targets are still in tension

The current code now exposes that tension more honestly than before, but it is not solved yet.

## 2. Post-op prediction workflow is not yet finalized in the shared branch

We discussed a post-op prediction flow with:

- baseline scaling
- pre-op calibrated parameter import
- post-op warm-start step
- direct surgery-state transition without re-optimization

Conceptually this is in motion, but the shared patient profiles and validated post-op target package are not complete yet in the current branch.

## 3. Preschool scaling remains extrapolated

Both latest runs explicitly state:

- `AgeValidityRegime: preschool_extrapolated`
- `AgeValidityEvidence: low`

So even when the fit is good, we should still report uncertainty because the scaling prior is not directly literature-validated for this age window.

## Recommended Next Steps

## If the goal is to stabilize `Reyna`

1. Keep the new baseline branch separate from the Reyna finetune branch.
2. Reapply the Reyna shunt-flow governance improvements onto the baseline branch.
3. Decide explicitly whether `LAP_mean` and `LVEDP` remain treated as fit-worthy evidence or are downgraded.
4. Decide explicitly whether Teichholz-derived LV volumes remain only consistency checks or get partially reintroduced with broader uncertainty.
5. Compare `zhang` and `lundquist_bsa` again for Reyna only after the above governance decision is fixed.

## If the goal is to move the project forward fast

1. Treat `Razka` as the stronger example for the new integrated baseline.
2. Treat `Reyna` as the stress-test case for mixed-source consistency governance.
3. Preserve the previous Reyna finetune branch as the best shunt-flow-specific refinement branch.
4. Avoid claiming the integrated baseline is globally better until Reyna is reconciled.

## Bottom Line

As of 2026-05-20:

- the codebase is much better organized, more reproducible, and more scientifically explicit than it was a few sessions ago
- the baseline integration branch is a meaningful step forward
- `Razka` looks better on the new baseline
- `Reyna` is not yet acceptable on the new baseline under the current scientific governance, even though the optimizer can find a much lower-RMSE candidate

That is actually useful progress: we are no longer just getting lower numbers, we are now seeing exactly where lower numbers stop being trustworthy.
