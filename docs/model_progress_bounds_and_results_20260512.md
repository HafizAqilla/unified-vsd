# Model Progress, Bounds, and Run Results

Date: 2026-05-12
Project: Unified VSD pre-surgery calibration workflow
Focus: progress to date on clinical data mapping, model bounds (`lb` / `ub`), calibration cleanup, and current run results

## Purpose

This note records what has been done so far to make the model more
scientifically defensible, especially:

- where the clinical data came from
- how the calibration bounds are defined
- what changes were made to improve the model
- what the most important runs have shown so far

The goal is to keep one readable progress record that explains both the
scientific reasoning and the operational results.

## Data Sources Used So Far

### 1. Reyna

Primary sources:

- `Filled Protokol_VSD Pre_Reyna (2).pdf`
- `docs/clinical_data_dictionary.md`
- `docs/pak_dipo_co_objective_memo.md`
- `docs/reyna_pak_dipo_systemic_flow_revision.md`

Key decisions:

- Active anthropometry uses the raw protocol values:
  - `weight_kg = 13.4`
  - `height_cm = 95.0`
  - `BSA = 0.588 m^2`
- Clinical `CO_Lmin = 3.423` is compared to model `Qs_Lmin`, not `LVCO_Lmin`
- Pulmonary calibration is intentionally reduced to `PAP_mean` only
- `PAP_sys_mmHg` and `PAP_dia_mmHg` are set to `NaN` in the active Reyna file

Why:

- This keeps Reyna strictly patient-specific
- It avoids spending optimization leverage on pulmonary waveform shape that is
  less trustworthy than the mean pressure
- It follows Pak Dipo's advice that the main remaining issue is systemic
  preload-flow-volume consistency

### 2. Razka

Primary sources:

- `config/patient_profile_Razka.m`
- procedure-log style catheter values embedded in the profile comments

Key characteristics:

- Real patient data
- Good cath pressure and `QpQs` information
- No measured `CO_Lmin`
- No ventricular volume block

Why it matters:

- Razka is useful for testing the sparse real-data path
- Razka is not the right case to test the Reyna CO-volume cleanup directly

### 3. Profile A and Profile B

Primary sources:

- `config/patient_profile_A.m`
- `config/patient_profile_B.m`
- synthetic benchmark comments embedded in those files

Key characteristics:

- These are synthetic benchmark profiles, not real patient cases
- They are useful for stress-testing the pipeline
- They should not be treated as scientific validation equal to Reyna or Razka

Important caution:

- Profile A in particular contains internal tension between stated `MAP`, `Qs`,
  `SVR`, and volume targets
- Therefore a poor fit on A may reflect benchmark inconsistency rather than a
  true modeling failure

## What We Changed to Make the Model Better

### 1. Locked the CO definition to the correct comparator

We changed the interpretation of clinical cardiac output so that:

```text
clinical CO target -> model Qs_Lmin
not -> model LVCO_Lmin
```

Reason:

- In unrepaired VSD, LV output includes recirculated pulmonary/shunt flow
- Effective systemic delivery is represented by `Qs`
- Comparing Fick CO to LV outflow can reward the wrong physiology

### 2. Treated the systemic quartet as one coupled physiological target

The model now treats these as linked rather than fully independent:

```text
SAP_mean, RAP_mean, CO_Lmin, SVR
```

Reason:

- `SVR = (SAP_mean - RAP_mean) / Qs`
- If these are optimized independently, the solver can improve pressure by
  suppressing flow or inflate flow while creating unrealistic resistance

### 3. Reduced pulmonary waveform pressure matching for Reyna

For Reyna, the active target now keeps:

- `PAP_mean`
- `QpQs`

and intentionally de-emphasizes:

- `PAP_min`
- `PAP_max`

Reason:

- The remaining bottleneck was not pulmonary polishing
- The model was already better at matching shunt and pulmonary mean pressure
  than at producing correct systemic flow and chamber-volume balance

### 4. Added targeted systemic flow / volume guardrails for Reyna

The Reyna-specific systemic-flow profile now:

- strengthens `CO_Lmin`, `SAP_mean`, and `RAP_mean`
- adds stronger penalty against high `RVEDV`
- adds stronger penalty against low `LVESV`
- adds stronger penalty against overly high `LVEF`
- keeps `CO_Lmin` in the preferred primary set

This was implemented in:

- [build_case_calibration_profile.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/src/calibration/build_case_calibration_profile.m:174)
- [objective_calibration.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/src/calibration/objective_calibration.m:442)
- [select_primary_metrics.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/src/utils/select_primary_metrics.m:69)

### 5. Added a regression test for the Reyna profile

The Reyna-specific assumptions are now checked in:

- [test_reyna_systemic_flow_profile.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/tests/test_reyna_systemic_flow_profile.m:1)

The test verifies:

- raw protocol anthropometry remains active
- `PAP_sys` / `PAP_dia` stay omitted
- primary metrics include `CO_Lmin`
- the volume-flow guard is active

## How `lb` and `ub` Are Defined

The current calibration bounds are no longer arbitrary wide boxes. They are
defined centrally in:

- [build_parameter_registry.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/config/build_parameter_registry.m:200)

### General rule

Each calibratable parameter gets:

- a `baseline_scaled` value after pediatric scaling
- a `seeded_value` after clinical seeding
- a `bound_anchor`
- an `lb`
- an `ub`
- a `bound_source`
- a `confidence_level`

This makes every active calibration knob traceable.

### Bound families

#### Coupled vascular resistance scales

Examples:

- `group.R_sys_scale`
- `group.R_pul_scale`

Default form:

```text
lb = 0.25
ub = 2.80
```

Then profile-level tightening can make the lower bound stricter. For Reyna:

- `group.R_sys_scale`: `0.25` to `2.80`
- `group.R_pul_scale`: `0.45` to `2.80`

Source tag:

- `kung2013_bsa_rc_coupled_resistance_prior`

#### VSD orifice coefficient

Example:

- `vsd.Cd`

Form:

```text
lb = 0.20
ub = 1.20
```

Source tag:

- `orifice_discharge_coefficient_literature_plus_expert_box`

#### Arterial compliances

Examples:

- `C.SAR`
- `C.PAR`

These are anchored to pulse-pressure estimates, not broad multiplier boxes.

For Reyna active run:

- `C.SAR`: `0.3785` to `0.6813`

Source tag:

- `windkessel_stroke_volume_over_pulse_pressure_anchor`

#### Resistances

Example:

- `R.SVEN`

Form:

```text
lb = lower_mult * baseline_scaled
ub = upper_mult * baseline_scaled
```

With venous resistances given slightly wider upper room.

For Reyna active run:

- `R.SVEN`: `0.0588` to `0.4119`

Source tag:

- `kung2013_sharp2000_scaled_resistance_prior`

#### Elastances

Examples:

- `E.LV.EA`
- `E.LV.EB`
- `E.RV.EA`
- `E.RV.EB`

These are anchored to pediatric elastance priors, then widened if needed to
include the clinically seeded start.

For Reyna active run:

- `E.LV.EA`: `7.0026` to `30.7613`
- `E.LV.EB`: `0.0840` to `0.3501`
- `E.RV.EA`: `1.3531` to `5.9046`
- `E.RV.EB`: `0.1041` to `0.4920`

Source tag:

- `zhang2019_pediatric_elastance_anchor`

#### Unstressed volumes

Examples:

- `V0.LV`
- `V0.RV`

These are not free-floating broad search terms anymore. They are anchored to
blood-volume and preload consistency.

For Reyna active run:

- `V0.LV`: `3.8237` to `7.1376`
- `V0.RV`: `9.0034` to `18.3538`

Source tag:

- `blood_volume_preload_consistency_prior`

### Why this matters

The model now has two protections:

1. The optimizer cannot drift into physiologically absurd values too easily.
2. We can explain afterward whether a good fit was achieved cleanly or only by
   parking near the edge of the allowable space.

## Reyna Active Bounds and Best Candidate

Source artifact:

- `results/runs/20260511_223816_reyna_pre_surgery/tables/parameter_registry_pre_surgery.csv`
- `results/runs/20260511_223816_reyna_pre_surgery/tables/parameter_plausibility_best_candidate_pre_surgery.csv`

| Parameter | BaselineScaled | LB | UB | Best fit | Ratio to baseline | Plausibility |
|---|---:|---:|---:|---:|---:|---|
| `group.R_sys_scale` | 1.0000 | 0.2500 | 2.8000 | 0.4016 | 0.4016 | WARNING |
| `R.SVEN` | 0.1471 | 0.0588 | 0.4119 | 0.3255 | 2.2126 | OK |
| `group.R_pul_scale` | 1.0000 | 0.4500 | 2.8000 | 0.5447 | 0.5447 | WARNING |
| `C.SAR` | 0.4079 | 0.3785 | 0.6813 | 0.5067 | 1.2423 | OK |
| `E.LV.EA` | 13.9824 | 7.0026 | 30.7613 | 7.0807 | 0.5064 | WARNING |
| `E.LV.EB` | 0.1400 | 0.0840 | 0.3501 | 0.2415 | 1.7246 | OK |
| `E.RV.EA` | 2.4602 | 1.3531 | 5.9046 | 5.6594 | 2.3004 | WARNING |
| `E.RV.EB` | 0.1892 | 0.1041 | 0.4920 | 0.1547 | 0.8174 | OK |
| `V0.LV` | 5.0983 | 3.8237 | 7.1376 | 6.0467 | 1.1860 | OK |
| `V0.RV` | 13.5954 | 9.0034 | 18.3538 | 9.1245 | 0.6711 | WARNING |

Interpretation:

- The best Reyna scientific candidate stays inside all active bounds
- Several active parameters sit close to one edge of the allowed range
- That is why the run is scientifically informative but not yet clean enough
  to claim a finished patient-specific calibration

## Run Results So Far

### Reyna sequence

| Run folder | Main change | Baseline RMSE | Best scientific RMSE | Status |
|---|---|---:|---:|---|
| `20260510_225025_reyna_pre_surgery` | first strong systemic-flow scientific candidate | `0.2513` | `0.1195` | `REJECT`, rollback |
| `20260511_134027_reyna_pre_surgery` | corrected-BSA experiment (`14.0 kg`, `98 cm`, `0.6173 m^2`) | `0.2454` | `0.1087` | `REJECT`, rollback |
| `20260511_203930_reyna_pre_surgery` | raw PDF anthropometry + `PAP_mean` only | `0.2699` | `0.1702` | `REJECT`, rollback |
| `20260511_223816_reyna_pre_surgery` | targeted systemic-flow polish | `0.2699` | `0.1287` | `REJECT`, rollback |

Latest Reyna best-candidate bottleneck metrics:

| Metric | Clinical | Best candidate | Error % |
|---|---:|---:|---:|
| `RAP_mean` | 5.000 | 5.343 | 6.86 |
| `PAP_mean` | 15.000 | 16.049 | 6.99 |
| `SAP_mean` | 71.300 | 63.559 | -10.86 |
| `QpQs` | 1.194 | 1.140 | -4.52 |
| `SVR` | 19.370 | 20.784 | 7.30 |
| `CO_Lmin` | 3.423 | 2.801 | -18.17 |
| `RVEDV` | 30.500 | 38.079 | 24.85 |
| `LVESV` | 19.300 | 16.464 | -14.69 |
| `LVEF` | 0.528 | 0.610 | 15.45 |

Interpretation:

- The model is better than before at systemic consistency
- The remaining problem is still low effective systemic flow with RV volume
  overshoot and high LV ejection bias
- This is no longer mainly a pulmonary waveform problem

### Profile A

Run folder:

- `20260512_070804_patient_profile_a_pre_surgery`

Result summary:

| Quantity | Value |
|---|---:|
| Baseline RMSE | `0.2444` |
| Best candidate RMSE | `0.2951` |
| Status | `REJECT`, rollback |

Main best-candidate errors:

| Metric | Error % |
|---|---:|
| `QpQs` | -62.52 |
| `SVR` | -52.36 |
| `CO_Lmin` | +52.74 |
| `LVEDV` | -38.44 |
| `RVEDV` | -26.50 |
| `LVEF` | -22.29 |

Interpretation:

- Calibration made A worse, not better
- This benchmark is not a good judge of the Reyna improvement strategy
- Profile A likely reflects benchmark inconsistency more than a clean model
  failure

### Profile B

Run folder:

- `20260512_072450_patient_profile_b_pre_surgery`

Result summary:

| Quantity | Value |
|---|---:|
| Baseline RMSE | `0.2093` |
| Best candidate RMSE | `0.1880` |
| Status | `REJECT`, rollback |

Main best-candidate errors:

| Metric | Error % |
|---|---:|
| `RAP_mean` | +38.67 |
| `CO_Lmin` | +38.86 |
| `SVR` | -23.10 |
| `QpQs` | +18.92 |
| `PVR` | -13.70 |
| `SAP_mean` | +10.11 |

Interpretation:

- B improves slightly, but still rejects
- The residual pattern is different from Reyna
- This again argues against promoting the Reyna-specific polish blindly to all
  full-data cases

### Razka

Run folder:

- `20260512_073515_razka_pre_surgery`

Result summary:

| Quantity | Value |
|---|---:|
| Baseline RMSE | `0.2290` |
| Best candidate RMSE | `0.1346` |
| Status | `PROMISING_NEAR_MISS` |
| Rollback applied | `No` |

Main best-candidate errors:

| Metric | Error % |
|---|---:|
| `RAP_mean` | +3.64 |
| `PAP_min` | -14.92 |
| `PAP_mean` | -15.82 |
| `SAP_mean` | +11.63 |
| `QpQs` | -1.57 |

Interpretation:

- Razka is the strongest non-Reyna result so far
- It behaves like a real sparse-cath case, not like a synthetic benchmark
- Because Razka has no clinical CO or ventricular volume targets, it should be
  improved through the sparse-data path, not through Reyna-style systemic
  flow-volume tuning

## Current Scientific Position

The work so far supports the following conclusion:

### What is working

- The calibration framework is much more defensible than before
- `lb` / `ub` are now traceable to clinical seeding and literature-informed
  priors
- Reyna improved meaningfully once CO was mapped to `Qs` and systemic variables
  were treated as one coupled target
- Razka shows that the real-data pipeline can produce a near-miss result
  without rollback

### What is not solved yet

- Reyna still undershoots effective systemic flow
- Reyna still overshoots `RVEDV`
- Some good Reyna fits still lean on boundary-adjacent parameters
- Synthetic A/B are not reliable enough to validate the same fix pathway

### Best next interpretation

- `Reyna` remains the best case for systemic preload-flow-volume cleanup
- `Razka` should be used to improve the sparse real-data path
- `Profile A/B` should remain synthetic stress benchmarks, not core scientific
  validation cases

## Files Most Relevant to This Progress

- [config/patient_reyna.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/config/patient_reyna.m:37)
- [docs/clinical_data_dictionary.md](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/docs/clinical_data_dictionary.md:56)
- [docs/reyna_pak_dipo_systemic_flow_revision.md](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/docs/reyna_pak_dipo_systemic_flow_revision.md:1)
- [config/build_parameter_registry.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/config/build_parameter_registry.m:222)
- [src/calibration/build_case_calibration_profile.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/src/calibration/build_case_calibration_profile.m:174)
- [src/calibration/objective_calibration.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/src/calibration/objective_calibration.m:442)
- [tests/test_reyna_systemic_flow_profile.m](/D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd/tests/test_reyna_systemic_flow_profile.m:1)
