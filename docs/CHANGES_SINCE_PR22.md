# Changes Since PR #22

Scope: all work on `main` after commit `251fcb3` (merge of PR #22), including
PR #23 and the current branch `codex/reyna-statistical-calibration`.

Patient: Reyna, 3 years 2 months, female, restrictive VSD.
Model: lumped-parameter cardiovascular model, pre-surgery and post-surgery.

---

## 1. Summary

Work since PR #22 falls into five groups.

| Group | What it did |
|---|---|
| A. Target governance | Made the model fit the same metrics it is graded on |
| B. Statistical methods | Replaced percentage acceptance with uncertainty-based statistics |
| C. Clinical data correction | Found and fixed three wrong input values from the source record |
| D. New evidence | First out-of-sample prediction test; joint pre/post fit |
| E. Bug fixes | Six defects found, four of which produced silently wrong output |

The most important single result: **two model inputs were wrong**, so every
number published for this patient before this work was computed for a
mis-specified child. After correction, the model fits better and now has one
genuine prediction test.

---

## 2. Group A: Target governance

### 2.1 Problem

The model was graded on metrics it was never asked to fit.

`PAP_min` and `PAP_max` are directly measured catheter pressures. They were
marked `UseForCalibration = true` but fell into the `validation_only` tier.
That means they were excluded from the objective function but included in the
reported RMSE. The optimiser was scored on two pressures it never optimised.

The same defect affected `LAP_mean`, `RVEDV` (pre-surgery) and `RVEF`
(post-surgery).

### 2.2 Fix

Seven target-governance defects were closed. `PAP_min` and `PAP_max` now enter
the objective. A test asserts that no metric can be graded without being
fitted.

### 2.3 Second problem: two code paths disagreed

`build_target_tiers(clinical, scenario)` called directly returned 11 governed
rows. The production path returned 10, because it honours
`recipe.primary_rmse_holdout`. Two code paths, two different denominators for
the same RMSE, and nothing checking they agreed.

This was fixed and later caught again in the Phase 4 driver (see §6.2).

---

## 3. Group B: Statistical methods

Four methods were added. All are opt-in or report-only, so no existing result
changed silently.

### 3.1 Sigma-weighted objective

A 10% acceptance gate treats every metric as equally well known. It is not a
statistical criterion.

Example from the data. `PAP_min` was recorded three times as 10/10/10 mmHg, so
its uncertainty is about ±0.5 mmHg. `SAP_min` has an uncertainty of ±5.7 mmHg.
Under a percentage gate, an error of 1.0 mmHg in `PAP_min` passes and an error
of 9.9 mmHg in `SAP_min` fails, even though the first is 2.0 standard
deviations and the second is 1.7.

The objective now optionally normalises each residual by that metric's own
declared uncertainty instead of one global percentage.

- Setting: `calib.objectiveWeighting = 'legacy'` (default) or `'sigma'`.
- Environment override: `UNIFIED_VSD_OBJECTIVE_WEIGHTING`.
- The `legacy` path is proven byte-identical to the pre-change objective by a
  regression test.

### 3.2 Chi-squared reporting

Chi-squared answers whether a fit is as good as the data allows. It cannot be
improved by redefining which metrics count.

Every run now prints and exports:

- `chi2`, `N` (observations), `p` (free parameters), `dof = N - p`
- `chi2/N` and `chi2/dof`
- an interpretation label
- the worst metric by z-score

Chi-squared is recorded on the run status but does not decide ACCEPT or
REJECT. Changing the objective and the acceptance rule in the same release
would make neither attributable.

### 3.3 Parameter identifiability

Global sensitivity analysis ranks parameters one at a time. It cannot detect
that two retained parameters are redundant with each other.

A scaled sensitivity matrix is now built at the calibrated operating point:

```
S(i,j) = (dy_i / dtheta_j) * (theta_j / sigma_i)
```

computed by central finite differences with a 1% relative step. The report
gives the condition number, pairwise column correlations, and per-parameter
column norms. A pair with |rho| > 0.9 is flagged as collinear. A condition
number above 1e3 is flagged as near-dependent.

### 3.4 Validation holdout machinery

Built as generic infrastructure keyed on the `validation_holdout` tier.

The PRD asked for `SVR` to be labelled a holdout. This was not done, and the
reason is recorded. `SVR = (SAP_mean - RAP_mean) / CO`, and all three of those
are already fitted targets. Predicting `SVR` well shows only that three fitted
quantities are mutually consistent. Labelling it a holdout would have created
the appearance of an independent prediction where none existed.

A test records the honest state at that time: no metric in the recipe was a
genuine holdout. This changed later (see §7).

---

## 4. Group C: Clinical data correction

The source record was obtained: RSAB Harapan Kita procedure log, MRN 01008971,
CaseID HA000557, dated 06/04/2026.

All pre-surgery pressures in the config matched the log. Three other values did
not.

### 4.1 The three corrections

| Field | Old | New | Source | Effect |
|---|---:|---:|---|---|
| `HR` | 119 bpm | **136 bpm** | log 09.24.57 "Nadi 136 bpm" | Model input. Sets cycle length: 60/119 = 0.504 s becomes 60/136 = 0.441 s. Changes every stroke-volume derivation. |
| `weight` / `height` / `BSA` | 14.0 kg / 98.0 cm / 0.6173 | **13.4 kg / 95.0 cm / 0.588** | log 07.53 | Model input. BSA drives all demographic scaling of the parameter prior. |
| `SAP_mean` (pre) | 71.3 mmHg | **77 mmHg** | log 10.39.12 "RFA 100/57 (77)" | Hard-tier fitted target. |

### 4.2 Why each correction is justified

**HR.** The old value 119 is also the NIBP systolic pressure on the adjacent
log line (`NIBP 119/83 (95)`). The likely cause is reading the wrong field.

**Demographics.** The old values were labelled "Keisya 2026-05-11 revision".
That date is five weeks after the 06/04/2026 catheterisation. The child grew
between the two dates. Pairing May anthropometry with April haemodynamics
scales the model to a larger patient than the one measured.

This correction restores documented intent rather than changing it.
`docs/clinical_data_dictionary.md` already stated that the active patient file
uses 13.4 / 95.0 / 0.588, and described 14.0 / 98.0 / 0.6173 as "an alternate
baseline-scaling experiment". The experiment values had leaked into both active
config files, so the dictionary's own statement was false until this fix.

The two BSA values also use different formulas. 0.588 is DuBois
(0.007184 x 95^0.725 x 13.4^0.425 = 0.5879). 0.6173 is Mosteller
(sqrt(98 x 14 / 3600)). The hospital-stamped DuBois value is used.

**SAP_mean.** The catheter reports its own mean pressure directly: 77 mmHg.
The config was reconstructing it as `dia + (sys - dia)/3 = 71.3`. The
reconstruction under-reads by about 5 to 6 mmHg systematically. The same offset
appears after closure (formula gives 75, catheter stamps 79). This answers by
measurement the MAP form-factor question raised in the 2026-08-28 assessment.

### 4.3 Silent-override trap

`recipe.demographics` is merged over `clinical.common` by
`apply_calibration_recipe_to_clinical.m:24`.

Correcting only `config/patient_reyna.m` would have been reverted at runtime
for weight, height and BSA. HR would have survived because it is not in that
struct. The result would be a run that appeared corrected but was not, with no
warning.

Both files were corrected together. A test now asserts the effective
post-merge values, so the two files cannot drift apart silently.

### 4.4 New post-surgery data

The same log contains post-closure pressures. The closure device was placed at
11.50.19 and released at 12.06.23. All readings below are stamped 12.15 to
12.32, under the same anaesthesia and ventilator settings as the pre-closure
readings.

`clinical.post_surgery` went from all `NaN` to 7 finite targets:

| Target | Value | Source rows |
|---|---:|---|
| `PAP_max` / `PAP_min` / `PAP_mean` | 17 / 9 / 13 | PA 17/8 (13), 17/9 (13), 17/9 (13) |
| `SAP_max` / `SAP_min` / `SAP_mean` | 89 / 68 / 79 | RFA 91/68 (79), 89/68 (78), 89/68 (79) |
| `RAP_mean` | 5 | RA 8/5 (5), 8/5 (5), 7/5 (5) |

Direction check passes. PA mean pressure falls from 15 to 13 after closure
while RAP holds at 5. This is consistent with removing a left-to-right shunt.

Because both states come from one session, this is a paired dataset, not two
separate studies.

### 4.5 A withdrawn finding: the "60% inconsistency"

Earlier documents described a "60% internal stroke-volume inconsistency,
severity critical" in the H+1 echo block. This is withdrawn.

Inside the H+1 block:

- `SV_LV` = 41.0 - 19.3 = 21.7 mL/beat
- `SV_RV` = 30.5 - 12.0 = 18.5 mL/beat
- difference = 14.7%, normal for echo
- `LVEF` check: 21.7 / 41.0 = 0.529 against the stated 0.528

The block is internally consistent. The 60% figure comes from comparing
post-closure volumes against pre-closure flow-derived stroke volumes:

- `SV_Qp` = 3.423 x 1.194 x 1000 / 119 = 34.35 mL/beat
- against `SV_LV` = 21.7, difference = 58.3%

That is the shunt itself, which is the expected physiology, not a data defect.
With the block excluded the audit now returns severity `none`.

### 4.6 Two open data questions, closed by sensitivity arms

Neither could be answered from the record. Both were closed by showing the
conclusions do not depend on them.

**Question 1: was `Qp = 4.087` Fick-derived from a BSA-indexed VO2?**

If yes, the BSA correction should scale it and move the `CO_Lmin` target by
-4.75%, from 3.4229 to 3.2603 L/min.

| Candidate | Error vs 3.4229 | Error vs 3.2603 |
|---|---:|---:|
| Seed 20260828 | -2.59% | +2.27% |
| Seed 20260830 | -7.50% | -2.89% |

The calibrated CO lands between the two candidate targets and passes the 10%
gate against either. The governed gate counts are unchanged under both
hypotheses. The Fick hypothesis would improve the apparent fit for seed
20260830, so keeping the reported value is the conservative choice.

**Question 2: what do `PARI 1.9` and `FR 1.19` mean?**

These would matter only if post-closure flow entered the analysis. It does not.
Post-closure `QpQs`, `CO`, `PVR` and `SVR` are all `NaN` and were deliberately
never entered rather than guessed. The 7 governed post targets, the joint fit
and the out-of-sample prediction all rest on pressures alone. Confirming these
abbreviations could add constraints but cannot revise anything reported.

---

## 5. Group D: Calibration results

All runs: Zhang scaling, historical seeds disabled, GSA on (N = 128, active set
reduced 14 to 7), 300 function evaluations, 6 multi-starts, sigma-weighted
objective.

### 5.1 Effect of the data correction

Same seed (20260828), same settings. The only changed variable is the clinical
data.

| | Before correction | After correction |
|---|---:|---:|
| Governed gate | 7 / 9 | **8 / 9** |
| All clinical targets | 9 / 11 | **10 / 11** |
| Within 5% band | 5 / 11 | **8 / 11** |
| Best RMSE (6 starts) | 0.0780 | **0.0480** |
| RMSE improvement vs baseline | 70.8% | **83.5%** |
| chi2 | 14.371 | **7.059** |
| chi2/N | 1.597 | **0.784** |
| Gate failures | PAP_max, SAP_min | **PAP_min only** |

The fit improved on every measure. This supports the correction being right:
fitting the patient who was actually measured gives a better fit than fitting a
mis-specified one.

`SAP_min` moved from the worst failure (-13.85%) to a clear pass (+4.27%). This
follows directly from an `SAP_mean` target that had been 5.7 mmHg too low and a
systemic waveform fitted at the wrong heart rate.

### 5.2 Per-metric result, seed 20260828 (corrected data)

Governed set, sorted worst first by |z|.

| Metric | Tier | Clinical | Model | Error % | sigma | z | 10% gate |
|---|---|---:|---:|---:|---:|---:|---|
| PAP_min | soft | 10 | 11.07 | +10.72 | 0.50 | +2.14 | FAIL |
| RAP_mean | hard | 5 | 5.26 | +5.20 | 0.25 | +1.04 | pass |
| PAP_max | soft | 20 | 19.00 | -5.00 | 1.00 | -1.00 | pass |
| SAP_min | soft | 57 | 59.43 | +4.27 | 5.70 | +0.43 | pass |
| SAP_max | soft | 100 | 96.05 | -3.95 | 10.00 | -0.39 | pass |
| CO_Lmin | hard | 3.423 | 3.334 | -2.59 | 0.50 | -0.18 | pass |
| PAP_mean | hard | 15 | 14.94 | -0.40 | 0.75 | -0.08 | pass |
| SAP_mean | hard | 77 | 76.85 | -0.19 | 3.85 | -0.04 | pass |
| QpQs | hard | 1.194 | 1.1946 | +0.05 | 0.0597 | +0.01 | pass |

### 5.3 Two seeds: a better gate score that is worse science

A second seed (20260830) was run on corrected data with identical settings.

| | Seed 20260828 | Seed 20260830 |
|---|---:|---:|
| Governed gate | 8 / 9 | **9 / 9** |
| Best RMSE | 0.0480 | **0.0434** |
| chi2/N | 0.784 | **0.494** |
| chi2 interpretation | consistent | **overfit** |
| Q_shunt_Lmin (not graded) | -2.26% | **-22.90%** |

The seed with the perfect gate score is the worse result. Three independent
signals agree:

1. `chi2/N` fell to 0.494, below the 0.5 overfit threshold. The model is
   fitting below the declared noise floor.
2. `Q_shunt_Lmin` degraded ten-fold. This is the one metric deliberately
   excluded from the governed RMSE. The graded set improved while the ungraded
   metric collapsed.
3. The perfect gate score coincides with both of the above and with no
   improvement in physiology.

`Q_shunt_Lmin` is a sensitive detector because it is `CO x (QpQs - 1)`, and
`QpQs - 1 = 0.194` is a small difference of two near-equal numbers. A -2.71%
error in `QpQs` becomes about -16.6% in the difference. The amplification is
about 6x.

### 5.4 Run-to-run spread

| Statistic | Value |
|---|---|
| Best-of-6, seed 20260828 | 0.0480 |
| Best-of-6, seed 20260830 | 0.0434 |
| Range | 0.0434 to 0.0480, about 10% of the value |
| Pooled 12 starts | min 0.0434, median 0.0779, max 0.1053 |

This is not yet a 95% confidence interval, for two reasons. Two seeds give a
range, not a distribution. More importantly, the reported figure is a
best-of-6 selection, which is an extremum. A confidence interval computed over
the 12 pooled starts would describe the spread of attempts, not the uncertainty
of the reported number.

The defensible statement: two independent 6-start selections gave 0.0434 and
0.0480, and the difference between them changes the gate count from 9/9 to 8/9.
Run-to-run variation is large enough to move the headline claim.

### 5.5 Parameter identifiability

At the corrected operating point, seed 20260828:

- Condition number: **2.06e3**, flagged as near-dependent (threshold 1e3)
- `E.LV.EA` and `E.LV.EB`: rho = -0.921
- `E.LV.EA` and `vsd.Cd`: rho = -0.917

The condition number is worse than the superseded run's 232. The LV elastance
pair appears in both runs, so it is structural rather than an artefact of the
wrong data. The new `E.LV.EA` to `vsd.Cd` coupling is consistent with a better
shunt fit: as the fit improves, LV contractility and orifice discharge trade
off against each other more sharply.

A fit that improves while its parameters become less separable is a warning.
It supports reducing the number of free parameters.

---

## 6. Group D continued: joint pre/post inversion (Phase 4)

### 6.1 Purpose

Pre-only fitting has 9 observations and 12 free parameters, so `dof = -3`. With
more parameters than observations, a low chi-squared is guaranteed and proves
nothing.

Fitting both scenarios at once with one shared parameter set raises the
observation count without new data collection.

| | Pre-only | Joint |
|---|---:|---:|
| N | 9 | **16** (9 pre + 7 post) |
| p | 12 | 12 to 14 |
| dof | **-3** | **+2** at p = 14, **+4** at p = 12 |

Note PRD §7.5 asks for `dof >= 6`. That is not reachable at N = 16 with p in
the 12 to 14 range. Either p comes down or that criterion needs revising.

### 6.2 Two bugs found while building it

**Bug 1: closing the VSD by resistance alone does nothing in orifice mode.**

`vsd_shunt_model` dispatches on `params.vsd.mode`. Resistive modes close by
setting a large `R.vsd`. The `orifice_bidirectional` mode never reads `R.vsd`
at all and closes only when `vsd.area_mm2 = 0`.

Reyna runs in orifice mode. Closing by `R.vsd` alone leaves the shunt fully
open. Any "post-closure" simulation built that way would silently be a
pre-closure simulation.

This defect also exists at `main_run.m:839`, in the pre-to-post seed handoff.
It is recorded and not yet fixed.

**Bug 2: the joint driver governed the wrong observation set.**

Building per-scenario tiers with a bare `build_target_tiers(clinical, scenario)`
call does not honour `recipe.primary_rmse_holdout`. It governed 10 pre-surgery
rows instead of 9, readmitting `Q_shunt_Lmin`, which is the metric that turned
out to be the overfitting detector in §5.3.

`chi2_pre` would then have been computed over a different observation set than
the reported governed RMSE, making joint and single-scenario statistics
incomparable. Tiers now come from each scenario's case profile. A test pins
`n_pre` at exactly 9.

### 6.3 First run was a non-result

The first joint run moved `J` from 661.6533 to 661.5832, a change of 0.01%, and
stopped after one iteration. That is an optimiser returning its starting point,
not a calibration result.

Cause: the driver left fmincon's finite-difference settings at their defaults.
The default forward-difference step is about `sqrt(eps)`, roughly 1.5e-8
relative. That is far below the noise floor of an objective built on an ODE
steady-state solve, so the finite differences measured integration noise rather
than the gradient.

Evidence: reported first-order optimality of 2.1e6 alongside steps of 5e-8, and
a halt on `StepTolerance = 1e-10`.

Fix: `FiniteDifferenceStepSize = 1e-5` and `StepTolerance = 1e-6`, matching
`run_calibration.m:728-734`, which was already proven to converge on this model.
First-order optimality dropped to 4.8e3 and the step size to 1.3e-1.

An `OPTIMIZER_DID_NOT_MOVE` guard now flags relative improvement below 1e-3 or
relative parameter step below 1e-6.

### 6.4 Converged result

| | Baseline | Joint fit |
|---|---:|---:|
| J | 661.6533 | **96.1514** (-85.5%) |
| chi2_pre | 491.369 | **79.289** |
| chi2_post | 170.284 | **16.862** |
| N / p / dof | | 16 / 14 / **+2** |
| chi2/N | 41.35 | **6.01** |

Multi-start check, 4 starts at 1500 evaluations each:

| Start | Origin | J |
|---|---|---:|
| 1 | baseline x0 | **96.1514** |
| 2 | warm start from pre-only fit | 205.8968 |
| 3 | bound-interior perturbation | 515.8069 |
| 4 | bound-interior perturbation | 115.5764 |

No start beat 96.15. Under-convergence is ruled out at this budget.

**Interpretation.** `chi2/N = 6.01` is above the consistent band (0.5 to 2.0),
so the joint fit is underfit. With one shared parameter set, this model cannot
reproduce both haemodynamic states within declared measurement uncertainty.

This is a negative result and is reported as one. The pre-only fit looked
excellent (`chi2/N = 0.78`) largely because it was underdetermined. Adding real
constraints exposes that.

Two caveats. The warm start was degraded (see §6.5), so the sharpest form of
the test has not been run exactly as intended. And this is a single-stage
optimiser measured against a single-scenario result built from a 6-stage
pipeline, so structure still differs, not only budget.

### 6.5 The calibration vector does not fully determine the model

Warm-starting the joint fit from the pre-only solution exposed this.

Reading that solution's 14 calibration parameters out and applying them to a
fresh baseline does not reproduce it. `V0.SVEN` comes back as 564.02 against
the original 595.934, because it is a coupled quantity that lives outside the
calibration vector. That is why the warm start scored J = 4413 instead of
landing near the pre-only optimum.

The joint objective's parameter-sharing claim is unaffected, because both
scenario structs are built from the same base and receive the same vector, so
everything outside the vector is identical by construction.

The out-of-sample result in §7 is also unaffected, because it loads the full
calibrated parameter struct directly and never round-trips through the vector.

---

## 7. Out-of-sample prediction test

This is the strongest result in this work.

### 7.1 Why it is a real test

Before this, nothing in the recipe was a genuine holdout. Every finite target
was either fitted, or was algebra over fitted quantities (`SVR`,
`Q_shunt_Lmin`), so predicting it demonstrated nothing new.

The post-closure pressures are different. They are independent measurements, in
a different haemodynamic state, and were never seen by the pre-only
calibrations.

Method: take a parameter set calibrated on pre-closure data only, close the
defect, simulate, and compare against the post-closure pressures. No refitting.
Closure is the only intervention.

Script: `scripts/evaluate_post_closure_prediction.m`. It asserts the shunt is
actually shut before measuring anything, because closure is mode-dependent
(§6.2).

### 7.2 Result

| Metric | Measured | Predicted (seed 828) | Error % | Predicted (seed 830) | Error % |
|---|---:|---:|---:|---:|---:|
| SAP_mean | 79 | 81.75 | **+3.48** | 78.49 | **-0.64** |
| PAP_mean | 13 | 13.83 | +6.35 | 14.05 | +8.04 |
| RAP_mean | 5 | 5.36 | +7.29 | 5.36 | +7.24 |
| SAP_min | 68 | 63.29 | -6.92 | 59.46 | -12.56 |
| PAP_max | 17 | 17.75 | +4.39 | 18.78 | +10.46 |
| PAP_min | 9 | 10.33 | +14.80 | 10.01 | +11.19 |
| SAP_max | 89 | 101.40 | +13.94 | 99.59 | +11.90 |
| **Within 10%** | | **5 / 7** | | **3 / 7** | |
| **chi2/N** | | **2.311** | | **2.439** | |

### 7.3 Overfitting confirmed with held-out data

| | Seed 20260828 | Seed 20260830 |
|---|---:|---:|
| In-sample gate | 8 / 9 | **9 / 9** |
| Out-of-sample prediction | **5 / 7** | 3 / 7 |
| Out-of-sample chi2/N | **2.311** | 2.439 |

The arm that fits the training data better predicts held-out data worse. This
converts the overfitting argument in §5.3 from an inference into a
demonstration.

### 7.4 Honest reading

- All 7 pressures land within about 15%. Mean pressures, which are the most
  reliably measured, do best. `SAP_mean` is predicted to 0.6 to 3.5%.
- `chi2/N` of 2.3 to 2.4 is just above the consistent band, so this is
  marginally underfit rather than in full agreement.
- Errors are systematic, not random. Six of seven are positive in both seeds,
  and pulmonary pressures are over-predicted throughout. The model predicts
  less pulmonary unloading after closure than actually occurred. This is a
  physiological discrepancy and a concrete lead for model improvement.
- Unlike `chi2/N` on the fitted set, this number needs no degrees-of-freedom
  caveat. The model cannot have absorbed targets it never saw.

### 7.5 The holdout is proven uncontaminated

The post-closure block was added to `config/patient_reyna.m` before the
corrected pre-only calibrations were run. If any part of the pre-surgery path
read `clinical.post_surgery`, the holdout would be contaminated.

Proof: blank every `post_surgery` field and re-derive the pre-surgery pipeline.

| Artefact | Identical with post data blanked |
|---|---|
| `get_calibration_targets('pre_surgery', ...)` | yes |
| Case-profile tier table | yes |
| Allowed metric fields | yes |
| `params_from_clinical(..., 'pre_surgery', ...)` | yes |

`apply_post_surgery_warm_start` also returns immediately unless the scenario is
`post_surgery` (`src/utils/apply_post_surgery_warm_start.m:39-40`).

Note on method: the first version of this check used `isequal` and reported a
difference, which looked like contamination. The cause was
`isequal(NaN, NaN) = false` and these structs containing many legitimately
`NaN` fields. A field-by-field diff found zero real differences. The check uses
`isequaln`.

### 7.6 What can be claimed

> A lumped-parameter model calibrated solely on pre-closure catheterisation
> data predicted all seven independently measured post-closure pressures within
> 15% (five of seven within 10%), with mean arterial pressure predicted to
> within 3.5%, after applying defect closure as the only intervention.

---

## 8. Group E: bugs found

| # | Bug | Consequence if unfixed |
|---|---|---|
| 1 | `PAP_min`, `PAP_max` graded but not fitted | Optimiser scored on metrics it never optimised |
| 2 | Two code paths gave 10 vs 11 governed rows | RMSE denominator ambiguous |
| 3 | `dof` clamped to 0 in chi-squared report | Printed a false line; hid over-parameterisation |
| 4 | Chi-squared band label ignored `dof` | `chi2/N = 0.494` reported as "consistent" |
| 5 | Closing VSD by `R.vsd` is a no-op in orifice mode | "Post-closure" run would silently be pre-closure |
| 6 | Joint driver governed 10 pre rows, not 9 | Joint and single-scenario chi-squared incomparable |
| 7 | Joint optimiser used default finite-difference step | Optimiser returned its starting point |

Bugs 3, 4, 5, 6 and 7 all produce output that looks valid. None of them raise
an error.

### 8.1 Detail on bugs 3 and 4

`compute_chi_squared_report.m` computed `dof = max(n_obs - n_parameters, 0)`.
With N = 9 and p = 12 the true value is -3, but the console printed
`dof = N - p : 0`. That line is arithmetically wrong, and `dof = 0` reads as
"exactly determined" when the truth is "over-parameterised".

`classify_chi2_per_obs` assigned "consistent" from `chi2/N` alone. With more
parameters than observations, small residuals are guaranteed. Seed 20260830's
`chi2/N = 0.494` was being reported as residuals matching measurement noise.

After the fix:

```
dof = N - p               : -3  [dof <= 0: MORE FREE PARAMETERS THAN
                                 OBSERVATIONS -- a low chi2/N is guaranteed
                                 here and is not evidence of fit quality]
interpretation            : consistent_but_underdetermined
[QUALIFIED] ... Do not quote chi2/N alone as validation.
```

An `underfit` verdict is deliberately not qualified. Failing to match the data
despite excess freedom is a real signal.

---

## 9. What this work does not establish

- **n = 1.** One patient.
- **The pre-only fit is underdetermined.** N = 9 against p = 12 gives
  `dof = -3`. Any claim must rest on the gate count, the per-metric residuals
  and the out-of-sample test, and must state the parameter-to-observation ratio.
- **Joint fitting does not resolve this.** `dof` becomes +2 to +4, which is
  still low, and the joint fit is underfit at `chi2/N = 6.01`.
- **Two seeds are a range, not a confidence interval.** See §5.4.
- **Run-to-run variation moves the headline** between 8/9 and 9/9.
- **Identifiability is poor.** Condition number 2.06e3 with two collinear
  parameter pairs.
- **Zhang versus Lundquist remains uninterpretable** at these budgets.
- **No post-closure flow data.** The post state constrains pressures only.
- **Reducing p is untested.** Stages D to F use a 12-parameter mask, wider than
  the GSA screen's 7. Lowering p is the other lever on degrees of freedom and
  has not been tried.

---

## 10. Testing

| Suite | Result |
|---|---|
| `functiontests`-style (11 files) | 96 / 96 pass |
| `test_post_closure_prediction.m` | 6 / 6 pass |
| `test_reyna_systemic_flow_profile.m` | 6 pass / 3 fail, same as pre-change baseline |

The 3 failures in the last file are pre-existing. This was verified by stashing
all changes and re-running.

Note the test suite is mixed-style. Script-style tests must be run with
`run('tests/x.m')`. Using `runtests` on them splits each file into isolated
blocks and destroys the shared counter variables, producing false failures.

---

## 11. Files added

| File | Purpose |
|---|---|
| `src/utils/compute_chi_squared_report.m` | Chi-squared statistic and reporting |
| `src/calibration/analyse_parameter_identifiability.m` | Sensitivity matrix, condition number, collinearity |
| `src/calibration/identifiability_tables_from_matrix.m` | Table construction |
| `src/calibration/objective_joint_pre_post.m` | Joint pre/post objective |
| `scripts/run_joint_pre_post_calibration.m` | Joint fit driver, multi-start |
| `scripts/evaluate_post_closure_prediction.m` | Out-of-sample prediction test |
| `src/utils/print_validation_holdout.m` | Holdout reporting block |
| `src/utils/export_full_metric_gate.m` | Per-metric CSV with sigma and z-scores |

Evidence tables are now tracked in git. `.gitignore` was extended to un-ignore
`chi_squared_*.csv` and `parameter_identifiability_*.csv`, matching the
existing exception for `full_metric_gate_*.csv`.

---

## 12. Reproduction

Single-scenario calibration on corrected data:

```bash
matlab -batch "cd('D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd'); addpath(genpath(pwd)); setenv('UNIFIED_VSD_UQLAB_PATH', fullfile(pwd,'toolbox','UQLab_Rel2.2.0','core')); setenv('UNIFIED_VSD_SCALING_MODE','zhang'); setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','1'); setenv('UNIFIED_VSD_DO_GSA','1'); setenv('UNIFIED_VSD_GSA_PCE_N','128'); setenv('UNIFIED_VSD_DO_PLOTS','0'); setenv('UNIFIED_VSD_MAX_FUN_EVALS','300'); setenv('UNIFIED_VSD_MAX_ITERATIONS','40'); setenv('UNIFIED_VSD_NUM_STARTS','6'); setenv('UNIFIED_VSD_MULTISTART_SEED','20260828'); setenv('UNIFIED_VSD_OBJECTIVE_WEIGHTING','sigma'); main_run('pre_surgery', patient_reyna())"
```

Out-of-sample prediction test:

```bash
matlab -batch "cd('D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd'); addpath(genpath(pwd)); evaluate_post_closure_prediction('results/runs/20260829_225931_reyna_pre_surgery/mat/params_accepted_candidate_pre_surgery.mat')"
```

Run folders:

- `results/runs/20260829_225931_reyna_pre_surgery` (seed 20260828, corrected)
- `results/runs/20260830_033638_reyna_pre_surgery` (seed 20260830, corrected)

---

## 13. Known issue not yet fixed

`main_run.m:839` builds the pre-to-post seed by setting `R.vsd = 1e6` only.
For orifice-mode patients such as Reyna this is a no-op, so that seed is not a
closed-VSD model. This should be fixed before any post-surgery run uses it.
