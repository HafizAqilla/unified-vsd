# Codex Reyna Systemic Flow Handoff

Date: 2026-05-08
Branch: `codex/reyna-systemic-flow-consistency`
Repo: `D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd`

## Purpose

This note summarizes the current state of the MATLAB pediatric VSD lumped-parameter model work so it can be moved into a new analysis session. The main question is no longer whether the ODE solver works. The current bottleneck is calibration governance, target interpretation, identifiability, and especially the systemic output / resistance tradeoff in the Reyna pre-surgery case.

## Current Scientific Position

The model is not obviously broken or useless. It produces coherent VSD behavior, internally consistent flow accounting, and reproducible calibration tradeoffs. However, it is not yet sufficient to match all of Reyna's clinical targets simultaneously under the current model structure, objective function, and target definitions.

Best current wording:

> The model is physiologically meaningful and reproducible, but the remaining bottleneck is structural identifiability and systemic-output consistency, especially in matching `CO_Lmin` together with systemic pressure and SVR targets.

Avoid saying:

> The model is wrong.

Better:

> The current model formulation is partially valid, but still incomplete for simultaneous patient-specific matching of all key hemodynamic targets.

## Important Literature Interpretation

### Bozkurt Foundation

Bozkurt 2019 is a valid structural foundation for a lumped cardiovascular model with time-varying elastance and vascular compartments, but it is not a complete foundation for patient-specific pediatric VSD inverse calibration.

Important caveat:

- Bozkurt's child model is generic, roughly BSA around 1 m2 and age 8-12 years.
- It was not a patient-specific congenital VSD calibration framework.
- Some parameters were manually adjusted to generic healthy/DCM ranges.

Use Bozkurt as:

- structural model foundation
- chamber and vascular magnitude reference
- generic adult/child physiology reference

Do not use Bozkurt as:

- hard pediatric VSD patient-specific truth
- proof that preschool VSD cases should fit easily
- source for all parameter bounds

### Pediatric Scaling Papers

Adult baseline to pediatric scaling is defensible as a starting prior, not as exact truth.

Current preferred scaling mode:

- `lundquist_bsa`

Reason:

- it aligns with adult baseline to pediatric BSA scaling logic
- it performed better on Reyna and Patient A than Zhang in earlier comparison
- Zhang remains useful as a sensitivity comparator

Age caveat:

- Reyna is about 3.17 years old.
- The code now annotates this as `preschool_extrapolated`.
- This lies between infant maturation literature and generic school-age child anchors.

## Major Work Completed On Branch

### Batch A/B Governance And Identifiability

Added explicit case modes:

- `full_data`
- `sparse_cath`
- `synthetic_benchmark`

Key file:

- `src/calibration/build_case_calibration_profile.m`

Effect:

- Reyna gets full-data calibration freedom.
- Razka is restricted as sparse cath.
- Synthetic profile A is treated as stress-test benchmark, not equal real-patient truth.

### Parameter Registry And Plausibility

Added centralized parameter registry:

- `config/build_parameter_registry.m`

Added registry-driven bounds and plausibility evaluation:

- `src/calibration/build_calibration_vector.m`
- `src/calibration/validate_bounds.m`
- `src/calibration/evaluate_parameter_plausibility.m`

The model now reports:

- fitted vs scaled baseline
- ratio to baseline
- bounds
- `OK`, `WARNING`, `FAIL`

This was in response to Pak Dipo's concern:

> Good output fit is not enough if fitted parameters become biologically meaningless.

### Candidate Saving And Rollback Governance

The run pipeline now saves both:

- `best_candidate`
- `accepted_candidate`

Important files:

- `main_run.m`

New artifacts inside each run:

- `mat/params_best_candidate_pre_surgery.mat`
- `mat/params_accepted_candidate_pre_surgery.mat`
- `tables/validation_best_candidate_pre_surgery.csv`
- `tables/validation_accepted_candidate_pre_surgery.csv`
- `tables/parameter_plausibility_best_candidate_pre_surgery.csv`
- `tables/parameter_plausibility_accepted_candidate_pre_surgery.csv`

This fixed an earlier issue where promising calibration states were hidden by rollback.

### Age Validity And Baseline Provenance

Added age-validity annotation:

- `src/utils/build_age_validity_annotation.m`

Added baseline provenance audit:

- `config/build_baseline_provenance.m`
- `scripts/audit_baseline_provenance.m`

Current baseline is explicitly recognized as a hybrid foundation, not one clean literature table.

### Vascular R-C Coupling And Priors

Added grouped vascular resistance scales and R-C coupling:

- `src/calibration/get_vascular_rc_coupling_rules.m`
- `src/calibration/enforce_vascular_rc_coupling.m`
- `src/calibration/set_calibration_param_value.m`
- `src/calibration/get_calibration_param_value.m`

Scientific intent:

- avoid independent unrealistic drift in serial vascular terms
- tie vascular compliance behavior to resistance where defensible
- allow direct arterial compliance `C.SAR` and `C.PAR` residual tuning from stroke volume / pulse pressure

## Key Runs And Results

### Earlier Good Benchmark Before Major Governance Tightening

Reyna preserved run:

- `results/runs/20260503_203819_reyna_pre_surgery`

Result:

- RMSE about `0.1132`

Strong:

- `RAP_mean = -1.03%`
- `PAP_min = +2.19%`
- `SAP_mean = -2.73%`
- `LVEDV = +4.69%`

Weak:

- `PAP_mean = -8.46%`
- `QpQs = -5.30%`
- `SVR = +19.98%`
- `CO_Lmin = -19.03%`
- `RVEDV = +16.57%`
- `SAP_max = -15.72%`
- `LVEF = +15.12%`

This is still the best scalar RMSE benchmark, but it predates much of the newer plausibility/governance work.

### New Governance Best Candidate Before Systemic Load Patch

Run:

- `results/runs/20260508_083312_reyna_pre_surgery`

Manifest:

- `RMSE_Baseline = 0.247357`
- `RMSE_Calibrated = 0.141028`
- `CalibrationStatus = PROMISING_NEAR_MISS`
- `RollbackApplied = 0`

Primary gate:

- `RAP_mean = -6.65%` fail
- `PAP_min = +16.04%` fail
- `PAP_mean = +0.75%` pass
- `SAP_mean = -5.18%` fail by tiny margin
- `QpQs = -4.70%` pass

Strong improvements:

- `SAP_mean = -38.19% -> -5.18%`
- `QpQs = -7.90% -> -4.70%`
- `LVEDV = -22.19% -> -2.62%`
- `PAP_mean = +3.09% -> +0.75%`

Still weak:

- `PAP_min = +16.04%`
- `SVR = +34.78%`
- `CO_Lmin = -29.57%`
- `SAP_min = -12.55%`
- `SAP_max = -12.12%`

Parameter plausibility:

- `0 FAIL`
- `5 WARNING`

### Systemic Load Patch Experiment

Run:

- `results/runs/20260508_102935_reyna_pre_surgery`

Purpose:

- allow lower `group.R_sys_scale`
- strengthen penalty against low `CO_Lmin` plus high `SVR`

Best candidate:

- `RMSE = 0.135276`
- but rolled back because status remained `REJECT`

Compared to previous promising run:

- `CO_Lmin` improved: `2.4109 -> 2.6676 L/min`
- `SVR` improved: `26.11 -> 22.12 WU`
- RMSE improved: `0.1410 -> 0.1353`

But tradeoff worsened:

- `SAP_mean` worsened: `-5.18% -> -10.95%`
- `SAP_max` worsened: `-12.12% -> -22.91%`
- `LVEF` worsened: `+8.18% -> +17.60%`
- `LVESV` worsened: `-11.29% -> -19.50%`
- `RVEDV` worsened: `+7.32% -> +15.98%`

Conclusion:

> Reducing systemic load can improve `CO_Lmin` and `SVR`, but the objective then sacrifices systemic pressure waveform and chamber function.

## Uncertainty Ensemble Around Reyna

Script:

- `scripts/run_reyna_uncertainty_ensemble.m`
- `scripts/analyze_reyna_uncertainty_ensemble.m`

Latest analysis:

- `results/tables/reyna_uncertainty_analysis_20260508_102447.csv`

Results:

| Metric | Target | Ensemble Range | Min Abs Error | Assessment |
|---|---:|---:|---:|---|
| `CO_Lmin` | 3.423 | 2.0225 to 2.0416 | 40.36% | `likely_structural_or_strongly_constrained` |
| `SVR` | 19.37 | 17.4906 to 21.1929 | 0.19% | `prior_sensitive_recoverable` |
| `PAP_min` | 10 | 10.6239 to 13.0275 | 6.24% | `prior_sensitive` |

Interpretation:

- `SVR` is recoverable by prior/load tuning.
- `PAP_min` is prior-sensitive but not fully solved by current prior ranges.
- `CO_Lmin` is the main structural or strongly constrained bottleneck.

## Systemic Output Audit

Script:

- `scripts/audit_systemic_output_bottleneck.m`

Important output from `20260508_083312_reyna_pre_surgery`:

| Candidate | CO_Qs_Lmin | Qao_Lmin | LVCO_Lmin | Qvsd_Lmin | SVR | Q_if_target_SVR |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 2.0335 | 2.0338 | 2.2224 | 0.2027 | 19.369 | 2.0334 |
| best/accepted | 2.4109 | 2.4112 | 2.7138 | 0.3323 | 26.107 | 3.2494 |

Key finding:

- `Qao_Lmin` and `Qs_Lmin` agree almost exactly.
- Therefore low `CO_Lmin` is not a flow accounting bug.
- It is a real low systemic through-flow state in the model.

Important clue:

- Best candidate has `MAP - RAP = 62.942 mmHg`.
- Model `SVR = 26.107 WU`, yielding `CO = 2.4109 L/min`.
- If the same pressure gradient used target SVR, implied flow would be `3.2494 L/min`, close to target `3.423 L/min`.

Conclusion:

> The model can create a near-correct systemic pressure gradient, but only by ending up with too-high effective SVR, which keeps systemic flow too low.

## Most Important Suspected Workflow Mistakes

### 1. Clinical CO Target Meaning May Be Misinterpreted

The code defines:

```matlab
metrics.CO_Lmin = Qsys_Lmin;
```

So model `CO_Lmin` means effective systemic output, i.e. `Qs`.

But clinical "CO" may mean:

- Fick systemic flow `Qs`
- echo LV stroke volume times HR
- catheter-derived cardiac output
- another derived value

If Reyna's clinical `CO_Lmin = 3.423` is not truly systemic `Qs`, then we may be forcing an unfair target comparison.

This is the single most important data-definition audit before further model tuning.

### 2. Adult Baseline Plus Scaling Was Treated Too Much Like Truth

Adult baseline to pediatric scaling is valid as a prior, not a patient truth. Reyna is in a preschool age regime where current literature support is extrapolated.

### 3. GSA Was Initially Used Too Much Like Identifiability

Sobol sensitivity means parameter affects outputs. It does not mean the available clinical data can estimate that parameter uniquely.

### 4. Correlated Clinical Targets Were Over-Weighted Together

`SAP_mean`, `RAP_mean`, `CO_Lmin`, and `SVR` are mathematically linked. Treating all as independent can overconstrain the objective and create impossible tradeoffs.

### 5. Sparse Cases Can Look Better Than Full-Data Cases

Razka can get low RMSE because it has fewer available targets. That is not proof of better patient-specific calibration.

## Current Best Diagnosis

The model is not failing globally. It is repeatedly finding coherent tradeoffs:

- pressure targets can improve
- volume targets can improve
- shunt behavior is reasonable
- parameter values mostly stay plausible

But systemic output remains too low when pressures are reasonable.

Current likely bottleneck:

> systemic resistance/load representation and/or target definition of clinical `CO_Lmin`.

## Recommended Next Questions For A Fresh Model Review

1. Is Reyna's clinical `CO_Lmin` truly systemic flow `Qs`, or is it LV cardiac output / echo-derived CO?
2. Should validation compare clinical CO to model `Qs_Lmin`, `Qao_Lmin`, or `LVCO_Lmin`?
3. Should `SVR` be treated as a direct target, a derived consistency check, or both with reduced weight?
4. Should `SAP_mean`, `CO_Lmin`, and `SVR` be handled as one coupled measurement bundle rather than separate loss terms?
5. Does the systemic circuit topology impose too much effective resistance after calibration?
6. Are chamber elastance and V0 terms compensating for systemic load mismatch?
7. Does the VSD shunt model need geometry/orifice formulation for Reyna instead of linear resistance?
8. Is the target `PAP_min = 10` reliable as a cycle minimum, or is it a cath diastolic notation not equivalent to simulated minimum?
9. Should acceptance be based on direct measured targets first, with derived quantities secondary?
10. Should the thesis frame this as a physiologically constrained digital shadow rather than exact patient-specific parameter identification?

## Good Prompt For Another Model

Use this as the opening prompt:

```text
We are working on a MATLAB 0D lumped-parameter pediatric VSD model in branch codex/reyna-systemic-flow-consistency.

The model now has calibration governance, parameter plausibility checks, baseline provenance, age-validity annotation, candidate saving, vascular R-C coupling, and uncertainty/systemic-output audits.

Current main benchmark is Reyna pre-surgery, a full-data real case.

Best recent near-miss run:
- results/runs/20260508_083312_reyna_pre_surgery
- RMSE calibrated = 0.141028
- status = PROMISING_NEAR_MISS
- no rollback
- 0 parameter plausibility FAIL, 5 WARNING
- QpQs = -4.70%
- PAP_mean = +0.75%
- LVEDV = -2.62%
- SAP_mean = -5.18%
- but CO_Lmin = -29.57%, SVR = +34.78%, PAP_min = +16.04%

Systemic audit shows Qao and Qs agree, so low CO_Lmin is not a flow accounting bug. It is a real low-throughput state. If the model pressure gradient used target SVR, implied Q would be near clinical target.

Uncertainty ensemble says:
- SVR is prior-sensitive and recoverable
- PAP_min is prior-sensitive but still difficult
- CO_Lmin is likely structural or strongly constrained

Question:
Please review the workflow scientifically. Are we comparing clinical CO to the correct model output? Should clinical CO be compared to Qs_Lmin, Qao_Lmin, or LVCO_Lmin in a VSD case? Are SAP_mean, RAP_mean, CO_Lmin, and SVR being overconstrained as independent targets? What is the most defensible next modeling step before more optimization?
```

## Suggested Stop Point

Stop coding for now. Do not tune more until the clinical definition of `CO_Lmin` is verified.

The next step should be a data-definition and objective-structure review, especially around:

- `CO_Lmin`
- `Qs_Lmin`
- `Qao_Lmin`
- `LVCO_Lmin`
- `SVR`
- direct vs derived targets
