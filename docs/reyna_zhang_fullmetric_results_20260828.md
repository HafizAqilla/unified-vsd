# Reyna Zhang Full-Metric Campaign — Results, 2026-08-28

> ## ⚠ Superseded — do not quote these numbers
>
> This document's results were calibrated against clinical inputs that the
> publication-readiness protocol-data reconciliation (2026-09-05) has since
> corrected: heart rate 136→119 bpm, VSD diameter 3.025→3.665 mm, and
> pre-surgery chamber volumes NaN→consistency-only 32/23.6/30.5/12 mL. See
> `docs/CHANGES_SINCE_PR22.md` §14 for the full before/after and
> `docs/publication_readiness_prd.md` for the re-run this triggered. This
> document is retained for its process findings (the seven governance
> defects it found and fixed), not for its numbers.

Branch: `codex/reyna-zhang-fullmetric-10pct`
Scope: Reyna `pre_surgery`, Zhang scaling prior, fair-prior (historical seeds disabled)
PRD: [reyna_zhang_full_metric_prd.md](reyna_zhang_full_metric_prd.md)
Audit that motivated it: [reyna_zhang_scientific_assessment_20260828.md](reyna_zhang_scientific_assessment_20260828.md)

---

## 1. What this branch changes

Seven defects, each found by execution rather than inspection, and each closed
with an assertion so it cannot silently return.

| # | Defect | Fix |
|---|---|---|
| D1 | `PAP_min` / `PAP_max` declared `UseForCalibration`, but assigned no tier — graded inside the governed RMSE while absent from the objective | promoted to `soft`; weights `0.45` / `0.50` |
| D2 | Same class in `LAP_mean`, `RVEDV` (pre) and `RVEF` (post) | promoted to `soft`; `RVEF` also gets the `LVEF` anti-double-counting rule |
| D3 | No guard against the class recurring | `build_target_tiers:ungovernedCalibrationTarget` assertion, swept over all 13 profiles x 2 scenarios |
| D4 | `primary_rmse_holdout` removed metrics from the **objective**, not just the RMSE, so `Q_shunt_Lmin` was never fitted despite the recipe asking for it | `merge_allowed_metrics_with_tiers` now excludes only genuinely report-only tiers |
| D5 | Objective had no term for the acceptance band it is reported against | C¹ gate hinge, band read from `recipe.acceptance.primary_gate_pct` |
| D6 | `ACCEPT` computed over 5 selected primaries only | acceptance now requires the whole governed RMSE mask |
| D7 | `RVESV` fitted, and `override_IC` seeding chamber V0/E, both from the H+1 post-operative echo block | `RVESV` demoted to consistency-only; `override_IC = false`; `assert_evidence_timing_governance` enforces it |

Two supporting additions: `export_full_metric_gate` publishes the per-metric
acceptance table on every run (previously unavailable in any reviewable form),
and `build_multistart_starts` provides a deterministic Sobol multi-start so a
scaling-prior comparison can distinguish a worse prior from a solver that never
moved.

## 2. The measurement that started it

Baseline, fair-prior Zhang at the published controls
(`MaxFunctionEvaluations = 60`, `MaxIterations = 8`, GSA off, seeds disabled):

- primary RMSE `0.1040`, full RMSE `0.1590`
- **11 of 16** clinical targets within 10%
- **8 of 10** governed-RMSE targets within 10%
- classified **`ACCEPT`**, `primary10_fail=0`

That last line is the point. The run passed because the acceptance label only
looked at the five selected primaries, while `PAP_min` — a directly measured
catheter pressure, inside the governed RMSE — sat at **23.83%**.

| Metric | \|Error\| | Tier | In objective |
|---|---:|---|:--:|
| `LVEF` | 32.02% | consistency_check_only | no |
| `LVESV` | 30.68% | consistency_check_only | no |
| `RVEDV` | 28.83% | consistency_check_only | no |
| `PAP_min` | **23.83%** | validation_only | **no** |
| `SAP_min` | 13.98% | soft | yes |
| `RVESV` | 9.40% | soft | yes |
| `CO_Lmin` | 9.02% | hard | yes |
| `PAP_max` | 8.00% | validation_only | **no** |
| `Q_shunt_Lmin` | 7.86% | soft | yes |
| `LVEDV` | 7.71% | consistency_check_only | no |
| `SAP_max` | 6.63% | soft | yes |
| `SVR` | 6.42% | derived_validation | no |
| `RAP_mean` | 4.84% | hard | yes |
| `PAP_mean` | 3.34% | hard | yes |
| `SAP_mean` | 2.62% | hard | yes |
| `QpQs` | 0.20% | hard | yes |

Committed as
[docs/evidence/reyna_zhang_fullmetric_baseline_20260828.csv](evidence/reyna_zhang_fullmetric_baseline_20260828.csv).

The three worst rows overall (28–32%) are the H+1 post-operative echo block.
They are not valid pre-operative targets, which is why **"16 of 16" is the
wrong objective** and D7 removes the last two channels through which that block
was still influencing a pre-operative fit.

## 3. Result after the fixes

UQLab 2.2.0 was installed from `toolbox/UQLab_Rel2.2.0` (present but never
installed, which is why the GSA preflight had been failing). **This is the
first GSA the project has run.**

### 3.1 The parameter-count problem is resolved

Sobol screening at `N = 128` reduced the active set from 14 to **7**:

`group.R_sys_scale`, `R.SVEN`, `group.R_pul_scale`, `C.SAR`, `C.PAR`,
`E.RV.EB`, `vsd.Cd`

| | Fitted parameters | Independent measurements | Degrees of freedom |
|---|---:|---:|---:|
| Before | 14 | 8 | **−6** |
| After GSA | 7 | 8 | **+1** |

Negative degrees of freedom was the central objection in the audit: under it a
low RMSE measures model flexibility, not fidelity. It is now positive, which is
the precondition for the fit carrying any evidential weight at all. PRD Phase 4
targeted "7 or fewer"; the screen landed exactly there without being told to.

Note also that the screen is now driven by nine metrics rather than five,
because Phase 1 put the diastolic and systolic pressures into the calibration
set. The sensitivity analysis is answering a better-posed question than before.

### 3.2 Iteration history

All arms: fair-prior Zhang, GSA on at `N = 128`, chamber volumes removed.

| Arm | Budget | Governed gate | Primary RMSE | Full RMSE | Status |
|---|---|---:|---:|---:|---|
| Baseline (pre-branch, no GSA) | 60 | 7 / 9 | 0.1040 | 0.1590 | `ACCEPT` (wrongly) |
| Iter 1 — volumes removed, GSA on | 60 | 7 / 9 | **0.0629** | 0.0827 | `REJECT` |
| Iter 2 — physiological ranges, balanced pulse | 300 | **8 / 9** | 0.0873 | 0.1026 | `PROMISING_NEAR_MISS` |
| Iter 3 — MAP form factor 0.40 | 300 | 6 / 9 | 0.1273 | 0.1743 | `PHYSIOLOGICAL_BUT_POOR_FIT` |

The baseline row is the important one: it was labelled `ACCEPT` while
`PAP_min` sat at 23.83% and two of nine governed metrics failed, because the
label only inspected the five selected primaries. Under the corrected
classifier no arm here reaches `ACCEPT`, which is the honest outcome.

### 3.3 The physiological range penalty was necessary

Iteration 1 removed the chamber volumes as targets, and GSA independently
dropped every LV chamber parameter from the active set. Nothing then
constrained the left ventricle, and the optimiser exploited it:

| Quantity | Iter 1 | Iter 2 | Paediatric band |
|---|---:|---:|---|
| `LVEF` | **0.891** | 0.659 | 0.40 – 0.85 |
| `LVESV` | **4.0 mL** | 15.6 mL | 1 – 80 mL |
| `LVEDV` direction vs post-op | **INCONSISTENT** | consistent (+11.9%) | — |
| Within paediatric range | 5 / 6 | **6 / 6** | — |

An `LVEF` of 0.891 with a 4 mL end-systolic volume is a near-empty ventricle:
numerically excellent, physiologically absurd. `clinical_reference_ranges` had
existed since the previous PR but was consulted only for the *baseline*
plausibility report, never in the objective and never on the calibrated
candidate.

Removing a target removes the obligation to match a measurement. It does not
remove the obligation to remain physiological. Iteration 2 adds that as an
explicit penalty over outputs with no patient comparator, and every chamber
prediction is now inside the paediatric band and directionally consistent with
shunt physiology.

### 3.4 The MAP form factor sensitivity refuted its own hypothesis

Iteration 2's single remaining governed failure was `SAP_min` at 17.40%. The
apparent explanation was attractive:

| | sys / dia | pulse | mean |
|---|---|---:|---|
| Model | 90.67 / 47.08 | **43.59** | 67.15 (true time average) |
| Clinical | 100 / 57 | **43.00** | 71.30 (one-third rule) |

The model matched the pulse to 1.4%, but its own MAP form factor is **0.460**
while the clinical `SAP_mean` was computed with the resting-adult one-third
rule (**0.333**). At 119 bpm the one-third rule is known to underestimate,
because systole occupies more of the cycle. So the three systemic targets look
mutually inconsistent under the model's definition of mean pressure.

PRD Phase 6 exists to test that rather than assert it.
`run_map_form_factor_sensitivity(0.40)` re-derives `SAP_mean` as
`57 + 0.40 x 43 = 74.20` — a form factor justified from heart rate, and
deliberately *not* the model's own 0.460, so the arm could not be read as
tuning the target until it passed.

**It made everything worse: 6 of 9, primary RMSE 0.1273.** The hypothesis is
refuted. Raising `SAP_mean` also raises the derived SVR target that the
systemic bundle anchors on, and `CO_Lmin` degraded from 8.46% to 12.73% as the
coupled constraints fought each other.

### 3.5 What actually causes the residual `SAP_min` error

Decomposing the iteration 2 systemic result gives two independent causes.

**Cause 1 — cardiac output is low, and mean pressure follows it.**
`MAP = SVR x CO + RAP` reproduces the model's mean exactly:

```
model    19.846 x 3.133 + 4.964 = 67.15 mmHg   (reported 67.15)
clinical 19.369 x 3.423 + 5.000 = 71.30 mmHg   (recorded 71.30)
```

SVR is accurate (+2.46%); `CO_Lmin` is **−8.46%**. Had CO reached target at the
model's own SVR, MAP would be 72.90 mmHg (+2.24%). CO sits inside its
deliberately widened uncertainty (±0.50 L/min, set because of the documented
catheter-versus-echo conflict) and inside the 10% band, so neither the systemic
bundle nor the gate hinge pushes it further.

**Cause 2 — waveform shape, independent of level.** Even with a corrected mean,
the model's 0.460 form factor puts diastolic pressure at 52.8 mmHg against a
measured 57 — still **−7.3%**. This is a systemic RC time-constant property
(`C.SAR` x `R_sys`), not a pressure-level property, and it is why moving the
MAP target could not fix it.

**A governance consequence.** Of the three systemic pressure rows, only two are
independent measurements. `SAP_sys = 100` and `SAP_dia = 57` are catheter
readings; `SAP_mean = 71.3` is a formula applied to them. Fitting all three as
independent targets counts one measurement pair twice and imports the
one-third assumption into a model that does not share it. The consistent
treatment — already applied to `Q_shunt_Lmin`, which is likewise an algebraic
identity — is to keep `SAP_mean` in the objective as a consistency term and
hold it out of the governed RMSE. **This branch does not make that change**,
because it would shrink the reported denominator, and that decision belongs in
review rather than in a commit that also reports a pass count.

### 3.6 Multi-start: 8 of 9 is a structural limit, not a search failure

Six deterministic scrambled-Sobol starts over the GSA-reduced 7-parameter box
(`UNIFIED_VSD_NUM_STARTS=6`, seed `20260828`, 300 evaluations each):

| Start | Label | Primary RMSE | Primary fail | Governed gate fail |
|---:|---|---:|---:|---:|
| 1 | `seed` | **0.087312** | 0 | **1** |
| 3 | `sobol_2` | 0.087682 | 0 | 2 |
| 6 | `sobol_5` | 0.091356 | 0 | 3 |
| 4 | `sobol_3` | 0.091689 | 0 | 2 |
| 2 | `sobol_1` | 0.095634 | 0 | 2 |
| 5 | `sobol_4` | 0.104133 | 0 | 2 |

**RMSE across starts: min 0.0873, median 0.0915, max 0.1041, IQR 0.0080.**

This is the first uncertainty interval this project has produced. Two things
follow, and the second matters more than the first.

**The optimisation is stable.** An IQR of 0.008 across independent starts means
the reported RMSE is a property of the problem, not of where the search
happened to begin. That is what the original `0.221945 -> 0.221945` Lundquist
arm could never have shown.

**No start reached 9 of 9.** Every one of six independent starts failed at
least one governed metric, and every one that failed exactly one failed
`SAP_min`. Combined with the decomposition in §3.5, this establishes that
`SAP_min` is not reachable by better search under the current target
definitions. It is a structural consequence of an 8.5% cardiac-output
shortfall plus a waveform form factor of 0.460 against a target derived
assuming 0.333.

**So the campaign stops at 8 of 9 deliberately.** The remaining routes to 9 of 9
are all denominator changes rather than model improvements — holding `SAP_mean`
out as the derived quantity it is, or widening the acceptance band. Both may be
defensible, and the `SAP_mean` holdout is arguably required for consistency
with `Q_shunt_Lmin`, but neither belongs in the same change that reports a pass
count. Reporting 8 of 9 with the reason is worth more than reporting 9 of 9
with a redefinition.

### 3.7 Published candidate

| Quantity | Value |
|---|---|
| Governed gate | **8 of 9 within 10%** |
| All clinical targets | 9 of 11 within 10% |
| Primary governed RMSE | **0.0873** (min 0.0873, max 0.1041 across 6 starts) |
| Full transparent RMSE | 0.1026 |
| Active parameters | 7 of 14, Sobol-screened |
| Chamber predictions | 6 of 6 inside paediatric bands, all directions consistent |
| Calibration status | `PROMISING_NEAR_MISS` |
| Outstanding | `SAP_min` 17.40%; `Q_shunt_Lmin` 21.57% (RMSE holdout, derived identity) |

Per-metric evidence:
[full_metric_gate_pre_surgery.csv](../results/runs/20260829_000626_reyna_pre_surgery/tables/full_metric_gate_pre_surgery.csv).
Figures: [calibration_error_comparison_pre_surgery.png](evidence/calibration_error_comparison_pre_surgery.png),
[calibration_error_slopes_pre_surgery.png](evidence/calibration_error_slopes_pre_surgery.png).

The candidate is **not** `ACCEPT`, and that is the correct label. It carries
five parameter-plausibility warnings and one governed metric outside the band.
Under the pre-branch classifier this same class of run was labelled `ACCEPT`.

## 4. Reading the denominator honestly

D7 removes `RVESV` from the governed set, taking it from 10 rows to 9. A
denominator that shrinks alongside a rising pass count is exactly the move this
project's own audit warned against, so the comparison below is stated on the
**fixed 9-row haemodynamic set** — the rows that are pre-operative measurements
or direct transforms of them — for both baseline and candidate:

`RAP_mean`, `PAP_min`, `PAP_max`, `PAP_mean`, `SAP_min`, `SAP_max`,
`SAP_mean`, `QpQs`, `CO_Lmin`.

On that fixed set the baseline scored **7 of 9** (`PAP_min` 23.83%,
`SAP_min` 13.98%). Per-metric, on the identical denominator:

| Metric | Baseline | Iter 2 | |
|---|---:|---:|---|
| `PAP_mean` | 3.34% | **0.19%** | ✅ |
| `SAP_mean` | 2.62% | 5.82% | ✅ |
| `RAP_mean` | 4.84% | **0.73%** | ✅ |
| `QpQs` | 0.20% | 2.33% | ✅ |
| `CO_Lmin` | 9.02% | 8.46% | ✅ |
| `PAP_max` | 8.00% | 8.64% | ✅ |
| `SAP_max` | 6.63% | 9.33% | ✅ |
| `PAP_min` | **23.83%** ✗ | **9.95%** | ✅ fixed |
| `SAP_min` | **13.98%** ✗ | **17.40%** ✗ | still failing |
| **Within 10%** | **7 / 9** | **8 / 9** | |

`PAP_min` — the pulmonary diastolic pressure that was graded but never fitted —
improved from 23.83% to 9.95% once it entered the objective. That single change
is the branch's clearest result, and it is the same metric the thesis reviewer
flagged in comment [295] when asking for explicit discussion of diastole.

The honest counterpoint is that `SAP_min` moved the wrong way (13.98% to
17.40%). Section 3.5 decomposes why: it is not a fitting failure but the
combination of an 8.5% cardiac-output shortfall and a waveform-shape mismatch,
neither of which the acceptance gate currently penalises, because CO sits
inside both its widened uncertainty and the 10% band.

## 5. Recommendation: the echo block is misfiled, not merely unusable

`LVEF` (32%), `LVESV` (31%), and `RVEDV` (29%) are not hard to fit. They are
**measuring a different physiological state**. Reyna's chambers at H+1 after VSD
closure are not Reyna's chambers with an open VSD; a model that reproduced them
under pre-operative loading would be wrong, not accurate. No optimiser budget,
parameter set, or objective reweighting can or should close that gap.

But the block is not worthless — it is filed under the wrong scenario. Checking
`patient_reyna()`:

- `clinical.pre_surgery` receives the H+1 echo block through the recipe, where
  its timing is invalid.
- `clinical.post_surgery` has **every field set to `NaN`** — the post-surgery
  scenario currently has no clinical validation targets at all and cannot be
  validated against anything.

The same five numbers that are invalid pre-operatively are exactly the
comparators the post-operative scenario lacks:

| Field | Value | Currently in | Belongs in |
|---|---:|---|---|
| `LVEDV_mL` | 41.0 | `pre_surgery` (invalid timing) | `post_surgery` |
| `LVESV_mL` | 19.3 | `pre_surgery` (invalid timing) | `post_surgery` |
| `RVEDV_mL` | 30.5 | `pre_surgery` (invalid timing) | `post_surgery` |
| `RVESV_mL` | 12.0 | `pre_surgery` (invalid timing) | `post_surgery` |
| `LVEF` / `EF` | 0.528 | `pre_surgery` (invalid timing) | `post_surgery` |

**This branch does not make that move.** Relocating patient measurements is a
clinical data-governance decision for the study owners, not something a
refactor should do silently, and it needs the H+1 timing confirmed against the
source record first. It is recorded here as the recommended next step.

Two caveats if it is taken up:

1. The consistency audit rates the block's internal stroke-volume disagreement
   as **critical (60%)**. That inconsistency travels with the data — moving it
   to `post_surgery` relocates the problem, it does not resolve it. The block
   would need the same tiering discipline there (`LVEF` demoted when both LV
   volumes are present, and now `RVEF` likewise).
2. H+1 is an early post-operative timepoint. Whether it represents the
   converged post-closure state, or a still-adapting one, is a clinical
   judgement that should be stated as a limitation either way.

The upside is real: it would convert the pipeline's most persistent liability
into the post-surgery scenario's first validation evidence.

## 6. What this does not establish

Unchanged by this branch, and still true:

- **n = 1.** One patient. Nothing here generalises to a cohort.
- **Degrees of freedom are barely positive.** GSA cut the active set to 7 against 8 independent
  measurements (`RAP_mean`, `PAP_dia`, `PAP_sys`, `PAP_mean`, `SAP_dia`,
  `SAP_sys`, `Qp`, `QpQs`); everything else is an algebraic transform. A gate
  count is a statement about flexibility until the active set is reduced.
- **A 6-start interval is not a 16-start interval.** PRD Phase 3 specifies
  `NumStarts >= 16`. Six starts give a usable IQR and settle the structural
  question, but the reported min/max should be quoted as a range across six
  starts, not as a 95% interval.
- **The GSA is one screen at one design.** `N = 128` at a single Halton design
  with no convergence check on the Sobol indices themselves. PRD Phase 4 asks
  for `N = 512` as a stability check before the reduced set is treated as
  settled.
- **Zhang vs Lundquist remains uninterpretable** at these budgets. The
  `OPTIMIZER_DID_NOT_MOVE` flag now marks the failure mode that made the
  original comparison unreadable, but neither arm has been run to convergence.
- **The clinical data conflict is documented, not resolved.** The consistency
  audit still returns `critical` at a 60% stroke-volume disagreement, and the
  cuff-versus-catheter MAP conflict (95 vs 71.3 mmHg) has not been run as a
  sensitivity arm (PRD Phase 6).

## 7. Reproduction

```bash
matlab -batch "cd('<repo>'); addpath(genpath(pwd)); setenv('UNIFIED_VSD_SCALING_MODE','zhang'); setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','1'); setenv('UNIFIED_VSD_DO_GSA','0'); setenv('UNIFIED_VSD_MAX_FUN_EVALS','60'); setenv('UNIFIED_VSD_MAX_ITERATIONS','8'); main_run('pre_surgery', patient_reyna())"
```

Every run writes `tables/full_metric_gate_<scenario>.csv` into its run folder,
and `.gitignore` now keeps that one file reviewable while continuing to ignore
the rest of `results/`.
