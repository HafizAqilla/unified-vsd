# Scientific Assessment: Is the Reyna Zhang Result Publishable?

> **Update, 2026-08-29 — three of this assessment's open items are now closed
> by the source clinical record, and one finding below is withdrawn.**
> The study "reyna" procedure log (06/04/2026; full provenance kept
> locally, not in this tracked file — see
> `config/private/patient_provenance.local.m`) was retrieved. See
> `docs/reyna_statistical_calibration_results_20260829.md` §0.
>
> - **§2.3's "60% internal stroke-volume inconsistency" is withdrawn as
>   stated.** The H+1 block is internally coherent: `SV_LV` = 41.0 − 19.3 =
>   21.7 and `SV_RV` = 30.5 − 12.0 = 18.5 mL/beat, 14.7% apart, with `LVEF`
>   checking exactly (21.7/41.0 = 0.529 vs 0.528 stated). The 60% arises only
>   from comparing those post-closure volumes against *pre*-closure flow-derived
>   stroke volumes — i.e. it is the shunt, the expected physiology, not a data
>   defect. With the block excluded the audit now returns severity `none`.
> - **§6's MAP question ("a reviewer will ask why MAP=95 became 71.3") is
>   answered by measurement, not a sensitivity arm.** The catheter stamps its
>   own mean: `RFA 100/57 (77)`. The form-factor reconstruction was
>   under-reading ~5–6 mmHg systematically (the same offset recurs
>   post-closure: formula 75 vs stamped 79). `SAP_mean` is now 77.
> - **§2.1's degrees-of-freedom problem is materially improved.** The same log
>   carries post-closure pressures from the same session, giving
>   `clinical.post_surgery` 7 finite targets where it had none. A joint
>   pre/post inversion now has `N = 16` against `p = 12`, i.e. `dof = 4`
>   rather than 0. Still not comfortable, but no longer structurally
>   impossible — and it is a genuine paired dataset, not relocated echo rows.
> - **Separately, two model *inputs* were found to be wrong** (`HR` 119 → 136,
>   BSA 0.6173 → 0.588), so every calibration number cited in this assessment
>   and its companion results documents is superseded and must be regenerated.
>   The reasoning in this assessment stands; its numbers do not.

Date: 2026-08-28
Subject: PR #23 `codex/luna-reyna-publishability-20260827`, fair-prior Zhang arm
Assessed candidate: governed primary RMSE `0.060445`, full RMSE `0.063985`,
calibration status `ACCEPT`
Companion document: [reyna_zhang_full_metric_prd.md](reyna_zhang_full_metric_prd.md)

---

## Verdict

**As evidence that the model reproduces Reyna's pre-operative haemodynamics —
weak but salvageable. As evidence that Zhang scaling is superior to Lundquist —
not supportable as currently run.**

The `ACCEPT` status is real within its own governance rules, and the governance
machinery in this repository is unusually careful. But three problems sit
between the current result and a defensible publication claim, and only one of
them is about the number `0.060445`.

| Question | Answer |
|---|---|
| Is the primary RMSE genuinely low? | Yes, mechanically. |
| Does it demonstrate model validity? | **No** — there are more free parameters than independent measurements. |
| Does it demonstrate Zhang > Lundquist? | **No** — the Lundquist arm's optimiser never moved. |
| Is the 10% claim complete? | **No** — it covers 5 of 16 rows. |
| Can this be fixed? | Yes. Mostly without new patient data. |

---

## 1. What the acceptance claim actually covers

PR #23 reports: *"all five primary metrics within 10%"* and lists `RAP_mean`
1.37%, `PAP_mean` 3.45%, `SAP_mean` 3.50%, `QpQs` 1.52%, `CO_Lmin` 6.18%.

Executing the target and tier construction directly against `patient_reyna()`
with the Reyna recipe applied gives the real denominator:

- **16** metric rows carry a finite clinical comparator
- **11** are inside the governed primary RMSE mask
- **9** are actually inside the calibration objective
- **5** are covered by the published acceptance claim

So the headline claim covers 31% of the rows that have a clinical comparator and
45% of the rows that the reported RMSE is computed over. The remaining rows are
not disclosed per-metric anywhere in the PR, the results document, or the
committed evidence — the run folders under `results/luna_experiments/` and
`results/runs/` are gitignored and absent from the repository.

**The per-metric full table for the accepted candidate does not currently exist
in reviewable form.** That is the first thing to fix, and it costs no modelling
work.

> **Update, 2026-08-28 (later same day).** Phase 0 has since been implemented
> and run. The measured per-metric table is in §1.2 and it confirms the
> inference below: `PAP_min` sits at **23.83%**. The estimates in §1.1 are
> retained because they show what the published evidence did and did not
> permit. Corrections found by execution are marked in place.

### 1.1 What can be recovered from the published aggregates

The two RMSE figures plus the five per-metric errors are enough to bound the
rest by arithmetic. RMSE here is `sqrt(mean(pct_err^2))` over the mask, so
sums of squares are recoverable:

| Group | n | Sum of squares | Implied RMS abs error |
|---|---:|---:|---:|
| Published gated metrics | 5 | 0.00665 | **3.65%** |
| In primary RMSE, not published | 6 | 0.03354 | **7.48%** |
| Finite comparator, outside primary RMSE | 5 | 0.02532 | **7.12%** |
| All finite rows | 16 | 0.06551 | 6.40% |

*(Computed at `N_primary = 11`. At `N_primary = 10` — if the recipe's
`Q_shunt_Lmin` holdout was honoured in the production path — the figures are
7.73% and 6.95%. The conclusion is unchanged.)*

The metrics behind the acceptance claim fit about **twice as well** as the ones
behind it. That is not surprising, and it is not fraud — it is what happens when
five metrics are named primary and the rest are soft, unfitted, or excluded. But
it means the headline understates the model's true error profile, and the
selection of which five to gate was made by the same pipeline that reports the
result.

Can we say a metric is over 10%? **No, and it would be wrong to claim so.** The
sum-of-squares budget for 11 metrics all at exactly 10% is 0.1100; the observed
value is 0.0402. The aggregate is entirely compatible with every metric passing.
It is equally compatible with one metric at 14% and four at 5%. **The published
evidence cannot distinguish these, which is precisely the problem.**

### 1.2 The measured table

Running the pipeline with a per-metric exporter in place (fair-prior Zhang,
`MaxFunctionEvaluations = 60`, seeds disabled, GSA off) gives the distribution
directly:

| Metric | \|Error\| | Tier | In objective | In governed RMSE |
|---|---:|---|:--:|:--:|
| `LVEF` | 32.02% | consistency_check_only | no | no |
| `LVESV` | 30.68% | consistency_check_only | no | no |
| `RVEDV` | 28.83% | consistency_check_only | no | no |
| `PAP_min` | **23.83%** | validation_only | **no** | **yes** |
| `SAP_min` | 13.98% | soft | yes | yes |
| `RVESV` | 9.40% | soft | yes | yes |
| `CO_Lmin` | 9.02% | hard | yes | yes |
| `PAP_max` | 8.00% | validation_only | **no** | **yes** |
| `Q_shunt_Lmin` | 7.86% | soft | yes | no |
| `LVEDV` | 7.71% | consistency_check_only | no | no |
| `SAP_max` | 6.63% | soft | yes | yes |
| `SVR` | 6.42% | derived_validation | no | no |
| `RAP_mean` | 4.84% | hard | yes | yes |
| `PAP_mean` | 3.34% | hard | yes | yes |
| `SAP_mean` | 2.62% | hard | yes | yes |
| `QpQs` | 0.20% | hard | yes | yes |

**11 of 16 overall; 8 of 10 in the governed primary RMSE set.**

Two corrections to §1.1 follow from having executed rather than inferred:

- The governed mask on the **production** path is **10 rows**, not 11 —
  `make_recipe_target_tier_config` does thread the recipe's
  `primary_rmse_holdout`, so `Q_shunt_Lmin` is correctly excluded. The 11-row
  count appears only on a bare `build_target_tiers` call. The defect is a
  disagreement between code paths, not a wrong published denominator.
- This run reached primary RMSE `0.1040` and full RMSE `0.1590`, not the
  published `0.060445` / `0.063985`. It was invoked through `main_run`
  directly rather than through `run_reyna_scaling_experiment`, which sets
  additional environment controls, so this is **not** a failed reproduction —
  but it does mean the published figure has not been independently reproduced
  here, and the run-to-run spread at this budget is large.

And the decisive one:

> This run was classified **`ACCEPT`** with `primary10_fail=0`, while `PAP_min`
> stood at 23.83% and two of ten governed metrics were outside the band.

The acceptance label was computed over the five selected primaries only. That
is the governance hole, demonstrated rather than argued.

---

## 2. The three structural problems

### 2.1 Negative degrees of freedom

The calibration fits **14 parameters**
(`recipe.active_parameters` in [reyna_pre_surgery.m:54](config/calibration_recipes/reyna_pre_surgery.m:54)).

The independent measurement content of the Reyna pre-surgery record is **eight
numbers**: `RAP_mean` 5, `PAP_dia` 10, `PAP_sys` 20, `PAP_mean` 15, `SAP_dia`
57, `SAP_sys` 100, `Qp` 4.087, `QpQs` 1.194.

Everything else in the 16-row table is an algebraic transform of those eight, or
belongs to the contested echo block:

| Reported row | Actually is | Check |
|---|---|---|
| `SAP_mean` = 71.3 | `dia + (sys-dia)/3` | `57 + 43/3` = **71.33** |
| `CO_Lmin` = 3.423 | `Qp / QpQs` | `4.087/1.194` = **3.42295** |
| `Q_shunt_Lmin` = 0.664 | `CO x (QpQs - 1)` | `3.423 x 0.194` = **0.66406** |
| `SVR` = 19.369 | `(SAP_mean - RAP_mean)/CO` | `66.3/3.423` = **19.369** |
| `VSD_frac_pct` | `100(1 - 1/QpQs)` | **16.25%** |
| `LVEDV/LVESV/RVEDV/RVESV/LVEF` | H+1 post-operative echo | see §2.3 |

Fitting 14 parameters to 8 independent numbers is structurally
under-determined. Under that condition **a low RMSE is guaranteed to be
attainable and carries no information about whether the model is right.** It
measures flexibility, not fidelity. This is the most important single sentence
in this assessment.

The correct response is not to abandon the result — it is to (a) run the GSA and
reduce the active set, (b) report parameter confidence intervals rather than
point estimates, and (c) stop counting algebraic identities as separate
validation successes. `Q_shunt_Lmin` matching to within a few percent is not an
independent test; it is arithmetic that must hold if `CO` and `QpQs` hold.

### 2.2 The Zhang-versus-Lundquist comparison is not interpretable

The results document reports, for the fair Lundquist arm, baseline primary RMSE
`0.221945` and calibrated primary RMSE `0.221945`. Identical to six decimal
places.

That is not a physiological finding about the Lundquist prior. **That is an
optimiser that returned its starting point.** Three "deterministic repeats"
returning the same value confirm the same thing rather than establishing
repeatability — repeating a solve that does nothing reproduces nothing.

The optimiser budget makes this predictable rather than mysterious:
`MaxFunctionEvaluations = 60` and `MaxIterations = 8` against **14 free
parameters**. One forward-difference gradient in 14 dimensions costs 15
evaluations. The whole confirmation budget is approximately **four gradient
steps**. Neither arm converged; Zhang simply started closer.

So the honest reading of PR #23's headline comparison is:

> Under a budget too small to converge either arm, the Zhang-scaled starting
> point was closer to the clinical targets than the Lundquist-scaled starting
> point.

That is a legitimate and mildly interesting statement about **priors**, and it
is worth reporting as such. It is not a statement about which scaling law is
correct, and the PR's own framing — "does not prove Zhang is the generally
correct pediatric scaling law" — is appropriately cautious. The results
document's interpretation section is, to its credit, more careful than the PR
title.

### 2.3 Post-operative evidence is shaping a pre-operative fit

[patient_reyna.m:81](config/patient_reyna.m:81) is explicit: the LV/RV volume
block "was confirmed to be H+1 after surgery, so it is not a valid pre-surgery
calibration target." It sets all five volume/EF fields to `NaN` and
`override_IC = false`.

[reyna_pre_surgery.m:38](config/calibration_recipes/reyna_pre_surgery.m:38) then
puts them back: `LVEDV = 41.0`, `LVESV = 19.3`, `RVEDV = 30.5`, `RVESV = 12.0`,
`LVEF = 0.528`, `override_IC = true`.

Four are demoted to `consistency_check_only`, which is good governance. But
**`RVESV` remains tier `soft` and is actively fitted**, and `override_IC = true`
seeds the initial state vector from the same post-operative block. Two channels,
both carrying H+1 evidence into a pre-operative calibration.

Meanwhile the consistency audit returns:

> severity `critical`; max relative stroke volume difference = **60.0%**

A 60% internal inconsistency in the evidence block that is seeding initial
conditions is a serious limitation. The PR acknowledges it in one line
("unchanged and remains a data-governance limitation"). For publication it needs
quantification and a decision, not acknowledgement.

---

## 3. Smaller but real defects

**Two measured targets are graded but never fitted.** `PAP_min` and `PAP_max`
are declared `UseForCalibration = true`, `Reliability = 'High'` in
[get_calibration_targets.m:36](src/utils/get_calibration_targets.m:36), but
appear in neither the hard nor soft list in
[build_target_tiers.m:153](src/calibration/build_target_tiers.m:153). They fall
through to `validation_only` — excluded from the objective, **included in the
reported primary RMSE**. The optimiser is scored on two directly measured
catheter pressures it is never asked to fit. They influence `J` only through
`PAP_pulse` at weight 0.20. This is the highest-yield fix available and costs
two lines.

**The tier policy differs between code paths.** The recipe declares
`primary_rmse_holdout = {'Q_shunt_Lmin'}`; `default_target_tier_config`
declares `{}`. The production path honours the recipe (10 rows); a bare
`build_target_tiers(clinical, scenario)` call does not (11 rows). Two code
paths, two denominators, and nothing asserting they agree.

**The same ungoverned-target defect affects three more metrics.** Sweeping all
13 patient profiles across both scenarios with an assertion in place shows
`LAP_mean` and `RVEDV` (pre-surgery) and `RVEF` (post-surgery) are likewise
declared `UseForCalibration` yet land in `validation_only`. `LVEF` already had
a rule demoting it to consistency-only when both LV volumes are present;
`RVEF` had no equivalent, so the right heart was governed differently from the
left for no stated reason.

**Primary-metric selection disagrees between paths.**
`select_primary_metrics(clinical, [], 'pre_surgery')` returns `RAP_mean,
PAP_min, PAP_mean, SAP_mean, QpQs` — substituting `PAP_min` for `CO_Lmin`. The
recipe and the published run use `CO_Lmin`. Nothing asserts the two agree.

**The objective does not minimise the reported metric.** The report is an
unweighted RMS of percent errors; `J` is a weighted sum of squared normalised
errors plus z-score bundles, waveform penalties, guards, plausibility terms and
regularisation. Nothing in `J` specifically rewards pulling a worst-case metric
under 10%, and the guard thresholds (`0.12` for pressures, `0.08` for mean
pressures) straddle the 10% gate inconsistently.

**No uncertainty, no holdout.** There are no confidence intervals on
`0.060445`. Every finite-comparator row is either fitted or reported; none is
held out. The document correctly states the three repeats "are a deterministic
repeatability check, not an uncertainty interval" — which means the result
currently has no uncertainty interval at all.

**Five parameter-plausibility warnings** on the accepted candidate. Inside
registry bounds, but flagged. With negative degrees of freedom, plausibility
warnings are the main remaining signal that a fit is being achieved by
physiologically awkward parameter values.

---

## 4. What is genuinely good here

This should not read as a demolition. Several things in this codebase are better
than typical practice in the lumped-parameter cardiovascular literature:

- **Explicit target tiering** with `hard` / `soft` / `consistency_check_only` /
  `derived_validation` / `validation_holdout`, and machine-readable reasons for
  every exclusion. Most published models do not disclose which targets were
  fitted at all.
- **`LVEF` auto-demoted** when `LVEDV` and `LVESV` are both present, explicitly
  to avoid double-counting one echo block. That is exactly right, and it is
  reasoning most papers never perform.
- **Dual RMSE reporting** (`primary_governed` and `full_transparent`), with the
  full figure never hidden. The full RMSE being *worse* than the primary is
  disclosed rather than suppressed.
- **The GSA preflight failed loudly** rather than silently degrading. The runner
  refused to present a calibration-only result as GSA-complete. That is the
  correct behaviour and it was reported honestly in the PR.
- **The clinical consistency audit runs and returns `critical`** rather than
  being switched off once it became inconvenient.
- **The results document's interpretation section** is materially more cautious
  than the PR title, and explicitly refuses the "Zhang is universally superior"
  reading.

The governance scaffolding is strong. The problems are that the scaffolding is
not yet fully wired to the objective, and that the optimiser budget is too small
to make any comparison meaningful.

---

## 5. Recommended path to a defensible claim

In dependency order. The full specification is in
[reyna_zhang_full_metric_prd.md](reyna_zhang_full_metric_prd.md).

| Step | Effort | Cost if skipped |
|---|---|---|
| **0.** Export and commit the per-metric full table | hours | Every claim below is unverifiable |
| **1.** Put `PAP_min`/`PAP_max` in the objective; assert no ungoverned target | hours | Graded on unfitted metrics; RMSE means two different things |
| **2.** Add a 10% gate hinge to `J`; align guard thresholds to the recipe gate | 1 day | Optimiser is not minimising the reported quantity |
| **3.** Raise budget to ≥3000 evals, ≥16 multi-starts; flag `OPTIMIZER_DID_NOT_MOVE` | 1 day + compute | **The Zhang-vs-Lundquist comparison stays uninterpretable** |
| **4.** Install UQLab, run GSA at N=128, reduce to ≤7 parameters, report CIs | 2–3 days | Negative DOF stands; no validity claim is possible |
| **5.** Resolve the `patient_reyna` / recipe contradiction; drop `RVESV`, `override_IC=false` | 1 day | Post-operative data silently shapes a pre-operative fit |
| **6.** Sensitivity arms on MAP form factor (0.40) and the cuff/catheter conflict | 1 day | A reviewer will ask why `MAP=95` became `71.3` |

Steps 0–3 need no new data and no new software. **Step 3 is the one that
converts the headline from an artefact into a finding.**

### 5.1 What to claim once this is done

If steps 0–6 complete and the gate count holds, the defensible claim is
approximately:

> For a single 3-year-old patient with a restrictive VSD, a Zhang-scaled
> parameter prior combined with governed multi-start calibration reproduced
> *n* of 16 clinical comparators within 10%, with a governed primary RMSE of
> *x* (95% CI *a–b* across 16 starts). The active parameter set was reduced from
> 14 to *k* by total-order Sobol screening. Results were robust to the systemic
> mean-pressure form factor. Pre-operative chamber volumes were unavailable; the
> available echo block was post-operative (H+1) with a 60% internal
> stroke-volume inconsistency and was excluded from fitting.

That is a modest, honest, publishable claim. It is much stronger than the
current framing precisely because it states its own limits.

### 5.2 What not to claim

- That Zhang scaling is superior to Lundquist — **n=1, and neither arm
  converged.**
- That low RMSE validates the model — **not with 14 parameters and 8
  independent measurements.**
- That `Q_shunt_Lmin`, `SVR`, or `VSD_frac_pct` agreement are independent
  validation successes — they are algebraic identities.
- That three identical repeats establish repeatability — they establish that a
  solve which did nothing did nothing three times.

---

## 6. On the specific goal of "as many metrics under 10% as possible"

Worth stating plainly, because it cuts against the framing of the request:
**maximising the count of metrics under 10% is not by itself a scientific
objective, and optimising for it directly is a way to get a worse paper.**

With 14 parameters and 8 independent measurements, the count can very likely be
driven to 16 of 16. That would demonstrate nothing except that the model has
enough knobs. A reviewer who understands parameter counting will treat a perfect
score as evidence of overfitting rather than accuracy.

What makes the count meaningful is the surrounding structure:

1. **The count must be over a disclosed denominator.** `n of 16`, with the full
   table published — not `5 of 5` against a set selected by the same pipeline.
2. **Metrics that are graded must be fitted.** Otherwise the count is measuring
   luck (Phase 1).
3. **The parameter set must be smaller than the independent data**, or the count
   must be accompanied by confidence intervals wide enough to be honest
   (Phase 4).
4. **Algebraic identities must not be counted as passes** (Phase 1 decision on
   `Q_shunt_Lmin`).
5. **The result must survive the sensitivity arms** (Phase 6). A gate count
   stable across the MAP form factor and the cuff/catheter conflict is worth far
   more than a higher count at one arbitrary data choice.

So: pursue the goal, via the PRD — but report it as *n of 16 with k active
parameters and stated confidence intervals*, never as a bare count. Under those
conditions a result of 13 or 14 of 16 with 7 parameters is a considerably
stronger paper than 16 of 16 with 14 parameters.

---

## 7. Immediate next action

Run Phase 0. It changes no modelling behaviour, needs no UQLab, and produces the
one artefact every statement in this assessment is currently forced to
approximate: **the per-metric error table for the accepted `0.060445`
candidate.**

Until that table exists, neither the optimistic reading nor the pessimistic one
can be settled.
