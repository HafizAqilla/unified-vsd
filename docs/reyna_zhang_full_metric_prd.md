# PRD: Reyna Zhang Full-Metric 10% Campaign

Status: implementation-ready
Version: `reyna_zhang_fullmetric_v1`
Scope: Reyna `pre_surgery`, Zhang scaling prior
Supersedes for this goal: nothing. Extends `reyna_p1_scaling_v1`.

---

## 1. Objective

Maximise the number of **full** validation metrics whose absolute error is
within 10%, for the Reyna pre-surgery Zhang candidate, **without** loosening
target governance, hiding rejected candidates, or expanding the free-parameter
count.

The current published claim is that "all five primary metrics are within 10%".
That claim covers 5 of 16 rows carrying a finite clinical comparator. This PRD
governs the other 11.

### 1.1 Non-goals

- Lowering the reported RMSE by removing hard metrics from the RMSE mask.
- Re-injecting H+1 post-operative echo volumes as valid pre-surgery evidence.
- Declaring Zhang the generally correct paediatric scaling law.
- Promoting a candidate whose parameters sit on registry bounds.

### 1.2 Success definition

The campaign succeeds when a single Zhang candidate satisfies **all** of:

| # | Criterion | Threshold |
|---|---|---|
| S1 | Metrics in the governed primary RMSE mask within 10% | 11 of 11 |
| S2 | Governed primary RMSE | ≤ 0.090 (recipe `primary_rmse_max`) |
| S3 | Full transparent RMSE | ≤ 0.100 |
| S4 | Gated primary metrics within 10% | 5 of 5, unchanged |
| S5 | Parameters strictly inside registry bounds | all, no bound-parking |
| S6 | Parameter-plausibility warnings | ≤ 2, each with written justification |
| S7 | Per-metric full table exported and committed | required |
| S8 | Result reproducible from a recorded seed and contract | required |

A candidate meeting S1–S5 but not S6 is a **sensitivity result**, not the
publication candidate.

---

## 2. Ground truth this PRD is built on

Verified by executing `get_calibration_targets` + `build_target_tiers` against
`patient_reyna()` with `config/calibration_recipes/reyna_pre_surgery.m` applied
(MATLAB R2025a, 2026-08-28). **16 rows carry a finite clinical comparator:**

| Metric | Unit | Clinical | Tier | In objective | In primary RMSE |
|---|---|---:|---|:--:|:--:|
| `RAP_mean` | mmHg | 5.000 | hard | yes | yes |
| `PAP_min` | mmHg | 10.000 | **validation_only** | **no** | **yes** |
| `PAP_max` | mmHg | 20.000 | **validation_only** | **no** | **yes** |
| `PAP_mean` | mmHg | 15.000 | hard | yes | yes |
| `SAP_min` | mmHg | 57.000 | soft | yes | yes |
| `SAP_max` | mmHg | 100.000 | soft | yes | yes |
| `SAP_mean` | mmHg | 71.300 | hard | yes | yes |
| `QpQs` | – | 1.194 | hard | yes | yes |
| `Q_shunt_Lmin` | L/min | 0.664 | soft | yes | yes |
| `SVR` | WU | 19.369 | derived_validation | no | no |
| `CO_Lmin` | L/min | 3.423 | hard | yes | yes |
| `LVEDV` | mL | 41.000 | consistency_check_only | no | no |
| `LVESV` | mL | 19.300 | consistency_check_only | no | no |
| `RVEDV` | mL | 30.500 | consistency_check_only | no | no |
| `RVESV` | mL | 12.000 | soft | yes | yes |
| `LVEF` | – | 0.528 | consistency_check_only | no | no |

Counts on the bare default path: `N_FINITE=16`, `N_IN_PRIMARY_RMSE=11`,
`N_IN_CALIBRATION=9`. On the production path the recipe's
`primary_rmse_holdout` is honoured, so `Q_shunt_Lmin` leaves the mask and the
governed set is **10 rows**.

Clinical consistency audit returns **severity `critical`, max relative stroke
volume difference 60.0%**.

### 2.0 Measured baseline (2026-08-28, fair-prior Zhang, 60 evals)

A run at the published controls, with the Phase 0 exporter in place, gives the
per-metric distribution the original report never published:

| Metric | \|Error\| | Tier | In objective | In primary RMSE |
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

**Headline: 11 of 16 overall, 8 of 10 in the governed primary RMSE set.**

Two conclusions follow directly:

1. `PAP_min` at 23.83% is the single largest governed failure, and it is
   exactly the metric the objective never sees. This is Defect A, measured.
2. The three worst rows overall (`LVEF`, `LVESV`, `RVEDV`, 28–32%) are the H+1
   post-operative echo block. They are already excluded from the governed
   RMSE, and they are not valid pre-operative targets, so "16 of 16" is the
   wrong goal. The correct goal is **all governed haemodynamic rows within
   10%**, with the echo block reported as excluded cross-timing evidence.

### 2.1 Defect A — measured targets are graded but never fitted

`PAP_min` and `PAP_max` are declared `UseForCalibration = true` and
`Reliability = 'High'` in [get_calibration_targets.m:36](src/utils/get_calibration_targets.m:36).
They appear in neither `config.hard` nor `config.soft` in
[build_target_tiers.m:153](src/calibration/build_target_tiers.m:153), so they fall
through to `validation_only`: excluded from the objective, but **included in the
governed primary RMSE**. The optimiser is scored on two catheter pressures it is
never asked to fit. They reach the objective only indirectly, through
`PAP_pulse` at weight 0.20 in `pressure_waveform_penalty`.

This is the single highest-yield defect for this PRD's goal.

**The defect is not limited to Reyna.** Sweeping all 13 patient profiles across
both scenarios with an assertion in place found the same class of failure in
three further targets, all declared `UseForCalibration` and all landing in
`validation_only`:

| Target | Scenario | Why it matters |
|---|---|---|
| `PAP_min`, `PAP_max` | pre_surgery | directly measured catheter pressures, High reliability |
| `LAP_mean` | pre_surgery | measured atrial pressure wherever a profile carries it |
| `RVEDV` | pre_surgery | measured RV volume, the right-heart counterpart of `LVEDV` |
| `RVEF` | post_surgery | derived from RVEDV/RVESV, so it also needs the `LVEF` anti-double-counting rule |

`LVEF` already had a rule demoting it to consistency-only when both LV volumes
are present. `RVEF` had no equivalent, so the right heart was governed
differently from the left for no stated reason.

### 2.2 Defect B — the tier policy differs between code paths

`reyna_pre_surgery.m` declares `primary_rmse_holdout = {'Q_shunt_Lmin'}`.
`make_recipe_target_tier_config` in `build_case_calibration_profile.m` does
thread that through, so **the production path honours it** and the governed set
is 10 rows. But a bare `build_target_tiers(clinical, scenario)` call falls back
to `default_target_tier_config()`, where `primary_rmse_holdout = {}`, giving 11
rows. Two code paths, two denominators for the same headline number, and no
assertion that they agree.

### 2.2b Defect B2 — `primary_rmse_holdout` also removes a metric from the objective

`merge_allowed_metrics_with_tiers` in `build_case_calibration_profile.m`
subtracts `target_tiers.excluded_from_primary_rmse` from
`profile.allowedMetricFields`, and that set includes `primary_rmse_holdout`.
`allowedMetricFields` gates `calib.metricFields`, which is what the objective
iterates over.

So declaring `Q_shunt_Lmin` an RMSE holdout silently removed it from the fit as
well. The recipe's intent — *inform the objective as a consistency term, but do
not count as an independent RMSE row* — was not achievable. The tier table
reported `IncludedInCalibration = true` for a metric the objective never saw,
so the governance table was describing a fit that was not happening.

The exclusion must cover only genuinely report-only tiers
(`consistency_only`, `derived_validation`, `validation_holdout`), never
`primary_rmse_holdout`.

### 2.3 Defect C — primary-metric selection disagrees across paths

`select_primary_metrics(clinical, [], 'pre_surgery')` without a `caseProfile`
returns `RAP_mean, PAP_min, PAP_mean, SAP_mean, QpQs` — it substitutes `PAP_min`
for `CO_Lmin`. The recipe declares `RAP_mean, PAP_mean, SAP_mean, QpQs,
CO_Lmin`, and the published results used the recipe list. Two code paths, two
answers, no assertion guarding the difference.

### 2.4 Why the aggregates could not settle this

Before Phase 0, the only available evidence was the published primary RMSE
`0.060445`, full RMSE `0.063985`, and five per-metric errors. Those decompose
into group RMS values of roughly 3.65% (published), 7.48% (governed but
undisclosed), and 7.12% (outside the governed mask) — enough to show the
undisclosed rows fit about twice as badly, but **not** enough to name a failing
metric: the sum-of-squares budget for 11 metrics at 10% is 0.1100 against an
observed 0.0402, so "all pass" and "one at 14%" were equally consistent with
the published numbers.

The measured table in §2.0 settles it: `PAP_min` is at 23.83%. The aggregate
was compatible with that all along, which is precisely why aggregates are not
an acceptable substitute for the per-metric table.

---

## 3. Phased implementation

Each phase is independently mergeable and independently revertable. Do not start
a phase before its predecessor's acceptance tests pass.

### Phase 0 — Publish the per-metric truth (blocking, no modelling change)

**Goal:** make the distribution visible. No optimiser, objective, or target
change is permitted in this phase.

**Changes**

1. New `src/utils/export_full_metric_gate.m`. Input: the `report` struct from
   `validation_report`. Output: a CSV with one row per finite-comparator metric
   and columns
   `Metric, Unit, Tier, Clinical, Model, Error_pct, AbsError_pct, InObjective,
   InPrimaryRMSE, Within10pct, Within5pct, ObjectiveWeight`.
2. Call it from `main_run.m` immediately after `validation_report`, writing to
   `<run_dir>/full_metric_gate_<scenario>.csv`.
3. Add a console block `--- FULL METRIC 10% GATE (n of N) ---` printing the
   count within 10% over all finite rows, and naming every row outside it.
4. Extend `run_reyna_scaling_experiment.m` to copy that CSV into the experiment
   folder and to add `MetricsWithin10pct` and `WorstMetric` /
   `WorstAbsError_pct` columns to the typed summary CSV.
5. Un-ignore the gate CSV specifically, so results are reviewable in git:
   add `!results/**/full_metric_gate_*.csv` to `.gitignore`.

**Acceptance tests**

- `tests/test_full_metric_gate_export.m`: given a synthetic `report` with known
  errors, the CSV has the expected row count, `Within10pct` flags, and worst-row
  identification.
- Re-run the fair Zhang confirmation at the *published* controls
  (`MaxFunctionEvaluations=60`, `MaxIterations=8`, polish on, GSA off, seeds
  disabled). Commit the resulting gate CSV as
  `docs/evidence/reyna_zhang_fullmetric_baseline_<date>.csv`.

**Exit criterion:** the per-metric table for the published `0.060445` candidate
is committed and reviewable. **Do not proceed until this exists.**

---

### Phase 1 — Fix target governance (Defects A, B, C)

**Goal:** every metric that is graded is also fitted, and every code path agrees
on the target set.

**Changes**

1. `build_target_tiers.m` — add `PAP_min` and `PAP_max` to `config.soft` in
   `default_target_tier_config()`. Add weight multipliers
   `PAP_max = 0.50`, `PAP_min = 0.45` to `config.metric_weight_multipliers`
   (mirroring the existing `SAP_max` / `SAP_min` treatment).
2. `reyna_pre_surgery.m` — add `'PAP_max','PAP_min'` to `recipe.soft_metrics`
   so the recipe is explicit rather than relying on the default.
3. Add an **invariant assertion** in `build_target_tiers.m`: any target with
   `UseForCalibration = true` and a finite clinical value must resolve to
   `hard`, `soft`, or an explicitly named exclusion tier. A silent fall-through
   to `validation_only` raises
   `build_target_tiers:ungovernedCalibrationTarget`.
4. Thread the recipe's tier policy explicitly into `build_target_tiers` from
   `run_calibration.m`, so `primary_rmse_holdout` is honoured on every path.
   Record the resolved policy in the run manifest as `TargetTierPolicySource`.
5. Add an assertion that `select_primary_metrics` with the case profile returns
   exactly `recipe.primary_metrics`, and fail the run loudly if not.

**Decision required — `Q_shunt_Lmin`:** it is algebraically
`CO_Lmin x (QpQs - 1)` = `3.423 x 0.194` = `0.66406`, which reproduces the
recorded `0.664` exactly. It is not an independent measurement. **Recommendation:
keep it in the objective as a soft consistency term, exclude it from the primary
RMSE** (honour the recipe holdout), and state the reason in the report. Do not
count a derived identity as an independent validation success.

**Acceptance tests**

- Extend `tests/test_validation_target_governance.m`: assert `PAP_min` and
  `PAP_max` resolve to `soft` and are `IncludedInCalibration`.
- New `tests/test_no_ungoverned_calibration_targets.m`: for every patient
  profile in `config/`, no `UseForCalibration` target with a finite value lands
  in `validation_only`.
- New `tests/test_primary_metric_path_agreement.m`: recipe list and
  `select_primary_metrics` output agree for Reyna pre-surgery.

**Expected effect:** the primary-RMSE mask becomes 10 rows (Q_shunt removed),
all 10 of which are now in the objective. Re-run and record the new baseline.
The primary RMSE may *rise* — that is correct and expected, because two
previously unfitted metrics are now honestly scored against an objective that
finally sees them.

---

### Phase 2 — Align the objective with the 10% gate

**Goal:** make the optimiser minimise the thing that is being reported.

The reported RMSE is an unweighted RMS of percent errors. The objective is a
weighted sum of `(err/target)^2` plus z-score bundles, waveform penalties,
guards, and regularisation. Minimising `J` does not minimise the reported RMSE,
and nothing in `J` specifically rewards pulling a worst-case metric under 10%.

**Changes**

1. `objective_calibration.m` — add an explicit **gate hinge** over every metric
   in the primary-RMSE mask:

   ```
   J_gate = lambda_gate * sum_m w_m * max(0, |err_rel_m| - gate)^2 / gate^2
   ```

   with `gate = 0.10` from `recipe.acceptance.primary_gate_pct` (never
   hard-coded) and `lambda_gate` exposed as `calib.gateLambda`, default `8.0`.
   This is zero for compliant metrics and grows quadratically past the gate, so
   it reshapes the basin only where it matters.
2. Align the existing `soft_excess_penalty` thresholds in
   `clinical_guard_penalty` to the governed gate. `pressure_rel_thresh = 0.12`
   and `mean_pressure_rel_thresh = 0.08` currently straddle the 10% gate
   inconsistently. Drive both from `recipe.acceptance`, defaulting to `0.10`.
3. Add an optional **minimax polish stage**: after the main solve, run a short
   local refinement on `max_m |err_rel_m|` restricted to
   `recipe.systemic_polish_parameters`. Accept its result only if it improves
   the worst metric **and** does not worsen primary RMSE by more than 5%
   relative. Record both candidates.
4. Log per-term objective contributions (`J_primary`, `J_secondary`,
   `J_systemic`, `J_pressure`, `J_shunt`, `J_gate`, `J_reg`,
   `J_param_plausibility`, `J_boundary`) to
   `<run_dir>/objective_decomposition.csv` at the accepted point.

**Acceptance tests**

- `tests/test_gate_hinge_penalty.m`: hinge is exactly 0 at 9.99% error, strictly
  positive at 10.01%, and continuous and C¹ at the knot.
- `tests/test_objective_decomposition_export.m`: terms sum to the reported `J`
  within `1e-9`.
- Regression: the Phase 1 candidate re-scored under the Phase 2 objective must
  not change its metric values (objective change only, no model change).

---

### Phase 3 — Make the optimisation actually converge

**Goal:** remove the confound that currently makes every scaling-prior
comparison uninterpretable.

The published controls are `MaxFunctionEvaluations = 60` and
`MaxIterations = 8` for a **14-parameter** calibration. A single forward-
difference gradient in 14 dimensions costs 15 evaluations. The entire budget is
therefore about four gradient steps. Independently, the fair Lundquist arm
reported baseline primary RMSE `0.221945` and calibrated primary RMSE
`0.221945` — identical to six decimals, the signature of an optimiser that never
left its start point.

**Changes**

1. Raise the confirmation budget to `MaxFunctionEvaluations >= 3000`,
   `MaxIterations >= 200`. Keep a separate, clearly labelled screening budget for
   smoke runs; never publish a screening number as a confirmation.
2. Add multi-start: `NumStarts >= 16` sampled by scrambled Sobol or Latin
   hypercube over the registry bounds, each with a recorded seed. Report the
   full distribution of final RMSE across starts, not only the best.
3. Add a global stage before local polish — particle swarm or CMA-ES — with the
   local solver used only for refinement.
4. Add a **convergence assertion**: if `|RMSE_calibrated - RMSE_baseline| < 1e-6`
   the run is marked `OPTIMIZER_DID_NOT_MOVE` and is never eligible to be a
   publication or comparison candidate.
5. Record first-order optimality, exit flag, and evaluation count in the
   manifest.

**Acceptance tests**

- `tests/test_optimizer_did_not_move_detection.m`: a stubbed no-op solve is
  classified `OPTIMIZER_DID_NOT_MOVE`.
- `tests/test_multistart_determinism.m`: same seed, same starts, same result.
- Empirical: across 16 starts, the interquartile range of final primary RMSE is
  reported. A spread wider than 0.03 is itself a publishable finding about
  identifiability and must not be averaged away.

---

### Phase 4 — Identifiability and parameter reduction (requires UQLab)

**Goal:** stop fitting 14 parameters to 8 independent measurements.

Of the 16 finite rows, the independent measurement content is only **eight**
numbers: `RAP_mean`, `PAP_dia`, `PAP_sys`, `PAP_mean`, `SAP_dia`, `SAP_sys`,
`Qp`, `QpQs`. The rest are algebraic transforms —
`SAP_mean = 57 + (100-57)/3 = 71.33`, `CO_Lmin = Qp/QpQs = 3.42295`,
`Q_shunt = CO x (QpQs-1) = 0.66406`,
`SVR = (SAP_mean - RAP_mean)/CO = 19.369`,
`VSD_frac_pct = 100(1 - 1/QpQs) = 16.25` — or belong to the H+1 echo block whose
internal stroke-volume inconsistency the audit rates at 60%.

Fourteen free parameters against eight independent numbers is **negative degrees
of freedom**. Under that condition a low RMSE is a statement about model
flexibility, not about model validity.

**Changes**

1. Install UQLab and set `UNIFIED_VSD_UQLAB_PATH`. Run the contracted GSA at
   `N = 128` (and, budget permitting, `N = 512` for a convergence check on the
   Sobol indices themselves).
2. Rank parameters by total-order Sobol index aggregated over the primary
   metrics. Fix every parameter with `ST < 0.02` across all primary metrics at
   its scaled-baseline value and record the fixing decision with its index.
3. Re-run Phase 3 on the reduced set. Target an active set of **7 or fewer**
   parameters.
4. Compute and report a practical-identifiability diagnostic on the reduced set
   — profile likelihood or the Fisher information condition number — and publish
   parameter confidence intervals alongside every point estimate.

**Acceptance tests**

- `tests/test_gsa_registry_bounds.m` extended to assert the reduced active set
  is a strict subset of `recipe.active_parameters`.
- Reduced-set primary RMSE within 15% relative of the full-set value. A large
  degradation means the dropped parameters were doing real work and the
  reduction is wrong.

---

### Phase 5 — Resolve the clinical data contradiction

**Goal:** stop two files disagreeing about what the patient's data is.

[patient_reyna.m:83](config/patient_reyna.m:83) states the LV/RV volume block
"was confirmed to be H+1 after surgery, so it is not a valid pre-surgery
calibration target" and sets `LVEDV_mL = NaN`, `override_IC = false`.
[reyna_pre_surgery.m:38](config/calibration_recipes/reyna_pre_surgery.m:38) then
re-injects `LVEDV_mL = 41.0`, `LVESV_mL = 19.3`, `RVEDV_mL = 30.5`,
`RVESV_mL = 12.0`, `LVEF = 0.528`, and sets `override_IC = true`.

Four of those become `consistency_check_only`, but **`RVESV` stays `soft` and is
actively fitted**, and `override_IC = true` seeds the initial state from the same
block. Post-operative evidence is therefore shaping a pre-operative fit through
two channels.

**Changes**

1. Decide and document, in the recipe, one of:
   - **(a) Recommended.** Demote `RVESV` to `consistency_check_only` and set
     `override_IC = false` for fair-prior arms. The pre-surgery fit becomes
     purely haemodynamic, matching the stated intent of `patient_reyna.m`.
   - **(b)** Keep the block, and state in the paper that pre-operative
     calibration is partly informed by H+1 post-operative echo, with the 60%
     stroke-volume inconsistency quantified in the limitations.
2. Add a provenance field `evidence_timing` (`pre_operative` /
   `post_operative_H1` / `derived`) to every recipe override, and export it in
   the validation report so the timing mismatch is visible in every run.
3. Add a guard: no override whose `evidence_timing` is `post_operative_H1` may
   hold a tier of `hard` or `soft` in a `pre_surgery` scenario without an
   explicit `allow_cross_timing_evidence = true` flag in the recipe.

**Acceptance tests**

- `tests/test_evidence_timing_governance.m`: a `post_operative_H1` override at
  `soft` tier without the flag raises an error.
- The exported validation report shows `evidence_timing` for all 16 rows.

---

### Phase 6 — Sensitivity to the two contested clinical numbers

**Goal:** show the result is not an artefact of two defensible-but-arguable data
choices.

1. **MAP form factor.** `SAP_mean = 71.3` uses `dia + (sys-dia)/3`, a resting-
   adult approximation. At HR 119 in a 3-year-old the systolic time fraction is
   larger. Re-run with form factor `0.40` (`SAP_mean = 74.2`) as a sensitivity
   arm and report both.
2. **Cuff-versus-catheter conflict.** The same visit recorded NIBP `119/83`,
   `MAP 95` — about 33% above the catheter-derived `71.3` used as the fit
   target. Run an explicit arm at `SAP_mean = 95` and report how the accepted
   parameter set and the 10% gate count change. This is a limitation to
   quantify, not to bury.
3. Report all three systemic-pressure arms in one table. If the gate count is
   stable across them, that is a genuine robustness result and strengthens the
   paper considerably.

---

## 4. Execution contract

```bash
setx UNIFIED_VSD_UQLAB_PATH "<path-to-UQLab-core>"
```

```matlab
addpath(genpath(pwd));
result = run_reyna_scaling_experiment( ...
    'config/experiments/reyna_p1_scaling_v1.m', ...
    'Arms', {'fair_prior_zhang'}, ...
    'Repeats', 1, ...
    'DoGSA', true, ...
    'GsaN', 128, ...
    'ScreeningMode', false, ...
    'NumStarts', 16, ...
    'MaxFunctionEvaluations', 3000, ...
    'MaxIterations', 200);
```

Every phase must record: git commit SHA, scaling mode, seed policy, RNG seeds
for all starts, optimiser exit flag, evaluation count, target-tier policy
source, and the full-metric gate CSV.

## 5. Reporting requirement

The dated results report must publish the **full 16-row table**, not a 5-row
excerpt, and must state the count within 10% as `n of 16` alongside `n of 11`
for the governed mask. A claim of the form "all primary metrics pass" without
the accompanying full table is not acceptable output of this PRD.

## 6. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Phase 1 raises primary RMSE above the 0.09 acceptance ceiling | high | medium | Expected and correct; renegotiate the ceiling against the *honest* target set rather than reverting governance |
| Larger budget finds a lower-RMSE but less plausible basin | medium | high | Phase 2 plausibility terms retained; S6 caps warnings at 2 |
| UQLab remains unavailable | medium | high | Phases 0–3 and 5–6 do not need it; Phase 4 blocks and must be declared as an open limitation |
| Multi-start reveals wide RMSE spread | medium | medium | Report it — it is a finding about identifiability, not a failure |
| Gate hinge over-weights one metric and destabilises the fit | low | medium | `gateLambda` is tunable; Phase 2 regression compares decomposition tables |

## 7. Ordering

Phase 0 blocks all. Phase 1 blocks 2. Phase 2 blocks 3. Phase 4 requires 3 and
UQLab. Phases 5 and 6 may run in parallel with 2–4 and should be complete before
any publication claim.
