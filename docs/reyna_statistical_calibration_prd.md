# PRD: Statistically Governed Calibration (σ-weighting, χ², Joint Pre/Post Inversion)

Status: implementation-ready
Version: `reyna_statistical_calibration_v1`
Predecessor: [reyna_zhang_full_metric_prd.md](reyna_zhang_full_metric_prd.md) (Phases 0–3, 5 complete; PR #24)
Target executor: Claude Sonnet, working autonomously
Scope: Reyna `pre_surgery` + `post_surgery`

---

## 0. Read this first

This PRD is written to be executed without re-deriving context. Everything in
§1 is **verified fact** as of 2026-08-29 on branch
`codex/reyna-zhang-fullmetric-10pct`. Do not re-discover it; do verify anything
you are about to depend on, because code moves.

**The single most important instruction:** this project has a documented
history of results that looked good because the scorecard was wrong, not
because the model was right. Every phase below must be evaluated on whether it
makes the *evidence* better, never on whether it makes the *number* better.
§3 lists the specific ways that can go wrong here. Read it before writing code.

---

## 1. Verified current state

### 1.1 Where the pipeline stands

| Fact | Value |
|---|---|
| Branch | `codex/reyna-zhang-fullmetric-10pct` (PR #24) |
| UQLab | 2.2.0, installed, working (`toolbox/UQLab_Rel2.2.0/core`) |
| GSA | Runs at `N=128`; reduces active set 14 → **7** |
| Active parameters | `group.R_sys_scale`, `R.SVEN`, `group.R_pul_scale`, `C.SAR`, `C.PAR`, `E.RV.EB`, `vsd.Cd` |
| Governed metric set | **9 rows** (chamber volumes removed) |
| Best candidate | primary RMSE `0.0873`, governed gate **8/9** |
| Multi-start spread | min 0.0873, median 0.0915, max 0.1041 (6 Sobol starts) |
| Calibration status | `PROMISING_NEAR_MISS` |
| Chamber predictions | 6/6 inside paediatric bands, directions consistent |

### 1.2 The current candidate, per metric

σ below is from `get_calibration_targets.m`: `UncertaintyAbs` when finite,
otherwise `UncertaintyFraction × |ClinicalValue|`, where the fraction comes
from the reliability grade (High 5%, Moderate 10%, otherwise 20%).

| Metric | Clinical | Model | σ | \|z\| | 10% gate |
|---|---:|---:|---:|---:|---|
| `RAP_mean` | 5.000 | 4.964 | 0.250 | 0.15 | PASS |
| `PAP_min` | 10.000 | 9.005 | 0.500 | **1.99** | PASS |
| `PAP_max` | 20.000 | 21.728 | 1.000 | 1.73 | PASS |
| `PAP_mean` | 15.000 | 14.497 | 0.750 | 0.67 | PASS |
| `SAP_min` | 57.000 | 47.082 | 5.700 | 1.74 | **FAIL** |
| `SAP_max` | 100.000 | 90.668 | 10.000 | 0.93 | PASS |
| `SAP_mean` | 71.300 | 67.148 | 3.565 | 1.16 | PASS |
| `QpQs` | 1.194 | 1.166 | 0.060 | 0.47 | PASS |
| `CO_Lmin` | 3.423 | 3.133 | 0.500 | 0.58 | PASS |

`χ² = 13.22`, `N = 9`, `p = 7`, `χ²/N = 1.47`, `χ²/(N−p) = 6.61`.

### 1.3 The three findings that motivate this PRD

**F1 — the percentage gate ranks metrics wrongly.** `PAP_min` is the worst
residual in the set at 1.99σ and passes the 10% gate; `SAP_min` is better at
1.74σ and fails it. A percentage band treats a pressure recorded as
`10/10/10` across three repeats identically to one known to ±10%. It is not a
statistical criterion.

**F2 — the objective does not use the σ values that already exist.**
`get_calibration_targets` assigns every target an uncertainty. Only
`systemic_bundle_penalty` consumes it. The main loop in
`objective_calibration.m` normalises by `calib.primaryTarget` (a single global
percentage) times a tier weight, so declared measurement precision has no
effect on the fit.

**F3 — degrees of freedom are 2.** 9 observations, 7 parameters. This is why
`χ²/(N−p)` reads 6.61 despite an average residual of only 1.21σ. No amount of
optimiser effort fixes this; it needs more independent data.

---

## 2. Objectives

| # | Objective | Measured by |
|---|---|---|
| O1 | Residuals weighted by declared measurement uncertainty | objective consumes σ for every fitted target |
| O2 | A defensible stopping criterion that cannot be gamed by redefining the denominator | `χ²/N` reported every run |
| O3 | The retained parameter subset is jointly identifiable, not merely individually sensitive | correlation/condition-number report on the active set |
| O4 | More independent observations without collecting new data | joint pre/post inversion; DOF 2 → ≥6 |
| O5 | At least one target never used in fitting | holdout metric reported as a prediction |

### 2.1 Non-goals

- Reaching 9/9 on the percentage gate. If σ-weighting is correct, the
  percentage gate becomes the *secondary* criterion and its count may go down
  while the fit improves. **That is an acceptable and expected outcome.**
- Removing metrics from the reported denominator to raise a pass count.
- Re-tuning `gateLambda`, tier weights, or bounds to recover a number.
- Changing any clinical measurement value.

---

## 3. Guardrails

These are hard constraints. Violating one invalidates the phase.

**G1 — Never change a target value to improve a fit.** The one legitimate
exception already ran and was reported as refuted (MAP form factor, PRD-1
Phase 6). Clinical numbers in `config/patient_reyna.m` are read-only.

**G2 — Never shrink the reported denominator in the same change that reports a
pass count.** If a metric should leave the governed set on principle, land that
in its own commit, with the before/after count on both denominators.

**G3 — Report every arm you run, including the ones that got worse.** The
refuted MAP hypothesis is in the results doc precisely because it failed.

**G4 — Do not invent literature values.** Reference ranges, σ values, and
physiological constants must come from the repo, from the patient record, or
from a citation the user supplies. If a number is needed and unavailable,
**stop and ask** rather than choosing a plausible one.

**G5 — σ values are governance, not tuning knobs.** You may change how σ is
*used*. You may not adjust a σ to make a metric fit. If a σ looks wrong,
report it as a finding.

**G6 — Every new penalty term is opt-out-able and reported separately.** Add it
behind a named lambda in `calibration_param_sets.m`, and surface its
contribution in the objective breakdown.

**G7 — A phase that makes the fit worse is still a completed phase** if the
diagnosis is sound. Write it up and move on; do not iterate until the number
improves.

---

## 4. Phase 1 — σ-weighted (χ²) residuals

**Goal (O1).** Every fitted residual expressed in units of that metric's own
measurement uncertainty, so "close" means the same thing for a pressure known
to ±5% and a flow known to ±15%.

### 4.1 Formula

Replace the current data term

```
J_m = w_m * ( |y_model - y_obs| / max(|y_obs|,1e-6) / calib.primaryTarget )^2
```

with

```
J_m = w_m * ( (y_model - y_obs) / sigma_m )^2
```

where `sigma_m` resolves as:

1. `targets(i).UncertaintyAbs` if finite and > 0;
2. otherwise `targets(i).UncertaintyFraction * |targets(i).ClinicalValue|`;
3. otherwise `0.10 * |ClinicalValue|`, and **emit a warning naming the metric**
   (a target with no declared uncertainty is a governance gap, not a default).

Floor `sigma_m` at `1e-9` to avoid division blow-up.

`w_m` remains the tier weight, default 1.0. Under σ-weighting the tier weights
should be *closer to uniform* than they are now, because precision is handled
by σ. **Do not retune them in this phase** — that is Phase 1b, gated on the A/B
result.

### 4.2 Files

| File | Change |
|---|---|
| `src/calibration/build_case_calibration_profile.m` | Build `profile.targetSigma` (struct: metric → σ) from `get_calibration_targets(scenario, clinical)`. Mirror the existing `profile.referenceRanges` pattern added in PR #24. |
| `src/calibration/calibration_param_sets.m` | Add `calib.objectiveWeighting` (`'sigma'` \| `'legacy'`), read from `UNIFIED_VSD_OBJECTIVE_WEIGHTING`, **default `'legacy'`** until §4.5 passes. Copy `calib.targetSigma` from the case profile. |
| `src/calibration/objective_calibration.m` | Branch the data term on `calib.objectiveWeighting`. Keep the legacy path byte-for-byte reachable. |

`systemic_bundle_penalty` already works in σ units — **leave it alone**. Verify
you have not double-counted: metrics routed to the bundle must still `continue`
before the data term, exactly as now.

### 4.3 Why this should move `SAP_min`

Under legacy weighting the objective sees `SAP_min` and `CO_Lmin` as similar
relative errors. Under σ-weighting `SAP_min` is 1.74σ and `CO_Lmin` is 0.58σ,
so `SAP_min` contributes ~9× more to `J`. The expected trade is `SAP_min`
improves and `CO_Lmin` degrades toward its (deliberately wide) tolerance.

**That is the intended behaviour, not a regression.** Reyna's pressures are
known far better than her Fick-derived flow.

### 4.4 Tests — `tests/test_sigma_weighted_objective.m`

1. σ resolution prefers `UncertaintyAbs` over `UncertaintyFraction`
   (`CO_Lmin` must resolve to 0.50, not 0.15 × 3.423).
2. σ resolution falls back to the fraction when abs is NaN
   (`SAP_min` → 5.70).
3. A target with neither declared warns and defaults to 10%.
4. σ is floored: a zero/negative declared σ does not produce Inf or NaN in `J`.
5. `objectiveWeighting = 'legacy'` reproduces the pre-change objective exactly
   (`AbsTol 1e-12`) for a fixed parameter vector. **This is the regression
   guard — write it first and confirm it fails if you break the legacy path.**
6. Under `'sigma'`, a metric at 1σ contributes `w_m` exactly.
7. Under `'sigma'`, two metrics at equal z contribute equally regardless of
   their absolute magnitudes (use `PAP_min` at 10 mmHg vs `SAP_max` at 100).
8. Systemic-bundle metrics are not double-counted under either mode.

### 4.5 A/B acceptance

Run both arms at identical settings (Zhang, fair prior, GSA on `N=128`,
`MAX_FUN_EVALS=300`, `NUM_STARTS=6`, seed 20260828):

| | legacy | sigma |
|---|---|---|
| `χ²/N` | 1.47 (known) | record |
| Governed gate | 8/9 (known) | record |
| Primary RMSE | 0.0873 (known) | record |
| Per-metric \|z\| | §1.2 | record |

**Promote σ-weighting to default only if `χ²/N` moves toward 1 and no metric
exceeds 3σ.** A drop in the percentage-gate count is acceptable (§2.1) provided
`χ²/N` improved — say so explicitly in the write-up rather than burying it.

If `χ²/N` moves *away* from 1, keep `'legacy'` as default, keep the code, and
report the negative result. That is a complete phase under G7.

---

## 5. Phase 2 — χ² reporting and the discrepancy criterion

**Goal (O2).** Every run reports a statistic that cannot be improved by
redefining the denominator.

### 5.1 Changes

**`src/utils/export_full_metric_gate.m`** — add columns `Sigma`, `ZScore`,
`ZScoreSquared` to the exported table. Sort remains worst-first, but sort by
`|z|`, not by percentage, and add `AbsError_pct` as a secondary display column
so both views survive.

**New `src/utils/compute_chi_squared_report.m`:**

```
report = compute_chi_squared_report(gate_table, n_active_parameters)
  .chi2              sum of z^2 over governed rows
  .n_obs             governed row count
  .n_parameters      active parameter count
  .dof               max(n_obs - n_parameters, 0)
  .chi2_per_obs      chi2 / n_obs
  .chi2_reduced      chi2 / dof, NaN when dof <= 0
  .interpretation    'underfit' | 'consistent' | 'overfit' | 'insufficient_dof'
  .worst_metric, .worst_z
```

Interpretation bands, applied to `chi2_per_obs`:

| Band | Label | Meaning |
|---|---|---|
| > 2.0 | `underfit` | model or data inconsistent |
| 0.5 – 2.0 | `consistent` | residuals match measurement noise |
| < 0.5 | `overfit` | fitting below the noise floor, or σ too generous |

When `dof <= 2`, set `.interpretation` regardless but **also** emit
`insufficient_dof` in the console line — with `p = 7` and `N = 9` the reduced
χ² is not stable and must not be quoted alone.

**`src/calibration/classify_calibration_run.m`** — record `chi2_per_obs` and
its label in `status`. **Do not gate `ACCEPT` on it in this phase.** Observe it
across several runs first; changing the acceptance rule and the objective in
the same change makes neither attributable.

**`src/utils/validation_report.m`** — print a `--- CHI-SQUARED ---` block after
the gate block: `χ²`, `N`, `p`, `χ²/N`, `χ²/(N−p)`, interpretation, worst
metric by `|z|`.

### 5.2 Tests — `tests/test_chi_squared_report.m`

1. Synthetic table where every metric is exactly 1σ → `chi2 == n_obs`,
   `chi2_per_obs == 1`, label `consistent`.
2. Every metric at 3σ → `chi2_per_obs == 9`, label `underfit`.
3. Every metric at 0.1σ → label `overfit`.
4. `dof <= 0` → `chi2_reduced` is NaN, no error thrown.
5. Only governed rows count: a row with `InPrimaryRMSE = false` is excluded
   from `chi2` and `n_obs`.
6. Rows with non-finite σ or model value are skipped and counted in a
   `.n_skipped` field.
7. Reproduces §1.2 exactly: feeding that table with `p = 7` gives
   `chi2 = 13.22 ± 0.01`, `chi2_per_obs = 1.47 ± 0.01`.

---

## 6. Phase 3 — joint identifiability of the retained subset

**Goal (O3).** Sobol ranks parameters one at a time. It cannot detect that two
retained parameters are collinear — e.g. `C.SAR` and `group.R_sys_scale` trading
off along the systemic RC time constant, which §3.5 of the results doc
identified as governing the waveform form factor.

### 6.1 Method

At the calibrated operating point, build the scaled sensitivity matrix

```
S(i,j) = (partial y_i / partial theta_j) * (theta_j / sigma_i)
```

by central finite differences with a 1% relative step, over governed metrics
`i` and active parameters `j`. Then report:

- `cond(S)` — condition number. `> 1e3` indicates near-dependence.
- Column correlation matrix. Any `|rho| > 0.9` pair is a collinear pair.
- Per-parameter column norm — a near-zero column is an inactive parameter that
  GSA should have dropped.

### 6.2 Files

New `src/calibration/analyse_parameter_identifiability.m`, called from
`main_run.m` after calibration, exporting
`<run_dir>/tables/parameter_identifiability_<scenario>.csv`.

**Report only. Do not automatically drop parameters in this phase.** A
collinear pair is a finding that needs physiological judgement about which
member to fix, and that judgement is the user's.

### 6.3 Tests — `tests/test_parameter_identifiability.m`

1. Analytic case: two exactly-duplicated columns → `|rho| = 1`, flagged.
2. Orthogonal columns → correlations ≈ 0, `cond` near 1.
3. A zero column is flagged as inactive without dividing by zero.
4. Output CSV has one row per active parameter plus a pair table.

---

## 7. Phase 4 — joint pre/post inversion

> **Status update, 2026-08-29 — §7.2's blocking question is ANSWERED and this
> phase's premise has improved. Read this before implementing §7.**
>
> The §7.2 questions were put to the study owner and answered, and the source
> catheterisation record (RSAB Harapan Kita, MRN 01008971, 06/04/2026) was
> then retrieved. Two things changed:
>
> 1. **Real post-closure haemodynamics exist**, from the *same* session as the
>    pre-closure readings (device released 12.06.23; readings 12.15–12.32 under
>    unchanged anaesthesia). `clinical.post_surgery` now holds 7 finite
>    pressure targets — `PAP` 17/9 (13), `SAP` 89/68 (79), `RAP` mean 5. This
>    is a genuine paired dataset, which is a **stronger** basis for §7.1's
>    shared-parameter assumption than the volume relocation this section was
>    written around. `N = 16` against `p = 12` gives `dof = 4`, not 2.
> 2. **The "60% internal stroke-volume inconsistency" in §7.2 item 3 is
>    withdrawn.** It was never internal to the H+1 block (`SV_LV` 21.7 vs
>    `SV_RV` 18.5, 14.7% apart, `LVEF` self-consistent). The 60% came from
>    comparing post-closure volumes against pre-closure flows — the shunt
>    itself. With the block excluded the audit returns severity `none`.
>
> Two amendments to §7 follow. **Do not relocate the H+1 echo volumes into
> `post_surgery`**: they are a ward echo at a different timepoint from the
> in-lab catheter pressures now stored there, and merging them would recreate
> the timing mismatch this project has repeatedly been bitten by. And note
> there is still **no post-closure flow** (see the open `PARI`/`FR` question),
> so the post state constrains pressures only.
>
> Full detail: `docs/reyna_statistical_calibration_results_20260829.md` §0, §6.

**Goal (O4).** The largest available win. Raises observation count without new
data collection and resolves the volume problem by putting the volumes where
their timing is valid.

### 7.1 The idea

One patient means one set of structural parameters. Between pre-op and post-op
essentially only the defect closes. So fit **both simulations simultaneously**
from a single shared parameter vector:

| | Pre-op simulation | Post-op simulation |
|---|---|---|
| Shared | chamber elastances, `V0`, vascular R and C, `R_sys`, `R_pul` | same values |
| Differs | `R.vsd` open (from geometry), `vsd.Cd` active | `R.vsd` closed |
| Constrained by | catheter pressures, `QpQs`, `CO` | H+1 echo volumes, `LVEF` |

Objective: `χ²_pre + χ²_post`, both σ-weighted per Phase 1.

| | Current | Joint |
|---|---:|---:|
| Observations | 9 | ~14 |
| Parameters | 7 (+1 for the closed shunt) | ~8 |
| **DOF** | **2** | **~6** |

### 7.2 BLOCKING PREREQUISITE — data governance

The H+1 chamber block must move from `pre_surgery` to `post_surgery` in
`config/patient_reyna.m`. Values are preserved in
`recipe.excluded_evidence` (`config/calibration_recipes/reyna_pre_surgery.m`):
`LVEDV 41.0`, `LVESV 19.3`, `RVEDV 30.5`, `RVESV 12.0`, `LVEF 0.528`.

**Do not make this move autonomously.** It relocates patient measurements and
requires the user to confirm:

1. The H+1 timing against the source record.
2. That H+1 is an acceptable proxy for the converged post-closure state, or
   that the limitation is stated.
3. How the block's **60% internal stroke-volume inconsistency** (severity
   `critical` from `audit_clinical_consistency`) is handled — it travels with
   the data and is not resolved by moving it.

**Stop and ask. Phase 4 does not start until this is answered.**

### 7.3 Implementation once unblocked

**New `src/calibration/objective_joint_pre_post.m`:**

```
J = objective_joint_pre_post(x, params0_pre, params0_post, clinical, calib, pce)
```

- Applies `x` to both parameter structs via `set_calibration_param_value`.
- Overrides `R.vsd` per scenario: pre from geometry, post from the closed value
  already used by `configure_vsd`.
- Integrates both, computes both metric sets.
- Returns `χ²_pre + χ²_post` plus the existing plausibility/boundary terms,
  applied once to the shared vector (**not** twice — that would double the
  regularisation).

**New `scripts/run_joint_pre_post_calibration.m`** — entry point, modelled on
the existing `scripts/run_pre_post_protocol.m` and
`config/run_protocols/reyna_pre_post_publishable.m`, both already present.

**Target governance:** the post-op chamber rows become genuine `soft` targets
in `post_surgery`. `assert_evidence_timing_governance` must then record their
timing as matching the scenario, so no violation is raised. Verify the existing
`RVEF` anti-double-counting rule fires when `RVEDV`/`RVESV`/`RVEF` are all
present post-op.

### 7.4 Tests — `tests/test_joint_pre_post_objective.m`

1. Shared parameters are identical in both structs after application; only
   `R.vsd` differs.
2. Post-op `R.vsd` is the closed value; pre-op is not.
3. `J` equals `χ²_pre + χ²_post` plus regularisation counted **once**
   (`AbsTol 1e-12`).
4. A parameter set matching pre-op perfectly but post-op badly scores worse
   than one splitting the difference — the joint fit must actually couple.
5. Missing post-op targets degrade gracefully: joint reduces to pre-only, with
   a warning, not an error.
6. `assert_evidence_timing_governance` raises no violation for post-op chamber
   rows in the `post_surgery` scenario.

### 7.5 Acceptance

- DOF ≥ 6.
- `χ²/N` over the combined set within the `consistent` band.
- Pre-op governed gate no worse than 7/9 (joint fitting trades some pre-op
  accuracy for identifiability; a small loss is expected and acceptable, a
  large one means the shared-parameter assumption is wrong and should be
  reported as such).
- Post-op scenario has ≥ 4 clinical targets where it currently has zero.

---

## 8. Phase 5 — validation holdout

**Goal (O5).** Nothing is currently held out; every finite target is either
fitted or reported. The AJP-Heart guidelines call this out directly.

Designate `SVR` as `validation_holdout` in the Reyna recipe. It is already
`derived_validation`, is excluded from the governed RMSE, and currently lands
at 2.46% — a genuine prediction test rather than a fitted result.

Report it in a `--- VALIDATION HOLDOUT ---` block, never in `χ²` or the gate
count. Add to `tests/test_no_ungoverned_calibration_targets.m`: `SVR` must be
absent from `included_in_calibration` **and** from `included_in_primary_rmse`.

---

## 9. Execution environment

Learned the hard way during PR #24. Follow exactly.

**MATLAB invocation.** `-batch` does not start in the shell's working
directory. Always:

```bash
matlab -batch "cd('D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd'); addpath(genpath(pwd)); <code>"
```

**Never run two MATLAB processes at once.** A calibration run and a test run
compete; a run that normally takes 20 minutes stretched past 75 and had to be
killed. Wait for completion.

**Never `git stash` / `checkout` while a run is in flight.** MATLAB lazily
loads `.m` files, so a mid-run working-tree change silently mixes code
versions.

**Nested-function scoping.** `objective_calibration.m` uses nested functions
sharing the parent scope. A loop variable named `idx` in a new local function
collides with the parent's. Use distinctive names (`range_idx`, `target_ix`).

**Seed the RNG before `scramble()`.** `sobolset` scrambling draws from the
global stream; without seeding, "deterministic" multi-start is not reproducible.
`build_multistart_starts.m` does this correctly — copy the pattern (save state,
seed, restore via `onCleanup`).

**Timing.** GSA at `N=128` ≈ 8 min. Single calibration at 300 evals ≈ 15 min.
Six starts ≈ 100 min. Budget accordingly and background long runs.

**Reference commands.**

```bash
# full new-test sweep
matlab -batch "cd('<repo>'); addpath(genpath(pwd)); r=runtests({'tests/test_full_metric_gate_export.m','tests/test_no_ungoverned_calibration_targets.m','tests/test_gate_hinge_penalty.m','tests/test_multistart_starts.m','tests/test_evidence_timing_governance.m','tests/test_governed_gate_acceptance.m'}); fprintf('P=%d F=%d\n',nnz([r.Passed]),nnz([r.Failed])); exit(0)"
```

```bash
# reference calibration arm
matlab -batch "cd('<repo>'); addpath(genpath(pwd)); setenv('UNIFIED_VSD_UQLAB_PATH',fullfile(pwd,'toolbox','UQLab_Rel2.2.0','core')); setenv('UNIFIED_VSD_SCALING_MODE','zhang'); setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','1'); setenv('UNIFIED_VSD_DO_GSA','1'); setenv('UNIFIED_VSD_GSA_PCE_N','128'); setenv('UNIFIED_VSD_DO_PLOTS','0'); setenv('UNIFIED_VSD_MAX_FUN_EVALS','300'); setenv('UNIFIED_VSD_MAX_ITERATIONS','40'); setenv('UNIFIED_VSD_NUM_STARTS','6'); main_run('pre_surgery', patient_reyna()); exit(0)"
```

**Known pre-existing failure**, present on `main`, unrelated to this work:
`tests/test_clinical_consistency_target_tiers.m` (a PVR tier assertion). Do not
"fix" it as part of these phases.

---

## 10. Ordering and dependencies

```
Phase 1 (sigma weighting) ──┬──> Phase 2 (chi-squared reporting)
                            │
                            └──> Phase 4 (joint inversion)  [BLOCKED on §7.2]

Phase 3 (identifiability) ── independent, may run any time after Phase 1

Phase 5 (holdout) ── independent, run last so it does not perturb the A/B
```

Phase 2 depends on Phase 1 only because χ² is most meaningful once σ is what
the objective actually minimises. It can be implemented first as
reporting-only if that is more convenient — say which you did.

---

## 11. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| σ-weighting lowers the percentage-gate count | **high** | low | Expected; §2.1 declares it acceptable. Report both criteria. |
| σ-weighting degrades `CO_Lmin` past 10% | medium | medium | Intended direction, but if `CO` exceeds 3σ the widened uncertainty is wrong and that is a finding to report, not a thing to tune away. |
| Reduced χ² unstable at dof = 2 | **certain** | medium | Always quote `χ²/N` alongside; never quote `χ²/(N−p)` alone until Phase 4 raises dof. |
| Joint inversion blocked on data governance | high | high | §7.2 is a hard stop. Phases 1, 2, 3, 5 all deliver without it. |
| Shared-parameter assumption is wrong (remodelling between pre and post) | medium | high | Phase 4 acceptance allows a small pre-op degradation; a large one is evidence against the assumption and must be reported. |
| Legacy objective path silently broken | medium | **high** | Test 4.4.5 is the regression guard. Write it first. |

---

## 12. Deliverables

- Code and tests on a branch off `codex/reyna-zhang-fullmetric-10pct`.
- `docs/reyna_statistical_calibration_results_<date>.md` reporting **every**
  arm run, including regressions, with per-metric `|z|` tables.
- Updated `docs/supervisor_comments_response_20260828.md` — Phase 2 gives a
  sharper answer to comment [221] ("apa itu primary governed? full
  transparent?"), because `χ²/N` is the honest single number that discussion
  was reaching for.
- A PR whose body states, for each phase: what changed, what the number did,
  and whether the phase is being promoted or retained as a negative result.

## 13. Definition of done

- [ ] Phase 1 A/B run and reported, default set by the §4.5 rule
- [ ] `χ²/N` printed and exported on every run
- [ ] Identifiability report exported; collinear pairs named or absence stated
- [ ] Phase 4 either complete, or blocked with the §7.2 question put to the user
- [ ] `SVR` held out and reported as a prediction
- [ ] All new tests passing; the pre-existing failure in §9 still the only one
- [ ] Every arm in the results doc, including any that got worse
