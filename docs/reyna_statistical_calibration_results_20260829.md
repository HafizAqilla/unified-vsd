# Reyna Statistical Calibration — Results, 2026-08-29

Branch: `codex/reyna-statistical-calibration` (off `codex/reyna-zhang-fullmetric-10pct`, PR #24)
PRD: [reyna_statistical_calibration_prd.md](reyna_statistical_calibration_prd.md)
Scope: Reyna `pre_surgery`, Zhang scaling prior, fair-prior (historical seeds disabled)

---

## 1. What this branch changes

Four of the PRD's five phases are complete. Phase 4 (joint pre/post inversion)
is deliberately blocked — see §6.

| Phase | Status | Summary |
|---|---|---|
| 1 — σ-weighted objective | **Complete** | Residuals normalised by declared measurement uncertainty instead of one global percentage; opt-in via `objectiveWeighting`, default remains `legacy` |
| 2 — χ² reporting | **Complete** | Discrepancy-principle goodness-of-fit statistic, printed and exported every run; does not gate `ACCEPT` |
| 3 — parameter identifiability | **Complete** | Scaled sensitivity matrix, condition number, pairwise correlation; report-only |
| 4 — joint pre/post inversion | **Blocked** | Requires a clinical data-governance decision — see §6 |
| 5 — validation holdout | **Deviated, with reasoning** | `SVR` was not relabelled `validation_holdout` — see §5 |

## 2. The measurement that motivated this branch

A percentage acceptance gate treats every metric as equally well known. It
isn't a statistical criterion. Checking the PR #24 final candidate:

| Metric | Absolute error | Declared σ | \|z\| (error in σ units) | 10% gate |
|---|---:|---:|---:|---|
| `PAP_min` | 1.0 mmHg | ±0.5 (repeated 10/10/10) | **1.99** ← worst residual | PASS |
| `SAP_min` | 9.9 mmHg | ±5.7 | 1.74 | **FAIL** |

`PAP_min` is the worst-fitting metric in the whole set by any statistical
measure, and the percentage gate passes it while failing a metric that fits
*better* in σ units. This is not a fitting failure — it is the gate measuring
the wrong thing.

By the proper statistical measure, the PR #24 candidate is already close to
appropriately fit: `χ² = 13.22`, `N = 9`, `p = 7`, `χ²/N = 1.47` — inside the
`consistent` band (0.5–2.0), meaning the residuals are broadly consistent with
declared measurement noise, not badly wrong.

## 3. Phase 1 — σ-weighted objective

### 3.1 What changed

`objective_calibration.m` gained a second data-term mode. Residuals are now,
optionally, normalised by each metric's own declared σ (from
`get_calibration_targets`: `UncertaintyAbs` when finite, else
`UncertaintyFraction × |ClinicalValue|`) instead of one global
`calib.primaryTarget`/`calib.secondaryTarget` percentage.

- `calib.objectiveWeighting`: `'legacy'` (default) or `'sigma'`, settable via
  `UNIFIED_VSD_OBJECTIVE_WEIGHTING`.
- The `'legacy'` path is byte-identical to the pre-change objective — this is
  asserted by a dedicated regression test
  (`test_sigma_weighted_objective.m::test_legacy_mode_is_byte_identical_to_pre_change`),
  which must never be allowed to fail.
- The gate hinge (10% acceptance band) is unchanged and still applied in
  addition to whichever data-term mode is active.

### 3.2 A/B comparison — first attempt invalidated, redo pending

Both arms: fair-prior Zhang, GSA on (`N=128`, active set 14→7), 300 evals,
6 Sobol starts, seed 20260828 — identical to the PR #24 final candidate.

**A first attempt at the sigma-side run was invalidated by a process error and
is not reported here.** The sigma-weighted calibration arm was launched
before Phases 2, 3, and 5 were implemented, and this branch's own files
(`main_run.m`, `validation_report.m`, `classify_calibration_run.m`,
`export_full_metric_gate.m`) were then edited *while that run was still
executing in the background* — a violation of this same PRD's §9 guidance
("never edit files a running calibration depends on"). The consequence:

- `main_run.m` is the top-level function and stays on the call stack for the
  full run duration, so its own control flow (which arguments it passes to
  `validation_report`) was frozen at the moment the run started — before the
  `NumActiveParameters` argument and the Phase 3 identifiability call existed.
- `validation_report.m` and the functions it calls are invoked fresh, so they
  picked up the finished Phase 2 code when Step 8 finally ran, near the end of
  the ~3.4-hour run.
- The result was a hybrid: Phase 2's chi-squared code ran, but received no
  active-parameter count (defaulting to `p=0`), and Phase 3's identifiability
  analysis never ran at all, because that code did not exist yet in the
  frozen `main_run.m` body.
- A genuine bug in the chi-squared CSV export (`compute_chi_squared_report.m`
  defaulted `dof_note` to `''`, an empty char; `struct2table` requires every
  field of a scalar struct to represent one row, and an empty-char field
  breaks that) then crashed the entire process at the very last reporting
  step — **after** the full 6-start optimisation had completed successfully,
  but **before** any artefact could be saved to disk.

The sigma-weighted RMSE values that were visible in the console before the
crash (`min 0.0780, median 0.1040, max 0.1144, IQR 0.0143` across 6 starts)
are plausible — Phase 1's objective code was complete and frozen before this
run started, so that part is not contaminated — but were not saved to any
file and are reported here only as an unverified data point, not as the
branch's result.

**Both bugs are fixed:**

- `compute_chi_squared_report.m`: `dof_note` and `worst_metric` now default
  to non-empty placeholders (`'sufficient_dof'`, `'none'`) instead of `''`,
  and a direct regression test (`test_report_is_always_struct2table_compatible`)
  asserts `struct2table` never throws on any code path.
- `validation_report.m`: the chi-squared CSV export is now wrapped in
  `try`/`catch`, matching the existing defensive pattern already used for the
  chamber-state and identifiability exports — a reporting bug can no longer
  discard a completed calibration.

**Re-run deferred, not scheduled.** All Phases 1–5 code is complete, tested
(85/85), and committed on this branch. The clean 6-start A/B run is deferred
to a later session at the user's request — see §10 for the exact command and
the rule that must hold for the duration of that run: nothing else touches
this branch's calibration-path files while it executes.

| | legacy (PR #24, cited) | sigma (pending re-run) |
|---|---:|---:|
| Governed gate | 8 / 9 | <!-- SIGMA_GATE --> |
| Primary RMSE (best of 6 starts) | 0.0873 | <!-- SIGMA_RMSE --> |
| RMSE spread (min–max, 6 starts) | 0.0873 – 0.1041 | <!-- SIGMA_SPREAD --> |
| χ²/N | 1.47 | <!-- SIGMA_CHI2 --> |
| Interpretation | consistent | <!-- SIGMA_INTERP --> |

<!-- SIGMA_PER_METRIC_TABLE -->

### 3.3 Decision

<!-- SIGMA_DECISION -->

### 3.4 Process lesson

This project's own PRD (§9) already warned against exactly this failure mode,
and it happened anyway because the warning was written for future executors
and not re-checked against live behaviour in the moment. The concrete rule
going forward: **once a calibration run is launched, the files it depends on
are read-only until it completes** — new phases queue behind the run, they do
not land mid-flight. This is now also the operating rule for the scheduled
re-run in §10: nothing else touches this branch's calibration-path files
while it is executing.

## 4. Phase 2 — χ² reporting

Every run now exports:

- `full_metric_gate_<scenario>.csv` with three new columns: `Sigma`, `ZScore`,
  `ZScoreSquared`. The table is now sorted worst-first by `|z|` rather than by
  percentage error, with `AbsError_pct` retained as a secondary column.
- `chi_squared_<scenario>.csv`: a one-row summary (`chi2`, `N`, `p`, `dof`,
  `chi2/N`, `chi2/dof`, `interpretation`, `dof_note`, `worst_metric`,
  `worst_z`).
- A `--- CHI-SQUARED ---` console block after the metric gate block.

`chi2_per_obs` and its interpretation are recorded on
`classification_status` (visible in the console summary line as
`chi2_per_obs=X.XX (label)`) but **do not** influence the `ACCEPT`/`REJECT`
label in this phase — changing the objective and the acceptance rule in the
same change would make neither attributable. This is locked in by
`test_governed_gate_acceptance.m::test_chi2_is_recorded_but_does_not_gate_accept`.

## 5. Phase 5 — validation holdout: deviation from the written PRD

The PRD's Phase 5 instruction was to designate `SVR` as `validation_holdout`,
on the stated premise that it is "a genuine prediction test rather than a
fitted result."

**That premise does not survive checking the derivation.**
`build_systemic_bundle` in `objective_calibration.m` computes:

```
SVR_target_WU = (SAP_mean - RAP_mean) / CO_Lmin
```

`SAP_mean`, `RAP_mean`, and `CO_Lmin` are all **hard-tier targets already in
the objective**. SVR is not an independent measurement; it is algebra over
three quantities the model is already being fitted to. Predicting it well
(currently 2.46% error) demonstrates that `SAP_mean`, `RAP_mean`, and `CO_Lmin`
are mutually consistent — it does not demonstrate anything the model was not
already told. This is the identical reasoning that already correctly excludes
`Q_shunt_Lmin` from the governed primary RMSE (`primary_rmse_holdout`, PR #24).

Relabelling SVR `validation_holdout` would have created the appearance of an
independent prediction success where none exists — precisely the kind of
overclaiming this whole campaign exists to eliminate.

**What was implemented instead:**

- The holdout reporting machinery (`print_validation_holdout.m`, a
  `--- VALIDATION HOLDOUT ---` console block) was built generically, keyed off
  the `validation_holdout` tier, so it is ready to use once a genuinely
  independent holdout candidate is identified.
- `SVR` stays `derived_validation`. A new test
  (`test_no_ungoverned_calibration_targets.m::test_svr_is_excluded_from_calibration_and_primary_rmse`)
  locks in the exclusion outcome (excluded from fitting and from the governed
  RMSE, which both tiers share) while documenting why the *label* is wrong.
- A second test
  (`test_no_current_metric_is_designated_a_genuine_validation_holdout`) locks
  in the honest current state: nothing in the Reyna pre-surgery recipe is
  presently a genuine independent holdout, with an explicit note on what would
  need to be true (independent measurement, deliberately excluded from
  fitting for that reason) before this assertion should be updated.

Of the 9 governed metrics, at least two others are also partially derived
rather than raw independent readings — `SAP_mean` (from `SAP_sys`/`SAP_dia` via
a form-factor formula) and `CO_Lmin` (`Qp / QpQs`) — documented in
`docs/reyna_zhang_scientific_assessment_20260828.md` §2. A genuinely
independent holdout candidate would need to come from the raw catheter
readings (`RAP_mean`, `PAP_min`, `PAP_max`, `SAP_min`, `SAP_max`, `QpQs`)
deliberately excluded from the objective and refit — a materially larger
change than a recipe relabel, closer in scope to another multi-start
calibration arm. Not attempted in this branch; recorded as a candidate for a
future PRD.

## 6. Phase 4 — blocked

Per PRD §7.2, joint pre/post inversion requires relocating the H+1
post-operative chamber-volume block from `clinical.pre_surgery` to
`clinical.post_surgery` in `config/patient_reyna.m`. This is a clinical
data-governance decision, not a code change, and the PRD explicitly instructs:
**stop and ask rather than make it autonomously.**

Three questions remain open for the user:

1. Confirm the H+1 timing against the source clinical record.
2. Confirm H+1 is an acceptable proxy for the converged post-closure state, or
   state the limitation explicitly if it is not.
3. Decide how to handle the block's **60% internal stroke-volume
   inconsistency** (severity `critical` per `audit_clinical_consistency`) once
   it becomes a `post_surgery` target — moving it relocates the
   inconsistency, it does not resolve it.

Phases 1, 2, 3, and 5 do not depend on this and are complete without it.

## 7. Parameter identifiability (Phase 3) — preliminary finding

A 2-parameter, 2-metric smoke check during test development (`R.SVEN` vs
`C.SAR`, evaluated at `RAP_mean`/`SAP_mean`) showed **ρ = −1.000** — perfect
collinearity, though on a toy subset and at the demographic-scaled baseline
rather than the calibrated operating point. This is consistent with the
systemic RC time-constant coupling identified in
`docs/reyna_zhang_fullmetric_results_20260828.md` §3.5 (the waveform
form-factor mismatch driving `SAP_min`'s residual).

<!-- IDENTIFIABILITY_FULL_REPORT -->

## 8. Reproduction

```bash
matlab -batch "cd('<repo>'); addpath(genpath(pwd)); setenv('UNIFIED_VSD_UQLAB_PATH', fullfile(pwd,'toolbox','UQLab_Rel2.2.0','core')); setenv('UNIFIED_VSD_SCALING_MODE','zhang'); setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','1'); setenv('UNIFIED_VSD_DO_GSA','1'); setenv('UNIFIED_VSD_GSA_PCE_N','128'); setenv('UNIFIED_VSD_MAX_FUN_EVALS','300'); setenv('UNIFIED_VSD_MAX_ITERATIONS','40'); setenv('UNIFIED_VSD_NUM_STARTS','6'); setenv('UNIFIED_VSD_MULTISTART_SEED','20260828'); setenv('UNIFIED_VSD_OBJECTIVE_WEIGHTING','sigma'); main_run('pre_surgery', patient_reyna())"
```

Every run now writes, into its run folder's `tables/`:
- `full_metric_gate_<scenario>.csv` (with `Sigma`, `ZScore`, `ZScoreSquared`)
- `chi_squared_<scenario>.csv`
- `parameter_identifiability_<scenario>.csv` and
  `parameter_identifiability_pairs_<scenario>.csv`

## 9. What this does not establish

- **n = 1.** Still one patient.
- **Phase 4 is unresolved.** DOF stays at 2 until it lands.
- **The identifiability report in §7 is a smoke test, not a governed-set
  analysis.** The full 9-metric × 7-parameter report from this branch's actual
  calibrated candidate is in §7's placeholder above once available.
- **Zhang vs Lundquist remains uninterpretable** at these budgets, unchanged
  from PR #24.

## 10. Resuming this work in a new session

State as of 2026-08-29, branch `codex/reyna-statistical-calibration` (off
`codex/reyna-zhang-fullmetric-10pct`, PR #24): **all code for Phases 1, 2, 3,
and 5 is complete, tested, and committed.** 85/85 new tests pass; the only
pre-existing failure is `test_clinical_consistency_target_tiers.m` (a PVR
tier assertion, present on `main`, unrelated to this branch). Nothing further
needs to be built before the final A/B run — this section exists so a fresh
session (or a fresh model instance) can execute it without re-deriving
anything above.

### 10.1 The one rule for this run

**Do not edit any file this run depends on while it is executing.** §3.2
documents exactly what goes wrong: `main_run.m` freezes its own control flow
at launch (it stays on the call stack for the whole run), while functions it
calls are reloaded fresh on each invocation — so a mid-run edit produces a
silent hybrid of old and new code. If a genuine bug is found while a run is
in flight, let the run finish (or kill it) before touching the file.

### 10.2 The command

```bash
matlab -batch "cd('D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd'); addpath(genpath(pwd)); setenv('UNIFIED_VSD_UQLAB_PATH', fullfile(pwd,'toolbox','UQLab_Rel2.2.0','core')); setenv('UNIFIED_VSD_SCALING_MODE','zhang'); setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','1'); setenv('UNIFIED_VSD_DO_GSA','1'); setenv('UNIFIED_VSD_GSA_PCE_N','128'); setenv('UNIFIED_VSD_DO_PLOTS','0'); setenv('UNIFIED_VSD_MAX_FUN_EVALS','300'); setenv('UNIFIED_VSD_MAX_ITERATIONS','40'); setenv('UNIFIED_VSD_NUM_STARTS','6'); setenv('UNIFIED_VSD_MULTISTART_SEED','20260828'); setenv('UNIFIED_VSD_OBJECTIVE_WEIGHTING','sigma'); main_run('pre_surgery', patient_reyna())"
```

Run this in the background and redirect output to a log file; do not run any
other MATLAB process concurrently (a prior 20-minute run stretched past 75
minutes when something else was competing for the same machine). Expected
duration: roughly 70–80 minutes per start × 6 starts ≈ **4–5 hours** based on
the timing actually observed on 2026-08-29 (slower than this PRD's original
~100-minute estimate — the six-stage-per-start pipeline, each stage its own
`fmincon` call with a full ODE steady-state solve per function evaluation, is
the reason; see the "why does it take so long" exchange in this session for
the full breakdown if useful).

### 10.3 What "done" looks like

The run succeeds when the console log contains a `FULL METRIC 10% GATE` block
followed by a `CHI-SQUARED` block with **`p (active parameters)` equal to
`numel(calib_out.names)` for the winning start — a small positive integer,
not 0** (0 was the symptom of the contamination bug this section exists to
prevent a repeat of). A `PARAMETER IDENTIFIABILITY` block should also appear
(Phase 3, previously missing entirely from the contaminated run).

On completion, the run folder's `tables/` directory should contain
`full_metric_gate_pre_surgery.csv`, `chi_squared_pre_surgery.csv`,
`parameter_identifiability_pre_surgery.csv`, and
`parameter_identifiability_pairs_pre_surgery.csv`. Fill in this document's
`<!-- SIGMA_* -->` placeholders in §3.2 and the report in §7 from these
files, then compare against the legacy PR #24 numbers already tabulated.

### 10.4 After a good result

If the sigma-weighted candidate is an improvement worth keeping, it can be
baked into `config/calibration_recipes/reyna_pre_surgery.m` as a new
`initial_parameter_values` seed (the existing seed there came from a 2026-05-21
run via the same pattern) — future runs would then take the
"Accepted explicit recipe seed" fast path and skip the multi-stage
optimisation entirely, rather than re-paying this multi-hour cost. This is a
deliberate follow-up decision, not something to do automatically.
