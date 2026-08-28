# Reyna Zhang Full-Metric Campaign — Results, 2026-08-28

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

<!-- RESULTS_TABLE_PLACEHOLDER -->

## 4. Reading the denominator honestly

D7 removes `RVESV` from the governed set, taking it from 10 rows to 9. A
denominator that shrinks alongside a rising pass count is exactly the move this
project's own audit warned against, so the comparison below is stated on the
**fixed 9-row haemodynamic set** — the rows that are pre-operative measurements
or direct transforms of them — for both baseline and candidate:

`RAP_mean`, `PAP_min`, `PAP_max`, `PAP_mean`, `SAP_min`, `SAP_max`,
`SAP_mean`, `QpQs`, `CO_Lmin`.

On that fixed set the baseline scored **7 of 9** (`PAP_min` 23.83%,
`SAP_min` 13.98%).

<!-- COMPARISON_PLACEHOLDER -->

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
- **Negative degrees of freedom.** 14 fitted parameters against 8 independent
  measurements (`RAP_mean`, `PAP_dia`, `PAP_sys`, `PAP_mean`, `SAP_dia`,
  `SAP_sys`, `Qp`, `QpQs`); everything else is an algebraic transform. A gate
  count is a statement about flexibility until the active set is reduced.
- **No GSA.** UQLab is not installed on this host, so PRD Phase 4 — Sobol
  screening and parameter reduction — has not run. `UNIFIED_VSD_UQLAB_PATH` is
  supported and the preflight still fails loudly rather than degrading.
- **No confidence intervals.** The multi-start machinery is in place and
  deterministic, but the published candidate here is a single start. Reporting
  an interval requires PRD Phase 3 at `NumStarts >= 16`.
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
