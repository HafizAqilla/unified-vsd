# Publication-Readiness Execution PRD

Status: **written 2026-09-05, executing same day.** This is the Phase 5
deliverable of the publication-readiness cleanup (see the approved plan in
this repo's PR history / `docs/CHANGES_SINCE_PR22.md` §14). It exists so the
compute-heavy Phase 6 steps below have a written contract to point back to,
independent of any one chat session's context window.

## Why this exists

Phases 0-4 (freeze, de-identification, dead-code removal, code-defect fixes,
protocol-data reconciliation, citations) are complete and committed on
`codex/publication-readiness`. §14 of `docs/CHANGES_SINCE_PR22.md` states
plainly that the protocol-data reconciliation (HR 136→119, VSD diameter
3.025→3.665, chamber volumes NaN→consistency-only 32/23.6/30.5/12)
**invalidates every existing calibration result in this repository.** No run
has yet been executed against the corrected inputs. This document is the
contract for that run.

## Guardrails (carried forward from the approved plan, G1-G7, plus G8-G9)

- **G1** Never change a target to improve a fit.
- **G2** Never shrink the denominator in a change that reports a pass count.
- **G3** Report every arm, including failures — a worse result is still a
  result and must be stated, not minimized or omitted.
- **G4** Never invent values. If a number is needed and not available, stop
  and say so rather than approximate it silently.
- **G5** (implicit throughout Phases 0-4) No commit that is not itself a
  calibration re-run may change an exported number.
- **G8** Only `fair_prior_zhang` / `fair_prior_lundquist` arms (historical
  seeds disabled) may be quoted as *the* scaling-method comparison. The
  `operational_*` arms measure warm-start vs cold-start, not Zhang vs
  Lundquist, because only the Lundquist operational path carries a
  historical accepted-seed warm start (`config/calibration_recipes/
  reyna_pre_surgery.m`). Reporting operational RMSE deltas as if they were a
  fair scaling comparison repeats the mistake found in PR #23.
- **G9** No file any currently-running calibration depends on
  (`config/patient_reyna.m`, `config/calibration_recipes/
  reyna_pre_surgery.m`, `main_run.m`, anything under `src/calibration/`,
  `src/utils/`, `src/validation/`, `config/experiments/
  reyna_p1_scaling_v1.m`) may be edited while a run is in flight. A prior
  mid-run edit here cost roughly 3.4 hours of compute and silently mixed
  code versions inside one run (2026-08-29 incident, see
  `docs/reyna_statistical_calibration_results_20260829.md` §0). Editing
  *unrelated* files (this PRD, a results-summary doc not yet written) while
  a run is in flight is fine.

## Environment

- MATLAB R2025a at `C:\Program Files\MATLAB\R2025a\bin\matlab.exe`.
- `-batch` needs an explicit `cd` + `addpath(genpath(pwd))` — it does not
  inherit the shell's working directory.
- UQLab is installed locally at `toolbox/UQLab_Rel2.2.0`; every GSA-enabled
  run needs `UNIFIED_VSD_UQLAB_PATH` set to that absolute path before
  MATLAB starts (env vars set with `setenv` *inside* the same `-batch`
  invocation, since each `-batch` call is a fresh process).
- Never run two MATLAB processes against this repo concurrently — both
  would write to `results/runs/` and could race on the same clinical-lock
  verification.
- Wall-clock budget observed previously: roughly 70-80 minutes for a single
  full multi-start calibration run at publication budget. A 6-start
  multi-start run at that budget is proportionally longer. This PRD's step
  sequencing (below) is deliberately built around that cost.

## Step sequence and per-step success criteria

### Step 1 — Scaling head-to-head (screening budget)

Command:
```matlab
cd('D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd');
addpath(genpath(pwd));
setenv('UNIFIED_VSD_UQLAB_PATH', fullfile(pwd, 'toolbox', 'UQLab_Rel2.2.0'));
result = run_reyna_scaling_experiment('config/experiments/reyna_p1_scaling_v1.m', ...
    'Repeats', 1, 'DoGSA', true, 'GsaN', 128, ...
    'MaxFunctionEvaluations', 300, 'MaxIterations', 40);
```

**Deviation from the original plan text, stated explicitly (G4):** the
original plan phrase "GSA N=128, 300 evals, 6 starts" is used here for GSA
sample count and function-evaluation budget, but **not** for multi-start
count: `run_reyna_scaling_experiment.m` does not expose a multi-start
parameter, and `src/calibration/run_calibration.m:173` defaults
`UNIFIED_VSD_NUM_STARTS` to **1** when unset. Running all 4 arms x 6 starts
x repeats at this budget was estimated (from the 70-80 min/full-multistart
figure above) to run into many hours to low tens of hours — impractical for
a single screening pass whose purpose is to pick a winner, not to be the
final reported number. Step 1 therefore runs single-start (`NUM_STARTS=1`,
the tool's own default) as a **screening** pass, and reserves the full
6-start multi-start budget for Step 2, applied only to the winning arm. This
is stated as a limitation of Step 1's numbers, not hidden.

Repeats is 1, not the contract's default of 3, for the same time-budget
reason — again a screening-pass deviation, stated explicitly. If time
allows after Step 2-5 complete, additional repeats of the fair-prior arms
can be added to strengthen the head-to-head; not doing so does not block
the rest of this PRD.

**Success criteria:**
- All 4 arms produce a row in `reyna_scaling_experiment_summary.csv`, pass
  or fail — a failed arm is still reported (G3).
- `decision_summary.txt` names a `FairWinner` distinct from
  `OperationalWinner` reporting (G8): only the fair-prior comparison is
  used to answer "does Zhang or Lundquist scale Reyna's data better".
- Every arm's manifest confirms `ClinicalProfilePolicy: frozen_input_profile`
  and the fair-prior arms confirm `HistoricalSeedPolicy: disabled_fair_prior`
  (enforced by the runner itself, `verify_clinical_profile_lock` /
  `enrich_record_from_run` — a run that fails this check errors loudly
  rather than silently reporting a compromised comparison).

### Step 2 — Full multi-start on the winning fair-prior arm

Take `decision.fair_winner` from Step 1 (`fair_prior_zhang` or
`fair_prior_lundquist`). Run once more, same GSA/eval budget, with
`UNIFIED_VSD_NUM_STARTS=6`:

```matlab
setenv('UNIFIED_VSD_NUM_STARTS', '6');
setenv('UNIFIED_VSD_MULTISTART_SEED', '20260905');
result2 = run_reyna_scaling_experiment('config/experiments/reyna_p1_scaling_v1.m', ...
    'Arms', {<winning arm id>}, 'Repeats', 1, 'DoGSA', true, 'GsaN', 128, ...
    'MaxFunctionEvaluations', 300, 'MaxIterations', 40);
```

**Success criteria:** ACCEPT or a stated reason it did not reach ACCEPT
(G3); multi-start report in the run folder shows more than one distinct
converged basin was actually explored (guards against "the optimizer never
moved" masquerading as convergence).

### Step 3 — Parameter reduction (p=12 → p=7)

`scripts/analyse_parameter_reduction.m` exists and was previously analysed
but never run against real data (`docs/CHANGES_SINCE_PR22.md` §8.2:
`cond(S)` 2060 → 21.7, `dof` -3 → +2, but "not yet run"). Run it against the
Step 2 candidate's identifiability report.

**Success criteria:** a concrete reduced parameter list is produced, with
`dof` positive at the reduced-p operating point (this is what makes χ²/N
interpretable at all, per the existing analysis).

### Step 4 — Out-of-sample post-closure prediction

```matlab
scripts/evaluate_post_closure_prediction.m
```
against the Step 2 (or Step 3, if the reduced-parameter set changes the
accepted candidate) fitted parameters, predicting post-surgery pressures
without fitting to them.

**Success criteria:** report exists whether the prediction is good or bad
(G3) — this is the genuine independent-data test in this repo's pipeline
(`SVR` was previously rejected as a holdout candidate because it is algebra
over already-fitted quantities, not new information; post-closure pressures
are the real thing).

### Step 5 — Identifiability report at the final operating point

`src/calibration/analyse_parameter_identifiability.m`, already wired into
`main_run.m`, run at the final (Step 3 or Step 2) operating point. Success:
condition number and pairwise-correlation table exported, `dof` positive.

## What happens after Step 5

Phase 7 (publication artifacts): a single authoritative results document,
inline-banner the two superseded pre-reconciliation results docs
(`reyna_zhang_fullmetric_results_20260828.md`,
`reyna_scaling_results_20260827.md`), rewrite `README.md`'s Current Status
and Latest Changes sections with the real numbers from Steps 1-5, and open
the PR with an honest summary — including whatever got worse, per G3.
