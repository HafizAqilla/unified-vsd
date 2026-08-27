# PRD: Luna Autonomous Reyna Scaling Experiment

Status: implementation-ready and executable
Version: `reyna_p1_scaling_v1`
Scope: Reyna `pre_surgery` only

## Objective

Determine whether the Zhang or Lundquist-BSA scaling prior produces the most
publishable Reyna result when the same clinical inputs, target governance,
calibration code, and acceptance gates are used. The experiment must separate
the effect of the scaling prior from the effect of inherited historical recipe
seeds.

The experiment is not allowed to silently change the Reyna clinical profile,
promote derived variables into the primary objective, or replace a rejected
candidate with a lower-RMSE candidate without reporting the rejection.

## Immutable scientific contract

- Patient: `patient_reyna()`.
- Scenario: `pre_surgery`.
- Primary metrics: `RAP_mean`, `PAP_mean`, `SAP_mean`, `QpQs`, `CO_Lmin`.
- Target tiers, uncertainty weights, consistency-only fields, and acceptance
  thresholds remain those of the Reyna pre-surgery recipe.
- Clinical profile is frozen during each experiment arm.
- Historical calibrated parameter and initial-condition seeds are disabled for
  fair-prior arms.
- The GSA contract remains `N=128`; reducing it must be reported as a failed
  or preliminary run, not as a completed publication GSA.
- Existing historical run folders are read-only evidence and are never edited.

## Four-arm design

| Arm | Scaling | Historical seeds | Purpose |
|---|---|---:|---|
| `fair_prior_zhang` | Zhang | no | Isolate Zhang from the inherited recipe seed |
| `fair_prior_lundquist` | Lundquist BSA | no | Isolate Lundquist from the inherited recipe seed |
| `operational_zhang` | Zhang | yes | Reproduce the current operational route |
| `operational_lundquist` | Lundquist BSA | yes | Reproduce the current operational route |

## Autonomous execution behavior

Luna should execute the following loop without changing the contract:

1. Run a dry contract check.
2. Run one fast screening repeat per selected arm with GSA disabled only when
   the external UQLab dependency is unavailable.
3. Repeat the fair-prior arm(s) at the screening budget to check deterministic
   repeatability.
4. Run matched higher-budget confirmations with polish enabled.
5. Run the GSA path at exactly 128 samples when UQLab is available.
6. Validate the run manifest, frozen clinical snapshot, target-tier export,
   primary RMSE export, calibration status, and parameter plausibility.
7. Select a winner separately for the fair-prior and operational families.
   Prefer `ACCEPT`; use validated RMSE only as a labelled fallback when a
   family has no accepted candidate.
8. Stop method promotion if the best candidate has unresolved critical
   clinical-data inconsistency, out-of-bound parameters, or a failed primary
   acceptance gate.
9. Export the run contract, typed summary CSV, decision summary, and run paths.

The runner is `scripts/run_reyna_scaling_experiment.m` and the locked contract
is `config/experiments/reyna_p1_scaling_v1.m`.

## Acceptance criteria

A result is a publishability candidate only if all of the following are true:

- runner status is `validated`;
- calibration status is `ACCEPT`;
- all five primary metrics are within the 10% patient-acceptance gate;
- primary target membership is unchanged;
- clinical profile snapshot matches `patient_reyna()` exactly;
- no parameter is outside its registry bounds;
- plausibility warnings are reported and scientifically reviewed;
- critical clinical consistency findings remain visible rather than hidden.

An arm may still be retained as a sensitivity result when it is validated but
not accepted. It must not be described as the publication winner.

## Required execution command

From the repository root:

```matlab
addpath(genpath(pwd));
result = run_reyna_scaling_experiment( ...
    'config/experiments/reyna_p1_scaling_v1.m', ...
    'Arms', {'fair_prior_zhang','fair_prior_lundquist', ...
             'operational_zhang','operational_lundquist'}, ...
    'Repeats', 3, ...
    'DoGSA', true, ...
    'GsaN', 128, ...
    'ScreeningMode', false, ...
    'MaxFunctionEvaluations', 60, ...
    'MaxIterations', 8);
```

If UQLab is installed outside the repository, set
`UNIFIED_VSD_UQLAB_PATH` before the MATLAB run. SoBioS can be supplied through
`UNIFIED_VSD_SOBIOS_PATH`. If UQLab is absent, the runner must preserve the
failed GSA evidence and may run a separate `DoGSA=false` preliminary screen;
that screen is not a GSA-complete result.

## Deliverables

- code and tests on a dedicated `codex/` branch;
- MATLAB run folders under ignored `results/runs/`;
- experiment summaries under ignored `results/luna_experiments/`;
- this PRD;
- the execution plan and dated results report;
- a GitHub PR whose body states the exact branch, controls, validation status,
  and whether the result is accepted or preliminary.
