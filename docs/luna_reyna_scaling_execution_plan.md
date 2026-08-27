# Luna Execution Plan: Reyna Scaling and Calibration

## Phase 1 — Lock the comparison

1. Refresh `origin/main` and work in an isolated branch.
2. Load `config/experiments/reyna_p1_scaling_v1.m`.
3. Assert the scenario, patient, five primary metrics, four arms, and GSA
   sample count before any simulation starts.
4. Preserve the original checkout and all historical result folders.

## Phase 2 — Make provenance executable

For every `main_run` invocation, set and record:

- scaling mode;
- clinical-profile policy;
- historical-seed policy;
- experiment ID;
- optimizer evaluation and iteration budgets;
- GSA enabled state and requested sample count.

The run package must contain a clinical snapshot that is byte-for-value
equivalent to a fresh `patient_reyna()` snapshot for the `common` and
`pre_surgery` sections.

## Phase 3 — Screening loop

Run a small, deterministic screening budget to catch configuration problems
and gross method differences. Screening results are diagnostic only. They must
not replace the higher-budget confirmation because optional polish stages can
change the result substantially.

For fair-prior repeats, expect deterministic output unless a future optimizer
or model change introduces a controlled random seed. Report repeated values and
their spread; do not imply statistical uncertainty from identical deterministic
runs.

## Phase 4 — Matched confirmation loop

Run Zhang and Lundquist with the same higher-budget settings and polish enabled.
For each run:

1. retain the baseline candidate;
2. retain the best raw optimizer candidate;
3. apply the pipeline's gate and plausibility guards;
4. report the accepted/scientific/rolled-back candidate separately;
5. read the exported RMSE and target-tier tables rather than parsing console
   output alone.

The comparison must include both:

- fair prior: scaling effect without historical seeds;
- operational route: what the current production recipe does with its seed.

## Phase 5 — GSA path

When UQLab is available:

1. run initial PCE GSA with `N=128`;
2. build the optimization mask from the registry-backed bounds and governed
   metrics;
3. calibrate;
4. run final PCE GSA at `N=128`;
5. verify checkpoint signature, metric membership, nominal point, bounds, and
   scenario before reusing any checkpoint.

When UQLab is unavailable, the code must fail explicitly with an actionable
preflight message. A calibration-only run can continue as a preliminary
diagnostic, but the PR must state that the GSA requirement is pending.

## Phase 6 — Decision logic

Use this order of precedence within each family:

1. validated `ACCEPT` candidates;
2. among them, lowest mean governed primary RMSE;
3. if no accepted candidate exists, lowest validated RMSE as a fallback,
   labelled `validated_rmse_fallback_no_accept`.

Never select a method only because it has a low full RMSE if it fails a hard
primary target, violates a bound, or hides the clinical consistency audit.

## Phase 7 — PR and handoff

Before opening the PR:

- run `git diff --check`;
- run the contract and targeted MATLAB tests;
- confirm generated results remain ignored;
- record the exact run directories and summary CSVs;
- commit only source, tests, configuration, and Markdown documentation;
- push the dedicated `codex/` branch;
- open a PR against `main` with the evidence table and known limitations.

## Stop conditions

Stop automatic method promotion and report a blocker when:

- clinical values differ between the input snapshot and the run package;
- primary target membership changes;
- a derived target is promoted into governed RMSE;
- the calibration status is rejected or poor-fit;
- a parameter is out of registry bounds;
- UQLab is absent for a run declared GSA-complete.

