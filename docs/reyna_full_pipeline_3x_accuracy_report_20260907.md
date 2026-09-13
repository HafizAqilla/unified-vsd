# Reyna Full Pipeline — 3× Pre-Surgery and 3× Post-Surgery Accuracy Report

**Date:** 2026-09-07. **Scope:** every clinical metric the model reports, for
every one of 6 independent full-pipeline runs (3 pre-surgery + 3
post-surgery), with no metric omitted. This supersedes the single-run
numbers in `docs/reyna_publication_readiness_results_20260906.md` for the
pre-surgery hemodynamic fit; that document's Step 4 out-of-sample
post-closure *prediction* (forward-simulating the pre-surgery fit with the
shunt closed) is a different analysis from the post-surgery *calibration*
runs reported here, and both are valid, distinct things to know.

## 0. What was actually run

Six full pipeline executions, each independent end to end (scaling → GSA
screening → 6-start multi-start calibration → validation/export), on the
corrected clinical data (`docs/CHANGES_SINCE_PR22.md` §14):

- **3× pre-surgery**, Zhang scaling, fair prior (historical seeds
  disabled), frozen clinical profile, GSA N=128, 300 function evaluations,
  40 iterations, 6-start multi-start, multi-start seeds **20260906 /
  20260907 / 20260908** (one distinct seed per repeat — using the same
  seed three times would have produced bit-identical results, which was
  checked and deliberately avoided).
- **3× post-surgery**, same configuration, same three seeds. **No
  dedicated post-surgery calibration recipe exists yet for Reyna**, so
  these runs used the pipeline's default target-tier governance, which — unlike
  the pre-surgery recipe — fits the chamber volumes (LVEDV, LVESV, RVEDV,
  RVESV) as hard/soft targets rather than consistency-only. This is a real
  difference in what "post-surgery calibration" means here versus the
  pre-surgery recipe's philosophy, and it is why LVEF/RVEF appear
  consistency-only below while the four volumes are fitted.

All six runs are tracked in git (evidence CSVs only, per the retention
rule) under `results/runs/<timestamp>_reyna_<scenario>/`.

## 1. Top-level summary, all 6 runs

| Scenario | Repeat | Seed | Status | Primary RMSE (governed) | Governed gate | χ²/N | cond(S) | Wall time |
|---|---|---|---|---:|---:|---:|---:|---:|
| Pre-surgery | 1 | 20260906 | REJECT | 0.0627 | 8/9 | 1.29 (consistent) | 421 | 4817.3 s (80.3 min) |
| Pre-surgery | 2 | 20260907 | **ACCEPT** | 0.0546 | **9/9** | 0.85 (consistent) | 299 | 6117.9 s (102.0 min) |
| Pre-surgery | 3 | 20260908 | REJECT | 0.0822 | 7/9 | 1.88 (consistent) | 325 | 5367.4 s (89.5 min) |
| Post-surgery | 1 | 20260906 | PHYSIOLOGICAL_BUT_POOR_FIT | 0.1550 | 8/11 | 3.52 (underfit) | 12.8 | 2758.8 s (46.0 min) |
| Post-surgery | 2 | 20260907 | PHYSIOLOGICAL_BUT_POOR_FIT | 0.1477 | 8/11 | 2.55 (underfit) | 10.5 | 3356.5 s (55.9 min) |
| Post-surgery | 3 | 20260908 | PHYSIOLOGICAL_BUT_POOR_FIT | 0.1488 | 8/11 | 3.04 (underfit) | 10.1 | 2800.5 s (46.7 min) |

Total wall-clock across all 6 runs: 26,568 s ≈ 7.38 hours.

**Immediate honest finding: the pre-surgery calibration status is NOT
stable across seeds.** χ²/N is consistently good (0.85–1.88, always in the
"consistent" band) across all three, but only one of three seeds (repeat
2) reaches full ACCEPT with a perfect 9/9 governed gate — the other two
each miss the gate on one or two specific metrics (repeat 1: `PAP_max`
only; repeat 3: `PAP_min` and `SAP_min`). The earlier single-run result in
`docs/reyna_publication_readiness_results_20260906.md` (governed gate 7/9,
χ²/N=1.14, `PHYSIOLOGICAL_BUT_POOR_FIT`) used yet a fourth seed
(20260905) not repeated here and landed between these three. **Which
result gets quoted as "the" Reyna pre-surgery fit therefore depends on
which multi-start seed is picked**, not on a single settled number — this
is exactly the kind of instability an uncertainty analysis is supposed to
surface, not hide.

Post-surgery is comparatively **stable** across all three seeds: same
governed-gate count (8/11) every time, and — as detailed below — the same
two metrics (`LVESV`, `SAP_max`) fail in every single run, meaning that is
a structural miss, not seed noise. Post-surgery is also far better
statistically conditioned than pre-surgery (cond(S) ≈ 10–13 vs. 299–421),
because it has close to as many parameters as governed metrics (11 vs.
11) rather than pre-surgery's underdetermined 12-vs-9.

---

## 2. Pre-surgery: every metric, every repeat (11 metrics, none omitted)

Mean ± SD is the sample standard deviation across the 3 repeats (n=3).

| Metric | Unit | Tier | Clinical | Repeat 1 | Repeat 2 | Repeat 3 | Mean model ± SD | Mean \|error\| % ± SD | Gate pass |
|---|---|---|---:|---:|---:|---:|---:|---:|:---:|
| RAP_mean | mmHg | hard | 5.000 | 5.300 | 5.261 | 5.497 | 5.353 ± 0.127 | 7.05 ± 2.53 | 3/3 |
| PAP_mean | mmHg | hard | 15.000 | 15.733 | 15.718 | 14.928 | 15.460 ± 0.461 | 3.39 ± 2.52 | 3/3 |
| SAP_mean | mmHg | hard | 77.000 | 73.274 | 74.093 | 82.901 | 76.756 ± 5.337 | 5.43 ± 2.01 | 3/3 |
| CO_Lmin | L/min | hard | 3.423 | 3.269 | 3.213 | 3.189 | 3.224 ± 0.041 | 5.82 ± 1.20 | 3/3 |
| QpQs | – | hard | 1.194 | 1.200 | 1.196 | 1.141 | 1.179 ± 0.033 | 1.73 ± 2.37 | 3/3 |
| PAP_max | mmHg | soft | 20.000 | 22.752 | 21.963 | 18.801 | 21.172 ± 2.091 | 9.86 ± 3.88 | **2/3** (fails repeat 1) |
| PAP_min | mmHg | soft | 10.000 | 10.078 | 10.330 | 11.249 | 10.552 ± 0.616 | 5.52 ± 6.16 | **2/3** (fails repeat 3) |
| SAP_max | mmHg | soft | 100.000 | 92.520 | 92.352 | 103.171 | 96.014 ± 6.198 | 6.10 ± 2.54 | 3/3 |
| SAP_min | mmHg | soft | 57.000 | 55.879 | 56.993 | 64.683 | 59.185 ± 4.794 | 5.15 ± 7.28 | **2/3** (fails repeat 3) |
| Q_shunt_Lmin | L/min | soft (primary-RMSE holdout, not gated) | 0.664 | 0.655 | 0.630 | 0.449 | 0.578 ± 0.113 | **12.90 ± 16.95** | n/a |
| SVR | WU | derived_validation (not gated) | 21.034 | 20.791 | 21.426 | 24.269 | 22.162 ± 1.852 | 6.13 ± 8.01 | n/a |

**Notable within this table:**
- Every hard metric (RAP_mean, PAP_mean, SAP_mean, CO_Lmin, QpQs) passes
  the gate in all 3 repeats — the model's core hemodynamic fit is
  genuinely stable.
- The three soft-tier misses are seed-specific, not the same metric every
  time: repeat 1 fails `PAP_max` only; repeat 3 fails `PAP_min` and
  `SAP_min`; repeat 2 fails nothing.
- `Q_shunt_Lmin`, although never part of the governed gate, has by far the
  largest relative spread of any metric (error ranging 1.30%–32.35%
  across the three repeats, SD 16.95 percentage points) — this residual
  shunt-flow estimate is the least seed-stable quantity the model
  produces, worth flagging for anyone tempted to quote it as a fitted
  result.
- LVEF/LVEDV/LVESV/RVEDV/RVESV/RVEF **do not appear in this table at
  all**, because — per the 2026-09-06 timing correction — there is no
  confirmed pre-surgery chamber-volume measurement for Reyna; the
  protocol form's volumes belong to the post-surgery state (Section 3).
  This is a deliberate absence, not an omission.

---

## 3. Post-surgery: every metric, every repeat (13 metrics, none omitted)

| Metric | Unit | Tier | Clinical | Repeat 1 | Repeat 2 | Repeat 3 | Mean model ± SD | Mean \|error\| % ± SD | Gate pass |
|---|---|---|---:|---:|---:|---:|---:|---:|:---:|
| RAP_mean | mmHg | hard | 5.000 | 4.915 | 5.143 | 5.018 | 5.025 ± 0.114 | 1.64 ± 1.25 | 3/3 |
| PAP_mean | mmHg | hard | 13.000 | 13.216 | 12.831 | 13.383 | 13.143 ± 0.283 | 1.97 ± 0.87 | 3/3 |
| SAP_mean | mmHg | hard | 79.000 | 86.894 | 86.889 | 86.870 | 86.884 ± 0.012 | 9.98 ± 0.02 | 3/3 |
| LVEDV | mL | hard | 32.000 | 34.463 | 31.442 | 30.915 | 32.273 ± 1.914 | 4.28 ± 3.07 | 3/3 |
| LVESV | mL | hard | 23.600 | 13.872 | 13.612 | 13.563 | 13.682 ± 0.166 | **42.02 ± 0.70** | **0/3 — fails every repeat** |
| PAP_min | mmHg | soft | 9.000 | 8.550 | 9.438 | 8.750 | 8.913 ± 0.466 | 4.22 ± 1.25 | 3/3 |
| SAP_min | mmHg | soft | 68.000 | 67.377 | 69.514 | 69.870 | 68.920 ± 1.348 | 1.96 ± 0.95 | 3/3 |
| PAP_max | mmHg | soft | 17.000 | 19.847 | 16.746 | 19.316 | 18.636 ± 1.659 | 10.62 ± 8.06 | **1/3** (fails repeats 1 and 3) |
| SAP_max | mmHg | soft | 89.000 | 107.159 | 104.616 | 104.149 | 105.308 ± 1.620 | **18.32 ± 1.82** | **0/3 — fails every repeat** |
| RVEDV | mL | soft | 30.500 | 32.788 | 28.655 | 29.298 | 30.247 ± 2.224 | 5.83 ± 1.79 | 3/3 |
| RVESV | mL | soft | 12.000 | 12.080 | 10.686 | 11.817 | 11.528 ± 0.741 | 4.38 ± 5.71 | **2/3** (fails repeat 2) |
| LVEF | – | consistency_check_only (not gated) | 0.2625 | 0.5975 | 0.5671 | 0.5613 | 0.5753 ± 0.0194 | **119.16 ± 7.41** | n/a |
| RVEF | – | consistency_check_only (not gated) | 0.6066 | 0.6316 | 0.6271 | 0.5967 | 0.6184 ± 0.0190 | 3.04 ± 1.27 | n/a |

**Notable within this table:**
- **`LVESV` fails all three repeats by a nearly identical margin
  (41.2%, 42.3%, 42.5% — SD only 0.70 percentage points).** This is not
  seed noise; it is a structural, reproducible miss. The model
  consistently predicts a post-closure LV end-systolic volume around
  13.5–13.9 mL against a measured 23.6 mL — always in the same direction
  (under-prediction).
- **`SAP_max` also fails all three repeats (17.0–20.4% error, always
  over-predicted)** — the model consistently predicts a higher
  post-closure systolic systemic pressure than measured.
- `LVEF` (not gated, consistency-only) is thrown off by the same LVESV
  miss: since EF = (EDV−ESV)/EDV and ESV is under-predicted, LVEF comes
  out roughly 113–128% too high in every run — a direct, mechanical
  consequence of the LVESV finding above, not an independent second
  problem.
- `RVEDV`, `RVESV`, `RVEF`, `PAP_min`, `SAP_min`, `RAP_mean`, `PAP_mean`
  all pass in either 2 or 3 of 3 repeats — a genuinely reasonable fit on
  the right-heart and lower-pressure side.
- `PAP_max` is seed-sensitive (fails 1, passes 2, fails 3) similar to the
  pre-surgery soft-tier pattern.

---

## 4. Reading this alongside the earlier post-closure *prediction*

`docs/reyna_publication_readiness_results_20260906.md` §4 reports a
different thing: forward-simulating the **pre-surgery-fitted** parameters
with the shunt closed, and comparing against the same post-surgery
measurements — a genuine zero-shot, no-refit test. That analysis found
LVEF over-predicted by 115% (0.564 vs. 0.263) via the same mechanism (ESV
under-predicted). **The post-surgery calibration runs in this report
independently reproduce that same LVEF/ESV problem even when the model is
allowed to fit post-surgery data directly** (LVEF still 113–128% too high,
LVESV still 41–43% too low, in every one of 3 repeats) — which is a
stronger and more specific finding than either analysis alone: this is not
an artifact of using pre-surgery parameters unchanged, since directly
calibrating to post-surgery data does not fix it either. Something about
how the model represents post-closure LV end-systolic volume specifically
does not fit this patient's data, independent of which parameters are
allowed to move.

## 5. What this means, stated plainly

- The pre-surgery hemodynamic (pressure/flow) fit is genuinely good and
  stable across repeats; the pass/fail *label* (ACCEPT vs. REJECT) is
  not stable, because it hinges on one or two soft-tier metrics that
  land just outside or inside the 10% gate depending on which multi-start
  seed is used.
- The post-surgery chamber-volume fit has one specific, reproducible
  failure mode (LV end-systolic volume, and everything downstream of it)
  that appears in every repeat regardless of seed, and is corroborated
  independently by the separate pre-surgery-to-post-closure prediction
  test. This is the single most concrete, repeatable limitation this
  model currently has for Reyna.
- `Q_shunt_Lmin` (pre-surgery) is the least seed-stable quantity in the
  entire report and should not be quoted as a precise number without
  stating its spread.

## 6. Where the raw data lives

All 6 run folders (evidence CSVs tracked in git):
- `results/runs/20260906_213635_reyna_pre_surgery/`
- `results/runs/20260906_225713_reyna_pre_surgery/` (the ACCEPT run)
- `results/runs/20260907_003935_reyna_pre_surgery/`
- `results/runs/20260907_021019_reyna_post_surgery/`
- `results/runs/20260907_025625_reyna_post_surgery/`
- `results/runs/20260907_035228_reyna_post_surgery/`

Each folder's `tables/full_metric_gate_<scenario>.csv` is the exact source
of every number in Sections 2 and 3 above.
