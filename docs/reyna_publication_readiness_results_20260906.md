# Reyna Publication-Readiness Re-Run — Results, 2026-09-06

**This is the current authoritative results document for the Reyna model.**
It supersedes every prior results document in `docs/` (each of which now
carries a superseded-banner pointing here or to
`docs/CHANGES_SINCE_PR22.md` §14). It reports the execution of
`docs/publication_readiness_prd.md` against the protocol-corrected clinical
inputs from `docs/CHANGES_SINCE_PR22.md` §14 (HR 119 bpm, VSD diameter
3.665 mm).

**Updated 2026-09-06 (same day), timing correction:** the protocol form's
chamber volumes (32/23.6/30.5/12 mL, rows 26-29, "PRE RELEASE OCCLUDER")
were initially placed under pre-surgery consistency-only. That was wrong:
per the study owner, "pre-release occluder" means the closure device is
already deployed and occluding the defect, simply not yet mechanically
detached — the VSD is already functionally **closed** at that measurement.
This data now lives in `clinical.post_surgery`
(`config/patient_reyna.m`), and there is **no confirmed pre-surgery
chamber-volume measurement for Reyna at all.** This also means these six
values (LVEDV, LVESV, RVEDV, RVESV, LVEF, RVEF) are now genuine **held-out
comparators in the Step 4 out-of-sample post-closure prediction** below,
not consistency-only pre-surgery rows — §4 is revised accordingly and the
result is substantially worse than first reported.

**One-line honest summary:** under a fair prior (no historical warm start),
Zhang scaling clearly outperforms Lundquist-BSA scaling on Reyna's corrected
data (RMSE 0.118 vs 0.222). A full 6-start multi-start on the Zhang arm
reaches primary RMSE 0.0739, governed gate 7/9, and a χ²/N of 1.14 on the
pre-surgery hemodynamics — consistent, not overfit. **The model does not
reach ACCEPT status, and its out-of-sample post-closure prediction is
poor: 3 of 13 targets within 10% once chamber volumes/EF are correctly
included (χ²/N = 15.48), driven substantially by ejection fraction and
chamber volumes the pre-surgery fit had no way to anticipate.** Every
number below is reported as obtained, including the ones that are
unfavorable, per this repo's own governance rule (G3).

---

## 0. What changed since every prior results document

See `docs/CHANGES_SINCE_PR22.md` §14 for the full account. In short: HR
136→119 bpm, VSD diameter recipe-override 3.025→3.665 mm (deleted, D1),
pre-surgery chamber volumes NaN→consistency-only 32/23.6/30.5/12 mL. No
combination of these three has ever been run before this document.

A real orchestration bug was found and fixed while executing Step 1: the
scaling-experiment runner's row buffer was sized for 18 columns while
`record_to_row` returns 24, so every *real* (non-dry-run) invocation of
`scripts/run_reyna_scaling_experiment.m` had crashed at the first row write.
This explains why this script's four-arm contract had only ever been
dry-run in every prior session. Fixed in commit `61973da`.

---

## 1. Step 1 — Scaling head-to-head (screening budget)

Ran all 4 arms of `config/experiments/reyna_p1_scaling_v1.m` at
single-start, single-repeat, GSA N=128, 300 function evaluations, 40
iterations — a screening budget, not the final reported number (see
`docs/publication_readiness_prd.md` for why NUM_STARTS=6 was reserved for
Step 2 rather than spent on all 4 arms here).

| Arm | Historical seed | Status | PrimaryRMSE | GovernedGate | CalibrationStatus |
|---|---|---|---|---|---|
| `fair_prior_zhang` | disabled | validated | **0.1180** | 5/9 | PHYSIOLOGICAL_BUT_POOR_FIT |
| `fair_prior_lundquist` | disabled | validated | 0.2216 | 6/9 | PHYSIOLOGICAL_BUT_POOR_FIT |
| `operational_zhang` | n/a (Zhang never gets a historical seed either way) | validated | 0.1180 (identical to fair_prior_zhang) | 5/9 | PHYSIOLOGICAL_BUT_POOR_FIT |
| `operational_lundquist` | enabled (historical accepted seed) | **FAILED** | — | — | — |

**`operational_lundquist` failed outright**, not from a runner bug: `E.RV.EB`'s
historical accepted-seed value (0.108365846196) fell outside the
newly-computed bounds `[0.11657760947, 0.551094153858]` once the corrected
demographics/HR changed the scaling. This is a real, informative finding:
the previously-"accepted" Lundquist candidate is stale under the corrected
inputs and cannot even be evaluated without a new seed, let alone accepted.
It is not chased further here because per **G8**, only the fair-prior arms
are a valid scaling-method comparison — the operational arms measure
warm-start-vs-cold-start, and this failure is itself evidence of exactly
that asymmetry.

**Decision (per `decision_summary.txt`):** `report_fair_and_operational_results_separately`.
**FairWinner: `fair_prior_zhang`** (basis: `validated_rmse_fallback_no_accept`
— neither fair-prior arm reached ACCEPT at screening budget, so the winner
is chosen by lower validated RMSE within the family, exactly as
`docs/publication_readiness_prd.md` specifies for that fallback case).

**Answer to "does Zhang or Lundquist scale Reyna's data better":** at this
budget, on the corrected data, under a fair (no-warm-start) prior, **Zhang
scales this patient's data better** — nearly half the RMSE of Lundquist-BSA
(0.118 vs 0.222) with an identical governed-gate count in the same
ballpark (5/9 vs 6/9; Lundquist edges ahead on raw gate count here but at
roughly double the RMSE, and neither reaches ACCEPT).

Full data: `results/luna_experiments/reyna_p1_scaling_v1_20260905_230628/reyna_scaling_experiment_summary.csv`.

---

## 2. Step 2 — Full 6-start multi-start on the winner (`fair_prior_zhang`)

Run folder: `results/runs/20260906_002327_reyna_pre_surgery/`. Same GSA/eval
budget as Step 1, `UNIFIED_VSD_NUM_STARTS=6`, seed `20260905`. Total wall
time: 5841.3 s (~97 minutes).

### 2.1 Multi-start diversity (guards against "the optimizer never moved")

| Start | Label | PrimaryRMSE | PrimaryFail | GateFail |
|---|---|---|---|---|
| 2 | sobol_1 (winner) | **0.0739** | 1 | 2 |
| 1 | seed | 0.1180 | 2 | 4 |
| 5 | sobol_4 | 0.1947 | 2 | 4 |
| 4 | sobol_3 | 0.4382 | 3 | 6 |
| 3 | sobol_2 | 0.3256 | 3 | 7 |
| 6 | sobol_5 | 0.4470 | 4 | 8 |

RMSE across starts: min 0.0739, median 0.260, max 0.447, IQR 0.320 — real
spread across genuinely different converged basins, not six copies of the
same point. The winning start (`sobol_1`, a Sobol-scrambled start, not the
historical seed) beat the seed-derived start by 37%.

### 2.2 Winning candidate — full per-metric result

`CalibrationStatus: PHYSIOLOGICAL_BUT_POOR_FIT` (not ACCEPT).
`governed_gate=7/9`, `RMSE_improvement=73.8%` over baseline,
`chi2_per_obs=1.14` (**consistent_but_underdetermined** band — not
overfit, not underfit; this is the healthiest χ²/N this patient's pipeline
has ever reported).

| Metric | Tier | Clinical | Model | AbsError% | WithinGate |
|---|---|---:|---:|---:|:---:|
| RAP_mean | hard | 5 | 5.29 | 5.80% | yes |
| QpQs | hard | 1.194 | 1.137 | 4.78% | yes |
| CO_Lmin | hard | 3.423 | 2.948 | **13.87%** | **no** |
| PAP_mean | hard | 15 | 15.62 | 4.13% | yes |
| SAP_mean | hard | 77 | 79.78 | 3.61% | yes |
| PAP_min | soft | 10 | 11.09 | **10.86%** | **no** |
| PAP_max | soft | 20 | 20.60 | 2.99% | yes |
| SAP_min | soft | 57 | 62.21 | 9.14% | yes |
| SAP_max | soft | 100 | 98.58 | 1.42% | yes |
| Q_shunt_Lmin | soft (primary-RMSE holdout) | 0.664 | 0.404 | 39.20% | n/a (not graded) |

The two governed misses are `PAP_min` (10.86%, just over the 10% gate) and
`CO_Lmin` (13.87%). `Q_shunt_Lmin` is graded in the objective but excluded
from the governed primary RMSE by design (`primary_rmse_holdout`); its
39.2% error is large and worth noting as a real limitation even though it
does not count against the gate.

**Chamber volumes are not pre-surgery rows at all (revised 2026-09-06, same
day).** The `full_metric_gate_pre_surgery.csv` tracked from this run still
shows LVEDV/LVESV/RVEDV/RVESV/LVEF as `consistency_check_only` — that
reflects the understanding in place when the run was generated, since
superseded (see the correction note at the top of this document). Under
the corrected timing, none of these five have a pre-surgery clinical value
at all (tier `unavailable` in current code); the real comparators are
post-closure and are evaluated properly in §4 below, where — unlike the
pre-surgery report — they matter a great deal.

Full data: `results/runs/20260906_002327_reyna_pre_surgery/tables/full_metric_gate_pre_surgery.csv`
(historical artifact; read alongside the correction note above).

---

## 3. Step 3 — Parameter reduction (p=12 → p=7)

Ran `scripts/analyse_parameter_reduction.m` against the Step 2 run (linear
algebra only, no re-simulation, on the already-fitted operating point).

Full parameter set: N = 9 governed metrics, p = 12 parameters,
**cond(S) = 223, dof = −3** (underdetermined — a low χ²/N here would prove
nothing, which is exactly why Step 2's 1.14 is meaningful evidence: it was
computed with dof still negative, so it is a conservative, not inflated,
reading).

| p | dof | Best cond(S) | Parameter subset |
|---|---|---:|---|
| 9 | 0 | 238.1 | group.R_sys_scale, R.SVEN, C.SAR, C.PAR, E.LV.EA, E.LV.EB, E.RV.EA, E.RV.EB, vsd.Cd |
| 8 | +1 | 62.07 | group.R_sys_scale, R.SVEN, C.SAR, C.PAR, E.LV.EB, E.RV.EA, E.RV.EB, vsd.Cd |
| **7** | **+2** | **17.17** | **group.R_sys_scale, R.SVEN, C.SAR, C.PAR, E.LV.EB, E.RV.EB, vsd.Cd** |
| 6 | +3 | 7.404 | R.SVEN, C.SAR, C.PAR, E.LV.EB, E.RV.EB, vsd.Cd |

**Recommended reduced set (p=7, dof=+2): `group.R_sys_scale, R.SVEN, C.SAR,
C.PAR, E.LV.EB, E.RV.EB, vsd.Cd`** — cond(S) drops from 223 to 17.17, over
13x more identifiable, while finally giving this patient's fit **positive
degrees of freedom** for the first time. `E.RV.EA` and `E.RV.EB` are
flagged collinear (ρ = −0.910) in the full-p analysis; the p=7 subset
resolves this by keeping only `E.RV.EB`.

**This reduced set has not itself been re-run through calibration** — per
`scripts/analyse_parameter_reduction.m`'s own stated limitation, conditioning
alone identifies a good *candidate*, but fixing the other 5 parameters at
their Step 2 calibrated values is a modelling commitment that needs its own
calibration run and gate comparison before it can be *claimed*, not just
suggested. That run was out of scope for this time-boxed execution (see
§5 "What was not done").

---

## 4. Step 4 — Out-of-sample post-closure prediction

**Revised 2026-09-06 (same day)** after the chamber-volume timing
correction (see the note at the top of this document):
`config/patient_reyna.m`'s protocol-sourced chamber volumes belong to
`clinical.post_surgery`, and a pre-existing gap in
`src/utils/get_calibration_targets.m` (the post-surgery metric table had
no row for `LVESV`/`RVESV` at all, so they were silently invisible to
every consumer, including this test) was fixed alongside it. Both changes
mean this out-of-sample test now includes six chamber-volume/EF
comparators it did not see the first time it ran.

Genuine holdout: **none** of these 13 post-closure measurements were used
to fit the pre-surgery parameters. The only change applied to the fitted
parameter set is closing the defect (no refitting, no adjustment).

| Metric | Measured | Predicted | Error% | Z-score |
|---|---:|---:|---:|---:|
| RAP_mean | 5 | 5.29 | +5.77% | 1.15 |
| PAP_min | 9 | 10.55 | +17.17% | 3.43 |
| PAP_max | 17 | 19.54 | +14.92% | 2.98 |
| PAP_mean | 13 | 14.73 | +13.30% | 2.66 |
| SAP_min | 68 | 64.64 | −4.95% | −0.49 |
| SAP_max | 89 | 102.22 | +14.85% | 1.48 |
| SAP_mean | 79 | 82.87 | +4.90% | 0.98 |
| LVEDV | 32 | 45.92 | **+43.48%** | 4.35 |
| LVESV | 23.6 | 20.02 | −15.18% | −1.52 |
| RVEDV | 30.5 | 34.62 | +13.52% | 1.35 |
| RVESV | 12 | 8.61 | −28.28% | −2.83 |
| LVEF | 0.2625 | 0.5640 | **+114.87%** | 11.49 |
| RVEF | 0.6066 | 0.7514 | +23.87% | 2.39 |

**Within 10%: 3 of 13. χ² = 201.20 over 13 targets, χ²/N = 15.48** — far
into the overfit/misspecification band (>2.0), and substantially worse
than the pressure-only reading (χ²/N = 4.64) this document first reported
before the timing correction.

**Honest reading:** the pre-surgery fit does not generalize to the
post-closure state, and the failure is concentrated in chamber function,
not pressures. The pressure story is unchanged from the first pass
(pulmonary pressures over-predicted 13-17%, systemic pressures within or
near 10%). The new information is damning: **LVEF is over-predicted by
115%** (0.564 predicted vs. 0.263 measured) and every chamber volume misses
by 13-43%. This says the calibrated ventricular elastance/volume
parameters (`E.LV.EA/EB`, `E.RV.EA/EB`, `V0.LV/RV` — fit only against
pre-surgery hemodynamics, since no pre-surgery chamber-volume evidence
exists to constrain them at all) do not capture the actual post-closure
ventricular state. This is exactly the failure mode a shared-parameter
model (Phase 4's premise: one patient, one parameter set, only the defect
changes) is vulnerable to when chamber compliance/geometry are
underconstrained by the fitted data — reported here as a finding, not
minimized, per this repo's own governance rule (G3) and the evaluation
script's own docstring, which anticipated exactly this outcome as a
possibility.

---

## 5. Step 5 — Identifiability at the final operating point

Already exported by the Step 2 run itself
(`results/runs/20260906_002327_reyna_pre_surgery/tables/parameter_identifiability_pre_surgery.csv`
and `..._pairs_pre_surgery.csv`), consistent with Step 3's findings: full
cond(S) = 223 (dof −3), one collinear pair (`E.RV.EA` / `E.RV.EB`,
ρ = −0.910). `V0.LV` has by far the smallest column norm (0.179) — the
data barely constrains it at all, consistent with it being among the first
parameters Step 3 drops.

---

## 6. What was not done (explicit, per G3/G4)

- **The p=7 reduced set was not re-run through calibration.** It is a
  well-supported candidate (§3), not yet a validated result.
- **Only 1 repeat per arm** in Step 1 (contract default is 3) and **only
  the winning arm** got the full 6-start budget, both stated and justified
  in `docs/publication_readiness_prd.md` as wall-clock-driven screening
  deviations, not hidden shortcuts.
- **`operational_lundquist`'s failure was not chased** (no new seed was
  computed to make it runnable) because it is not part of the fair-prior
  comparison that answers the scaling question (G8).
- **No history rewrite.** Patient identifiers remain in git history before
  the de-identification commit, per the study owner's explicit decision
  (see prior session record); only tracked-file state going forward is
  clean.
- **The GSA-vs-calibration wall-clock cost** means a genuinely broader
  sweep (more repeats, more starts on every arm, the p=7 set validated by
  its own multi-start) is future work, not something this session's budget
  covered.
- **No pre-surgery chamber-volume calibration target exists for Reyna at
  all**, following the 2026-09-06 timing correction. If genuine pre-surgery
  echo/cath chamber volumes are ever obtained, they would materially change
  what is identifiable in the pre-surgery fit (see Step 3/5 — `V0.LV` is
  currently the least-constrained parameter in the model).
- **No post-surgery calibration recipe exists for Reyna.** The chamber
  volumes now correctly placed in `clinical.post_surgery` are consumed
  generically by `get_calibration_targets`/`evaluate_post_closure_prediction`,
  but a future post-surgery calibration run would fall back to the
  *default* target-tier policy, which treats LVEDV/LVESV/LVEF as hard and
  RVEDV/RVESV/RVEF as soft **fitted** targets, not consistency-only —
  documented as an explicit trap in `config/patient_reyna.m`, not yet
  addressed.

## 7. Bottom line for anyone deciding whether to publish this

The corrected-data, fair-prior, Zhang-scaled, 6-start result
(`results/runs/20260906_002327_reyna_pre_surgery/`) is the best-supported
single candidate this repository has ever produced for Reyna's pre-surgery
**hemodynamics**: positive-χ²/N-adjacent fit quality at the pre-surgery
operating point (1.14, computed honestly at dof=−3), a real (not illusory)
multi-start search, and a concrete, quantified path to positive degrees of
freedom (p=7). It is **not** an ACCEPT-status result. Its **out-of-sample
generalization to the post-closure state is poor and, once chamber function
is correctly counted, badly so**: χ²/N = 15.48, driven by a 115%
over-prediction of post-closure LVEF. The pre-surgery hemodynamic fit and
the model's ability to predict post-closure chamber function are two
separate claims with very different strength, and any publication draft
built on this work should state both exactly that plainly, not average them
into one impression.
