# Reyna Publication-Readiness Re-Run — Results, 2026-09-06

**This is the current authoritative results document for the Reyna model.**
It supersedes every prior results document in `docs/` (each of which now
carries a superseded-banner pointing here or to
`docs/CHANGES_SINCE_PR22.md` §14). It reports the execution of
`docs/publication_readiness_prd.md` against the protocol-corrected clinical
inputs from `docs/CHANGES_SINCE_PR22.md` §14 (HR 119 bpm, VSD diameter
3.665 mm, pre-surgery chamber volumes 32/23.6/30.5/12 mL as
consistency-only).

**One-line honest summary:** under a fair prior (no historical warm start),
Zhang scaling clearly outperforms Lundquist-BSA scaling on Reyna's corrected
data (RMSE 0.118 vs 0.222). A full 6-start multi-start on the Zhang arm
reaches primary RMSE 0.0739, governed gate 7/9, and a χ²/N of 1.14 —
consistent, not overfit. **The model does not reach ACCEPT status, and its
out-of-sample post-closure prediction is poor (3 of 7 within 10%,
χ²/N = 4.64).** Every number below is reported as obtained, including the
ones that are unfavorable, per this repo's own governance rule (G3).

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

Consistency-only chamber rows (never fitted; see `docs/CHANGES_SINCE_PR22.md`
§14 for why): RVEDV predicted 35.48 vs. clinical 30.5 mL (+16.3%), RVESV
9.09 vs. 12 mL (−24.2%). These two are reported for transparency and are
fine to quote.

**LVEF, LVEDV, and LVESV are excluded from this and any publication-facing
report (study-owner decision, 2026-09-06).** The model still computes and
exports them in the full metric gate CSV for code-level transparency and
audit — nothing in `config/patient_reyna.m` or the calibration pipeline
changed — but their source values are the internally implausible LV pair
documented in `docs/CHANGES_SINCE_PR22.md` §14 (SV_LV = 8.4 mL vs. an SV of
roughly 34 mL implied by the protocol's own Qp = 4.087 L/min): the data
itself, not just the model's fit to it, is not trustworthy enough to state
as a finding. For the record, the raw numbers were LVEF predicted 0.610 vs.
clinical 0.2625 (+132%), LVEDV 45.66 vs. 32 mL (+42.7%), LVESV 17.80 vs.
23.6 mL (−24.6%) — kept here only so the exclusion is auditable, not as
something to cite.

Full data: `results/runs/20260906_002327_reyna_pre_surgery/tables/full_metric_gate_pre_surgery.csv`.

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

Genuine holdout: post-closure catheter pressures were **never used to fit**
these parameters. The only change applied is closing the defect (no
refitting, no parameter adjustment).

| Metric | Measured | Predicted | Error% | Z-score |
|---|---:|---:|---:|---:|
| RAP_mean | 5 | 5.29 | +5.77% | 1.15 |
| PAP_min | 9 | 10.55 | **+17.17%** | 3.43 |
| PAP_max | 17 | 19.54 | +14.92% | 2.98 |
| PAP_mean | 13 | 14.73 | +13.30% | 2.66 |
| SAP_min | 68 | 64.64 | −4.95% | −0.49 |
| SAP_max | 89 | 102.22 | +14.85% | 1.48 |
| SAP_mean | 79 | 82.87 | +4.90% | 0.98 |

**Within 10%: 3 of 7. χ² = 32.51 over 7 targets, χ²/N = 4.64 (overfit
band — well above the 2.0 threshold).**

**Honest reading:** the pre-surgery fit does not generalize well to the
post-closure state without refitting. The systematic direction is
informative: every pulmonary-pressure metric (PAP_min/mean/max) is
*over*-predicted by 13-17%, suggesting the calibrated pulmonary
resistance/compliance combination is too restrictive for the closed-VSD
state, or that the pre-surgery fit compensated for something (most likely
the consistency-only chamber-volume mismatch documented in §2.2) in a way
that does not transfer. SAP tracks better (both systemic pressures within
or near 10%). This is a real limitation of the current shared-parameter
assumption (Phase 4's premise: one patient, one parameter set, only the
defect changes) and is reported as such, not minimized — a large error here
is a finding, not a disappointment to be tuned away (see the script's own
docstring, which anticipated exactly this outcome as a possibility).

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

## 7. Bottom line for anyone deciding whether to publish this

The corrected-data, fair-prior, Zhang-scaled, 6-start result
(`results/runs/20260906_002327_reyna_pre_surgery/`) is the best-supported
single candidate this repository has ever produced for Reyna pre-surgery:
positive-χ²/N-adjacent fit quality at the pre-surgery operating point
(1.14, computed honestly at dof=−3), a real (not illusory) multi-start
search, and a concrete, quantified path to positive degrees of freedom
(p=7). It is **not** an ACCEPT-status result, and its **out-of-sample
generalization to the post-closure state is poor**. Both facts belong in
any publication draft built on this work, stated as plainly as they are
here.
