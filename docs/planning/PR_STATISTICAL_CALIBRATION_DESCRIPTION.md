# Statistical calibration, and a correction to the clinical inputs

Implements PRD `reyna_statistical_calibration_v1` (σ-weighted objective, χ²
reporting, parameter identifiability, validation-holdout machinery) — and then,
in the course of validating it, finds and fixes three errors in the clinical
inputs that supersede the branch's own first set of results.

## Why this is worth reviewing carefully

The headline is not the σ-weighted objective. It is that **two model inputs
were wrong**, so every calibration number this project has published for this
patient was fitting a mis-specified child:

| Field | Was | Now | Why it matters |
|---|---:|---:|---|
| `HR` | 119 bpm | **136** | Sets cycle length (0.504 s → 0.441 s) and every stroke-volume derivation |
| `BSA` (+ weight/height) | 0.6173 (14.0 kg / 98 cm) | **0.588** (13.4 kg / 95 cm) | Drives the entire Zhang demographic scaling of the parameter prior |
| `SAP_mean` (pre) | 71.3 mmHg | **77** | Hard-tier *fitted* target |

Source: the study "reyna" procedure log, 06/04/2026 (full provenance kept
locally, not in this tracked file — see
`config/private/patient_provenance.local.m`) — the same catheterisation that
produced the pressures being fitted. The pre-surgery
pressures in the config all cross-checked correctly against it; these three
fields did not.

Two details that make the case rather than just asserting it:

- **`HR = 119` coincides exactly with the NIBP systolic** on the adjacent log
  line (`NIBP 119/83 (95)`), which is the likely origin of the error.
- **The demographics correction is a restoration, not a new decision.**
  `docs/clinical_data_dictionary.md` *already* documented 13.4/95/0.588 as the
  active anthropometry, and 14.0/98/0.6173 as "an alternate baseline-scaling
  experiment". The experiment's values had leaked into both active config
  files, making the dictionary's own claim false. The 2026-05-11 revision that
  introduced them post-dates the catheterisation by five weeks — the child had
  grown.

## New clinical data: `post_surgery` goes from empty to 7 targets

The same log carries post-closure pressures **from the same session** (device
released 12.06.23; readings 12.15–12.32, unchanged anaesthesia and ventilator
settings): `PAP` 17/9 (13), `SAP` 89/68 (79), `RAP` mean 5. Direction check
passes — PA pressure falls while RAP holds, consistent with removing the
left-to-right shunt.

This is a genuine **paired** pre/post dataset rather than two separate studies,
which changes the degrees-of-freedom outlook that
`reyna_zhang_scientific_assessment_20260828.md` §2.1 identified as the central
obstacle to publication:

| | Before | Now |
|---|---:|---:|
| Observations `N` | 9 | **16** |
| Parameters `p` | 12 | 12 |
| `dof = N − p` | **0** (`insufficient_dof`) | **4** |

## A silent-override trap, caught and guarded

`recipe.demographics` is merged **over** `clinical.common` by
`apply_calibration_recipe_to_clinical.m:24`. Correcting only
`config/patient_reyna.m` would have been silently reverted at runtime for
weight/height/BSA — a run that *looked* corrected but was not, with no warning.
Both files are corrected together, and
`test_reyna_systemic_flow_profile.m::Test2` now asserts the **effective
post-merge** values so the pair cannot drift apart silently again.

## Two findings withdrawn from earlier documents

- **The "60% critical internal stroke-volume inconsistency" was never
  internal.** The H+1 block is coherent: `SV_LV` = 21.7, `SV_RV` = 18.5
  (14.7% apart), `LVEF` self-consistent to 3 decimal places. The 60% came from
  comparing post-closure volumes against *pre*-closure flows — i.e. the shunt,
  the expected physiology. With the block excluded the audit returns severity
  `none`. `reyna_zhang_scientific_assessment_20260828.md` §2.3 and the PRD
  §7.2 are annotated accordingly.
- **The MAP form-factor question is answered by measurement, not a sensitivity
  arm.** The catheter stamps its own mean (`RFA 100/57 (77)`); the config was
  reconstructing it as `dia + (sys−dia)/3`, under-reading ~5–6 mmHg
  *systematically* (same offset post-closure: formula 75 vs stamped 79).

## Result on corrected data

Re-run with the identical seed and settings, so the data correction is the
only changed variable:

| | Before correction | **After** |
|---|---:|---:|
| Governed gate | 7 / 9 | **8 / 9** |
| All clinical targets | 9 / 11 | **10 / 11** |
| Within 5% excellent | 5 / 11 | **8 / 11** |
| Best RMSE (6 starts) | 0.0780 | **0.0480** |
| RMSE improvement | 70.8% | **83.5%** |
| Gate failures | `PAP_max`, `SAP_min` | **`PAP_min` only** |

The fit improving on corrected data is itself evidence the correction was
right. `SAP_min` moved from the worst failure (−13.85%) to a comfortable pass
(+4.27%) — a direct consequence of an `SAP_mean` target that had been
5.7 mmHg too low and a systemic waveform fitted at the wrong heart rate.

### Two χ² reporting defects found and fixed in review

Both in the statistic the publication claim rests on; the arithmetic was
correct, the reporting was not.

- **`dof` was clamped to 0.** At `N = 9`, `p = 12` the console printed
  `dof = N - p : 0` — arithmetically false, and it reads as *exactly
  determined* when the truth is **over-parameterised**, which is this
  project's central scientific criticism. Now reports `−3` with an explicit
  warning.
- **The `consistent` label ignored `dof`.** χ²/N = 0.784 was being labelled
  "residuals match measurement noise". With more parameters than
  observations, residuals that small are *guaranteed*. Now
  `consistent_but_underdetermined`, printing "do not quote chi2/N alone as
  validation". An `underfit` verdict is deliberately **not** softened —
  failing despite excess freedom is a real signal.

**The write-up refuses the flattering reading**: χ²/N sits inside the nominal
band and is reported as carrying no evidential weight. The gate count and
per-metric residuals are the defensible results.

### The second seed scored 9/9 — and that is the bad arm

Two independent 6-start selections on corrected data:

| | Seed `20260828` | Seed `20260830` |
|---|---:|---:|
| Governed gate | 8 / 9 | **9 / 9** |
| Best RMSE | 0.0480 | **0.0434** |
| χ²/N | 0.784 | **0.494** ← below overfit threshold |
| `Q_shunt_Lmin` (**not** graded) | −2.26% | **−22.90%** |

The arm with the perfect gate score is the one that crossed into overfitting
*and* saw the one deliberately-ungraded metric degrade ten-fold. The graded
set improved while the ungraded metric collapsed — the signature of fitting to
the scoreboard.

This is the empirical demonstration of
`reyna_zhang_scientific_assessment_20260828.md` §6's warning that maximising
the sub-10% count is not itself a scientific objective. It is no longer
hypothetical.

`Q_shunt_Lmin` is a sharp detector because it is `CO × (QpQs − 1)` and
`QpQs − 1 = 0.194` is a small difference of near-equal quantities — roughly
**6× amplification** of `QpQs` error. Excluded from grading for sound reasons;
worth reporting alongside the gate count *precisely because* it is not fitted.

**Run-to-run variation moves the headline claim** (9/9 vs 8/9), which is
itself a reportable result.

### Identifiability got worse, and is reported as such

Condition number **232 → 2.06 × 10³** (now flagged near-dependent). The LV
elastance collinearity recurs (ρ = −0.921), confirming it as structural rather
than an artefact of the bad data, and a new `E.LV.EA` ↔ `vsd.Cd` coupling
(−0.917) appears. A fit that improves while its parameters become less
separable is a warning, not a success — it strengthens the case that `p` must
come down.

## Phases delivered

| Phase | Status |
|---|---|
| 1 — σ-weighted objective | Complete; `legacy` default proven byte-identical by regression test |
| 2 — χ² reporting | Complete; recorded on status but deliberately **not** gating `ACCEPT` |
| 3 — parameter identifiability | Complete; governed-set report, found 3 collinear parameter pairs |
| 4 — joint pre/post inversion | **Objective + driver built and tested** (§ below); calibration run in progress |
| 5 — validation holdout | Machinery built; deviates from PRD with reasoning (see below) |

**Phase 5 deviates deliberately.** The PRD said to relabel `SVR` as
`validation_holdout` on the premise it is "a genuine prediction test". It is
not — `SVR = (SAP_mean − RAP_mean)/CO`, all three of which are already fitted
targets. Relabelling it would have manufactured the appearance of an
independent prediction success. `SVR` stays `derived_validation`, the holdout
machinery is built generically for a future genuine candidate, and a test locks
in the honest current state: **nothing in this recipe is presently a real
holdout.**

## Evidence is now reviewable

`.gitignore` un-ignored `full_metric_gate_*.csv` so gate-count claims stayed
verifiable. This PR extends the same reasoning to the new statistical
artefacts (`chi_squared_*.csv`, `parameter_identifiability_*.csv`) — a χ²/N or
condition-number claim is no more verifiable than a gate count if its table is
untracked.

## The headline: a genuine validation holdout, and overfitting confirmed

The post-closure pressures were **never used in fitting**. Taking a
pre-calibrated parameter set, closing the defect, and comparing is therefore
the first true out-of-sample test this model has had — no refitting, closure
is the only intervention.

| Metric | Measured | Predicted | Error |
|---|---:|---:|---:|
| `SAP_mean` | 79 | 81.7 | **+3.5%** |
| `PAP_mean` | 13 | 13.8 | +6.4% |
| `RAP_mean` | 5 | 5.36 | +7.3% |
| `SAP_min` | 68 | 63.3 | −6.9% |
| `PAP_max` | 17 | 17.7 | +4.4% |
| `SAP_max` | 89 | 101.4 | +13.9% |
| `PAP_min` | 9 | 10.3 | +14.8% |

**5 of 7 within 10%**, χ²/N = 2.31. This number needs no degrees-of-freedom
caveat — the model cannot have absorbed targets it never saw.

**And it settles the overfitting question with held-out data:**

| | Seed `20260828` | Seed `20260830` |
|---|---:|---:|
| In-sample gate | 8 / 9 | **9 / 9** ← looks better |
| **Out-of-sample** | **5 / 7** ← actually better | 3 / 7 |

The arm that fit the training data better predicts unseen data **worse**.
Overfitting demonstrated, not inferred — and the clearest possible argument
for why a 9/9 headline would have been the wrong thing to publish.

Errors are systematic (6 of 7 positive, pulmonary pressures over-predicted
throughout): the model predicts **less pulmonary unloading after closure than
actually occurred** — an interpretable lead, not just a residual.

## Phase 4: degrees of freedom are now positive

`src/calibration/objective_joint_pre_post.m` + `scripts/run_joint_pre_post_calibration.m`
fit both haemodynamic states from **one shared parameter vector**, so the
observation count rises without adding parameters:

| | `N` | `p` | `dof` |
|---|---:|---:|---:|
| Pre-only (previous state) | 9 | 12 | **−3** |
| Joint, full recipe set | 16 | 14 | **+2** |
| Joint, GSA-masked set | 16 | 12 | **+4** |

Pinned by test: the vector is genuinely shared (every parameter but the shunt
bit-identical across both structs), regularisation is applied **once** not per
scenario, σ resolution matches Phase 1 so `χ²_pre` stays comparable with the
reported statistic, and absent post targets degrade to pre-only *with a
warning that the DOF benefit no longer applies*.

### Two latent bugs this surfaced

**1. "Closed VSD" was mode-dependent — and wrong for this patient.**
`vsd_shunt_model` dispatches on `vsd.mode`. Resistive modes close via a large
`R.vsd`; but `orifice_bidirectional` — **the mode Reyna uses** — never reads
`R.vsd` at all and closes only when `vsd.area_mm2 = 0`. Closing by `R.vsd`
alone is a **no-op** here: the "post-closure" simulation would keep shunting
at full strength. The test asserts closure *behaviourally* (zero flow at a
70 mmHg gradient), because a field-based assertion would have passed while the
physics was wrong.

> This reaches beyond Phase 4: `main_run.m:839` builds the pre-to-post seed
> with `R.vsd = 1e6` and nothing else, so **that seed is not a closed-VSD
> model** for an orifice-mode patient. Left unfixed (outside this PRD's
> scope) but flagged — it must be fixed before any post-closure result is
> published from that seed.

**2. The two-code-path tier disagreement, hit in practice.** A bare
`build_target_tiers(clinical, scenario)` ignores
`recipe.primary_rmse_holdout` and governs **10** pre-surgery rows where the
production path governs **9** — silently readmitting `Q_shunt_Lmin`, the very
metric that acts as the overfitting detector above. `χ²_pre` would then have
been computed over a different set than the reported governed RMSE. Tiers now
come from each scenario's case profile, and a test pins `n_pre = 9`.

## Parameter reduction: analysed, not yet run

`dof = −3` is the binding limitation. GSA screens the set to 7, but stages D–F
expand it back to 12, undoing most of that.

Scoring subsets by `cond(S)` at the calibrated point (no recalibration needed —
fixed columns, so it is linear algebra):

| p | dof | best cond(S) |
|---:|---:|---:|
| 12 (current) | −3 | 2060 |
| **7** | **+2** | **21.7** |
| 6 | +3 | 6.82 |

Cutting 12 → 7 makes `dof` positive and improves conditioning **95×**. The
resulting set shares **6 of 7 members** with the Sobol-screened set — two
independent criteria converging, which is the best available evidence the
subset reflects the data rather than the method.

**Not validated.** A dropped parameter is fixed at its calibrated value, which
is a modelling commitment. A run at `p = 7` is the highest-value next step.

## What this does NOT establish

- **n = 1.** One patient.
- **The fit is underdetermined — this is the binding limitation.** `N = 9`
  against `p = 12` gives `dof = −3`. Any claim must rest on the gate count and
  per-metric residuals, and must state the parameter/observation ratio
  alongside. Phase 4 would give `N = 16`, `dof = 4` — an improvement, not a
  resolution; reducing `p` below 12 is likely needed too.
- **Zhang vs Lundquist stays uninterpretable** at these budgets.
- **No post-closure flow.** The log records `PARI 1.9` / `FR 1.19` whose
  meaning is unconfirmed; `FR 1.19` is near-identical to the pre-op `QpQs`
  1.194 while post-dating device release, so it was **not** entered rather
  than guessed.
- **`CO_Lmin` = 3.423 unchanged**, derived from a protocol `Qp` that may carry
  a BSA-indexed Fick derivation the BSA correction should propagate into. Not
  altered without confirming.

## Testing

- 85/85 `functiontests`-style tests pass.
- `test_reyna_systemic_flow_profile.m`: 6 passed / 3 failed — **identical to
  the pre-change baseline**, verified by stashing the edits and re-running.
  Those 3 failures are pre-existing on this branch and unrelated to this work.
- Note the suite is mixed-style: script-style tests must be run via
  `run('tests/x.m')`, not `runtests`, which splits them and destroys their
  shared counter scope.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
