# Reyna Statistical Calibration — Results, 2026-08-29

> ## ⚠ Read §0 before quoting any number in this document
>
> On 2026-08-29 the source catheterisation record (RSAB Harapan Kita
> PROCEDURE LOG, MRN 01008971, 06/04/2026) was obtained and revealed three
> errors in the clinical inputs — including **HR and BSA, which are model
> inputs, not just report labels**.
>
> **The authoritative result is §3.5** (corrected data). §3.2 and §7's
> original numbers were produced against the uncorrected inputs and are
> **superseded** — they are retained because the process findings around them
> matter, not because the numbers do.
>
> **The one-line honest summary:** on corrected data the model reproduces
> **8 of 9 governed clinical targets within 10%** (10 of 11 across all
> targets), best primary RMSE **0.0480** across 6 starts. χ²/N = 0.784 sits
> inside the nominal "consistent" band but **carries no evidential weight**,
> because at `p = 12` against `N = 9` the fit has `dof = −3` — more free
> parameters than observations. The gate count and per-metric table are
> defensible; χ²/N is not. See §3.5 and §4.

## 0. Clinical data correction, 2026-08-29 (supersedes all results below)

The procedure log for the catheterisation that produced this patient's
haemodynamics was retrieved after the runs in §3 completed. Cross-checking
the config against it confirmed the pre-surgery *pressures* were all correct
(`PA 20/10 (15)` ✓, `RFA 100/57` ✓, `RA` mean 5 ✓) but exposed three errors:

| Field | Was | Corrected to | Source row | Why it matters |
|---|---:|---:|---|---|
| `HR` | 119 bpm | **136 bpm** | log 09.24.57 `Nadi 136 bpm` | **Model input.** Sets cycle length: 60/119 = 0.504 s → 60/136 = 0.441 s. Also changes every stroke-volume derivation (`SV_Qs` 28.8 → 25.2 mL/beat). |
| `weight` / `height` / `BSA` | 14.0 kg / 98.0 cm / 0.6173 | **13.4 kg / 95.0 cm / 0.588** | log 07.53.01–07.53.20 | **Model input.** BSA drives the entire Zhang demographic scaling of the parameter prior. |
| `SAP_mean` (pre) | 71.3 mmHg | **77 mmHg** | log 10.39.12 `RFA 100/57 (77)` | Hard-tier **fitted target**. |

**On HR.** The prior value of 119 coincides exactly with the NIBP *systolic*
on the adjacent log line (`NIBP 119/83 (95)`), which is the likely origin of
the error.

**On demographics.** The prior values were labelled *"Keisya 2026-05-11
revision"* — dated five weeks **after** this 06/04/2026 catheterisation. The
child had grown; pairing May anthropometry with April haemodynamics scaled
the model to a larger patient than the one measured. The measurement-day
values are the correct ones. (The stamped BSA 0.588 also uses DuBois, where
the old config value used Mosteller; the stamped figure is kept so the model
matches the source record exactly.)

**On MAP.** Three candidates existed: NIBP cuff 95 (different method),
form-factor reconstruction 71.3 (`dia + (sys−dia)/3`), and the transducer's
own stamped mean 77. The catheter's direct reading is now used. The
form-factor assumption was under-reading by ~5–6 mmHg *systematically* — the
same offset recurs post-closure (formula 75 vs stamped 79), so it was a
method error, not noise. This resolves the MAP form-factor question the
2026-08-28 assessment raised (§5, step 6) with a measurement rather than a
sensitivity arm.

### 0.0 The correction restores documented intent — it does not change it

`docs/clinical_data_dictionary.md` already specified the corrected values as
the active ones, stating that *"the active Reyna patient file follows the raw
protocol anthropometry"* and listing `13.4 kg / 95.0 cm / BSA 0.588`. It
explicitly described `14.0 / 98.0 / 0.6173` as *"an alternate
baseline-scaling experiment"*.

The experiment's values had nonetheless leaked into both active config files,
making the dictionary's own claim false and silently scaling every Reyna run
to a larger child than the one measured. **The correction below is therefore
a restoration of the documented design, not a new decision** — which is the
strongest possible footing for it. The drift is now recorded in the
dictionary itself, and the new assertion in
`test_reyna_systemic_flow_profile.m::Test2` prevents a silent recurrence.

### 0.1 A silent-override trap that was caught

`recipe.demographics` in `config/calibration_recipes/reyna_pre_surgery.m` is
merged **over** `clinical.common` by
`apply_calibration_recipe_to_clinical.m:24`. Correcting only
`config/patient_reyna.m` would therefore have been silently reverted at
runtime for weight/height/BSA (HR would have survived, as it is not in that
struct), producing a run that looked corrected but was not. Both files were
updated together, and
`tests/test_reyna_systemic_flow_profile.m::Test2` — which asserts the
*effective* post-merge values — was updated to match and now guards the
pairing.

### 0.2 Post-surgery data now exists

The same log carries a full set of **post-closure** pressures (device placed
11.50.19, released 12.06.23; all readings below stamped 12.15–12.32, same
anaesthesia and ventilator settings). `clinical.post_surgery` went from
entirely `NaN` to **7 finite targets**:

| Target | Value | Source rows |
|---|---:|---|
| `PAP_max` / `PAP_min` / `PAP_mean` | 17 / 9 / 13 | `PA 17/8 (13)`, `17/9 (13)`, `17/9 (13)` @12.26–12.27 |
| `SAP_max` / `SAP_min` / `SAP_mean` | 89 / 68 / 79 | `RFA 91/68 (79)`, `89/68 (78)`, `89/68 (79)` @12.31–12.32 |
| `RAP_mean` | 5 | `RA 8/5 (5)`, `8/5 (5)`, `7/5 (5)` @12.28 |

Direction check passes: PA pressure falls (mean 15 → 13) while `RAP` holds at
5 — consistent with removal of the left-to-right shunt.

**This is a genuine paired pre/post dataset from a single session**, not two
separate studies, which makes the two states directly comparable and changes
the outlook for Phase 4 substantially (§6).

### 0.3 What is still open

- **Post-closure flow is not established.** The log notes `no CO`
  pre-closure, and records `Hasil PARI 1.9` / `Hasil FR 1.19` at 12.10–12.11.
  If `FR` is a flow ratio, 1.19 is near-identical to the pre-operative
  `QpQs` 1.194 and is timestamped *after* device release, which would be
  unexpected for a complete closure. **Not entered pending clarification of
  what PARI and FR denote.** `post.QpQs`, `post.CO_Lmin`, `post.PVR_WU` and
  `post.SVR_WU` remain `NaN`.
- **`CO_Lmin` = 3.423 is unchanged**, derived as `Qp/QpQs` = 4.087/1.194 from
  the *protocol* document, not this log. If `Qp` was Fick-derived using a
  BSA-indexed VO₂, the BSA correction above would propagate into it. Not
  altered without confirming that derivation.
- **The H+1 echo volumes were not moved into `post_surgery`.** They are a
  different timepoint from these in-lab pressures (ward echo vs catheter
  table under anaesthesia); combining them into one "post" state would repeat
  the timing-mismatch error this correction exists to fix.

### 0.4 Verification status of the correction

- `patient_reyna()` and the recipe parse; effective post-merge values confirmed
  as `HR=136 BSA=0.588 W=13.4 H=95 SAP_mean=77`.
- The 10 `functiontests`-style suites: **85/85 pass**.
- `test_reyna_systemic_flow_profile.m`: **6 passed / 3 failed**, identical to
  the pre-change baseline measured by stashing the edits — the 3 failures
  (PAP waveform tiering, H+1 volume governance, recipe seed/manifest drift)
  are pre-existing on this branch and unrelated to this correction.


Branch: `codex/reyna-statistical-calibration` (off `codex/reyna-zhang-fullmetric-10pct`, PR #24)
PRD: [reyna_statistical_calibration_prd.md](reyna_statistical_calibration_prd.md)
Scope: Reyna `pre_surgery`, Zhang scaling prior, fair-prior (historical seeds disabled)

---

## 1. What this branch changes

Four of the PRD's five phases are complete. Phase 4 (joint pre/post inversion)
is deliberately blocked — see §6.

| Phase | Status | Summary |
|---|---|---|
| 1 — σ-weighted objective | **Complete; re-run on corrected data** | Residuals normalised by declared measurement uncertainty instead of one global percentage; opt-in via `objectiveWeighting`, default remains `legacy`. Authoritative result in **§3.5**: governed gate **8/9**, best RMSE **0.0480** across 6 starts. (§3.2's 7/9 is superseded — wrong clinical inputs.) |
| 2 — χ² reporting | **Complete** | Discrepancy-principle goodness-of-fit statistic, printed and exported every run; does not gate `ACCEPT` |
| 3 — parameter identifiability | **Complete** | Scaled sensitivity matrix, condition number, pairwise correlation; report-only |
| 4 — joint pre/post inversion | **Governance resolved, awaiting data** | The three governance questions are answered (§6); `clinical.post_surgery` is currently all-`NaN` and the study owner is retrieving the post-operative record. Value depends on which rows it yields — see §6.0 |
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

By the proper statistical measure, the PR #24 candidate had
`χ² = 13.22`, `N = 9`, `p = 7`, `χ²/N = 1.47` — inside the nominal
`consistent` band (0.5–2.0).

> **Correction, 2026-08-29.** The original wording here read that this meant
> "the residuals are broadly consistent with declared measurement noise, not
> badly wrong". That is precisely the inference §4 now forbids, and it is
> withdrawn. Even at that candidate's `p = 7`, `dof = 2` — inside the
> `insufficient_dof` band, where the reduced statistic is not stable. Being
> inside the χ²/N band is a *necessary* condition for a good fit, not
> evidence of one, and this document should not have implied otherwise while
> arguing for more statistical rigour. The point §2 actually establishes
> stands unaffected: the **percentage gate ranks the wrong metric worst**,
> which is what motivated Phase 1.
>
> Note also that these PR #24 figures were computed against the uncorrected
> clinical inputs (§0), so the specific σ values and errors in the table above
> would differ if recomputed today.

## 3. Phase 1 — σ-weighted objective

> **⚠ Every numeric result in §3 is SUPERSEDED by the §0 data correction.**
> These runs used `HR = 119` (should be 136) and `BSA = 0.6173` (should be
> 0.588) — both model inputs, so the fits are to a mis-specified patient. The
> *code*, the *method*, and the *process findings* in this section stand and
> are not affected; the numbers must be regenerated. A re-run on corrected
> data with the identical seed (`20260828`) is what §3.5 will report.

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

**Re-run complete, 2026-08-29.** All Phases 1–5 code is complete, tested
(85/85), and committed on this branch. The clean 6-start sigma-weighted run
was executed successfully on this session's second attempt (the first attempt
was independently killed by an unrelated tooling issue — an externally
terminated process, not a code or data bug — partway through start 3/6; no
partial artefacts existed to salvage, so it was simply relaunched from
scratch). Total wall time: 9743 s (~2.7 h), run folder
`results/runs/20260829_161536_reyna_pre_surgery`.

Per-start results (`UNIFIED_VSD_MULTISTART_SEED=20260828`):

| Start | Label | RMSE | primary_fail | gate_fail |
|---|---|---:|---:|---:|
| 1/6 | seed | 0.0780288 | 0 | 2 |
| 2/6 | sobol_1 | 0.104133 | 0 | 2 |
| 3/6 | sobol_2 | 0.114366 | 0 | 3 |
| 4/6 | sobol_3 | 0.0897911 | 0 | 2 |
| 5/6 | sobol_4 | 0.104133 | 0 | 2 |
| 6/6 | sobol_5 | 0.103836 | 1 | 2 |

RMSE across starts: min 0.0780288, median 0.103984, max 0.114366,
IQR 0.0143417 — this matches the unverified console figures glimpsed before
the earlier crashed attempt (§3.2 above), confirming Phase 1's objective code
was not itself in question. The winning candidate is start 1/6 (seed).

| | legacy (PR #24, cited) | sigma (this run) |
|---|---:|---:|
| Governed gate | 8 / 9 | **7 / 9** |
| Primary RMSE (best of 6 starts) | 0.0873 | **0.0780** |
| RMSE spread (min–max, 6 starts) | 0.0873 – 0.1041 | **0.0780 – 0.1144** |
| χ²/N | 1.47 | **1.60** |
| Interpretation | consistent | **consistent** |
| Active parameters (p) | 7 | **12** |
| Condition number | — | **232** |

`governed_gate_failures = PAP_max, SAP_min`. `p` came back as 12 (a genuine
positive integer, not 0), confirming the §3.2 contamination bug did not
recur on this run.

### Per-metric table (governed set, sorted worst-first by \|z\|)

| Metric | Tier | Clinical | Calibrated | Error % | σ | z | z² | 10% gate |
|---|---|---:|---:|---:|---:|---:|---:|---|
| PAP_max | soft | 20 | 22.609 | +13.05 | 1.00 | +2.609 | 6.808 | FAIL |
| PAP_min | soft | 10 | 10.857 | +8.57 | 0.50 | +1.714 | 2.939 | PASS |
| PAP_mean | hard | 15 | 16.053 | +7.02 | 0.75 | +1.404 | 1.971 | PASS |
| SAP_min | soft | 57 | 49.104 | -13.85 | 5.70 | -1.385 | 1.919 | FAIL |
| RAP_mean | hard | 5 | 5.129 | +2.59 | 0.25 | +0.518 | 0.268 | PASS |
| SAP_max | soft | 100 | 95.601 | -4.40 | 10.00 | -0.440 | 0.194 | PASS |
| CO_Lmin | hard | 3.423 | 3.223 | -5.85 | 0.50 | -0.401 | 0.161 | PASS |
| SAP_mean | hard | 71.3 | 70.157 | -1.60 | 3.565 | -0.321 | 0.103 | PASS |
| QpQs | hard | 1.194 | 1.188 | -0.47 | 0.0597 | -0.094 | 0.009 | PASS |

Worst by |z|: `PAP_max` (z = +2.61), same metric flagged worst by the
CHI-SQUARED block. Under the sigma-normalised view, `PAP_max` and `PAP_min`
switch places relative to the legacy percentage view — `PAP_max` is now the
single largest statistical outlier despite passing the old percentage read
comfortably at earlier stages, which is the exact effect §2 predicted a
sigma-aware view would surface.

(`Q_shunt_Lmin` and `SVR` excluded from this table — both are
`primary_rmse_holdout`/`derived_validation`, not governed primary metrics;
see §5.)

### 3.3 Decision

The sigma-weighted arm trades one governed-gate pass for a better primary
RMSE (0.0780 vs 0.0873) and a materially tighter best-start result, at the
cost of `PAP_max` flipping from PASS to FAIL under the 10% gate (barely:
13.05% vs the legacy candidate's presumed pass) while `SAP_min` remains the
same persistent failure in both arms. χ²/N moved from 1.47 to 1.60 — both
values sit inside the `[0.5, 2.0]` "consistent" band, so this is not a
regression by the statistical criterion the whole PRD exists to introduce;
it is the same overall fit quality, redistributed across metrics differently
because the objective now weights by declared measurement uncertainty
instead of a flat percentage.

This is not an unambiguous win. Recommendation: do **not** promote the
sigma-weighted result to the recipe's `initial_parameter_values` seed
(§10.4) on this evidence alone — a single n=1 A/B comparison with one flipped
gate metric is not a strong enough signal, and `objectiveWeighting` stays
opt-in (`legacy` default) pending a second run or a decision from the study
owners on which regime should be authoritative. What this run does establish
is that the sigma-weighted code path is real, executes end-to-end, and
produces defensible, differently-weighted results — not that it is strictly
better.

### 3.5 Corrected-data re-run (authoritative result)

Launched 2026-08-29 after the §0 correction, with **the identical seed
(`20260828`) and identical settings** as the superseded run above, so that the
clinical data correction is the *only* changed variable and its effect is
cleanly attributable.

Log: `results/runs/_logs/sigma_correcteddata_20260828seed.log`
Run folder: `results/runs/20260829_225931_reyna_pre_surgery`
Wall time: 16 265 s (~4.5 h).

**The corrected data fits substantially better on every measure.** This is
itself evidence that the correction was right: fitting the patient who was
actually measured produces a better fit than fitting a mis-specified one.

| | Superseded (§3.2) | **Corrected (authoritative)** |
|---|---:|---:|
| Governed gate | 7 / 9 | **8 / 9** |
| All clinical targets | 9 / 11 | **10 / 11** |
| Within 5% excellent band | 5 / 11 | **8 / 11** |
| Best RMSE (of 6 starts) | 0.0780 | **0.0480** |
| RMSE spread (min–max) | 0.0780 – 0.1144 | **0.0480 – 0.1053** |
| RMSE improvement vs baseline | 70.8% | **83.5%** |
| χ² | 14.371 | **7.059** |
| χ²/N | 1.597 | **0.784** |
| Gate failures | `PAP_max`, `SAP_min` | **`PAP_min` only** |
| Worst by \|z\| | `PAP_max` (+2.61) | **`PAP_min` (+2.14)** |

Per-start (seed `20260828`, identical to the superseded run):

| Start | Label | RMSE | primary_fail | gate_fail |
|---|---|---:|---:|---:|
| 1/6 | seed | 0.0905017 | 0 | 3 |
| 2/6 | sobol_1 | 0.0778811 | 0 | 2 |
| **3/6** | **sobol_2** | **0.0480448** | **0** | **1** ← winner |
| 4/6 | sobol_3 | 0.0778811 | 0 | 2 |
| 5/6 | sobol_4 | 0.0778811 | 0 | 2 |
| 6/6 | sobol_5 | 0.105326 | 1 | 2 |

min 0.0480448, median 0.0778811, max 0.105326, IQR 0.0126206.

### Per-metric, governed set, worst-first by \|z\|

| Metric | Tier | Clinical | Model | Error % | σ | z | 10% gate |
|---|---|---:|---:|---:|---:|---:|---|
| `PAP_min` | soft | 10 | 11.07 | +10.72 | 0.50 | **+2.14** | **FAIL** |
| `RAP_mean` | hard | 5 | 5.26 | +5.20 | 0.25 | +1.04 | PASS |
| `PAP_max` | soft | 20 | 19.00 | −5.00 | 1.00 | −1.00 | PASS |
| `SAP_min` | soft | 57 | 59.43 | +4.27 | 5.70 | +0.43 | PASS |
| `SAP_max` | soft | 100 | 96.05 | −3.95 | 10.00 | −0.39 | PASS |
| `CO_Lmin` | hard | 3.423 | 3.334 | −2.59 | 0.50 | −0.18 | PASS |
| `PAP_mean` | hard | 15 | 14.94 | −0.40 | 0.75 | −0.08 | PASS |
| `SAP_mean` | hard | 77 | 76.85 | −0.19 | 3.85 | −0.04 | PASS |
| `QpQs` | — | 1.194 | 1.1946 | +0.05 | 0.0597 | +0.01 | PASS |

Notably **`SAP_min` moved from the worst failure (−13.85%) to a comfortable
pass (+4.27%)**, and `SAP_mean` now fits to 0.19%. Both are direct consequences
of the §0 corrections — `SAP_min` was being pulled by an `SAP_mean` target that
was 5.7 mmHg too low, and the whole systemic waveform was being fitted at the
wrong heart rate.

### χ² after the §4 reporting fix

```
  N (governed observations) : 9
  p (active parameters)     : 12
  dof = N - p               : -3  [dof <= 0: MORE FREE PARAMETERS THAN
                                   OBSERVATIONS -- a low chi2/N is guaranteed
                                   here and is not evidence of fit quality]
  chi2 / N                  : 0.784
  interpretation            : consistent_but_underdetermined
  [QUALIFIED] ... Do not quote chi2/N alone as validation.
```

**This is the honest reading and it must not be softened in the write-up.**
χ²/N = 0.784 sits inside the nominal `consistent` band, but with `dof = −3`
that band carries no evidential weight: the model has three more free
parameters than observations, so residuals this small are expected whether or
not the model is correct. The gate count (8/9) and the per-metric table above
are the defensible results; χ²/N is not.

### Identifiability at the corrected operating point

Condition number **2.06 × 10³** — flagged `near-dependence: cond > 1e3`, and
notably *worse* than the superseded run's 232. Two collinear pairs:

- `E.LV.EA` ↔ `E.LV.EB` (ρ = −0.921)
- `E.LV.EA` ↔ `vsd.Cd` (ρ = −0.917)

The LV elastance pair reappears, confirming it as structural rather than an
artefact of the wrong data. The new `E.LV.EA` ↔ `vsd.Cd` coupling is
consistent with a better-fitting shunt: as the fit improves, LV contractility
and orifice discharge trade off more sharply against each other. **This
strengthens rather than weakens the case that `p` must come down** — the
better fit is being bought partly with parameter redundancy.

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

> **Two defects found in review, 2026-08-29 — fix queued, not yet applied.**
> Both were found while preparing this branch for publication and both affect
> how the headline statistic reads to a reviewer. Neither could be fixed
> immediately because `compute_chi_squared_report.m` is on the calibration
> path and a run was in flight (§10.1); they are to be applied as soon as it
> completes, followed by a regeneration of the χ² numbers.
>
> **4.D1 — reported `dof` is clamped to 0 and prints a false statement.**
> `compute_chi_squared_report.m:83` computes
> `dof = max(n_obs - n_parameters, 0)`. On the actual governed set
> (`N = 9`, `p = 12`) the true value is **−3**, but the console prints
> `dof = N - p : 0`, which is arithmetically wrong as written. Worse, it
> conceals the finding that matters most: `dof = 0` reads as *exactly
> determined*, whereas `dof = −3` reads as **over-parameterised** — which is
> precisely the central criticism in
> `reyna_zhang_scientific_assessment_20260828.md` §2.1. The clamp suppresses
> the signal the reader most needs. Fix: keep the true (possibly negative)
> value for reporting; retain the clamp only where a non-negative divisor is
> required.
>
> **4.D2 — the `consistent` label does not account for `dof`.**
> `classify_chi2_per_obs` assigns `underfit` / `consistent` / `overfit` from
> `χ²/N` alone. When `p > N` the model has more freedom than data, so small
> residuals are guaranteed rather than earned; labelling `χ²/N = 1.60`
> "consistent — residuals match measurement noise" therefore overstates the
> evidence. The band is only meaningful with positive `dof`. Fix: qualify the
> label when `dof <= 0` so the statistic cannot be quoted as validation of a
> model that is over-parameterised.
>
> Note both defects are *reporting* faults, not errors in the χ² arithmetic
> itself: `chi2` and `chi2_per_obs` are computed correctly.


Every run now exports:

- `full_metric_gate_<scenario>.csv` with three new columns: `Sigma`, `ZScore`,
  `ZScoreSquared`. The table is now sorted worst-first by `|z|` rather than by
  percentage error, with `AbsError_pct` retained as a secondary column.
- `chi_squared_<scenario>.csv`: a one-row summary (`chi2`, `N`, `p`, `dof`,
  `chi2/N`, `chi2/dof`, `interpretation`, `dof_note`, `worst_metric`,
  `worst_z`).
- A `--- CHI-SQUARED ---` console block after the metric gate block.

**Note on `p`**: the Sobol GSA screen at Step 5 selects a 7-parameter subset
for the core mechanistic optimization (Stages A–C). Stages D–F (systemic
polish, plausibility polish, validation-gate polish) are deliberately built
on a broader active set — `calib.names_all(calib.mask)`, 12 parameters,
including cardiac elastances and unstressed volumes the GSA screen never
selected (`run_validation_gate_polish.m:145`). `p` in the chi-squared/
identifiability report is this 12-parameter set, since that is what actually
had freedom to move in producing the winning candidate (Stage F won in the
2026-08-29 run — see §3.2). This is correct and intentional, not a
double-counting bug; using the narrower 7 would understate real DOF.

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

## 6. Phase 4 — governance resolved, awaiting post-operative data

**Status as of 2026-08-29: the three governance questions are answered; Phase
4 is now blocked on data entry rather than on a decision.**

The PRD framed Phase 4 as a *governance* decision: relocate the H+1
post-operative chamber-volume block from `clinical.pre_surgery` to
`clinical.post_surgery`, then run a joint pre/post inversion to buy degrees
of freedom. All three of its open questions were put to the study owner and
answered:

1. **H+1 timing confirmed.** The chamber volumes and every quantity derived
   from them were taken from the medical record at H+1 post-operatively. The
   values, held as `recipe.excluded_evidence` in
   [reyna_pre_surgery.m:72](config/calibration_recipes/reyna_pre_surgery.m:72):
   `LVEDV` 41.0 mL, `LVESV` 19.3 mL, `RVEDV` 30.5 mL, `RVESV` 12.0 mL,
   `LVEF` 0.528.
2. **H+1 accepted as a valid post-closure state.**
3. **The "inconsistency" needs no resolution — see §6.1, it was never real.**
   The study owner's instruction was to accept it, on the grounds that a
   pre-closure target cannot be compared against H+1 volumes anyway, and to
   validate predicted volumes against published paediatric ranges instead.
   Both points are correct; §6.1 shows the arithmetic, and the range check is
   already implemented (see §6.3 on strengthening it).

### 6.0 RESOLVED — the post-operative data exists and is now encoded

The source record was retrieved on 2026-08-29 and **the high-value case
obtained**. `clinical.post_surgery`, previously entirely `NaN`, now carries
**7 finite haemodynamic targets** taken in the same catheterisation session
as the pre-closure readings, after device release. Full values, provenance
and direction check are in §0.2.

*(An earlier draft of this section concluded no post-operative data existed
and closed Phase 4 permanently. That was a misreading of the study owner's
answer and is withdrawn.)*

**Phase 4 is now genuinely worth building.** The arithmetic that previously
argued against it has reversed:

| | Before | With post-closure pressures |
|---|---:|---:|
| Governed observations `N` | 9 | **16** (9 pre + 7 post) |
| Active parameters `p` | 12 | 12 |
| `dof = N − p` | **0** (`insufficient_dof`) | **4** |

A joint inversion sharing patient parameters across both states — with the
VSD orifice as the difference between them — is exactly the fix the
2026-08-28 assessment's §2.1 demands, and it is now supported by real paired
measurements rather than relocated echo rows. This is the single largest
remaining step toward a defensible publication claim.

Two caveats before building it:

- **No post-closure flow yet** (§0.3). The post state constrains pressures
  only, so `SVR`/`PVR`/`QpQs` cannot be evaluated post-closure until the
  PARI/FR question is answered.
- **`dof = 4` is an improvement, not a resolution.** Four degrees of freedom
  still makes reduced χ² unstable. The complementary lever is reducing `p`
  (currently 12 because Stages D–F use the broader mask, §4) — genuinely
  resolving DOF likely needs both.

### 6.1 The 60% stroke-volume inconsistency was an artefact, and is already gone

The PRD and the 2026-08-28 assessment both treat a "60% internal
stroke-volume inconsistency, severity `critical`" as an open problem in the
H+1 block. Checking the arithmetic directly, it is not internal to that block:

| Quantity | Value | Source |
|---|---:|---|
| `SV_LV` = LVEDV − LVESV | 21.7 mL/beat | H+1 echo |
| `SV_RV` = RVEDV − RVESV | 18.5 mL/beat | H+1 echo |
| **SV_LV vs SV_RV** | **14.7% apart** | both H+1 |
| `LVEF` check: 21.7 / 41.0 | 0.529 vs stated 0.528 | internally consistent ✓ |

Internally the H+1 block is coherent to within ordinary echo measurement
error. The 60% figure appears only when those post-closure volumes are
compared against **pre-closure** flow-derived stroke volumes:

- `SV_Qs` = 3.423 × 1000 / 119 = 28.77 mL/beat
- `SV_Qp` = 3.423 × 1.194 × 1000 / 119 = 34.35 mL/beat
- `SV_Qp` vs `SV_LV` = (34.35 − 21.7) / 21.7 = **58.3% ≈ the reported 60%**

That is not an inconsistency in the data. It is the expected physiological
difference between a shunt-loaded pre-operative ventricle and the same
ventricle after closure — the very reason the block was excluded from
pre-surgery fitting in the first place.

**It is also already resolved in the current code.** With the block excluded,
the 2026-08-29 run's audit reports severity `none`, max relative SV difference
**17.7%** — and that 17.7% is simply `SV_Qs` vs `SV_Qp`, i.e. the Qp/Qs shunt
ratio itself, with `SV_LV` and `SV_RV` correctly `NaN`. Nothing further is
required. Any remaining text in the PRD or the 2026-08-28 assessment
describing a live `critical` inconsistency is stale and should be read
against this section.

### 6.2 Validating predicted volumes against the literature

Per the study owner's answer to question 3, predicted chamber volumes should
be judged against published paediatric ranges rather than against this
patient's H+1 block. **This is already implemented** — the
`PREDICTED CHAMBER STATE` console block screens all six predictions and the
2026-08-29 run passes 6 of 6:

| Metric | Predicted | Screening range |
|---|---:|---|
| LVEDV | 45.9 mL | 5 – 120 |
| LVESV | 14.5 mL | 1 – 80 |
| RVEDV | 41.7 mL | 5 – 140 |
| RVESV | 12.2 mL | 1 – 90 |
| LVEF | 0.684 | 0.40 – 0.85 |
| RVEF | 0.708 | 0.30 – 0.80 |

**But these ranges are too wide to constitute evidence.** An LVEDV band of
5–120 mL admits essentially any physiologically possible value for a child;
passing it demonstrates only that the model is not absurd. For the literature
comparison to carry weight in a publication it needs **BSA-indexed normative
values** for this patient's demographics (3.17 years, 0.617 m²) — e.g.
LVEDV/BSA in mL/m² against published paediatric echo normals with a stated
z-score or percentile, not a fixed absolute band. That is a materially
tighter test and one the model could actually fail, which is what makes it
worth reporting. Not implemented; recorded as a concrete follow-up.

### 6.3 Superseded — the original three open questions

*(Retained for provenance; all three are answered above.)*

1. Confirm the H+1 timing against the source clinical record.
2. Confirm H+1 is an acceptable proxy for the converged post-closure state, or
   state the limitation explicitly if it is not.
3. Decide how to handle the block's **60% internal stroke-volume
   inconsistency** (severity `critical` per `audit_clinical_consistency`) once
   it becomes a `post_surgery` target — moving it relocates the
   inconsistency, it does not resolve it.

Phases 1, 2, 3, and 5 do not depend on this and are complete without it.

## 7. Parameter identifiability (Phase 3)

> **⚠ The numbers in §7 are SUPERSEDED by the §0 data correction** (same
> reason as §3: `HR` and `BSA` were wrong, and both are model inputs). The
> *structural* finding — that the elastance and unstressed-volume parameters
> are strongly collinear, and that `p = 12` rather than the GSA screen's 7 —
> is a property of the model's parameterisation rather than of the clinical
> values, so it is expected to persist; but the specific condition number and
> correlations must be regenerated before being quoted.

A 2-parameter, 2-metric smoke check during test development (`R.SVEN` vs
`C.SAR`, evaluated at `RAP_mean`/`SAP_mean`) showed **ρ = −1.000** — perfect
collinearity, though on a toy subset and at the demographic-scaled baseline
rather than the calibrated operating point. This is consistent with the
systemic RC time-constant coupling identified in
`docs/reyna_zhang_fullmetric_results_20260828.md` §3.5 (the waveform
form-factor mismatch driving `SAP_min`'s residual).

### Full report (this run's winning candidate, 9 metrics × 12 parameters)

Scaled sensitivity matrix `S(i,j) = dy_i/dtheta_j * theta_j/sigma_i`.
Condition number: **232**.

| Parameter | Column norm | Inactive | Max \|corr\| | Most correlated with | Flagged |
|---|---:|---|---:|---|---|
| group.R_sys_scale | 15.192 | false | 0.666 | R.SVEN | false |
| R.SVEN | 5.833 | false | 0.666 | group.R_sys_scale | false |
| group.R_pul_scale | 28.551 | false | 0.665 | E.LV.EB | false |
| C.SAR | 6.430 | false | 0.599 | group.R_pul_scale | false |
| C.PAR | 8.965 | false | 0.380 | E.LV.EA | false |
| E.LV.EA | 5.553 | false | 0.903 | E.LV.EB | **true** |
| E.LV.EB | 9.029 | false | 0.903 | E.LV.EA | **true** |
| E.RV.EA | 2.583 | false | 0.811 | V0.RV | false |
| E.RV.EB | 6.391 | false | 0.944 | V0.LV | **true** |
| V0.LV | 0.631 | false | 0.944 | E.RV.EB | **true** |
| V0.RV | 0.501 | false | 0.904 | E.RV.EB | **true** |
| vsd.Cd | 3.987 | false | 0.716 | E.LV.EA | false |

Collinear pairs flagged (|ρ| threshold exceeded):
- `E.RV.EB` ↔ `V0.LV` (ρ = -0.944)
- `E.RV.EB` ↔ `V0.RV` (ρ = +0.904)
- `E.LV.EA` ↔ `E.LV.EB` (ρ = -0.903)

This confirms the §7 smoke-test finding at governed-set scale rather than a
toy 2-parameter subset: the RV end-systolic/diastolic elastance pair and both
ventricular unstressed volumes are collinear with each other, and the LV
elastance pair is collinear with itself. None of the 12 active parameters are
fully inactive (`Inactive = false` throughout), so the mask itself is not
retaining dead weight — the identifiability problem is redundancy between
retained parameters, not inclusion of irrelevant ones. Condition number 232
is high enough to warrant caution interpreting individual parameter values
from this fit as uniquely determined, though not so high as to indicate the
fit itself is numerically degenerate.

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
- **The fit is underdetermined, and this is the binding limitation.**
  `N = 9` governed observations against `p = 12` free parameters gives
  `dof = −3`. A low χ²/N is *guaranteed* in that regime and is not evidence
  the model is correct — the report now says so explicitly (§4). Any claim
  from this branch must rest on the gate count and per-metric residuals, and
  must state the parameter/observation ratio alongside. Two levers exist and
  probably both are needed: raise `N` via Phase 4 joint pre/post inversion
  (now viable, §6.0 — would give `N = 16`, `dof = 4`), and lower `p` below
  12 (stages D–F currently use a broader mask than the GSA screen's 7).
- **The better fit came partly at the cost of identifiability.** Condition
  number rose from 232 to 2.06 × 10³ between the superseded and corrected
  runs, and a new `E.LV.EA` ↔ `vsd.Cd` collinearity (ρ = −0.917) appeared.
  A fit that improves while its parameters become less separable is a
  warning, not a success.
- **Single seed.** §3.5 is one seed (`20260828`). A second (`20260830`) is
  running for a spread estimate; until it lands there is no confidence
  interval, and §5.1 of the 2026-08-28 assessment requires one.
- **The identifiability report in §7 is now a governed-set analysis** (9
  metrics × 12 parameters, at this run's actual calibrated operating point)
  — no longer just the earlier 2-parameter smoke test, but still a single
  run's sensitivity matrix, not a distribution over runs.
- **One A/B comparison, one seed.** `UNIFIED_VSD_MULTISTART_SEED=20260828` was
  reused from PR #24 for comparability; a different seed could shift which
  start wins and by how much. This is not evidence the sigma-weighted mode
  generalizes past this single patient/seed/scenario combination.
- **Zhang vs Lundquist remains uninterpretable** at these budgets, unchanged
  from PR #24.

## 10. This work, and how to reproduce it

State as of 2026-08-29, branch `codex/reyna-statistical-calibration` (off
`codex/reyna-zhang-fullmetric-10pct`, PR #24): **all code for Phases 1, 2, 3,
and 5 is complete, tested, and committed**, and the deferred 6-start
sigma-weighted A/B run (§3.2–3.3) has now executed successfully — results are
folded into §3 and §7 above. 85/85 new tests pass; the only pre-existing
failure is `test_clinical_consistency_target_tiers.m` (a PVR tier assertion,
present on `main`, unrelated to this branch). This section is kept so the run
can be reproduced or repeated with a different seed.

### 10.1 The one rule for this run

**Do not edit any file this run depends on while it is executing.** §3.2
documents exactly what goes wrong: `main_run.m` freezes its own control flow
at launch (it stays on the call stack for the whole run), while functions it
calls are reloaded fresh on each invocation — so a mid-run edit produces a
silent hybrid of old and new code. If a genuine bug is found while a run is
in flight, let the run finish (or kill it) before touching the file.

**A second, unrelated failure mode surfaced launching this run from an
agent session on Windows**, worth recording alongside the above: a manually
backgrounded process (`nohup ... & disown` from within a shell tool call) is
not reliably detached from the console/job object on Windows the way it is
on Linux, and can be killed outright when the launching shell session is
recycled between tool invocations — independent of anything in this
codebase. The run that eventually produced §3's results was launched twice:
the first attempt died partway through start 3/6 with `Exit Status:
0x40010004` (Windows' external-termination code) for exactly this reason,
losing ~2.25 hours of progress with nothing to salvage (artefacts are only
written at the very end). The second attempt used the calling tool's native
tracked-background execution instead of a manual `nohup`, and completed
without incident. If reproducing this from an agent/automation context on
Windows: use whatever backgrounding mechanism that context tracks natively,
not a manually detached shell process.

### 10.2 The command

```bash
matlab -batch "cd('D:/Kuliah/Skripsi/CollabHafizKeisya/unified_vsd'); addpath(genpath(pwd)); setenv('UNIFIED_VSD_UQLAB_PATH', fullfile(pwd,'toolbox','UQLab_Rel2.2.0','core')); setenv('UNIFIED_VSD_SCALING_MODE','zhang'); setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','1'); setenv('UNIFIED_VSD_DO_GSA','1'); setenv('UNIFIED_VSD_GSA_PCE_N','128'); setenv('UNIFIED_VSD_DO_PLOTS','0'); setenv('UNIFIED_VSD_MAX_FUN_EVALS','300'); setenv('UNIFIED_VSD_MAX_ITERATIONS','40'); setenv('UNIFIED_VSD_NUM_STARTS','6'); setenv('UNIFIED_VSD_MULTISTART_SEED','20260828'); setenv('UNIFIED_VSD_OBJECTIVE_WEIGHTING','sigma'); main_run('pre_surgery', patient_reyna())"
```

Run this in the background and redirect output to a log file; do not run any
other MATLAB process concurrently (a prior 20-minute run stretched past 75
minutes when something else was competing for the same machine). Budget
70–80 minutes per start × 6 starts as a planning estimate, though the
2026-08-29 run that produced §3's results completed in **9743 s (~2.7
hours)** end to end — faster than that estimate in practice, but treat
4–5 hours as the number to plan around for scheduling purposes.

### 10.3 What "done" looks like (confirmed, 2026-08-29)

The run succeeds when the console log contains a `FULL METRIC 10% GATE` block
followed by a `CHI-SQUARED` block with **`p (active parameters)` equal to
`numel(calib_out.names)` for the winning start — a small positive integer,
not 0** (0 was the symptom of the contamination bug this section exists to
prevent a repeat of). A `PARAMETER IDENTIFIABILITY` block should also appear
(Phase 3). On the confirmed run this landed as `p = 12`, `condition number =
232` — see §7.

The run folder's `tables/` directory contains `full_metric_gate_pre_surgery.csv`,
`chi_squared_pre_surgery.csv`, `parameter_identifiability_pre_surgery.csv`,
and `parameter_identifiability_pairs_pre_surgery.csv`, at
`results/runs/20260829_161536_reyna_pre_surgery/tables/`. These are the
source for §3's and §7's tables above.

### 10.4 After a good result

If the sigma-weighted candidate is an improvement worth keeping, it can be
baked into `config/calibration_recipes/reyna_pre_surgery.m` as a new
`initial_parameter_values` seed (the existing seed there came from a 2026-05-21
run via the same pattern) — future runs would then take the
"Accepted explicit recipe seed" fast path and skip the multi-stage
optimisation entirely, rather than re-paying this multi-hour cost. This is a
deliberate follow-up decision, not something to do automatically.
