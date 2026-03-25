# Unified VSD — Batched Improvement Plan (Revised)

Ordered **correctness-first, modularity-second**. Each batch is self-contained and committable independently.

> [!NOTE]
> Revision incorporates all 9 technical corrections from user review (rounds 1+2):
> 1. `epsilon_valve` → 0.1 mmHg; centralise Q_VSD in one helper (DRY)
> 2. Sobol default N → 512; explicit convergence criterion
> 3. SS flow check uses relative tolerance
> 4. Calibration barrier replaced with uncertainty-weighted loss
> 5. Thesis text sync added as Batch 3.5
> 6. GSA relaxed cycle count (`gsa_ss_max_cycles = 10`) for wall-clock feasibility
> 7. State matrix variable renamed `X_k` with unit verification comment
> 8. Sigma values require per-value citation in calibration_param_sets.m
> 9. Batch 3.5 verification step added to table

---

## Batch 1 — Correctness-Critical Fixes (Physics & Numerics)

Items that can produce **wrong results** if unfixed.

---

### 1.1  VSD Shunt Diode — Centralised, Unidirectional

**Problem:** `Q_VSD = (P_LV - P_RV) / R.vsd` is duplicated in 3 files and allows reverse shunt.

**Fix — two parts:**

**(A) Create a single helper** `src/models/vsd_shunt_model.m`:

```matlab
function Q_VSD = vsd_shunt_model(P_LV, P_RV, params)
% VSD_SHUNT_MODEL — Smooth-diode VSD shunt flow (L→R only)
%   Q_VSD = gate(ΔP) · ΔP / R_vsd     [mL/s]
%   gate  = 0.5 + 0.5·tanh(ΔP / ε)    sigmoid, ε = 0.1 mmHg
%
% epsilon_valve = 0.1 mmHg chosen to preserve >99% of diastolic
% L→R flow at ΔP ≈ 2 mmHg (restrictive pediatric VSD).
% At 0.5 mmHg the gate clips ~10% of diastolic shunt — too aggressive.

dP   = P_LV - P_RV;                                       % [mmHg]
gate = 0.5 + 0.5 * tanh(dP / params.epsilon_valve);       % smooth diode
Q_VSD = gate .* dP / params.R.vsd;                         % [mL/s]
end
```

**(B) Replace all 3 inline computations** with calls to the helper:

| File | Line(s) | Replace with |
|---|---|---|
| system_rhs.m | 110 | `Q_VSD = vsd_shunt_model(P_LV, P_RV, params);` |
| compute_clinical_indices.m | 185 | `Q.VSD = vsd_shunt_model(P_LV, P_RV, params);` |
| plotting_tools.m | 227 | `Q.VSD = vsd_shunt_model(P_LV, P_RV, params);` |

**(C) Add default** in default_parameters.m:

```matlab
params.epsilon_valve = 0.1;   % [mmHg] — VSD diode gate width
                               % 0.1 preserves >99% of diastolic L→R flow at ΔP=2 mmHg
                               % See vsd_shunt_model.m header for derivation.
```

---

### 1.2  Sobol Sample Size — N=512 with Convergence Criterion

**Changes:**

| File | Change |
|---|---|
| `gsa_sobol_setup.m:34` | `cfg.N = 512;` — delete the misleading "temporarily reduced" comment entirely |
| gsa_run_sobol.m | Add bootstrap CI (95%) for S1 and ST |

**Convergence criterion** (add as documented check in main_run.m or a helper):

$$\text{converged if } \max_i \lvert\hat{S}_{T_i}^{(N)} - \hat{S}_{T_i}^{(N/2)}\rvert < 0.05 \quad \forall\, i$$

```matlab
cfg.N_convergence = [128 256 512 1024];
% For each consecutive pair, compute max |ΔST|
% Report and save stability table to results/gsa/
```

**Wall-clock feasibility:** For d=19, N_total = N(d+2) = 512×21 = 10,752 evaluations.

> [!WARNING]
> The ~3s/eval estimate assumes **GSA-mode relaxed cycle count**, not the full 80-cycle calibration warmup. Without this, N=512 at full fidelity may take **days**, not hours.

Add to gsa_sobol_setup.m:
```matlab
% GSA evaluations use reduced warmup for speed (SS not critical per-sample).
cfg.gsa_sim_overrides.nCyclesSteady = 10;   % not 80
cfg.gsa_sim_overrides.ss_tol_P      = 1.0;  % [mmHg] relaxed
cfg.gsa_sim_overrides.ss_tol_V      = 1.0;  % [mL]   relaxed
```

These overrides must be applied inside gsa_run_sobol.m before each `integrate_system` call. Use N=256 for development; report N=512 in thesis.

---

### 1.3  Steady-State Convergence — Add Flow States with Relative Tolerance

| File | Change |
|---|---|
| `integrate_system.m:80-82` | Add flow state indices to convergence check |
| default_parameters.m | Add `params.sim.ss_rtol = 0.005` (0.5% relative) |

```matlab
flow_state_idx = [sidx.Q_SAR sidx.Q_SVEN sidx.Q_PAR sidx.Q_PVEN];

% Inside cycle loop, after computing delta_V and delta_P:
% NOTE: X_k is the full state matrix [n×14]. Columns at flow_state_idx
% hold Q states in [mL/s] — verified against system_rhs.m state layout.
mean_Q = mean(abs(X_k(:, flow_state_idx)), 1);   % cycle-mean |Q| [mL/s]
delta_Q_rel = abs(peak_now(flow_state_idx) - peak_prev(flow_state_idx)) ...
              ./ max(mean_Q, 1e-6);
ss_converged_Q = all(delta_Q_rel < params.sim.ss_rtol);

% Combined check:
if max(delta_V) < ss_tol_V && max(delta_P) < ss_tol_P && ss_converged_Q
    ss_reached = true;
    ...
end
```

> [!NOTE]
> The existing variable `V_k` in integrate_system.m holds the **full state vector** (volumes, flows, pressure), not just volumes. Rename it to `X_k` during this edit to avoid confusion with the volume-only naming convention used elsewhere.

---

### 1.4  Add `Refine = 4` to ODE Solver

| File | Change |
|---|---|
| `integrate_system.m:64` | Add `'Refine', 4` to `odeset` call |

---

## Batch 2 — Missing Clinical Outputs

No physics changes — required for thesis tables.

### 2.1  Add Metrics to compute_clinical_indices.m

| Metric | Formula | Unit |
|---|---|---|
| `CO_Lmin` | `LVSV * HR / 1000` | L/min |
| `Qp_mean_mLs` | store existing `Qpul_mLs` | mL/s |
| `Qs_mean_mLs` | store existing `Qsys_mLs` | mL/s |
| `VSD_frac_pct` | `100 * Q_shunt_mean_mLs / max(Qp_mean_mLs, 1e-6)` | % |
| `RAP_min` / `RAP_max` | `min/max(Pc.RA)` | mmHg |
| `LAP_min` / `LAP_max` | `min/max(Pc.LA)` | mmHg |

### 2.2  Add Rows to validation_report.m

Add CO, VSD_frac, RA/LA sys/dia to the report table.

---

## Batch 3 — Documentation ↔ Code Sync

### 3.1  Fix README

- Replace "Nelder-Mead / `fminsearch`" → "fmincon interior-point + L-BFGS"
- Remove `DO_CALIBRATION` / `DO_GSA` toggle block (these don't exist in main_run.m)
- Replace Torczon reference with Byrd/Lu/Nocedal (1995)

### 3.2  Centralise Heuristic Constants

Move from params_from_clinical.m into default_parameters.m:
```matlab
params.const.DP_ref_vsd_mmHg = 20;
params.const.rho_blood_kg_m3 = 1060;
params.const.Cd_orifice      = 0.7;
```

### 3.3  Document V0 Scaling Decision in `theory_notes.md`

Add a section explaining BV-ratio vs BSA-law with Guyton (1991) §15 citation.

---

## Batch 3.5 — Thesis Text Alignment

> [!IMPORTANT]
> These are mandatory changes to the thesis manuscript when the corresponding code changes land.

| Batch | Thesis Section | Change |
|---|---|---|
| 1.1 | Ch.3 §3.8.3 (VSD model) | Add: *"Shunt flow is constrained to L→R via a smooth sigmoid gate with ε = 0.1 mmHg"* |
| 1.2 | Ch.3 §2.8.2 (GSA) | State N=512, d=19, N_total=10752; cite convergence criterion |
| 1.3 | Ch.3 §3.9 (solver) | Add flow-state relative convergence to SS criterion list |
| 1.4 | Ch.3 §3.9 (solver) | Add `Refine = 4` to the odeset config table |
| 5.1 | Ch.3 §3.10 (calibration) | Describe uncertainty-weighted objective; cite inverse-σ² rationale |

---

## Batch 4 — Visualization & Export

### 4.1  Combined Valve Flow Panel

Add Figure 6 to plotting_tools.m: overlay Q.AV, Q.MV, Q.TV, Q.PVv, Q.VSD.

### 4.2  Pre-vs-Post Overlay Plot

New function `plot_comparison(sim_pre, params_pre, sim_post, params_post)`.

### 4.3  CSV Waveform Export

Export last-cycle pressure/flow table to `results/waveforms_<tag>_<scenario>.csv`.

---

## Batch 5 — Calibration Objective Improvement

### 5.1  Replace Hard Barrier with Uncertainty-Weighted Loss

**Root cause:** The 100× barrier is an ad-hoc substitute for proper measurement uncertainty weighting.

**Fix:** Weight each target by inverse measurement uncertainty σ_i⁻², estimated from RSAB inter-patient variance or clinical measurement repeatability literature:

```matlab
% In calibration_param_sets.m, define per-metric uncertainty:
%
% sigma values: measurement uncertainty (1 SD) — each MUST have a citation.
%   QpQs:     ±0.15 [-]    echo Doppler Qp/Qs repeatability  [cite: e.g. Valenti 2023 Table 3.6]
%   SAP_mean: ±5.0  [mmHg] oscillometric cuff BP SD           [cite: e.g. Stergiou 2018]
%   LVEF:     ±0.05 [-]    2D echo biplane Simpson EF         [cite: e.g. Lang 2015 ASE guidelines]
%   PAP_mean: ±3.0  [mmHg] right heart cath measurement       [cite: e.g. Hoeper 2006]
%
% Without citations, a thesis examiner WILL ask "where did these σ values come from?"
%
calib.sigma.QpQs     = 0.15;   % [-]    Source: [fill citation]
calib.sigma.SAP_mean = 5.0;    % [mmHg] Source: [fill citation]
calib.sigma.LVEF     = 0.05;   % [-]    Source: [fill citation]
calib.sigma.PAP_mean = 3.0;    % [mmHg] Source: [fill citation]
% ... etc for all targets

% In objective_calibration.m, replace barrier with:
w = 1 / max(calib.sigma.(mf), 1e-6)^2;   % inverse-variance weighting
J = J + w * (y_model - y_clin)^2;         % absolute, not normalised
```

**Fallback** (if σ estimates are not available): use the smooth ramp `w = w * (1 + 20 * max(0, err_rel - 0.05))` but document the magic number 20 as approximately `1/0.05` — the reciprocal of the acceptable relative error threshold.

---

## Batch 6 — Project Hygiene

### 6.1  Create `run_all.m`

Full regression: tests → virtual patients → both scenarios.

### 6.2  Consolidate Duplicate `results/` Directories

Remove `src/results/` and consolidate into top-level `results/`.

---

## Verification Plan

### Automated (per batch)

| Batch | Command | Pass Criterion |
|---|---|---|
| 1 | `test_baseline.m` + `test_valve_logic.m` | All pass; `Q_VSD ≥ 0` everywhere in pre-surgery |
| 1.3 | `test_baseline.m` | SS reached with new flow check active |
| 2 | `main_run('pre_surgery', patient_profile_A())` | `metrics.CO_Lmin`, `.VSD_frac_pct`, `.RAP_min`, `.LAP_max` exist and are finite |
| 3.5 | Manual cross-check | N=512, ε=0.1, ss_rtol=0.005, Refine=4 in thesis text match `default_parameters.m` / `gsa_sobol_setup.m` exactly |
| 5 | Compare `calib_out.fbest` before/after | Secondary targets improve without degrading primaries |
| 6 | `run_all()` completes without error | — |

### Manual (user-required)

1. **Batch 1.2:** Verify wall-clock time at N=512 is acceptable (with `gsa_ss_max_cycles=10`, expect ~3–9 hrs)
2. **Batch 3.5:** Update thesis LaTeX sections — cross-check every numeric constant against code defaults
3. **Batch 5:** Fill citation placeholders for each σ value before thesis submission
