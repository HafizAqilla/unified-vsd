# Unified VSD Lumped-Parameter Model

**A single, physically complete 0-D cardiovascular simulation for pre- and post-operative paediatric ventricular septal defect (VSD) analysis.**

---

## What this model does

This codebase simulates patient-specific paediatric cardiovascular haemodynamics using a **14-state Valenti RLC lumped-parameter model** (Valenti 2023, Eqs. 2.1–2.7).  It handles two clinical scenarios from one unified codebase:

| Scenario | Physics | `R.vsd` | Primary calibration targets |
|---|---|---|---|
| `pre_surgery` | Full VSD shunt active (finite `R.vsd`) | Computed from clinical gradient/orifice via Ohm's law | Qp/Qs, PAP_mean, PVR, LVEDV, LVESV, LVEF |
| `post_surgery` | Shunt closed (`R.vsd → ∞`) | Fixed at 10⁶ mmHg·s/mL | SAP, MAP, SVR, LVEF, RVEF, chamber volumes |

---

## How to run

```matlab
% 0. Add everything to path (run once per MATLAB session)
addpath(genpath('D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd'));

% 1. Load a patient profile (or edit patient_template.m with your own data)
clinical = patient_profile_A();   % Profile A: 1.6-month high-flow infant
clinical = patient_profile_B();   % Profile B: older infant, high-pressure
clinical = patient_template();    % Blank template to fill manually

% 2. Run the desired scenario
main_run('pre_surgery',  clinical)
main_run('post_surgery', clinical)

% 3. Run virtual patient sweep (both profiles, both scenarios)
run_virtual_patients()
```

Toggle `DO_CALIBRATION`, `DO_GSA`, and `DO_PLOTS` at the top of `main_run.m`.

---

## Project structure

```
unified_vsd/
│
├── main_run.m                    ← Entry point only. No physics, no plotting.
├── run_virtual_patients.m        ← Batch runner: both profiles × both scenarios.
├── README.md                     ← This file.
├── AGENTS.md                     ← World-class MATLAB coding standards.
├── .gitignore                    ← Excludes .venv, .vscode, PDFs, *.asv, results/*.mat
│
├── config/
│   ├── default_parameters.m      ← Adult Valenti reference (BSA = 1.73 m²).
│   │                                nCyclesSteady = 80 (adequate for high-flow VSD).
│   ├── patient_template.m        ← Unified clinical struct template (blank).
│   ├── patient_profile_A.m       ← Profile A: 1.6-month, 3.7 kg, VSD=5.7 mm,
│   │                                Qp/Qs=3.44, PAP_mean=28, MAP=52 mmHg.
│   └── patient_profile_B.m       ← Profile B: older infant, high PA pressure.
│
├── models/
│   ├── system_rhs.m              ← 14-state RLC ODE (Valenti Eqs. 2.1–2.7).
│   ├── elastance_model.m         ← Time-varying piecewise-cosine elastance.
│   └── valve_model.m             ← Non-ideal diode valve with tanh smoothing.
│
├── solvers/
│   └── integrate_system.m        ← ode15s cycle-by-cycle integration;
│                                    steady-state convergence check each cycle;
│                                    returns last nCyclesKeep cycles.
│
├── utils/
│   ├── apply_scaling.m           ← BSA allometric scaling + blood-volume IC.
│   ├── params_from_clinical.m    ← Deterministic Phase 0 mapping: HR, SVR, PVR,
│   │                                R.vsd, elastances from 28 clinical inputs.
│   ├── compute_clinical_indices.m← All haemodynamic metrics from ODE output.
│   ├── validation_report.m       ← Error% table vs clinical targets; RMSE summary.
│   ├── unit_conversion.m         ← All unit conversion factors (single source).
│   └── plotting_tools.m          ← Publication-ready PDF figures.
│
├── calibration/
│   ├── calibration_param_sets.m  ← Free-parameter lists, weights, bounds per
│   │                                scenario. pre_surgery: R.vsd + 4 elastances.
│   │                                regLambda = 0 (no Tikhonov regularisation).
│   ├── run_calibration.m         ← Nelder-Mead (fminsearch) optimiser with
│   │                                quadratic bound penalty. ~100–200 evals,
│   │                                ~5–15 min for d=5 parameters.
│   └── objective_calibration.m   ← Weighted normalised least-squares cost.
│                                    Calibration override: nCyclesSteady=40,
│                                    ss_tol=0.5 (speed). Validation uses 80/0.1.
│
├── gsa/
│   ├── gsa_sobol_setup.m         ← Saltelli sample matrices + metric config.
│   ├── gsa_run_sobol.m           ← Jansen S1/ST estimators.
│   └── make_gsa_summary_table.m  ← Ranked sensitivity summary table.
│
├── tests/
│   ├── test_baseline.m           ← Physiological range assertions.
│   │                                Must pass before any patient simulation.
│   ├── test_valve_logic.m        ← Valve flow direction and switching tests.
│   └── test_ic_perturbation.m    ← ±10% IC perturbation → same limit cycle check.
│
├── results/
│   ├── figures/                  ← *.pdf figures (tracked). *.png ignored.
│   └── tables/                   ← *.csv tables (tracked). *.mat ignored.
│
└── docs/
    ├── theory_notes.md           ← Governing equations, assumptions, references.
    └── clinical_data_dictionary.md ← Maps 28 clinical fields to MATLAB variables.
```

---

## Physics engine

The model is the **14-state Valenti RLC** system (`models/system_rhs.m`):

- 4 cardiac chambers: RA, RV, LA, LV (time-varying elastance)
- Systemic circuit:   SAR (RLC) → SC (RC) → SVEN (RLC)
- Pulmonary circuit:  PAR (RLC) → P_PC (combined capillary pressure state, Eq. 2.7) → PVEN (RLC)
- VSD shunt:         Q_VSD = (P_LV − P_RV) / R_vsd

When `R.vsd` is very large (post-surgery), Q_VSD ≈ 0 automatically — no new equation needed.

---

## Allometric scaling

`utils/apply_scaling.m` scales from the adult Valenti reference to a paediatric patient:

| Quantity | Exponent | Physical basis |
|---|---|---|
| Resistance R | s⁻¹ | Poiseuille (wider vessels → lower R) |
| Compliance C | s⁺¹ | Elastic volume ∝ body size |
| Inertance L | s⁻¹ | Same geometry as R |
| Unstressed V₀ | s⁺¹ | Compartment size ∝ body size |
| Heart rate HR | s⁻⁰·³³ | Metabolic / basal scaling |
| Elastance E_LV/LA | s⁻¹ | Systemic pressure load |
| Elastance E_RV/RA | s⁻¹·⁵ | Pulmonary (lower pressure) |

Blood volume initial conditions are corrected using:
- **Infant (< 1 year):** 82 mL/kg
- **Older patient (≥ 1 year):** 70 mL/kg

---

## Calibration architecture

Calibration uses a **two-phase** approach:

### Phase 0 — Deterministic parameter mapping (zero optimizer calls)

`utils/params_from_clinical.m` maps all 28 clinical inputs to model parameters analytically before any simulation, using Ohm's law and mass conservation:

```matlab
% Systemic resistance from Ohm's law
params.R.SAR = (SAP_mean - RAP_mean) / (Qs_mLs * R_SAR_fraction);

% VSD resistance from measured shunt gradient and flow
params.R.vsd = DeltaP_VSD_mmHg / Q_shunt_mLs;   % Option A (preferred)

% Elastance bounds from echo volumes + pressure (Valenti Eq. 2.1)
params.E.LV.EA = (SAP_mean / (LVEDV - V0_LV)) - params.E.LV.EB;
```

After Phase 0, the model baseline is typically within ±15% of all clinical targets without any optimisation.

### Phase 1 — Nelder-Mead polish (elastances + R.vsd)

`run_calibration` refines the pre-conditioned baseline using **Nelder-Mead** (`fminsearch`) with a quadratic bound penalty. No gradients are computed — immune to ODE solver noise.

Free parameters for `pre_surgery`:
```
R.vsd    — VSD shunt resistance  (bounds: 0.05×–20× analytical x0)
E.LV.EA  — LV active elastance   (bounds: 0.3×–3.0× echo-derived x0)
E.LV.EB  — LV passive elastance
E.RV.EA  — RV active elastance
E.RV.EB  — RV passive elastance
```

Objective (no regularisation — `regLambda = 0`):

$$J(x) = \sum_k w_k \left(\frac{y_k(x) - y_k^{\text{clin}}}{y_k^{\text{clin}}}\right)^2$$

| Metric | Weight | Rationale |
|---|---|---|
| LVEF | 4.0 | Primary elastance target |
| LVEDV, LVESV | 3.0 | LV volume overload signature |
| RVEDV, RVESV | 2.5 | RV pressure load |
| SAP_mean | 2.0 | MAP check |
| PAP_mean, PVR, QpQs | 1.5 | Pulmonary/shunt targets |
| All others | 1.0 | Equal weight |

**Calibration speed settings** (set inside `objective_calibration.m`):
- `nCyclesSteady = 40` per eval (validation/baseline uses 80)
- `ss_tol_P/V = 0.5 mmHg/mL` per eval (validation uses 0.1)
- Convergence: ~100–200 evaluations ≈ **5–15 minutes**

**Solver comparison:**

| Solver | Evals | Time (d=5) | Gradient needed | Used here |
|---|---|---|---|---|
| `fmincon` | 200–500 | 15–50 min | Yes — fails with ODE noise | ❌ Removed |
| `particleswarm` | 200–800 | 15–50 min | No | ❌ Removed (overkill) |
| `fminsearch` (Nelder-Mead) | 100–200 | **5–15 min** | No | ✅ Current |

---

## Global Sensitivity Analysis

`gsa_sobol_setup` + `gsa_run_sobol` implement the **Saltelli / Jansen Sobol** method.  The output metrics of interest are scenario-aware:

- **Pre-surgery:** Qp/Qs, PAP_mean, PVR
- **Post-surgery:** LVEF, SAP_mean, SVR

---

## Units

| Quantity | Internal unit | Clinical reporting |
|---|---|---|
| Pressure | mmHg | mmHg |
| Volume | mL | mL |
| Flow | mL/s (ODE) | L/min (metrics) |
| Resistance | mmHg·s/mL | Wood units (÷ 0.06) |
| Compliance | mL/mmHg | mL/mmHg |

---

## Dependencies

To run the Global Sensitivity Analysis (GSA) module, the following third-party toolboxes are required. **They are not tracked in this repository** due to size limits and licensing. You must download them and place them in the `Toolbox/` directory (or update the paths in `main_run.m`):

1. **UQLab** (Release 2.2.0 or later): Framework for uncertainty quantification.
2. **SoBioS**: Surrogate-based optimization and sensitivity analysis toolbox.

---

## Coding Standards

All code in this repository follows the guardrails defined in [AGENTS.md](AGENTS.md). Key rules:

- State vector accessed only via `params.idx` — never hardcoded indices
- All unit conversions via `params.conv.*` named constants (see `utils/unit_conversion.m`)
- Valve switching uses smooth `tanh` blending (`params.epsilon_valve`)
- Solver choice and tolerances justified in `solvers/integrate_system.m`
- All figures exported as vector PDF; publication-ready from day one
- `tests/test_baseline.m` must pass before any patient simulation

---

## Numerical settings summary

| Setting | Baseline / Validation | Calibration evals |
|---|---|---|
| `nCyclesSteady` | **80** (high-flow VSD needs ~60–80 cycles) | 40 (speed override) |
| `ss_tol_P` | **0.1 mmHg** (AGENTS.md §8.3) | 0.5 mmHg |
| `ss_tol_V` | **0.1 mL** | 0.5 mL |
| `rtol` / `atol` | 1e-7 / 1e-8 | same |
| Solver | `ode15s` (stiff BDF) | same |
| Wall time / run | ~30–40 s | ~3–4 s |

---

## Running tests

Always run before any patient simulation (AGENTS.md §10.1):

```matlab
addpath(genpath('D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd'));
run('tests/test_baseline.m')
run('tests/test_valve_logic.m')
run('tests/test_ic_perturbation.m')
```

---

## Authors

| Contributor | Primary contributions |
|---|---|
| **Hafiz** | Lundquist allometric scaling, patient blood-volume IC formula, 8-state RC reference model |
| **Keisya** | Valenti 14-state RLC core model, Sobol GSA pipeline |
| **Joint** | Calibration rewrite (Nelder-Mead), patient profiles A/B, test suite, `.gitignore`, documentation |

Initial merge: 2026-02-26 — Last updated: 2026-03-03

---

## References

1. Valenti (2023). *Thesis: Full-order 0D cardiovascular model*. Table 3.3, Eqs. 2.1–2.7.
2. Lundquist et al. (2025). *Patient-Specific Pediatric Cardiovascular Lumped Parameter Modeling*. ASAIO J. DOI: 10.1097/MAT.0000000000002528.
3. Saltelli et al. (2010). *Variance based sensitivity analysis of model output*. CPC 181:259–270.
4. Jansen (1999). *Analysis of variance designs for model output*. CPC 117:35–43.
5. Torczon V (1991). *On the convergence of pattern search algorithms*. SIAM J. Optimization 1(2):123–145. (Nelder-Mead convergence theory)
6. Gorlin R & Gorlin SG (1951). *Hydraulic formula for valve area*. Am Heart J 41(1):1–29.
