# Unified VSD Lumped-Parameter Model

**A single, physically complete 0-D cardiovascular simulation for pre- and post-operative paediatric ventricular septal defect analysis.**

---

## What this model does

This codebase simulates patient-specific paediatric cardiovascular haemodynamics using a **14-state Valenti RLC lumped-parameter model**.  It handles two clinical scenarios from one unified codebase:

| Scenario | Physics | R_VSD | Calibration targets |
|---|---|---|---|
| `pre_surgery` | Full VSD shunt active (finite R_VSD) | Computed from clinical gradient / orifice equation | Qp/Qs, PAP, PVR, shunt flow |
| `post_surgery` | Shunt closed (R_VSD → ∞) | Fixed at 10⁶ mmHg·s/mL | SAP, MAP, SVR, LVEF, chamber volumes |

---

## How to run

```matlab
% 1. Fill in your patient's clinical data
clinical = patient_template();      % edit this file with real measurements

% 2. Run the desired scenario
main_run('pre_surgery',  clinical)
main_run('post_surgery', clinical)
```

Toggle `DO_CALIBRATION`, `DO_GSA`, and `DO_PLOTS` inside `main_run.m`.

---

## Project structure

```
unified_vsd/
│
├── main_run.m                    ← Entry point. No physics. No plotting logic.
├── README.md                     ← This file.
│
├── config/
│   ├── default_parameters.m      ← Adult Valenti reference (BSA = 1.73 m²)
│   └── patient_template.m        ← Unified clinical struct template
│
├── models/
│   ├── system_rhs.m              ← 14-state RLC ODE (Valenti Eqs. 2.1–2.7)
│   ├── elastance_model.m         ← Time-varying chamber elastance (Eqs. 2.2–2.4)
│   └── valve_model.m             ← Non-ideal diode valve (Eq. 2.6)
│
├── solvers/
│   └── integrate_system.m        ← ode15s wrapper; returns last N cycles
│
├── utils/
│   ├── apply_scaling.m           ← BSA allometric scaling + blood-vol IC correction
│   ├── params_from_clinical.m    ← Map HR/SVR/PVR/R_VSD from clinical struct
│   ├── compute_clinical_indices.m← Derive all haemodynamic metrics from sim output
│   ├── validation_report.m       ← Scenario-aware comparison table + RMSE
│   └── plotting_tools.m          ← Publication figures (pressures, PV loops, flows)
│
├── calibration/
│   ├── calibration_param_sets.m  ← Scenario-specific free-parameter lists + bounds
│   ├── run_calibration.m         ← fmincon / MultiStart wrapper
│   └── objective_calibration.m   ← Weighted normalised least-squares cost
│
├── gsa/
│   ├── gsa_sobol_setup.m         ← Saltelli sample matrices + scenario metric config
│   ├── gsa_run_sobol.m           ← Jansen S1/ST estimators
│   └── make_gsa_summary_table.m  ← Pretty summary of top influential parameters
│
├── tests/
│   └── test_baseline.m           ← Must pass before any patient simulation
│
├── results/
│   ├── figures/                  ← Auto-generated; not version-controlled
│   └── tables/                   ← Saved .mat / .csv outputs
│
└── docs/
    ├── theory_notes.md           ← Governing equations, assumptions, references
    └── clinical_data_dictionary.md ← Maps clinical terms to MATLAB variables
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

## Calibration

`run_calibration` uses `fmincon` with MultiStart.  The free-parameter set is **scenario-specific** (see `calibration/calibration_param_sets.m`) but the objective function and solver are identical:

$$J(x) = \sum_k w_k \left(\frac{y_k(x) - y_k^{\text{clin}}}{y_k^{\text{clin}}}\right)^2 + \lambda \sum_i \left(\frac{x_i - x_0^i}{x_0^i}\right)^2$$

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

## Authors

Unified from:
- **Hafiz** — Lundquist allometric scaling, patient blood-volume formula, 8-state RC model (physics superseded)
- **Keisya** — Valenti 14-state RLC model (adopted as core), fmincon calibration, Sobol GSA

Merged and extended: 2026-02-26

---

## References

1. Valenti (2023). *Thesis: Full-order 0D cardiovascular model*. Table 3.3, Eqs. 2.1–2.7.
2. Lundquist et al. (2025). *Allometric Scaling in Pediatric Cardiology*. Eqs. 3.2, 4.1, 4.3.
3. Saltelli et al. (2010). *Variance based sensitivity analysis*. Computer Physics Communications 181:259–270.
4. Jansen (1999). *Analysis of variance designs for model output*. Computer Physics Communications 117:35–43.
