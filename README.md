# Unified VSD Lumped-Parameter Model

**A cohesive, physically complete 0-D cardiovascular simulation for pre- and post-operative paediatric Ventricular Septal Defect (VSD) analysis.**

---

## 📖 Overview

This codebase simulates patient-specific paediatric cardiovascular haemodynamics using a robust **14-state Valenti RLC lumped-parameter model** (Valenti 2023, Eqs. 2.1–2.7). Built for high-reliability medical physics, it manages pre- and post-surgical clinical scenarios from one unified codebase:

| Scenario | Physics | `R.vsd` | Primary Calibration Targets |
|---|---|---|---|
| `pre_surgery` | Full VSD shunt active | Computed from clinical gradient/orifice via Ohm's law | Qp/Qs, PAP_mean, PVR, LVEDV, LVESV, LVEF |
| `post_surgery` | Shunt closed (`R.vsd → ∞`) | Fixed rigidly at 10⁶ mmHg·s/mL | SAP, MAP, SVR, LVEF, RVEF, chamber volumes |

## ✨ Key Features

- **14-State RLC Physics Engine**: Covers time-varying elastances in 4 cardiac chambers, advanced systemic/pulmonary circuitry, and piecewise VSD modeling.
- **Pediatric Allometric Scaling**: Employs Lundquist-informed volume adjustments to effectively scale baseline references across variable infant body surface areas (BSA).
- **Expanded Parameter Calibration (Latest Upgrade)**: `pre_surgery` modeling now relies on a highly rigorous **11-parameter phase** bounds-constrained refinement.
- **Global Sensitivity Analysis (GSA)**: Includes comprehensive Saltelli / Jansen Sobol implementations mapped across parameters.
- **Strict Scientific Guardrails**: Ensures adherence to strict medical physics reproducibility standards (see [AGENTS.md](AGENTS.md)).

---

## 🚀 Getting Started

### 1. Prerequisites and Setup
The codebase runs seamlessly on modern **MATLAB** distributions. To utilize the Global Sensitivity Analysis (GSA) tools effectively, please manually install the following toolboxes:

- **UQLab** (Release 2.2.0 or later): Framework for uncertainty quantification.
- **SoBioS**: Surrogate-based optimization and sensitivity analysis toolbox.

*(Note: These dependencies are not tracked in the repository due to licensing and size footprint. Place them in your corresponding `Toolbox/` paths or adjust the mappings inside `main_run.m`)*

### 2. Execution Guide

To simulate scenarios, edit or leverage the main runner script:

```matlab
% 0. Add everything to path (run once per MATLAB session)
addpath(genpath('D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd'));

% 1. Load a patient profile 
clinical = patient_profile_A();   % Profile A: 1.6-month high-flow infant
clinical = patient_profile_B();   % Profile B: older infant, high-pressure
% clinical = patient_template();  % Utilize this blank template structure for novel patients

% 2. Run the desired scenario
main_run('pre_surgery',  clinical)
main_run('post_surgery', clinical)

% 3. Run virtual patient sweep (both profiles, both scenarios)
run_virtual_patients()
```

You can toggle simulation behaviors directly inside `main_run.m`:
```matlab
DO_CALIBRATION = true; % Optimized for recent 11-dimension refinements
DO_GSA         = false;
DO_PLOTS       = true;
```

---

## ⚙️ Under The Hood: The Calibration Engine

Calibration employs a powerful **two-phase approach** mitigating ODE solver noise and prioritizing physiological intuition.

### Phase 0: Deterministic Mapping
`utils/params_from_clinical.m` maps 28 detailed patient clinical inputs directly to fundamental model parameters before simulating. Utilizing Ohm's law and basic mass conservation, it secures a pre-conditioned starting state mathematically verified to generally sit right around ±15% of all clinical targets sans optimization.

### Phase 1: High-Dimensional Nelder-Mead Polish
`run_calibration.m` refines the Phase 0 mapping using a **bound-constrained Nelder-Mead** optimizer (`fminsearch`). Our most recent upgrade amplifies this dimensional loop inside `calibration_param_sets.m`, tracking **11 separate physical constants** concurrently in `pre_surgery`:

1. `R.SAR` & `R.SVEN`: Systemic arterial and venous resistance
2. `R.PAR` & `R.PVEN`: Pulmonary arterial and venous resistance
3. `C.SVEN`: Systemic venous compliance
4. `E.LV.EA` & `E.LV.EB`: Left ventricular active and passive elastance
5. `E.RV.EA`: Right ventricular active elastance
6. `V0.LV` & `V0.RV`: Left and right ventricular unstressed volume
7. `R.vsd`: VSD shunt resistance 

**Solver Profile:**
We utilize `MaxFunEvals = 500` to properly accommodate the heightened parameters, completing convergence comfortably in 5–15 minutes bypassing traditional gradient disruptions.

---

## 🧪 Global Sensitivity Analysis (GSA)

`gsa_sobol_setup.m` + `gsa_run_sobol.m` handles our Sobol (Saltelli/Jansen) variance analysis, outputting sensitivity matrices aligned differently with core targets per operation tier:
- **Pre-surgery:** Qp/Qs, PAP_mean, PVR
- **Post-surgery:** LVEF, SAP_mean, SVR

---

## 📂 Project Architecture

```
unified_vsd/
├── main_run.m                    ← Scenario entry point (configuration toggles inside).
├── run_virtual_patients.m        ← Batch runner: handles comprehensive profile sweeps.
├── README.md                     ← Standard project orientation.
├── AGENTS.md                     ← Strict modeling guardrails and standardization laws.
│
├── config/                       ← Clinical profiles, default references, and patient templates.
├── models/                       ← Base physical rules: 14-state ODE, elastance limits, valve modeling.
├── solvers/                      ← Main simulation drivers (`ode15s` integrations, convergence bounds).
├── utils/                        ← Allometric sizing, structural determinism, plotting wrappers.
├── calibration/                  ← Core parameter mapping sets, cost definitions, optimization logic.
├── gsa/                          ← Sobol execution environments and summary visualizations.
├── tests/                        ← Core health verifications against physiological limits.
├── results/                      ← Un-tracked raw matrix outputs safely separated from vector figures.
└── docs/                         ← Underlying mathematical schemas and nomenclature structures.
```

---

## 🔧 Numerical Reference Bounds

| Setting | Baseline / Validation | Calibration Context |
|---|---|---|
| `nCyclesSteady` | **80** (Crucial mapping limit for high-flow VSD) | 40 (Fast mapping override) |
| `ss_tol_P` | **0.1 mmHg** | 0.5 mmHg |
| `ss_tol_V` | **0.1 mL** | 0.5 mL |
| `rtol` / `atol` | 1e-7 / 1e-8 | same |
| Integrator | `ode15s` (Stiff setup, BDF) | same |
| Execution timing | ~30–40 s | ~3–4 s |

---

## ✅ Quality Standards & Testing

Always execute the testing suite prior to developing any patient's specific profile or injecting code (see AGENTS.md §10.1):

```matlab
addpath(genpath('D:\Kuliah\Skripsi\CollabHafizKeisya\unified_vsd'));
run('tests/test_baseline.m')
run('tests/test_valve_logic.m')
run('tests/test_ic_perturbation.m')
```

---

## 👥 Authors 

| Contributor | Focus Areas |
|---|---|
| **Hafiz** | Lundquist allometric scaling, patient IC formulas, baseline structural determinism. |
| **Keisya** | Valenti 14-state model mapping, Sobol GSA integrations. |
| **Joint** | 11-Parameter Multi-Dimensional Calibration Engine, architecture management, clinical abstractions. |

*Last Updated: March 2026*

---

## 📚 References

1. **Valenti (2023).** *Thesis: Full-order 0D cardiovascular model*. Table 3.3, Eqs. 2.1–2.7.
2. **Lundquist et al. (2025).** *Patient-Specific Pediatric Cardiovascular Lumped Parameter Modeling*. ASAIO J.
3. **Saltelli et al. (2010).** *Variance based sensitivity analysis of model output*. CPC 181:259–270.
4. **Jansen (1999).** *Analysis of variance designs for model output*. CPC 117:35–43.
5. **Torczon V (1991).** *On the convergence of pattern search algorithms*. SIAM J. Optimization 1(2):123–145. (Nelder-Mead convergence theory)
6. **Gorlin R & Gorlin SG (1951).** *Hydraulic formula for valve area*. Am Heart J 41(1):1–29.
