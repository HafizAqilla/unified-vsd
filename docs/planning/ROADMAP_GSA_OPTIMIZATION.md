# Research Team Roadmap: GSA-Informed L-BFGS-B Calibration Pipeline

Based on the team consensus, the codebase architecture must be refactored to support a strict **GSA-First $\rightarrow$ L-BFGS-B Optimization $\rightarrow$ Final GSA** pipeline. This document outlines the step-by-step implementation plan to achieve sub-5% error while strictly adhering to the agreed methodology.

---

## 📌 1. Core Agreements Summarized
1. **ODE Engine:** Strictly handles the physical derivation of Pressure, Volume, and Flow rates.
2. **Initial GSA (Sobol):** Must be run *first* to screen the 11 theoretical parameters and identify only the most "useful/influential" ones.
3. **L-BFGS-B Optimization:** Used by the main team (excluding Aida) to find the unknown physiological values, optimizing *only* the parameters selected by the initial GSA.
4. **Final GSA (Sobol):** Conducted at the very end of the pipeline on the optimized model to report final parameter sensitivities for publication.

---

## 🛠️ Phase 1: Initial GSA (Parameter Screening)
Before running any optimizer, the system must use Sobol variance analysis to map the landscape.

*   **Task 1.1: Decouple GSA from Calibration.** 
    Currently, GSA and Calibration are parallel toggles. Restructure `main_run.m` so that Phase 0 mapping feeds *directly* into the Initial GSA.
*   **Task 1.2: Define Sensitivity Thresholds.** 
    Establish a mathematical cutoff (e.g., Total Sensitivity Index $S_T > 0.05$). Parameters below this threshold are frozen at their deterministic Phase 0 values.
*   **Task 1.3: Dynamic Parameter Masking.**
    Write a bridging function (`utils/create_optimization_mask.m`) that takes the GSA output table and automatically constructs a reduced subset of free parameters (e.g., reducing the 11-parameter space to 4–6 highly influential ones).

## 🚀 Phase 2: L-BFGS-B Optimization Implementation
Migrating from Nelder-Mead (`fminsearch`) to L-BFGS-B (Limited-memory Broyden–Fletcher–Goldfarb–Shanno with Bounds).

*   **Task 2.1: MATLAB Solver Mapping.**
    *Note for Claude/Devs:* MATLAB's native optimization toolbox uses `fmincon` instead of a standalone L-BFGS-B function. We must configure `fmincon` with:
    `options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'HessianApproximation', 'lbfgs');`
    *(Alternatively, integrate a third-party L-BFGS-B MEX/wrapper if the team strictly demands the exact Fortran implementation).*
*   **Task 2.2: Manage ODE Gradient Noise.**
    L-BFGS-B is gradient-based. ODE solvers (`ode15s`) create numerical noise that breaks finite-difference gradients.
    *   *Action:* Enforce strict ODE tolerances (`RelTol = 1e-8`) during calibration.
    *   *Action:* Configure the optimizer's finite difference step size (`TypicalX` and `FiniteDifferenceStepSize`) to step *over* the numerical noise floor.
*   **Task 2.3: Objective Function Penalty Weights.**
    Implement heavily asymmetric weighting in `objective_calibration.m`. If primary targets (Qp/Qs, LVEF, MAP) drift beyond the target 5% error margin, multiply their cost penalty by 100x to force the L-BFGS-B algorithm within physiological bounds.

## 📊 Phase 3: Final GSA and Output Validation
Confirming the physiological integrity of the calibrated parameters.

*   **Task 3.1: Post-Calibration GSA.**
    Pass the optimal parameter vector ($X_{opt}$) back into `gsa_sobol_setup.m`. Run a final Sobol sweep to analyze how parameter dominance shifted after the bounds were tightened. 
*   **Task 3.2: 5% Error Assertion.**
    Expand `validation_report.m` to automatically flag any primary metric (Pressure, Volume, Flow) that exhibits $>5\%$ divergence from the clinical targets in the patient profile.

---

## 📝 Proposed Codebase Execution Flow (`main_run.m` Concept)

1. `params = params_from_clinical(patient);` *(Phase 0 determinism)*
2. `[S1, ST] = gsa_run_sobol(params, 'initial');` *(Initial screening)*
3. `active_params = filter_useful_parameters(ST, threshold);`
4. `optimal_params = run_lbfgsb_calibration(active_params, patient);` *(Targeting < 5% error)*
5. `final_metrics = simulate_ODE(optimal_params);` *(Resolving P, V, Q)*
6. `[S1_final, ST_final] = gsa_run_sobol(optimal_params, 'final');` *(Wrap-up)*
7. `generate_figures_and_tables();`
