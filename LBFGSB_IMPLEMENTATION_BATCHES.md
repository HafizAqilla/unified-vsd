# Implementation Batches: GSA-Informed L-BFGS-B Pipeline

This document outlines the sequential coding batches required to fully refactor the codebase to hit the < 5% error target using the GSA $\rightarrow$ L-BFGS-B $\rightarrow$ Final GSA pipeline.

---

## 📦 BATCH 1: Foundation & The Objective Function Barrier
**Goal:** Modify how the simulation calculates its penalty score so L-BFGS-B respects the 5% error boundary.
**Files to touch:** `calibration/objective_calibration.m`, `solvers/integrate_system.m`
* **Step 1A:** Rewrite `objective_calibration.m` to use Valenti’s normalized squared relative error: `L = sum(w * ((predicted - target)/target)^2)`.
* **Step 1B:** Build the "100x Barrier" into the objective function. Add logic: *if the absolute error of Qp/Qs, MAP, or LVEF > 0.05, multiply their specific weight by 100*.
* **Step 1C:** Tighten ODE solver tolerances inside `integrate_system.m` during calibration (`RelTol = 1e-8`, `AbsTol = 1e-10`) to eliminate numerical noise that will break L-BFGS-B's gradients.

## 📦 BATCH 2: The GSA "Masking" Bridge
**Goal:** Create the mathematical logic that links the initial GSA results to the optimizer by dynamically freezing useless parameters.
**Files to touch:** `utils/create_optimization_mask.m` (New File), `calibration/calibration_param_sets.m`
* **Step 2A:** Create the new `create_optimization_mask.m` function. It will take the Total Sobol Indices ($S_T$) array and a threshold (e.g., `0.1`), outputting a boolean array (e.g., `[1 1 0 0 1]`) indicating which parameters are "active".
* **Step 2B:** Modify `calibration_param_sets.m` to gracefully accept this boolean mask, passing only the *active* upper/lower bounds and initial guesses to the optimizer, rather than the full 11.

## 📦 BATCH 3: The L-BFGS-B Engine (MATLAB `fmincon`)
**Goal:** Gut the old Nelder-Mead algorithm and wire up the gradient-aware L-BFGS-B algorithm.
**Files to touch:** `calibration/run_calibration.m`
* **Step 3A:** Completely rip out `fminsearch`.
* **Step 3B:** Configure `optimoptions('fmincon')` with:
    * `'Algorithm', 'interior-point'`
    * `'HessianApproximation', 'lbfgs'`
    * `'FiniteDifferenceStepSize', 1e-5` *(Crucial to step over the ODE noise from Batch 1)*
* **Step 3C:** Rewrite the `fmincon` call to only optimize the array subset defined by the `mask` from Batch 2, then automatically re-insert those optimized variables back into the full 11-parameter array.

## 📦 BATCH 4: Re-Plumbing the Main Pipeline
**Goal:** Force `main_run.m` to execute sequentially: Phase 0 $\rightarrow$ Initial GSA $\rightarrow$ Optimization $\rightarrow$ Final GSA, rather than using isolated toggles.
**Files to touch:** `main_run.m`
* **Step 4A:** Remove the parallel `DO_CALIBRATION` and `DO_GSA` toggles.
* **Step 4B:** Wire the pipeline directly:
    1. Call `gsa_run_sobol` to get Initial $S_T$.
    2. Call `create_optimization_mask` using that $S_T$.
    3. Call `run_calibration` using the mask.
    4. Call `gsa_run_sobol` *again* using the newly optimized parameters for final reporting.

## 📦 BATCH 5: Validation & Final Output
**Goal:** Ensure the system explicitly proves it hit the < 5% target and saves the final sensitivity data.
**Files to touch:** `utils/validation_report.m`
* **Step 5A:** Update the validation printout table. If any primary metric (Qp/Qs, MAP, LVEF) is `> 5%`, flag it in red/warning text in the MATLAB console.
* **Step 5B:** Ensure the figures/tables export the *Final* GSA results overlaying the initialized Phase 0 parameters.
