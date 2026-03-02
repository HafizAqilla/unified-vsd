function main_run(scenario, clinical)
% MAIN_RUN
% -----------------------------------------------------------------------
% Single entry point for the Unified VSD Lumped-Parameter Model.
%
% USAGE:
%   clinical = patient_template();
%   % ... fill in clinical fields for your patient ...
%   main_run('pre_surgery',  clinical)   % simulate unrepaired VSD
%   main_run('post_surgery', clinical)   % simulate post-closure recovery
%
% INPUTS:
%   scenario  - string: 'pre_surgery' | 'post_surgery'
%   clinical  - unified clinical struct from patient_template.m
%               Must contain .common.weight_kg and .common.height_cm
%               (or .common.BSA) for allometric scaling.
%               All other fields default to NaN if not supplied.
%
% WORKFLOW (executed automatically by this function):
%   1. Apply allometric scaling  →  patient-specific params
%   2. Map clinical SVR/PVR/HR  →  params_from_clinical
%   3. Baseline simulation       →  integrate_system
%   4. Compute baseline metrics  →  compute_clinical_indices
%   5. Calibration (optional)    →  run_calibration
%   6. Post-calibration sim + metrics
%   7. Validation report         →  validation_report
%   8. Plots                     →  plotting_tools
%   9. GSA (optional)            →  gsa_sobol_setup + gsa_run_sobol
%
% SIMULATION TOGGLES (edit below):
%   DO_CALIBRATION  — run fmincon calibration (slow; set false for quick check)
%   DO_GSA          — run Sobol GSA (very slow; set false by default)
%   DO_PLOTS        — generate haemodynamic figures
%
% FILE STRUCTURE:
%   Entry point:    main_run.m            (this file — no physics)
%   Physics:        models/system_rhs.m
%   Solver:         solvers/integrate_system.m
%   Parameters:     config/default_parameters.m
%   Scaling:        utils/apply_scaling.m
%   Calibration:    calibration/run_calibration.m
%   GSA:            gsa/gsa_sobol_setup.m + gsa/gsa_run_sobol.m
%   Validation:     utils/validation_report.m
%
% REFERENCES:
%   [1] Valenti (2023). Full-order 0D cardiovascular model.
%   [2] Lundquist et al. (2025). Allometric scaling.
%   [3] docs/theory_notes.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

%% ---- housekeeping ------------------------------------------------------
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));

%% ---- user toggles (edit here) -----------------------------------------
DO_CALIBRATION = true;
DO_GSA         = false;
DO_PLOTS       = true;

%% ---- validate inputs ---------------------------------------------------
validatestring(scenario, {'pre_surgery', 'post_surgery'}, ...
    'main_run', 'scenario');

if nargin < 2 || isempty(clinical)
    clinical = patient_template();
    warning('main_run:noClinicalData', ...
        'No clinical struct provided; using patient_template() defaults (all NaN).');
end

fprintf('\n[main_run] Scenario: %s\n', scenario);
fprintf('[main_run] Patient: %.1f kg, %.1f cm, age %.2f yr\n', ...
    clinical.common.weight_kg, clinical.common.height_cm, clinical.common.age_years);

%% =====================================================================
%  STEP 1 — Allometric scaling
%% =====================================================================
params_ref  = default_parameters();
patient.age_years  = clinical.common.age_years;
patient.weight_kg  = clinical.common.weight_kg;
patient.height_cm  = clinical.common.height_cm;
patient.sex        = clinical.common.sex;
if isfield(clinical.common, 'BSA') && ~isnan(clinical.common.BSA)
    patient.BSA    = clinical.common.BSA;
end

params0 = apply_scaling(params_ref, patient);

%% =====================================================================
%  STEP 2 — Map clinical measurements (HR, SVR, PVR, R_VSD)
%% =====================================================================
params0 = params_from_clinical(params0, clinical, scenario);

%% =====================================================================
%  STEP 3 — Baseline simulation
%% =====================================================================
fprintf('\n[main_run] Running baseline simulation...\n');
sim_base     = integrate_system(params0);
metrics_base = compute_clinical_indices(sim_base, params0);
fprintf('[main_run] Baseline complete.\n');

%% =====================================================================
%  STEP 4 — Calibration
%% =====================================================================
params_cal  = params0;
metrics_cal = [];
calib_out   = [];

if DO_CALIBRATION
    fprintf('\n[main_run] Running calibration (fmincon)...\n');
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario);

    sim_cal     = integrate_system(params_cal);
    metrics_cal = compute_clinical_indices(sim_cal, params_cal);
    fprintf('[main_run] Calibration complete. Best J = %.6f\n', calib_out.fbest);

    %-- Calibration parameter summary table
    fprintf('\n--- Calibrated parameter changes ---\n');
    param_tbl = table(calib_out.names(:), calib_out.x0(:), calib_out.xbest(:), ...
        (calib_out.xbest(:) - calib_out.x0(:)) ./ abs(calib_out.x0(:)) * 100, ...
        'VariableNames', {'Parameter', 'Initial', 'Calibrated', 'Change_pct'});
    disp(param_tbl);
end

%% =====================================================================
%  STEP 5 — Validation report
%% =====================================================================
report = validation_report(clinical, metrics_base, metrics_cal, scenario);

%% =====================================================================
%  STEP 6 — Plots
%% =====================================================================
if DO_PLOTS
    plotting_tools(sim_base, params0, 'Baseline', scenario);
    if DO_CALIBRATION && ~isempty(metrics_cal)
        plotting_tools(sim_cal, params_cal, 'Calibrated', scenario);
    end
end

%% =====================================================================
%  STEP 7 — Global Sensitivity Analysis (Sobol)
%% =====================================================================
if DO_GSA
    fprintf('\n[main_run] Setting up Sobol GSA (scenario: %s)...\n', scenario);
    params_for_gsa = params_cal;   % use calibrated params as GSA baseline
    gsa_cfg = gsa_sobol_setup(params_for_gsa, scenario);

    fprintf('[main_run] Running Sobol GSA (%d × %d evaluations)...\n', ...
            gsa_cfg.N, numel(gsa_cfg.names) + 2);
    gsa_out = gsa_run_sobol(gsa_cfg, params_for_gsa);

    gsa_summary = make_gsa_summary_table(gsa_out);
    disp(gsa_summary);

    save(fullfile(root, 'results', 'tables', ...
         sprintf('gsa_results_%s.mat', scenario)), 'gsa_out', 'gsa_summary');
    fprintf('[main_run] GSA results saved.\n');
end

%% =====================================================================
%  STEP 8 — Save results
%% =====================================================================
results_dir = fullfile(root, 'results', 'tables');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

if DO_CALIBRATION
    save(fullfile(results_dir, sprintf('params_calibrated_%s.mat', scenario)), ...
         'params_cal', 'calib_out', 'report');
    fprintf('\n[main_run] Calibrated parameters saved to results/tables/\n');
end

fprintf('\n[main_run] Done. Scenario: %s\n', scenario);

end  % main_run
