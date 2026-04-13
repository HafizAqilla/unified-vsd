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
%   9. GSA (optional)            →  gsa_pce_setup   + gsa_run_pce  (PCE via UQLab)
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
%   GSA:            gsa/gsa_pce_setup.m   + gsa/gsa_run_pce.m   (PCE via UQLab/SoBioS)
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
% DO_CALIBRATION: set true only when you want fmincon/MultiStart tuning.
% Leave false for a quick baseline simulation (calibration takes ~5–30 min).
DO_CALIBRATION = true;
DO_GSA         = false;
DO_PLOTS       = true;

%% ---- validate inputs ---------------------------------------------------
if nargin < 1 || isempty(scenario)
    scenario = 'pre_surgery';
    fprintf('[main_run] No scenario specified — defaulting to "%s".\n', scenario);
end

if nargin < 2 || isempty(clinical)
    clinical = patient_template();
    warning('main_run:noClinicalData', ...
        'No clinical struct provided; using patient_template() defaults (all NaN).');
end

validatestring(scenario, {'pre_surgery', 'post_surgery'}, ...
    'main_run', 'scenario');

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
%  STEP 7 — Global Sensitivity Analysis (PCE via UQLab / SoBioS)
%% =====================================================================
if DO_GSA
    % ---- GSA Paths (update these to your actual install locations) ----
    UQLAB_PATH  = 'C:\Users\asus\Documents\MATLAB\VSD\Toolbox\UQLab_Rel2.2.0\core';
    SOBIOS_PATH = 'C:\Users\asus\Documents\MATLAB\VSD\Toolbox\americocunhajr-SoBioS-426756c';
    % -------------------------------------------------------------------

    fprintf('\n[main_run] Setting up PCE GSA (scenario: %s)...\n', scenario);
    params_for_gsa = params_cal;   % use calibrated params as PCE training baseline
    gsa_cfg = gsa_pce_setup(params_for_gsa, scenario, UQLAB_PATH, SOBIOS_PATH);

    fprintf('[main_run] Running PCE GSA (%d training evaluations)...\n', gsa_cfg.PCEOpts.ExpDesign.NSamples);
    gsa_out = gsa_run_pce(gsa_cfg, params_for_gsa);

    %% -----------------------------------------------------------------
    %  SAVE ALL GSA DATA immediately after computation
    %  (saved BEFORE display so data is never lost if display crashes)
    %  File: results/gsa/gsa_pce_<scenario>_<YYYYMMDD_HHMMSS>.mat
    %% -----------------------------------------------------------------
    gsa_results_dir = fullfile(root, 'results', 'gsa');
    if ~exist(gsa_results_dir, 'dir'), mkdir(gsa_results_dir); end

    timestamp  = datestr(now, 'yyyymmdd_HHMMSS');
    gsa_fname  = fullfile(gsa_results_dir, ...
                          sprintf('gsa_pce_%s_%s.mat', scenario, timestamp));

    % Package everything into a single struct for clean file management
    gsa_save.scenario        = scenario;
    gsa_save.timestamp       = timestamp;
    gsa_save.clinical        = clinical;         % patient clinical inputs
    gsa_save.params_baseline = params0;          % allometric-scaled params
    gsa_save.params_gsa      = params_for_gsa;   % params used as PCE centre
    gsa_save.metrics_base    = metrics_base;     % baseline haemodynamic metrics
    gsa_save.gsa_cfg         = gsa_cfg;          % PCE setup config
    gsa_save.gsa_out         = gsa_out;          % per-metric S1 / ST / tables

    % Optionally include calibration artefacts if calibration was run
    if DO_CALIBRATION && ~isempty(metrics_cal)
        gsa_save.calib_out   = calib_out;
        gsa_save.metrics_cal = metrics_cal;
    end

    save(gsa_fname, '-struct', 'gsa_save');
    fprintf('[main_run] All GSA data saved to:\n          %s\n', gsa_fname);

    %% -----------------------------------------------------------------
    %  DISPLAY summary table (non-critical — crash here won't lose data)
    %% -----------------------------------------------------------------
    gsa_summary = make_gsa_summary_table(gsa_out);
    disp(gsa_summary);

    % Matrix-style table: parameters × metrics with Sobol ST values
    make_gsa_matrix_table(gsa_out, 0.1, true);   % yellow threshold = 0.1, save PNG

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