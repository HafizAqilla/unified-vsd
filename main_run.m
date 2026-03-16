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
%   4. Initial Sobol GSA         →  gsa_sobol_setup + gsa_run_sobol
%   5. Build optimisation mask   →  create_optimization_mask
%   6. Calibration (masked)      →  run_calibration
%   7. Post-calibration sim + metrics
%   8. Final Sobol GSA           →  gsa_sobol_setup + gsa_run_sobol
%   9. Validation report         →  validation_report
%   10. Plots                    →  plotting_tools

% USER TOGGLE:
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
% DATE:     2026-03-16
% VERSION:  2.0
% -----------------------------------------------------------------------

%% ---- housekeeping ------------------------------------------------------
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));

%% ---- user toggles (edit here) -----------------------------------------
% Plotting only: the core pipeline is always executed sequentially.
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
%  STEP 4 — Initial GSA (Sobol)
%% =====================================================================
fprintf('\n[main_run] Running initial Sobol GSA (pre-calibration)...\n');
gsa_init_cfg = gsa_sobol_setup(params0, scenario);
gsa_init_out = gsa_run_sobol(gsa_init_cfg, params0);

%% =====================================================================
%  STEP 5 — Create optimisation mask from initial Sobol ST
%% =====================================================================
fprintf('\n[main_run] Building optimisation mask from initial Sobol ST...\n');

% Threshold for active parameter screening from Sobol ST.
sobol_ST_threshold = 0.10;   % [dimensionless]

calib_seed = calibration_param_sets(scenario, params0);
calib_names_all = calib_seed.names_all;

primary_metrics = gsa_init_cfg.primary_metrics;
n_calib_param = numel(calib_names_all);
n_primary = numel(primary_metrics);
ST_calib = zeros(n_calib_param, n_primary);

[isFound, idx_in_gsa] = ismember(calib_names_all, gsa_init_cfg.names);
if ~all(isFound)
    missing_names = calib_names_all(~isFound);
    error('main_run:missingGsaMapping', ...
          'Calibration parameters missing in GSA names: %s', strjoin(missing_names, ', '));
end

for m = 1:n_primary
    mf = primary_metrics{m};
    if ~isfield(gsa_init_out, mf) || ~isfield(gsa_init_out.(mf), 'ST')
        error('main_run:missingPrimaryMetricST', ...
              'Primary metric %s is missing ST values in initial GSA output.', mf);
    end
    ST_vec = gsa_init_out.(mf).ST(:);
    ST_calib(:, m) = ST_vec(idx_in_gsa);
end

optMask = create_optimization_mask(ST_calib, sobol_ST_threshold);
fprintf('[main_run] Active parameters after Sobol mask: %d/%d (threshold=%.2f)\n', ...
        nnz(optMask), n_calib_param, sobol_ST_threshold);

active_names = calib_names_all(optMask);
fprintf('[main_run] Active set: %s\n', strjoin(active_names, ', '));

%% =====================================================================
%  STEP 6 — Calibration (masked)
%% =====================================================================
fprintf('\n[main_run] Running masked calibration (fmincon)...\n');
[params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask);

sim_cal     = integrate_system(params_cal);
metrics_cal = compute_clinical_indices(sim_cal, params_cal);
fprintf('[main_run] Calibration complete. Best J = %.6f\n', calib_out.fbest);

% Calibration parameter summary table (active subset only)
fprintf('\n--- Calibrated parameter changes (active subset) ---\n');
param_tbl = table(calib_out.names(:), calib_out.x0(:), calib_out.xbest(:), ...
    (calib_out.xbest(:) - calib_out.x0(:)) ./ abs(calib_out.x0(:)) * 100, ...
    'VariableNames', {'Parameter', 'Initial', 'Calibrated', 'Change_pct'});
disp(param_tbl);

%% =====================================================================
%  STEP 7 — Final GSA (Sobol) with calibrated parameters
%% =====================================================================
fprintf('\n[main_run] Running final Sobol GSA (post-calibration)...\n');
gsa_final_cfg = gsa_sobol_setup(params_cal, scenario);
gsa_final_out = gsa_run_sobol(gsa_final_cfg, params_cal);

%% =====================================================================
%  STEP 8 — Validation report
%% =====================================================================
report = validation_report(clinical, metrics_base, metrics_cal, scenario);

%% =====================================================================
%  STEP 9 — Plots
%% =====================================================================
if DO_PLOTS
    plotting_tools(sim_base, params0, 'Baseline', scenario);
    if ~isempty(metrics_cal)
        plotting_tools(sim_cal, params_cal, 'Calibrated', scenario);
    end
end

%% =====================================================================
%  STEP 10 — Save GSA and calibration artefacts
%% =====================================================================
gsa_results_dir = fullfile(root, 'results', 'gsa');
if ~exist(gsa_results_dir, 'dir'), mkdir(gsa_results_dir); end

timestamp  = datestr(now, 'yyyymmdd_HHMMSS');
gsa_fname  = fullfile(gsa_results_dir, ...
                      sprintf('gsa_sobol_pipeline_%s_%s.mat', scenario, timestamp));

gsa_save.scenario        = scenario;
gsa_save.timestamp       = timestamp;
gsa_save.clinical        = clinical;
gsa_save.params_baseline = params0;
gsa_save.params_cal      = params_cal;
gsa_save.metrics_base    = metrics_base;
gsa_save.metrics_cal     = metrics_cal;
gsa_save.sobol_threshold = sobol_ST_threshold;
gsa_save.optMask         = optMask;
gsa_save.calib_names_all = calib_names_all;
gsa_save.gsa_init_cfg    = gsa_init_cfg;
gsa_save.gsa_init_out    = gsa_init_out;
gsa_save.gsa_final_cfg   = gsa_final_cfg;
gsa_save.gsa_final_out   = gsa_final_out;
gsa_save.calib_out       = calib_out;

save(gsa_fname, '-struct', 'gsa_save');
fprintf('[main_run] Pipeline artefacts saved to:\n          %s\n', gsa_fname);

% Display summary table for initial and final GSA
fprintf('\n[main_run] Initial Sobol summary table:\n');
gsa_summary_init = make_gsa_summary_table(gsa_init_out);
disp(gsa_summary_init);

fprintf('\n[main_run] Final Sobol summary table:\n');
gsa_summary_final = make_gsa_summary_table(gsa_final_out);
disp(gsa_summary_final);

% Matrix-style table for final sensitivities
make_gsa_matrix_table(gsa_final_out, sobol_ST_threshold, true);

%% =====================================================================
%  STEP 11 — Save calibrated parameter result package
%% =====================================================================
results_dir = fullfile(root, 'results', 'tables');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

save(fullfile(results_dir, sprintf('params_calibrated_%s.mat', scenario)), ...
     'params_cal', 'calib_out', 'report', 'optMask', 'sobol_ST_threshold');
fprintf('\n[main_run] Calibrated parameters saved to results/tables/\n');

fprintf('\n[main_run] Done. Scenario: %s\n', scenario);

end  % main_run
