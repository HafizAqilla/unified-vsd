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
%   4. Initial GSA (PCE surrogate) →  gsa_pce_setup + gsa_run_pce
%   5. Build optimisation mask   →  create_optimization_mask
%   6. Calibration (masked)      →  run_calibration
%   7. Final GSA (direct Sobol)     →  gsa_sobol_setup + gsa_run_sobol
%   8. Validation report            →  validation_report
%   9. Plots                        →  plotting_tools
%   10. Save artefacts              →  results/*.mat outputs

% USER TOGGLE:
%   DO_PLOTS        — generate haemodynamic figures
%   DO_OVERLAY      — baseline vs calibrated overlap curves in one canvas
%   DO_GSA          — run Sobol screening + post-calibration Sobol
%   GSA_SOBOL_N     — optional override for Sobol base sample N
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
root_paths = strsplit(genpath(root), pathsep);
is_shadow = contains(root_paths, [filesep '.claude' filesep]);
addpath(strjoin(root_paths(~is_shadow), pathsep));

%% ---- user toggles (edit here) -----------------------------------------
% Plotting only: the core pipeline is always executed sequentially.
DO_PLOTS       = true;
DO_OVERLAY     = true;
DO_GSA         = false;   % Skip GSA by default; focus on calibration
DO_FAST_CALIBRATION = true; % Fast test run — set false for final publication run
GSA_SOBOL_N    = 16;  % Fast end-to-end runtime for iterative runs; [] uses default/ENV

% Optional runtime override from environment variable:
%   UNIFIED_VSD_DO_PLOTS=0|false|off -> disable plotting
%   UNIFIED_VSD_DO_PLOTS=1|true|on   -> enable plotting
do_plots_env = getenv('UNIFIED_VSD_DO_PLOTS');
if ~isempty(do_plots_env)
    DO_PLOTS = any(strcmpi(strtrim(do_plots_env), {'1', 'true', 'yes', 'on'}));
end

% Optional runtime override from environment variable:
%   UNIFIED_VSD_DO_OVERLAY=0|false|off -> disable overlay
%   UNIFIED_VSD_DO_OVERLAY=1|true|on   -> enable overlay
do_overlay_env = getenv('UNIFIED_VSD_DO_OVERLAY');
if ~isempty(do_overlay_env)
    DO_OVERLAY = any(strcmpi(strtrim(do_overlay_env), {'1', 'true', 'yes', 'on'}));
end

% Optional runtime override from environment variable:
%   UNIFIED_VSD_DO_GSA=0|false|off  -> disable GSA
%   UNIFIED_VSD_DO_GSA=1|true|on    -> enable GSA
do_gsa_env = getenv('UNIFIED_VSD_DO_GSA');
if ~isempty(do_gsa_env)
    DO_GSA = any(strcmpi(strtrim(do_gsa_env), {'1', 'true', 'yes', 'on'}));
end

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

run_timer = tic;   % wall-clock timer for the full pipeline

%% =====================================================================
%  STEP 1 — Allometric scaling
%% =====================================================================
fprintf('\n=== [Step 1/10] Allometric scaling ===\n');
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
fprintf('\n=== [Step 2/10] Mapping clinical measurements (%.1fs elapsed) ===\n', toc(run_timer));
params0 = params_from_clinical(params0, clinical, scenario);

%% =====================================================================
%  STEP 3 — Baseline simulation
%% =====================================================================
fprintf('\n=== [Step 3/10] Baseline simulation (%.1fs elapsed) ===\n', toc(run_timer));
sim_base     = integrate_system(params0);
metrics_base = compute_clinical_indices(sim_base, params0);
fprintf('[main_run] Baseline complete.\n');

%% =====================================================================
%  STEP 4 — Initial PCE GSA
%% =====================================================================
fprintf('\n=== [Step 4/10] Initial PCE GSA (%.1fs elapsed) ===\n', toc(run_timer));
gsa_init_cfg = [];
gsa_init_out = [];
gsa_final_cfg = [];
gsa_final_out = [];

calib_seed = calibration_param_sets(scenario, params0);
calib_names_all = calib_seed.names_all;
n_calib_param  = numel(calib_names_all);

if DO_GSA
    fprintf('\n[main_run] Running initial PCE GSA (pre-calibration)...\n');
    gsa_init_cfg = gsa_pce_setup(params0, scenario);
    gsa_init_out = gsa_run_pce(gsa_init_cfg, params0);
    gsa_pce_out = gsa_init_out;
end

%% =====================================================================
%  STEP 5 — Create optimisation mask from initial Sobol ST
%% =====================================================================
fprintf('\n=== [Step 5/10] Building optimisation mask (%.1fs elapsed) ===\n', toc(run_timer));

if DO_GSA
    % Threshold for active parameter screening from Sobol ST.
    sobol_ST_threshold = 0.10;   % [dimensionless]

    % Build mask from primary AND secondary metrics (conservative union):
    %   Primary-only masking screened out V0.LV / V0.RV because their ST for
    %   QpQs/PAP/PVR is low, even though they directly control LVEDV/RVEDV.
    %   Including secondary metrics ensures volume-relevant parameters remain
    %   active without lowering the threshold globally.
    mask_metrics = unique([gsa_init_cfg.primary_metrics, gsa_init_cfg.secondary_metrics], 'stable');
    n_mask_metrics = numel(mask_metrics);
    ST_calib = zeros(n_calib_param, n_mask_metrics);

    [isFound, idx_in_gsa] = ismember(calib_names_all, gsa_init_cfg.names);
    if ~all(isFound)
        missing_names = calib_names_all(~isFound);
        error('main_run:missingGsaMapping', ...
              'Calibration parameters missing in GSA names: %s', strjoin(missing_names, ', '));
    end

    for m = 1:n_mask_metrics
        mf = mask_metrics{m};
        if ~isfield(gsa_init_out, mf) || ~isfield(gsa_init_out.(mf), 'ST')
            warning('main_run:missingMaskMetricST', ...
                    'Mask metric %s missing ST values in GSA output; column skipped.', mf);
            continue;
        end
        ST_vec = gsa_init_out.(mf).ST(:);
        ST_calib(:, m) = ST_vec(idx_in_gsa);
    end

    optMask = create_optimization_mask(ST_calib, sobol_ST_threshold);
    fprintf('[main_run] Active parameters after Sobol mask: %d/%d (threshold=%.2f, metrics=%s)\n', ...
            nnz(optMask), n_calib_param, sobol_ST_threshold, strjoin(mask_metrics, '+'));
else
    sobol_ST_threshold = NaN;
    optMask = true(n_calib_param, 1);
    mask_metrics = {'ALL_ACTIVE_NO_GSA'};
    fprintf('[main_run] GSA disabled. Using full calibration set: %d/%d active parameters.\n', ...
            nnz(optMask), n_calib_param);
end

active_names = calib_names_all(optMask);
fprintf('[main_run] Active set: %s\n', strjoin(active_names, ', '));

%% =====================================================================
%  STEP 6 — Calibration (masked)
%% =====================================================================
fprintf('\n=== [Step 6/10] Masked calibration — fmincon (%.1fs elapsed) ===\n', toc(run_timer));
if exist('gsa_pce_out', 'var') && isfield(gsa_pce_out, 'QpQs')
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask, DO_FAST_CALIBRATION, gsa_pce_out);
else
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask, DO_FAST_CALIBRATION, []);
end

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
if DO_GSA
    fprintf('\n=== [Step 7/10] Final Sobol GSA (%.1fs elapsed) ===\n', toc(run_timer));
    gsa_final_cfg = gsa_sobol_setup(params_cal, scenario, GSA_SOBOL_N);
    gsa_final_out = gsa_run_sobol(gsa_final_cfg, params_cal);
else
    fprintf('\n[main_run] Final Sobol GSA skipped (DO_GSA=false).\n');
end

%% =====================================================================
%  STEP 8 — Validation report
%% =====================================================================
fprintf('\n=== [Step 8/10] Validation report (%.1fs elapsed) ===\n', toc(run_timer));
results_dir = fullfile(root, 'results', 'tables');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

report = validation_report( ...
    clinical, metrics_base, metrics_cal, scenario, ...
    'ResultsDir', results_dir, ...
    'GsaInitOut', gsa_init_out, ...
    'GsaFinalOut', gsa_final_out);

%% =====================================================================
%  STEP 9 — Plots
%% =====================================================================
fprintf('\n=== [Step 9/10] Plots (%.1fs elapsed) ===\n', toc(run_timer));
if DO_PLOTS
    plotting_tools(sim_base, params0, 'Baseline', scenario);
    if ~isempty(metrics_cal)
        plotting_tools(sim_cal, params_cal, 'Calibrated', scenario);
        if DO_OVERLAY
            plot_overlay_comparison(sim_base, params0, 'Baseline', ...
                sim_cal, params_cal, 'Calibrated', scenario);
        end
    end
end

%% =====================================================================
%  STEP 10 — Save artefacts
%% =====================================================================
fprintf('\n=== [Step 10/10] Saving artefacts (%.1fs elapsed) ===\n', toc(run_timer));
timestamp  = datestr(now, 'yyyymmdd_HHMMSS');
if DO_GSA
    gsa_results_dir = fullfile(root, 'results', 'gsa');
    if ~exist(gsa_results_dir, 'dir'), mkdir(gsa_results_dir); end

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
else
    fprintf('[main_run] GSA pipeline artefact skipped (DO_GSA=false).\n');
end

%% Compact calibration diagnostics artifact
% A small, fast-loading file for run-to-run comparison and root-cause analysis.
% Contains objective trace, multi-start summary, metric errors, and parameter changes.
calib_diag = struct();
calib_diag.scenario      = scenario;
calib_diag.timestamp     = timestamp;

% Objective trace
calib_diag.J0            = calib_out.J0;
calib_diag.fbest         = calib_out.fbest;
calib_diag.improvement   = calib_out.improvement;
calib_diag.best_stage    = calib_out.best_stage;
calib_diag.best_restart  = calib_out.best_restart;
if isfield(calib_out, 'restart_J'); calib_diag.restart_J = calib_out.restart_J; end
if isfield(calib_out, 'restart_flag'); calib_diag.restart_flag = calib_out.restart_flag; end

% Active parameter changes
calib_diag.active_params = calib_out.names(:);
calib_diag.x0            = calib_out.x0;
calib_diag.xbest         = calib_out.xbest;
calib_diag.change_pct    = (calib_out.xbest - calib_out.x0) ./ ...
                            max(abs(calib_out.x0), 1e-9) * 100;

% Validation summary
calib_diag.rmse_baseline = report.rmse_baseline;
calib_diag.rmse_cal      = report.rmse_cal;
calib_diag.table_delta   = report.table_delta;
calib_diag.sorted_errors = report.sorted_errors;
calib_diag.primary_gate  = report.primary_gate;

% Mask configuration used
calib_diag.sobol_threshold = sobol_ST_threshold;
calib_diag.mask_metrics    = mask_metrics;
calib_diag.optMask         = optMask;

diag_fname = fullfile(results_dir, ...
    sprintf('calib_diagnostics_%s_%s.mat', scenario, timestamp));
save(diag_fname, 'calib_diag');
fprintf('[main_run] Calibration diagnostics saved to:\n          %s\n', diag_fname);

if DO_GSA
    % Display summary table for initial and final GSA
    fprintf('\n[main_run] Initial Sobol summary table:\n');
    gsa_summary_init = make_gsa_summary_table(gsa_init_out);
    disp(gsa_summary_init);

    fprintf('\n[main_run] Final Sobol summary table:\n');
    gsa_summary_final = make_gsa_summary_table(gsa_final_out);
    disp(gsa_summary_final);

    % Matrix-style table for final sensitivities
    make_gsa_matrix_table(gsa_final_out, sobol_ST_threshold, true);
end

%% =====================================================================
%  STEP 11 — Save calibrated parameter result package
%% =====================================================================
save(fullfile(results_dir, sprintf('params_calibrated_%s.mat', scenario)), ...
     'params_cal', 'calib_out', 'report', 'optMask', 'sobol_ST_threshold');
fprintf('\n[main_run] Calibrated parameters saved to results/tables/\n');

fprintf('\n[main_run] Done. Scenario: %s  |  Total time: %.1f s\n', scenario, toc(run_timer));

end  % main_run


