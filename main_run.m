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
%   7. Final GSA (PCE, post-calib)  →  gsa_pce_setup + gsa_run_pce
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
DO_GSA         = true;    % Run pre- and post-calibration GSA for paper
DO_FAST_CALIBRATION = false; % Full run: Stage 1 + Stage 2 restarts
DO_PARALLEL_FMINCON = true;  % Use parallel finite-difference in fmincon
USE_PCE_IN_CALIBRATION = false; % recovery default: direct ODE calibration (PCE optional)

% MASK_FILE: optional path to a previously saved params_calibrated_*.mat
% that contains an optMask computed from a prior patient's PCE GSA.
% Use this to skip Step 4 (initial PCE GSA) when running a second patient
% of the same scenario and similar physiology — the Sobol mask is assumed
% stable across the cohort and is reused directly.
%
% Set MASK_FILE = '' (empty) to run a full PCE GSA for the new patient.
% Set MASK_FILE = 'results/tables/params_calibrated_pre_surgery.mat' to
% reuse the mask from the most recent saved calibration run.
%
% WARNING: Only reuse a mask from patients of the SAME scenario and roughly
% the same defect physiology. Do NOT reuse a post_surgery mask for pre_surgery.
MASK_FILE = '';   % '' = full GSA; path string = reuse saved mask (irrelevant when DO_GSA=false)

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

do_fmincon_par_env = getenv('UNIFIED_VSD_FMINCON_PARALLEL');
if ~isempty(do_fmincon_par_env)
    DO_PARALLEL_FMINCON = any(strcmpi(strtrim(do_fmincon_par_env), {'1', 'true', 'yes', 'on'}));
end

use_pce_calib_env = getenv('UNIFIED_VSD_USE_PCE_IN_CALIBRATION');
if ~isempty(use_pce_calib_env)
    USE_PCE_IN_CALIBRATION = any(strcmpi(strtrim(use_pce_calib_env), {'1', 'true', 'yes', 'on'}));
end

% Optional parpool control:
%   UNIFIED_VSD_USE_PARPOOL=1 to auto-start parpool if none exists.
%   UNIFIED_VSD_PARPOOL_WORKERS=<N> to request N workers (optional).
use_parpool = false;
parpool_env = getenv('UNIFIED_VSD_USE_PARPOOL');
if ~isempty(parpool_env)
    use_parpool = any(strcmpi(strtrim(parpool_env), {'1', 'true', 'yes', 'on'}));
end

if use_parpool
    pool = gcp('nocreate');
    if isempty(pool)
        workers_env = getenv('UNIFIED_VSD_PARPOOL_WORKERS');
        if ~isempty(workers_env)
            n_workers = str2double(workers_env);
            if ~isnan(n_workers) && isfinite(n_workers) && n_workers > 0
                parpool('local', round(n_workers));
            else
                parpool('local');
            end
        else
            parpool('local');
        end
    end
end

% Keep run_calibration.m parallel setting aligned with this run.
if DO_PARALLEL_FMINCON
    setenv('UNIFIED_VSD_FMINCON_PARALLEL', '1');
else
    setenv('UNIFIED_VSD_FMINCON_PARALLEL', '0');
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
patient.age_days   = clinical.common.age_years * 365.25;
patient.weight_kg  = clinical.common.weight_kg;
patient.height_cm  = clinical.common.height_cm;
patient.sex        = clinical.common.sex;
patient.maturation_mode = 'normal';
if isfield(clinical.common, 'maturation_mode') && ~isempty(clinical.common.maturation_mode)
    patient.maturation_mode = clinical.common.maturation_mode;
end
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
validity_base = evaluate_simulation_validity(sim_base, params0, metrics_base, scenario);
fprintf('[main_run] Baseline complete.\n');

%% =====================================================================
%  STEP 4 — Initial PCE GSA
%% =====================================================================
fprintf('\n=== [Step 4/10] Initial PCE GSA (%.1fs elapsed) ===\n', toc(run_timer));
gsa_init_cfg   = [];
gsa_init_out   = [];
gsa_final_cfg  = [];
gsa_final_out  = [];
use_saved_mask = false;   % default; overridden inside DO_GSA block if MASK_FILE is set

calib_seed = calibration_param_sets(scenario, params0);
calib_names_all = calib_seed.names_all;
n_calib_param  = numel(calib_names_all);

if DO_GSA
    use_saved_mask = ~isempty(MASK_FILE) && isfile(MASK_FILE);
    if use_saved_mask
        fprintf('\n[main_run] MASK_FILE set — loading saved optMask from:\n          %s\n', MASK_FILE);
        saved = load(MASK_FILE, 'optMask', 'sobol_ST_threshold');
        optMask_saved = saved.optMask;
        fprintf('[main_run] Skipping initial PCE GSA (Step 4). Will use saved mask (%d/%d active).\n', ...
                nnz(optMask_saved), numel(optMask_saved));
        gsa_init_cfg = [];
        gsa_init_out = [];
        gsa_pce_out  = [];
    else
        fprintf('\n[main_run] Running initial PCE GSA (pre-calibration)...\n');
        gsa_init_cfg = gsa_pce_setup(params0, scenario);
        gsa_init_out = gsa_run_pce(gsa_init_cfg, params0);
        gsa_pce_out  = gsa_init_out;
    end
end

%% =====================================================================
%  STEP 5 — Create optimisation mask from initial Sobol ST
%% =====================================================================
fprintf('\n=== [Step 5/10] Building optimisation mask (%.1fs elapsed) ===\n', toc(run_timer));

if DO_GSA && use_saved_mask
    % Fast path: reuse mask from a previous patient's PCE run.
    sobol_ST_threshold = saved.sobol_ST_threshold;
    optMask = optMask_saved;
    mask_metrics = {'REUSED_FROM_SAVED_MASK'};
    fprintf('[main_run] Active parameters (reused mask): %d/%d (threshold=%.2f)\n', ...
            nnz(optMask), n_calib_param, sobol_ST_threshold);

elseif DO_GSA
    % Full PCE GSA path: build mask from Sobol ST.
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
    use_saved_mask     = false;
    optMask = true(n_calib_param, 1);
    mask_metrics = {'ALL_ACTIVE_NO_GSA'};
    fprintf('[main_run] GSA disabled. Using full calibration set: %d/%d active parameters.\n', ...
            nnz(optMask), n_calib_param);
end

active_names = calib_names_all(optMask);
fprintf('[main_run] Active set: %s\n', strjoin(active_names, ', '));

[primary_metrics, primary_selection_table] = select_primary_metrics(clinical, gsa_init_out, scenario);
fprintf('[main_run] GSA-guided primary metrics: %s\n', strjoin(primary_metrics, ', '));
disp(primary_selection_table(primary_selection_table.Selected, :));

%% =====================================================================
%  STEP 6 — Calibration (masked)
%% =====================================================================
fprintf('\n=== [Step 6/10] Masked calibration — fmincon (%.1fs elapsed) ===\n', toc(run_timer));
if USE_PCE_IN_CALIBRATION && exist('gsa_pce_out', 'var') && isfield(gsa_pce_out, 'QpQs')
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask, DO_FAST_CALIBRATION, gsa_pce_out, primary_metrics);
else
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask, DO_FAST_CALIBRATION, [], primary_metrics);
end

sim_cal     = integrate_system(params_cal);
metrics_cal = compute_clinical_indices(sim_cal, params_cal);
validity_cal = evaluate_simulation_validity(sim_cal, params_cal, metrics_cal, scenario);
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
    % PCE-based post-calibration GSA (same method as Step 4, applied to
    % calibrated parameters).  Sobol indices are extracted analytically
    % from the trained PCE — no additional Monte Carlo sampling required.
    % Runtime: ~200 ODE runs (same as pre-calibration), not N*(d+2)=5376.
    fprintf('\n=== [Step 7/10] Final PCE GSA on calibrated params (%.1fs elapsed) ===\n', toc(run_timer));
    gsa_final_cfg = gsa_pce_setup(params_cal, scenario);
    gsa_final_out = gsa_run_pce(gsa_final_cfg, params_cal);
else
    fprintf('\n[main_run] Final PCE GSA skipped (DO_GSA=false).\n');
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
    'GsaFinalOut', gsa_final_out, ...
    'PrimaryMetrics', primary_metrics);

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
                          sprintf('gsa_pce_pipeline_%s_%s.mat', scenario, timestamp));

    gsa_save.scenario        = scenario;
    gsa_save.timestamp       = timestamp;
    gsa_save.clinical        = clinical;
    gsa_save.params_baseline = params0;
    gsa_save.params_cal      = params_cal;
    gsa_save.metrics_base    = metrics_base;
    gsa_save.metrics_cal     = metrics_cal;
    gsa_save.validity_base   = validity_base;
    gsa_save.validity_cal    = validity_cal;
    gsa_save.sobol_threshold = sobol_ST_threshold;
    gsa_save.optMask         = optMask;
    gsa_save.primary_metrics = primary_metrics;
    gsa_save.primary_selection_table = primary_selection_table;
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
calib_diag.primary_metrics = primary_metrics;
calib_diag.primary_selection_table = primary_selection_table;
calib_diag.validity_base = validity_base;
calib_diag.validity_cal = validity_cal;
if isfield(calib_out, 'objective_breakdown')
    calib_diag.objective_breakdown = calib_out.objective_breakdown;
end

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
     'params_cal', 'calib_out', 'report', 'optMask', 'sobol_ST_threshold', ...
     'primary_metrics', 'primary_selection_table', 'validity_base', 'validity_cal');
fprintf('\n[main_run] Calibrated parameters saved to results/tables/\n');

if strcmp(scenario, 'pre_surgery')
    % Export a dedicated handoff package so the calibrated pre-op state can
    % be reused later as a seed for post-op experiments.
    pre_to_post_seed = struct();
    pre_to_post_seed.source_scenario = scenario;
    pre_to_post_seed.timestamp = timestamp;
    pre_to_post_seed.clinical = clinical;
    pre_to_post_seed.params_baseline = params0;
    pre_to_post_seed.params_calibrated = params_cal;
    pre_to_post_seed.metrics_baseline = metrics_base;
    pre_to_post_seed.metrics_calibrated = metrics_cal;
    pre_to_post_seed.validity_baseline = validity_base;
    pre_to_post_seed.validity_calibrated = validity_cal;
    pre_to_post_seed.report = report;
    pre_to_post_seed.calib_out = calib_out;
    pre_to_post_seed.optMask = optMask;
    pre_to_post_seed.sobol_ST_threshold = sobol_ST_threshold;
    pre_to_post_seed.primary_metrics = primary_metrics;
    pre_to_post_seed.primary_selection_table = primary_selection_table;

    % Explicit matrix/vector payload for downstream scripts.
    pre_to_post_seed.parameter_names_all = calib_out.names_all(:);
    pre_to_post_seed.parameter_names_active = calib_out.names(:);
    pre_to_post_seed.parameter_vector_initial_all = calib_out.x0_all(:);
    pre_to_post_seed.parameter_vector_calibrated_all = calib_out.xbest_all(:);
    pre_to_post_seed.parameter_vector_initial_active = calib_out.x0(:);
    pre_to_post_seed.parameter_vector_calibrated_active = calib_out.xbest(:);
    pre_to_post_seed.parameter_change_pct_active = ...
        (calib_out.xbest(:) - calib_out.x0(:)) ./ max(abs(calib_out.x0(:)), 1e-9) * 100;

    % Convenience copy: same calibrated pre-op parameters, ready to be
    % loaded and modified for a post-op closure experiment.
    post_seed_params = params_cal;
    if isfield(post_seed_params, 'R') && isfield(post_seed_params.R, 'vsd')
        post_seed_params.R.vsd = 1e6;
    end
    pre_to_post_seed.params_post_seed_closed_vsd = post_seed_params;

    seed_fname_timestamped = fullfile(results_dir, ...
        sprintf('pre_to_post_seed_%s.mat', timestamp));
    seed_fname_latest = fullfile(results_dir, 'pre_to_post_seed_latest.mat');
    save(seed_fname_timestamped, 'pre_to_post_seed');
    save(seed_fname_latest, 'pre_to_post_seed');
    fprintf('[main_run] Pre-to-post seed package saved to:\n          %s\n          %s\n', ...
        seed_fname_timestamped, seed_fname_latest);
end

fprintf('\n[main_run] Done. Scenario: %s  |  Total time: %.1f s\n', scenario, toc(run_timer));

end  % main_run
