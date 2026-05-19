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
canonical_root = resolve_canonical_project_root(root);
restoredefaultpath();
remove_shadow_project_paths(canonical_root);
cd(canonical_root);
addpath(build_clean_project_path(canonical_root));

%% ---- user toggles (edit here) -----------------------------------------
% Plotting only: the core pipeline is always executed sequentially.
DO_PLOTS       = true;
DO_OVERLAY     = true;
DO_GSA         = true;    % Run pre- and post-calibration GSA for paper
DO_FAST_CALIBRATION = false; % Full run: Stage 1 + Stage 2 restarts
DO_PARALLEL_FMINCON = false; % Serial by default; parallel remains opt-in via env
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

fast_calibration_env = getenv('UNIFIED_VSD_FAST_CALIBRATION');
if ~isempty(fast_calibration_env)
    DO_FAST_CALIBRATION = any(strcmpi(strtrim(fast_calibration_env), {'1', 'true', 'yes', 'on'}));
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
case_profile = build_case_calibration_profile(clinical, scenario);
fprintf('[main_run] Calibration case mode: %s (%s)\n', ...
    case_profile.mode, case_profile.description);
if isfield(case_profile, 'clinicalConsistencyAudit') && ...
        isfield(case_profile.clinicalConsistencyAudit, 'summary')
    fprintf('[main_run] %s\n', case_profile.clinicalConsistencyAudit.summary);
end

run_ctx = init_run_output(canonical_root, scenario, clinical);
run_cleanup = onCleanup(@() cleanup_run_output(run_ctx)); %#ok<NASGU>
fprintf('[main_run] Run folder: %s\n', run_ctx.root);

run_timer = tic;   % wall-clock timer for the full pipeline

%% =====================================================================
%  STEP 1 — Allometric scaling
%% =====================================================================
fprintf('\n=== [Step 1/10] Allometric scaling ===\n');
params_ref  = default_parameters();
baseline_provenance = build_baseline_provenance(params_ref);
patient.age_years  = clinical.common.age_years;
patient.age_days   = clinical.common.age_years * 365.25;
patient.weight_kg  = clinical.common.weight_kg;
patient.height_cm  = clinical.common.height_cm;
patient.sex        = clinical.common.sex;
patient.maturation_mode = 'normal';
if isfield(clinical.common, 'maturation_mode') && ~isempty(clinical.common.maturation_mode)
    patient.maturation_mode = clinical.common.maturation_mode;
end
patient.scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
if isfield(clinical.common, 'scaling_mode') && ~isempty(clinical.common.scaling_mode)
    patient.scaling_mode = clinical.common.scaling_mode;
end
if isempty(patient.scaling_mode)
    if isfield(case_profile, 'preferredScalingMode') && ~isempty(case_profile.preferredScalingMode)
        patient.scaling_mode = case_profile.preferredScalingMode;
    else
        patient.scaling_mode = 'lundquist_bsa';
    end
end
if isfield(clinical.common, 'BSA') && ~isempty(clinical.common.BSA) && isfinite(clinical.common.BSA)
    patient.BSA    = clinical.common.BSA;
end

params_scaled = apply_scaling(params_ref, patient);
scaling_comparison = compare_scaling_methods(params_ref, patient);
age_validity = build_age_validity_annotation( ...
    patient.age_years, patient.scaling_mode, patient.maturation_mode);
registry_context = struct('params_adult', params_ref, 'params_scaled', params_scaled);
params0 = params_scaled;

%% =====================================================================
%  STEP 1b - Warm-start post-surgery runs from calibrated pre-op physiology
%% =====================================================================
[params0, post_warm_start] = apply_post_surgery_warm_start(params0, clinical, scenario);
registry_context.post_warm_start = post_warm_start;
params_reference_for_clinical = params_scaled;
if post_warm_start.applied
    params_reference_for_clinical = params0;
    fprintf('[main_run] Post-surgery warm-start applied from %s.\n', ...
        post_warm_start.source);
    fprintf('[main_run] Warm-start copied %d parameter fields; R.vsd remains scenario-controlled.\n', ...
        numel(post_warm_start.copiedFields));
elseif strcmp(scenario, 'post_surgery')
    fprintf('[main_run] No pre-surgery CalibParams found; post-surgery starts from demographic scaling.\n');
end

%% =====================================================================
%  STEP 2 — Map clinical measurements (HR, SVR, PVR, R_VSD)
%% =====================================================================
fprintf('\n=== [Step 2/10] Mapping clinical measurements (%.1fs elapsed) ===\n', toc(run_timer));
params0 = params_from_clinical(params0, clinical, scenario, params_reference_for_clinical, case_profile);

%% =====================================================================
%  STEP 3 — Baseline simulation
%% =====================================================================
fprintf('\n=== [Step 3/10] Baseline simulation (%.1fs elapsed) ===\n', toc(run_timer));
sim_base     = integrate_system(params0);
metrics_base = compute_clinical_indices(sim_base, params0);
validity_base = evaluate_simulation_validity(sim_base, params0, metrics_base, scenario, clinical);
fprintf('[main_run] Baseline complete.\n');
print_systemic_consistency_summary(metrics_base, clinical, scenario, 'Baseline');

%% =====================================================================
%  STEP 4 — Initial PCE GSA
%% =====================================================================
fprintf('\n=== [Step 4/10] Initial PCE GSA (%.1fs elapsed) ===\n', toc(run_timer));
gsa_init_cfg   = [];
gsa_init_out   = [];
gsa_final_cfg  = [];
gsa_final_out  = [];
use_saved_mask = false;   % default; overridden inside DO_GSA block if MASK_FILE is set

calib_seed = calibration_param_sets(scenario, params0, [], {}, case_profile, registry_context);
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
    optMask = create_optimization_mask(double(optMask_saved(:)), 0.50, calib_names_all, case_profile);
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
    mask_metrics = restrict_metrics_by_case_profile(mask_metrics, case_profile);
    n_mask_metrics = numel(mask_metrics);
    ST_calib = zeros(n_calib_param, n_mask_metrics);

    for m = 1:n_mask_metrics
        mf = mask_metrics{m};
        if ~isfield(gsa_init_out, mf) || ~isfield(gsa_init_out.(mf), 'ST')
            warning('main_run:missingMaskMetricST', ...
                    'Mask metric %s missing ST values in GSA output; column skipped.', mf);
            continue;
        end
        ST_vec = gsa_init_out.(mf).ST(:);
        ST_calib(:, m) = map_gsa_sensitivities_to_calibration_names( ...
            ST_vec, gsa_init_cfg.names, calib_names_all, case_profile);
    end

    optMask = create_optimization_mask(ST_calib, sobol_ST_threshold, calib_names_all, case_profile);
    fprintf('[main_run] Active parameters after Sobol mask: %d/%d (threshold=%.2f, metrics=%s)\n', ...
            nnz(optMask), n_calib_param, sobol_ST_threshold, strjoin(mask_metrics, '+'));

else
    sobol_ST_threshold = NaN;
    use_saved_mask     = false;
    optMask = create_optimization_mask(ones(n_calib_param, 1), 0.50, calib_names_all, case_profile);
    mask_metrics = {'ALL_ACTIVE_NO_GSA'};
    fprintf('[main_run] GSA disabled. Using full calibration set: %d/%d active parameters.\n', ...
            nnz(optMask), n_calib_param);
end

active_names = calib_names_all(optMask);
fprintf('[main_run] Active set: %s\n', strjoin(active_names, ', '));

[primary_metrics, primary_selection_table] = select_primary_metrics(clinical, gsa_init_out, scenario, case_profile);
fprintf('[main_run] Case-profile primary metrics: %s\n', strjoin(primary_metrics, ', '));
disp(primary_selection_table(primary_selection_table.Selected, :));

%% =====================================================================
%  STEP 6 — Calibration (masked)
%% =====================================================================
fprintf('\n=== [Step 6/10] Masked calibration — fmincon (%.1fs elapsed) ===\n', toc(run_timer));
if USE_PCE_IN_CALIBRATION && exist('gsa_pce_out', 'var') && isfield(gsa_pce_out, 'QpQs')
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask, DO_FAST_CALIBRATION, gsa_pce_out, primary_metrics, case_profile, registry_context);
else
    [params_cal, calib_out] = run_calibration(params0, clinical, scenario, optMask, DO_FAST_CALIBRATION, [], primary_metrics, case_profile, registry_context);
end

sim_cal     = integrate_system(params_cal);
metrics_cal = compute_clinical_indices(sim_cal, params_cal);
validity_cal = evaluate_simulation_validity(sim_cal, params_cal, metrics_cal, scenario, clinical);
fprintf('[main_run] Calibration complete. Best J = %.6f\n', calib_out.fbest);
print_systemic_consistency_summary(metrics_cal, clinical, scenario, 'Calibrated');

% Calibration parameter summary table (active subset only)
fprintf('\n--- Calibrated parameter changes (active subset) ---\n');
param_tbl = table(calib_out.names(:), calib_out.x0(:), calib_out.xbest(:), ...
    (calib_out.xbest(:) - calib_out.x0(:)) ./ abs(calib_out.x0(:)) * 100, ...
    'VariableNames', {'Parameter', 'Initial', 'Calibrated', 'Change_pct'});
disp(param_tbl);
if isfield(calib_out, 'parameterPlausibility') && isfield(calib_out.parameterPlausibility, 'table')
    plaus_tbl = calib_out.parameterPlausibility.table;
    fprintf('\n--- Parameter plausibility (active subset) ---\n');
    disp(plaus_tbl(:, {'Parameter','BaselineScaled','FittedValue', ...
        'RatioToBaseline','WithinBounds','PlausibilityFlag'}));
    fprintf('[main_run] Plausibility summary: %d OK | %d WARNING | %d FAIL\n', ...
        calib_out.parameterPlausibility.n_ok, ...
        calib_out.parameterPlausibility.n_warning, ...
        calib_out.parameterPlausibility.n_fail);
end

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
results_dir = run_ctx.tables_dir;

report = validation_report( ...
    clinical, metrics_base, metrics_cal, scenario, ...
    'ResultsDir', results_dir, ...
    'GsaInitOut', gsa_init_out, ...
    'GsaFinalOut', gsa_final_out, ...
    'PrimaryMetrics', primary_metrics, ...
    'ValidationHoldoutMetrics', profile_holdout_metrics(case_profile), ...
    'TargetTiers', profile_target_tiers(case_profile), ...
    'ClinicalConsistencyAudit', profile_clinical_audit(case_profile));
calibration_status = classify_calibration_run(report, calib_out.parameterPlausibility);
fprintf('[main_run] Calibration status: %s\n', calibration_status.label);
fprintf('[main_run] Status summary: %s\n', calibration_status.summary);
best_candidate = build_candidate_snapshot( ...
    'best_candidate', params_cal, sim_cal, metrics_cal, validity_cal, ...
    report, calibration_status, calib_out.parameterPlausibility, false, '');
scientific_candidate = build_candidate_snapshot( ...
    'scientific_candidate', params_cal, sim_cal, metrics_cal, validity_cal, ...
    report, calibration_status, calib_out.parameterPlausibility, false, '');
if should_rollback_calibration(report, calibration_status)
    rollback_reason = sprintf('status=%s, RMSE baseline=%.4f, calibrated=%.4f', ...
        calibration_status.label, report.rmse_baseline, report.rmse_cal);
    fprintf('[main_run] Calibration rollback applied: %s\n', rollback_reason);
    calib_out.rollback_applied = true;
    calib_out.rollback_reason = rollback_reason;
    accepted_report = validation_report( ...
        clinical, metrics_base, metrics_base, scenario, ...
        'ResultsDir', results_dir, ...
        'GsaInitOut', gsa_init_out, ...
        'GsaFinalOut', gsa_final_out, ...
        'PrimaryMetrics', primary_metrics, ...
        'ValidationHoldoutMetrics', profile_holdout_metrics(case_profile), ...
        'TargetTiers', profile_target_tiers(case_profile), ...
        'ClinicalConsistencyAudit', profile_clinical_audit(case_profile));
    accepted_plausibility = struct();
    if isfield(calib_out, 'x0_active') && isfield(calib_out, 'parameterRegistryActive')
        accepted_plausibility = evaluate_parameter_plausibility( ...
            calib_out.x0_active, calib_out.parameterRegistryActive);
    end
    accepted_status = classify_calibration_run(accepted_report, accepted_plausibility);
    accepted_status.rollback_applied = true;
    accepted_status.rollback_reason = rollback_reason;
    fprintf('[main_run] Post-rollback status: %s\n', accepted_status.label);
    accepted_candidate = build_candidate_snapshot( ...
        'accepted_candidate', params0, sim_base, metrics_base, validity_base, ...
        accepted_report, accepted_status, accepted_plausibility, true, rollback_reason);
else
    calib_out.rollback_applied = false;
    calib_out.rollback_reason = '';
    accepted_candidate = build_candidate_snapshot( ...
        'accepted_candidate', params_cal, sim_cal, metrics_cal, validity_cal, ...
        report, calibration_status, calib_out.parameterPlausibility, false, '');
end

params_cal = accepted_candidate.params;
sim_cal = accepted_candidate.sim;
metrics_cal = accepted_candidate.metrics;
validity_cal = accepted_candidate.validity;
report = accepted_candidate.report;
calibration_status = accepted_candidate.calibration_status;

%% =====================================================================
%  STEP 9 — Plots
%% =====================================================================
fprintf('\n=== [Step 9/10] Plots (%.1fs elapsed) ===\n', toc(run_timer));
if DO_PLOTS
    plotting_tools(sim_base, params0, 'Baseline', scenario, ...
        'ResultsDir', run_ctx.figures_dir);
    if ~isempty(metrics_cal)
        plotting_tools(sim_cal, params_cal, 'Calibrated', scenario, ...
            'ResultsDir', run_ctx.figures_dir);
        if DO_OVERLAY
            plot_overlay_comparison(sim_base, params0, 'Baseline', ...
                sim_cal, params_cal, 'Calibrated', scenario, ...
                'ResultsDir', run_ctx.figures_dir);
        end
    end
end

%% =====================================================================
%  STEP 10 — Save artefacts
%% =====================================================================
fprintf('\n=== [Step 10/10] Saving artefacts (%.1fs elapsed) ===\n', toc(run_timer));
timestamp  = run_ctx.timestamp;
if DO_GSA
    gsa_results_dir = run_ctx.gsa_dir;

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
    gsa_save.case_profile    = case_profile;
    gsa_save.age_validity    = age_validity;
    gsa_save.scaling_comparison = scaling_comparison;
    gsa_save.primary_selection_table = primary_selection_table;
    gsa_save.calib_names_all = calib_names_all;
    gsa_save.gsa_init_cfg    = gsa_init_cfg;
    gsa_save.gsa_init_out    = gsa_init_out;
    gsa_save.gsa_final_cfg   = gsa_final_cfg;
    gsa_save.gsa_final_out   = gsa_final_out;
    gsa_save.calib_out       = calib_out;
    gsa_save.best_candidate  = best_candidate;
    gsa_save.accepted_candidate = accepted_candidate;

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
calib_diag.case_profile  = case_profile;

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
if isfield(report, 'rmse_summary')
    calib_diag.rmse_summary = report.rmse_summary;
end
if isfield(report, 'target_tiers')
    calib_diag.target_tiers = report.target_tiers;
end
if isfield(report, 'clinical_consistency_audit')
    calib_diag.clinical_consistency_audit = report.clinical_consistency_audit;
end
calib_diag.table_delta   = report.table_delta;
calib_diag.sorted_errors = report.sorted_errors;
calib_diag.primary_gate  = report.primary_gate;
calib_diag.primary_metrics = primary_metrics;
calib_diag.primary_selection_table = primary_selection_table;
calib_diag.validity_base = validity_base;
calib_diag.validity_cal = validity_cal;
calib_diag.calibration_status = calibration_status;
calib_diag.best_candidate = summarize_candidate_snapshot(best_candidate);
calib_diag.accepted_candidate = summarize_candidate_snapshot(accepted_candidate);
if isfield(calib_out, 'parameterPlausibility')
    calib_diag.parameter_plausibility = calib_out.parameterPlausibility;
end
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
    make_gsa_matrix_table(gsa_final_out, sobol_ST_threshold, true, ...
        'ResultsDir', run_ctx.gsa_dir);
end

%% =====================================================================
%  STEP 11 — Save calibrated parameter result package
%% =====================================================================
params_package_file = fullfile(run_ctx.mat_dir, sprintf('params_calibrated_%s.mat', scenario));
save(params_package_file, 'params_cal', 'calib_out', 'report', 'optMask', ...
     'sobol_ST_threshold', 'primary_metrics', 'primary_selection_table', ...
     'calibration_status', ...
    'case_profile', 'validity_base', 'validity_cal', 'best_candidate', ...
     'scientific_candidate', ...
     'accepted_candidate', 'age_validity', 'baseline_provenance', ...
     'scaling_comparison');
fprintf('\n[main_run] Calibrated parameters saved to:\n          %s\n', params_package_file);

best_candidate_file = fullfile(run_ctx.mat_dir, sprintf('params_best_candidate_%s.mat', scenario));
save(best_candidate_file, 'best_candidate', 'calib_out', 'optMask', ...
    'sobol_ST_threshold', 'primary_metrics', 'primary_selection_table', ...
    'case_profile', 'age_validity', 'baseline_provenance', 'scaling_comparison');
fprintf('[main_run] Best candidate package saved to:\n          %s\n', best_candidate_file);

scientific_candidate_file = fullfile(run_ctx.mat_dir, sprintf('params_scientific_candidate_%s.mat', scenario));
save(scientific_candidate_file, 'scientific_candidate', 'calib_out', ...
    'optMask', 'sobol_ST_threshold', 'primary_metrics', ...
    'primary_selection_table', 'case_profile', 'age_validity', ...
    'baseline_provenance', 'scaling_comparison');
fprintf('[main_run] Scientific candidate package saved to:\n          %s\n', scientific_candidate_file);

accepted_candidate_file = fullfile(run_ctx.mat_dir, sprintf('params_accepted_candidate_%s.mat', scenario));
save(accepted_candidate_file, 'accepted_candidate', 'report', 'calibration_status', ...
    'optMask', 'sobol_ST_threshold', 'primary_metrics', 'primary_selection_table', ...
    'case_profile', 'age_validity', 'baseline_provenance', 'scaling_comparison');
fprintf('[main_run] Accepted candidate package saved to:\n          %s\n', accepted_candidate_file);

write_validation_exports(report, run_ctx.tables_dir, scenario);
plausibility_export_source = calib_out;
if isfield(accepted_candidate, 'parameter_plausibility') && ~isempty(accepted_candidate.parameter_plausibility)
    plausibility_export_source = struct('parameterPlausibility', accepted_candidate.parameter_plausibility);
end
write_parameter_plausibility_exports(plausibility_export_source, run_ctx.tables_dir, scenario);
write_parameter_registry_exports(calib_out, run_ctx.tables_dir, scenario);
best_candidate_tables = write_candidate_exports(best_candidate, run_ctx.tables_dir, scenario);
scientific_candidate_tables = write_candidate_exports(scientific_candidate, run_ctx.tables_dir, scenario);
accepted_candidate_tables = write_candidate_exports(accepted_candidate, run_ctx.tables_dir, scenario);
[baseline_provenance_file, baseline_provenance_summary_file] = ...
    write_baseline_provenance_exports(baseline_provenance, run_ctx.tables_dir, scenario);
scaling_comparison_file = fullfile(run_ctx.tables_dir, ...
    sprintf('scaling_method_comparison_%s.csv', scenario));
writetable(scaling_comparison, scaling_comparison_file);
co_definition_audit_file = write_co_definition_audit_exports( ...
    metrics_base, scientific_candidate, accepted_candidate, clinical, scenario, run_ctx.tables_dir);

run_package = struct();
run_package.scenario = scenario;
run_package.timestamp = timestamp;
run_package.run_folder = run_ctx.root;
run_package.patient_label = run_ctx.patient_label;
run_package.clinical = clinical;
run_package.params_baseline = params0;
run_package.params_calibrated = params_cal;
run_package.sim_baseline = sim_base;
run_package.sim_calibrated = sim_cal;
run_package.metrics_baseline = metrics_base;
run_package.metrics_calibrated = metrics_cal;
run_package.validity_baseline = validity_base;
run_package.validity_calibrated = validity_cal;
run_package.report = report;
run_package.calib_out = calib_out;
run_package.calibration_status = calibration_status;
run_package.best_candidate = best_candidate;
run_package.scientific_candidate = scientific_candidate;
run_package.accepted_candidate = accepted_candidate;
run_package.optMask = optMask;
run_package.sobol_ST_threshold = sobol_ST_threshold;
run_package.primary_metrics = primary_metrics;
run_package.primary_selection_table = primary_selection_table;
run_package.case_profile = case_profile;
run_package.age_validity = age_validity;
run_package.baseline_provenance = baseline_provenance;
run_package.baseline_provenance_file = baseline_provenance_file;
run_package.scaling_comparison = scaling_comparison;
run_package.scaling_comparison_file = scaling_comparison_file;
run_package.co_definition_audit_file = co_definition_audit_file;
run_package.candidate_table_exports = struct( ...
    'best_candidate', best_candidate_tables, ...
    'scientific_candidate', scientific_candidate_tables, ...
    'accepted_candidate', accepted_candidate_tables);
run_package.mask_metrics = mask_metrics;
run_package.gsa_init_cfg = gsa_init_cfg;
run_package.gsa_init_out = gsa_init_out;
run_package.gsa_final_cfg = gsa_final_cfg;
run_package.gsa_final_out = gsa_final_out;
run_package.paths = run_ctx;
run_package_file = fullfile(run_ctx.mat_dir, ...
    sprintf('run_package_%s_%s.mat', scenario, timestamp));
save(run_package_file, 'run_package', '-v7.3');
fprintf('[main_run] Full run package saved to:\n          %s\n', run_package_file);

if strcmp(scenario, 'pre_surgery')
    % Export a dedicated handoff package so the calibrated pre-op state can
    % be reused later as a seed for post-op experiments.
    pre_to_post_seed = struct();
    pre_to_post_seed.source_scenario = scenario;
    pre_to_post_seed.timestamp = timestamp;
    pre_to_post_seed.clinical = clinical;
    pre_to_post_seed.params_baseline = params0;
    pre_to_post_seed.params_calibrated = params_cal;
    pre_to_post_seed.best_candidate = best_candidate;
    pre_to_post_seed.scientific_candidate = scientific_candidate;
    pre_to_post_seed.accepted_candidate = accepted_candidate;
    pre_to_post_seed.metrics_baseline = metrics_base;
    pre_to_post_seed.metrics_calibrated = metrics_cal;
    pre_to_post_seed.validity_baseline = validity_base;
    pre_to_post_seed.validity_calibrated = validity_cal;
    pre_to_post_seed.report = report;
    pre_to_post_seed.calib_out = calib_out;
    pre_to_post_seed.calibration_status = calibration_status;
    pre_to_post_seed.optMask = optMask;
    pre_to_post_seed.sobol_ST_threshold = sobol_ST_threshold;
    pre_to_post_seed.primary_metrics = primary_metrics;
    pre_to_post_seed.primary_selection_table = primary_selection_table;
    pre_to_post_seed.case_profile = case_profile;

    % Explicit matrix/vector payload for downstream scripts.
    accepted_x_all = calib_out.xbest_all(:);
    accepted_x_active = calib_out.xbest(:);
    if isfield(accepted_candidate, 'rollback_applied') && accepted_candidate.rollback_applied
        accepted_x_all = calib_out.x0_all(:);
        accepted_x_active = calib_out.x0(:);
    end
    pre_to_post_seed.parameter_names_all = calib_out.names_all(:);
    pre_to_post_seed.parameter_names_active = calib_out.names(:);
    pre_to_post_seed.parameter_vector_initial_all = calib_out.x0_all(:);
    pre_to_post_seed.parameter_vector_calibrated_all = accepted_x_all;
    pre_to_post_seed.parameter_vector_initial_active = calib_out.x0(:);
    pre_to_post_seed.parameter_vector_calibrated_active = accepted_x_active;
    pre_to_post_seed.parameter_change_pct_active = ...
        (accepted_x_active - calib_out.x0(:)) ./ max(abs(calib_out.x0(:)), 1e-9) * 100;
    pre_to_post_seed.parameter_vector_best_candidate_all = calib_out.xbest_all(:);
    pre_to_post_seed.parameter_vector_best_candidate_active = calib_out.xbest(:);
    pre_to_post_seed.parameter_vector_scientific_candidate_all = calib_out.xbest_all(:);
    pre_to_post_seed.parameter_vector_scientific_candidate_active = calib_out.xbest(:);

    % Convenience copy: same calibrated pre-op parameters, ready to be
    % loaded and modified for a post-op closure experiment.
    post_seed_params = params_cal;
    if isfield(post_seed_params, 'R') && isfield(post_seed_params.R, 'vsd')
        post_seed_params.R.vsd = 1e6;
    end
    pre_to_post_seed.params_post_seed_closed_vsd = post_seed_params;

    seed_fname_timestamped = fullfile(run_ctx.mat_dir, ...
        sprintf('pre_to_post_seed_%s.mat', timestamp));
    seed_fname_latest = fullfile(run_ctx.mat_dir, 'pre_to_post_seed_latest.mat');
    save(seed_fname_timestamped, 'pre_to_post_seed');
    save(seed_fname_latest, 'pre_to_post_seed');
    fprintf('[main_run] Pre-to-post seed package saved to:\n          %s\n          %s\n', ...
        seed_fname_timestamped, seed_fname_latest);
end

write_run_manifest(run_ctx, scenario, timestamp, clinical, patient.scaling_mode, case_profile, calibration_status, DO_GSA, DO_PLOTS, ...
    report, params_package_file, run_package_file, best_candidate_file, scientific_candidate_file, ...
    accepted_candidate_file, best_candidate, scientific_candidate, accepted_candidate, age_validity, ...
    baseline_provenance_file, baseline_provenance_summary_file);

fprintf('\n[main_run] Done. Scenario: %s  |  Total time: %.1f s\n', scenario, toc(run_timer));

end  % main_run

function project_path = build_clean_project_path(root)
root_paths = strsplit(genpath(root), pathsep);
root_paths = root_paths(~cellfun('isempty', root_paths));
is_shadow = contains(root_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, root_paths);
project_path = strjoin(root_paths(~is_shadow & is_existing), pathsep);
end

function canonical_root = resolve_canonical_project_root(root)
% RESOLVE_CANONICAL_PROJECT_ROOT - avoid running from shadow worktree copies.
canonical_root = root;
shadow_marker = [filesep '.claude' filesep 'worktrees' filesep];
marker_idx = strfind(lower(root), lower(shadow_marker));
if isempty(marker_idx)
    return;
end
canonical_root = root(1:marker_idx(1)-1);
end

function remove_shadow_project_paths(project_root)
% REMOVE_SHADOW_PROJECT_PATHS - remove shadow worktree copies from MATLAB path.
path_entries = strsplit(path, pathsep);
path_entries = path_entries(~cellfun('isempty', path_entries));
shadow_prefixes = { ...
    fullfile(project_root, '.claude'), ...
    fullfile(project_root, '.clone')};

for idx = 1:numel(path_entries)
    entry = path_entries{idx};
    for shadow_idx = 1:numel(shadow_prefixes)
        shadow_prefix = shadow_prefixes{shadow_idx};
        if startsWith(lower(entry), lower(shadow_prefix))
            rmpath(entry);
            break;
        end
    end
end
end

function metrics = restrict_metrics_by_case_profile(metrics, case_profile)
% RESTRICT_METRICS_BY_CASE_PROFILE — keep GSA screening aligned with evidence.
if ~isfield(case_profile, 'allowedMetricFields') || isempty(case_profile.allowedMetricFields)
    return;
end
metrics = intersect(metrics, case_profile.allowedMetricFields, 'stable');
if isempty(metrics)
    metrics = case_profile.allowedMetricFields(:)';
end
end

function metrics = profile_holdout_metrics(case_profile)
% PROFILE_HOLDOUT_METRICS - metrics reported but excluded from fit gates.
metrics = {};
if isstruct(case_profile) && isfield(case_profile, 'validationHoldoutMetrics') && ...
        ~isempty(case_profile.validationHoldoutMetrics)
    metrics = case_profile.validationHoldoutMetrics(:)';
end
end

function target_tiers = profile_target_tiers(case_profile)
% PROFILE_TARGET_TIERS - target governance passed into validation reports.
target_tiers = [];
if isstruct(case_profile) && isfield(case_profile, 'targetTiers') && ...
        isstruct(case_profile.targetTiers)
    target_tiers = case_profile.targetTiers;
end
end

function audit = profile_clinical_audit(case_profile)
% PROFILE_CLINICAL_AUDIT - clinical consistency audit for validation output.
audit = [];
if isstruct(case_profile) && isfield(case_profile, 'clinicalConsistencyAudit') && ...
        isstruct(case_profile.clinicalConsistencyAudit)
    audit = case_profile.clinicalConsistencyAudit;
end
end

function tf = should_rollback_calibration(report, calibration_status)
% SHOULD_ROLLBACK_CALIBRATION - protect final artifacts from rejected candidates.
rmse_worse = isfield(report, 'rmse_baseline') && isfield(report, 'rmse_cal') && ...
    isfinite(report.rmse_baseline) && isfinite(report.rmse_cal) && ...
    report.rmse_cal > report.rmse_baseline;
status_reject = isfield(calibration_status, 'label') && ...
    any(strcmp(calibration_status.label, {'REJECT','OUTPUT_FIT_ONLY'}));
tf = rmse_worse || status_reject;
end

function run_ctx = init_run_output(root, scenario, clinical)
% INIT_RUN_OUTPUT â€” create a single dated run folder and route sub-artifacts into it.
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
patient_label = resolve_patient_label(clinical);
run_name = sprintf('%s_%s_%s', timestamp, patient_label, scenario);

run_ctx = struct();
run_ctx.timestamp = timestamp;
run_ctx.patient_label = patient_label;
run_ctx.root = fullfile(root, 'results', 'runs', run_name);
run_ctx.tables_dir = fullfile(run_ctx.root, 'tables');
run_ctx.figures_dir = fullfile(run_ctx.root, 'figures');
run_ctx.gsa_dir = fullfile(run_ctx.root, 'gsa');
run_ctx.mat_dir = fullfile(run_ctx.root, 'mat');
run_ctx.logs_dir = fullfile(run_ctx.root, 'logs');
run_ctx.console_log = fullfile(run_ctx.logs_dir, sprintf('console_%s.log', run_name));
run_ctx.gsa_checkpoint_file = fullfile(run_ctx.gsa_dir, ...
    sprintf('gsa_pce_checkpoint_%s_%s.mat', patient_label, scenario));
run_ctx.previous_env.UNIFIED_VSD_RUN_ROOT = getenv('UNIFIED_VSD_RUN_ROOT');
run_ctx.previous_env.UNIFIED_VSD_GSA_DIR = getenv('UNIFIED_VSD_GSA_DIR');
run_ctx.previous_env.UNIFIED_VSD_GSA_CHECKPOINT_FILE = getenv('UNIFIED_VSD_GSA_CHECKPOINT_FILE');

mkdir_if_missing(run_ctx.root);
mkdir_if_missing(run_ctx.tables_dir);
mkdir_if_missing(run_ctx.figures_dir);
mkdir_if_missing(run_ctx.gsa_dir);
mkdir_if_missing(run_ctx.mat_dir);
mkdir_if_missing(run_ctx.logs_dir);

setenv('UNIFIED_VSD_RUN_ROOT', run_ctx.root);
setenv('UNIFIED_VSD_GSA_DIR', run_ctx.gsa_dir);
setenv('UNIFIED_VSD_GSA_CHECKPOINT_FILE', run_ctx.gsa_checkpoint_file);

if strcmpi(get(0, 'Diary'), 'on')
    diary off;
end
diary(run_ctx.console_log);
diary on;
fprintf('[main_run] Console diary started: %s\n', run_ctx.console_log);
end

function cleanup_run_output(run_ctx)
% CLEANUP_RUN_OUTPUT â€” stop diary and restore environment routing.
if strcmpi(get(0, 'Diary'), 'on')
    diary off;
end
restore_env_field(run_ctx.previous_env, 'UNIFIED_VSD_RUN_ROOT');
restore_env_field(run_ctx.previous_env, 'UNIFIED_VSD_GSA_DIR');
restore_env_field(run_ctx.previous_env, 'UNIFIED_VSD_GSA_CHECKPOINT_FILE');
end

function restore_env_field(previous_env, field_name)
% RESTORE_ENV_FIELD â€” restore previous environment value or clear it.
if ~isfield(previous_env, field_name)
    return;
end
if isempty(previous_env.(field_name))
    setenv(field_name, '');
else
    setenv(field_name, previous_env.(field_name));
end
end

function patient_label = resolve_patient_label(clinical)
% RESOLVE_PATIENT_LABEL â€” derive a run-folder-safe patient label.
patient_label = '';
if isfield(clinical, 'common')
    if isfield(clinical.common, 'patient_name') && ~isempty(clinical.common.patient_name)
        patient_label = clinical.common.patient_name;
    elseif isfield(clinical.common, 'patient_id') && ~isempty(clinical.common.patient_id)
        patient_label = clinical.common.patient_id;
    end
end

if isempty(patient_label)
    patient_label = sprintf('patient_%0.1fkg_%0.0fcm', ...
        clinical.common.weight_kg, clinical.common.height_cm);
end

patient_label = lower(char(patient_label));
patient_label = strrep(patient_label, ' ', '_');
patient_label = strrep(patient_label, '.', 'p');
patient_label = regexprep(patient_label, '[^a-z0-9_]+', '');
if isempty(patient_label)
    patient_label = 'patient_unknown';
end
end

function mkdir_if_missing(folder_path)
% MKDIR_IF_MISSING â€” create output directory when absent.
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end
end

function write_validation_exports(report, tables_dir, scenario)
% WRITE_VALIDATION_EXPORTS â€” persist key validation tables as CSV artifacts.
writetable(report.table_baseline, fullfile(tables_dir, ...
    sprintf('validation_baseline_%s.csv', scenario)));
if ~isempty(report.table_cal)
    writetable(report.table_cal, fullfile(tables_dir, ...
        sprintf('validation_calibrated_%s.csv', scenario)));
end
writetable(report.table_delta, fullfile(tables_dir, ...
    sprintf('validation_delta_%s.csv', scenario)));
writetable(report.sorted_errors, fullfile(tables_dir, ...
    sprintf('validation_sorted_errors_%s.csv', scenario)));
writetable(report.primary_gate, fullfile(tables_dir, ...
    sprintf('validation_primary_gate_%s.csv', scenario)));
if isfield(report, 'rmse_summary') && ~isempty(report.rmse_summary)
    writetable(report.rmse_summary, fullfile(tables_dir, ...
        sprintf('validation_rmse_summary_%s.csv', scenario)));
end
if isfield(report, 'target_tier_table') && ~isempty(report.target_tier_table)
    writetable(report.target_tier_table, fullfile(tables_dir, ...
        sprintf('validation_target_tiers_%s.csv', scenario)));
end
if isfield(report, 'clinical_data_rank') && ~isempty(report.clinical_data_rank)
    writetable(report.clinical_data_rank, fullfile(tables_dir, ...
        sprintf('clinical_data_availability_rank_%s.csv', scenario)));
    rank_summary = groupsummary(report.clinical_data_rank, 'EvidenceClass');
    writetable(rank_summary, fullfile(tables_dir, ...
        sprintf('clinical_data_availability_rank_summary_%s.csv', scenario)));
end
if isfield(report, 'clinical_consistency_audit') && ...
        isstruct(report.clinical_consistency_audit)
    write_clinical_consistency_exports(report.clinical_consistency_audit, tables_dir, scenario);
end
end

function write_clinical_consistency_exports(audit, tables_dir, scenario)
if isfield(audit, 'stroke_volume_table') && ~isempty(audit.stroke_volume_table)
    writetable(audit.stroke_volume_table, fullfile(tables_dir, ...
        sprintf('clinical_consistency_stroke_volumes_%s.csv', scenario)));
end
if isfield(audit, 'pairwise_table') && ~isempty(audit.pairwise_table)
    writetable(audit.pairwise_table, fullfile(tables_dir, ...
        sprintf('clinical_consistency_pairwise_%s.csv', scenario)));
end
end

function write_parameter_plausibility_exports(calib_out, tables_dir, scenario)
if ~isfield(calib_out, 'parameterPlausibility') || ...
        ~isfield(calib_out.parameterPlausibility, 'table') || ...
        isempty(calib_out.parameterPlausibility.table)
    return;
end

writetable(calib_out.parameterPlausibility.table, fullfile(tables_dir, ...
    sprintf('parameter_plausibility_%s.csv', scenario)));
end

function write_parameter_registry_exports(calib_out, tables_dir, scenario)
if isfield(calib_out, 'parameterRegistry') && ~isempty(calib_out.parameterRegistry)
    writetable(calib_out.parameterRegistry, fullfile(tables_dir, ...
        sprintf('parameter_registry_%s.csv', scenario)));
end
if isfield(calib_out, 'parameterRegistryActive') && ~isempty(calib_out.parameterRegistryActive)
    writetable(calib_out.parameterRegistryActive, fullfile(tables_dir, ...
        sprintf('parameter_registry_active_%s.csv', scenario)));
end
end

function candidate_paths = write_candidate_exports(candidate, tables_dir, scenario)
% WRITE_CANDIDATE_EXPORTS â€” persist explicit best/accepted candidate tables.
candidate_paths = struct('validation', '', 'plausibility', '');
if isempty(candidate) || ~isstruct(candidate) || ~isfield(candidate, 'label')
    return;
end

label = char(candidate.label);
if isfield(candidate, 'report') && isstruct(candidate.report)
    validation_path = fullfile(tables_dir, sprintf('validation_%s_%s.csv', label, scenario));
    if isfield(candidate.report, 'table_cal') && ~isempty(candidate.report.table_cal)
        writetable(candidate.report.table_cal, validation_path);
    else
        writetable(candidate.report.table_baseline, validation_path);
    end
    candidate_paths.validation = validation_path;
end

if isfield(candidate, 'parameter_plausibility') && isstruct(candidate.parameter_plausibility) && ...
        isfield(candidate.parameter_plausibility, 'table') && ~isempty(candidate.parameter_plausibility.table)
    plausibility_path = fullfile(tables_dir, sprintf('parameter_plausibility_%s_%s.csv', label, scenario));
    writetable(candidate.parameter_plausibility.table, plausibility_path);
    candidate_paths.plausibility = plausibility_path;
end
end

function csv_path = write_co_definition_audit_exports(metrics_base, scientific_candidate, accepted_candidate, clinical, scenario, tables_dir)
% WRITE_CO_DEFINITION_AUDIT_EXPORTS - persist Qs/Qp/LVCO/SVR definitions.
csv_path = fullfile(tables_dir, sprintf('co_definition_audit_%s.csv', scenario));
if ~isstruct(clinical) || ~isfield(clinical, scenario)
    return;
end

src = clinical.(scenario);
target_co = clinical_field_or_nan(src, 'CO_Lmin');          % [L/min]
target_qpqs = clinical_field_or_nan(src, 'QpQs');           % [-]
target_svr = clinical_field_or_nan(src, 'SVR_WU');          % [WU]
target_qp = target_co * target_qpqs;                        % [L/min]
co_target_sigma = clinical_field_or_nan(src, 'CO_uncertainty_Lmin'); % [L/min]
co_comparator = '';
if isfield(src, 'CO_comparator') && ~isempty(src.CO_comparator)
    co_comparator = char(src.CO_comparator);
end

candidate_names = {'baseline','scientific_candidate','accepted_candidate'};
metrics_list = {metrics_base, scientific_candidate.metrics, accepted_candidate.metrics};
rows = cell(numel(candidate_names), 23);
for idx = 1:numel(candidate_names)
    metrics = metrics_list{idx};
    q_sys = field_or_nan(metrics, 'CO_Lmin');               % [L/min]
    q_pul = field_or_nan(metrics, 'Qp_Lmin');               % [L/min]
    q_ao = field_or_nan(metrics, 'Qao_Lmin');               % [L/min]
    lv_co = field_or_nan(metrics, 'LVCO_Lmin');             % [L/min]
    q_vsd = field_or_nan(metrics, 'Qvsd_Lmin');             % [L/min]
    sap_mean = field_or_nan(metrics, 'SAP_mean');           % [mmHg]
    rap_mean = field_or_nan(metrics, 'RAP_mean');           % [mmHg]
    map_minus_rap = sap_mean - rap_mean;                    % [mmHg]
    svr_qs = field_or_nan(metrics, 'SVR');                  % [WU]
    svr_qao = map_minus_rap / max(q_ao, 1e-9);              % [WU]
    svr_lvco = map_minus_rap / max(lv_co, 1e-9);            % [WU]
    q_from_model_svr = map_minus_rap / max(svr_qs, 1e-9);   % [L/min]
    q_if_target_svr = map_minus_rap / max(target_svr, 1e-9); % [L/min]

    rows(idx, :) = { ...
        candidate_names{idx}, q_sys, q_pul, q_ao, lv_co, q_vsd, ...
        field_or_nan(metrics, 'QpQs'), field_or_nan(metrics, 'Qs_ao_gap_Lmin'), ...
        sap_mean, rap_mean, map_minus_rap, svr_qs, svr_qao, svr_lvco, ...
        q_from_model_svr, q_if_target_svr, target_qp, lv_co - q_pul, ...
        lv_co - target_qp, 100 * (q_sys - target_co) / max(abs(target_co), 1e-9), ...
        target_co, co_target_sigma, co_comparator};
end

audit_tbl = cell2table(rows, 'VariableNames', { ...
    'Candidate','Qs_Lmin','Qp_Lmin','Qao_Lmin','LVCO_Lmin','Qvsd_Lmin', ...
    'QpQs','QaoMinusQs_Lmin','SAP_mean','RAP_mean','MAPMinusRAP', ...
    'SVR_from_Qs','SVR_from_Qao','SVR_from_LVCO','Q_from_model_SVR', ...
    'Q_if_target_SVR','Clinical_Qp_Target_Lmin','LVCO_minus_Qp_Lmin', ...
    'LVCO_minus_clinical_Qp_Lmin','CO_ErrorPct','Clinical_Qs_Target_Lmin', ...
    'CO_TargetUncertainty_Lmin','CO_Comparator'});
writetable(audit_tbl, csv_path);
end

function [csv_path, txt_path] = write_baseline_provenance_exports(baseline_provenance, tables_dir, scenario)
% WRITE_BASELINE_PROVENANCE_EXPORTS â€” persist baseline source lineage tables.
csv_path = fullfile(tables_dir, sprintf('baseline_provenance_%s.csv', scenario));
txt_path = fullfile(tables_dir, sprintf('baseline_provenance_summary_%s.txt', scenario));
writetable(baseline_provenance, csv_path);

summary_tbl = groupsummary(baseline_provenance, 'FoundationSource');
fid = fopen(txt_path, 'w');
if fid < 0
    warning('main_run:baselineProvenanceOpenFailed', ...
        'Unable to write baseline provenance summary.');
    txt_path = '';
    return;
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Unified VSD Baseline Provenance Summary\n');
fprintf(fid, '======================================\n');
fprintf(fid, 'Scenario: %s\n', scenario);
fprintf(fid, 'Rows: %d\n\n', height(baseline_provenance));
for idx = 1:height(summary_tbl)
    fprintf(fid, '%s : %d\n', ...
        summary_tbl.FoundationSource{idx}, summary_tbl.GroupCount(idx));
end
end

function write_run_manifest(run_ctx, scenario, timestamp, clinical, scaling_mode, case_profile, calibration_status, do_gsa, do_plots, report, params_package_file, run_package_file, best_candidate_file, scientific_candidate_file, accepted_candidate_file, best_candidate, scientific_candidate, accepted_candidate, age_validity, baseline_provenance_file, baseline_provenance_summary_file)
% WRITE_RUN_MANIFEST â€” create a plain-text run summary inside the archive.
manifest_file = fullfile(run_ctx.root, 'run_manifest.txt');
fid = fopen(manifest_file, 'w');
if fid < 0
    warning('main_run:manifestOpenFailed', 'Unable to write run manifest: %s', manifest_file);
    return;
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

clinical_bsa = NaN;
if isfield(clinical, 'common') && isfield(clinical.common, 'BSA') && ...
        ~isempty(clinical.common.BSA) && isfinite(clinical.common.BSA)
    clinical_bsa = clinical.common.BSA;
end

fprintf(fid, 'Unified VSD Run Manifest\n');
fprintf(fid, '========================\n');
fprintf(fid, 'Timestamp: %s\n', timestamp);
fprintf(fid, 'Scenario: %s\n', scenario);
fprintf(fid, 'PatientLabel: %s\n', run_ctx.patient_label);
fprintf(fid, 'ScalingMode: %s\n', scaling_mode);
fprintf(fid, 'CalibrationCaseMode: %s\n', case_profile.mode);
fprintf(fid, 'CalibrationCaseDescription: %s\n', case_profile.description);
if isfield(case_profile, 'targetGovernance')
    fprintf(fid, 'TargetGovernance: %s\n', case_profile.targetGovernance);
end
if isfield(case_profile, 'validationHoldoutMetrics') && ~isempty(case_profile.validationHoldoutMetrics)
    fprintf(fid, 'ValidationHoldoutMetrics: %s\n', ...
        strjoin(case_profile.validationHoldoutMetrics(:)', ','));
end
if isfield(case_profile, 'ageValidityRegime')
    fprintf(fid, 'CalibrationAgeValidityRegime: %s\n', case_profile.ageValidityRegime);
end
if isfield(case_profile, 'preschoolPriorAdjustmentApplied')
    fprintf(fid, 'PreschoolPriorAdjustmentApplied: %d\n', ...
        case_profile.preschoolPriorAdjustmentApplied);
end
if isfield(case_profile, 'paramPlausibilityLambdaScale')
    fprintf(fid, 'ParamPlausibilityLambdaScale: %.4f\n', ...
        case_profile.paramPlausibilityLambdaScale);
end
fprintf(fid, 'Weight_kg: %.4f\n', clinical.common.weight_kg);
fprintf(fid, 'Height_cm: %.4f\n', clinical.common.height_cm);
fprintf(fid, 'Age_years: %.4f\n', clinical.common.age_years);
fprintf(fid, 'BSA_m2: %.4f\n', clinical_bsa);
fprintf(fid, 'AgeValidityRegime: %s\n', age_validity.regime);
fprintf(fid, 'AgeValidityEvidence: %s\n', age_validity.evidence_level);
fprintf(fid, 'AgeValiditySummary: %s\n', age_validity.summary);
fprintf(fid, 'AgeValidityAction: %s\n', age_validity.action);
fprintf(fid, 'MaturationMode: %s\n', age_validity.maturation_mode);
fprintf(fid, 'PlotsEnabled: %d\n', do_plots);
fprintf(fid, 'GSAEnabled: %d\n', do_gsa);
fprintf(fid, 'RMSE_Baseline: %.6f\n', report.rmse_baseline);
fprintf(fid, 'RMSE_Calibrated: %.6f\n', report.rmse_cal);
if isfield(report, 'rmse_full_baseline')
    fprintf(fid, 'RMSE_Full_Baseline: %.6f\n', report.rmse_full_baseline);
    fprintf(fid, 'RMSE_Full_Calibrated: %.6f\n', report.rmse_full_cal);
    fprintf(fid, 'RMSE_Primary_Baseline: %.6f\n', report.rmse_primary_baseline);
    fprintf(fid, 'RMSE_Primary_Calibrated: %.6f\n', report.rmse_primary_cal);
    fprintf(fid, 'RMSE_Hard_Calibrated: %.6f\n', report.rmse_hard_cal);
    fprintf(fid, 'RMSE_Soft_Calibrated: %.6f\n', report.rmse_soft_cal);
end
if isfield(report, 'clinical_consistency_audit') && ...
        isfield(report.clinical_consistency_audit, 'summary')
    fprintf(fid, 'ClinicalConsistencyAudit: %s\n', report.clinical_consistency_audit.summary);
end
if isfield(report, 'target_tiers') && isfield(report.target_tiers, 'consistency_only')
    fprintf(fid, 'ConsistencyOnlyTargets: %s\n', ...
        manifest_join_or_none(report.target_tiers.consistency_only));
end
fprintf(fid, 'CalibrationStatus: %s\n', calibration_status.label);
fprintf(fid, 'CalibrationStatusSummary: %s\n', calibration_status.summary);
if isfield(best_candidate, 'report') && isfield(best_candidate.report, 'rmse_cal')
    fprintf(fid, 'BestCandidateRMSE: %.6f\n', best_candidate.report.rmse_cal);
end
if isfield(best_candidate, 'calibration_status') && isfield(best_candidate.calibration_status, 'label')
    fprintf(fid, 'BestCandidateStatus: %s\n', best_candidate.calibration_status.label);
end
if isfield(scientific_candidate, 'report') && isfield(scientific_candidate.report, 'rmse_cal')
    fprintf(fid, 'ScientificCandidateRMSE: %.6f\n', scientific_candidate.report.rmse_cal);
end
if isfield(scientific_candidate, 'calibration_status') && isfield(scientific_candidate.calibration_status, 'label')
    fprintf(fid, 'ScientificCandidateStatus: %s\n', scientific_candidate.calibration_status.label);
end
if isfield(accepted_candidate, 'rollback_applied')
    fprintf(fid, 'RollbackApplied: %d\n', accepted_candidate.rollback_applied);
end
if isfield(accepted_candidate, 'rollback_reason') && ~isempty(accepted_candidate.rollback_reason)
    fprintf(fid, 'RollbackReason: %s\n', accepted_candidate.rollback_reason);
end
fprintf(fid, 'RunFolder: %s\n', run_ctx.root);
fprintf(fid, 'ConsoleLog: %s\n', run_ctx.console_log);
fprintf(fid, 'FiguresDir: %s\n', run_ctx.figures_dir);
fprintf(fid, 'TablesDir: %s\n', run_ctx.tables_dir);
fprintf(fid, 'GsaDir: %s\n', run_ctx.gsa_dir);
fprintf(fid, 'MatDir: %s\n', run_ctx.mat_dir);
fprintf(fid, 'ParamsPackage: %s\n', params_package_file);
fprintf(fid, 'BestCandidatePackage: %s\n', best_candidate_file);
fprintf(fid, 'ScientificCandidatePackage: %s\n', scientific_candidate_file);
fprintf(fid, 'AcceptedCandidatePackage: %s\n', accepted_candidate_file);
fprintf(fid, 'RunPackage: %s\n', run_package_file);
fprintf(fid, 'BaselineProvenanceFile: %s\n', baseline_provenance_file);
fprintf(fid, 'BaselineProvenanceSummary: %s\n', baseline_provenance_summary_file);
end

function text = manifest_join_or_none(values)
if isempty(values)
    text = 'none';
else
    text = strjoin(values(:)', ',');
end
end

function candidate = build_candidate_snapshot(label, params, sim, metrics, validity, report, calibration_status, parameter_plausibility, rollback_applied, rollback_reason)
% BUILD_CANDIDATE_SNAPSHOT â€” capture one candidate state without overwriting another.
candidate = struct();
candidate.label = label;
candidate.params = params;
candidate.sim = sim;
candidate.metrics = metrics;
candidate.validity = validity;
candidate.report = report;
candidate.calibration_status = calibration_status;
candidate.parameter_plausibility = parameter_plausibility;
candidate.rollback_applied = rollback_applied;
candidate.rollback_reason = rollback_reason;
end

function summary = summarize_candidate_snapshot(candidate)
% SUMMARIZE_CANDIDATE_SNAPSHOT â€” compact metadata for diagnostics MAT files.
summary = struct();
if isempty(candidate) || ~isstruct(candidate)
    return;
end
summary.label = candidate.label;
summary.rollback_applied = candidate.rollback_applied;
summary.rollback_reason = candidate.rollback_reason;
if isfield(candidate, 'report') && isstruct(candidate.report)
    if isfield(candidate.report, 'rmse_baseline')
        summary.rmse_baseline = candidate.report.rmse_baseline;
    end
    if isfield(candidate.report, 'rmse_cal')
        summary.rmse_cal = candidate.report.rmse_cal;
    end
    if isfield(candidate.report, 'rmse_full_cal')
        summary.rmse_full_cal = candidate.report.rmse_full_cal;
    end
    if isfield(candidate.report, 'rmse_primary_cal')
        summary.rmse_primary_cal = candidate.report.rmse_primary_cal;
    end
end
if isfield(candidate, 'calibration_status') && isstruct(candidate.calibration_status)
    if isfield(candidate.calibration_status, 'label')
        summary.status = candidate.calibration_status.label;
    end
    if isfield(candidate.calibration_status, 'summary')
        summary.status_summary = candidate.calibration_status.summary;
    end
end
end

function print_systemic_consistency_summary(metrics, clinical, scenario, tag)
% PRINT_SYSTEMIC_CONSISTENCY_SUMMARY — concise systemic flow/resistance diagnostics.
if ~isstruct(clinical) || ~isfield(clinical, scenario)
    return;
end

src = clinical.(scenario);
fprintf('\n[main_run] %s systemic consistency summary\n', tag);
fprintf('           Qs(model)=%.4f L/min', field_or_nan(metrics, 'Qs_Lmin'));
if isfield(src, 'CO_Lmin') && ~isnan(src.CO_Lmin)
    fprintf(' | Qs(target)=%.4f L/min', src.CO_Lmin);
end
fprintf('\n');

fprintf('           Qao(model)=%.4f L/min | Qp(model)=%.4f L/min | Qpv(model)=%.4f L/min\n', ...
    field_or_nan(metrics, 'Qao_Lmin'), field_or_nan(metrics, 'Qp_Lmin'), field_or_nan(metrics, 'Qpv_Lmin'));
fprintf('           SVR(model)=%.4f WU | SVR(from Qao)=%.4f WU', ...
    field_or_nan(metrics, 'SVR'), field_or_nan(metrics, 'SVR_from_Qao'));
if isfield(src, 'SVR_WU') && ~isnan(src.SVR_WU)
    fprintf(' | SVR(target)=%.4f WU', src.SVR_WU);
end
fprintf('\n');

fprintf('           Qao-Qs gap=%.4f L/min | Qp-Qpv gap=%.4f L/min\n', ...
    field_or_nan(metrics, 'Qs_ao_gap_Lmin'), field_or_nan(metrics, 'Qp_pv_gap_Lmin'));
fprintf('           SAP min/max/pulse = %.3f / %.3f / %.3f mmHg\n', ...
    field_or_nan(metrics, 'SAP_min'), field_or_nan(metrics, 'SAP_max'), field_or_nan(metrics, 'SAP_pulse'));
fprintf('           PAP min/max/pulse = %.3f / %.3f / %.3f mmHg\n', ...
    field_or_nan(metrics, 'PAP_min'), field_or_nan(metrics, 'PAP_max'), field_or_nan(metrics, 'PAP_pulse'));
end

function value = field_or_nan(s, field_name)
% FIELD_OR_NAN — helper for optional diagnostic fields.
if isstruct(s) && isfield(s, field_name) && isfinite(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

function value = clinical_field_or_nan(s, field_name)
% CLINICAL_FIELD_OR_NAN - scalar clinical field lookup with NaN fallback.
if isstruct(s) && isfield(s, field_name) && isnumeric(s.(field_name)) && ...
        isscalar(s.(field_name)) && isfinite(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

function ST_mapped = map_gsa_sensitivities_to_calibration_names(ST_vec, gsa_names, calib_names, case_profile)
% MAP_GSA_SENSITIVITIES_TO_CALIBRATION_NAMES - aggregate ST onto grouped knobs.
ST_mapped = zeros(numel(calib_names), 1);
for calib_idx = 1:numel(calib_names)
    member_names = expand_calibration_parameter_names(calib_names{calib_idx}, case_profile);
    [is_found, gsa_idx] = ismember(member_names, gsa_names);
    if ~all(is_found)
        missing_names = member_names(~is_found);
        error('main_run:missingGsaMapping', ...
            'Calibration parameter members missing in GSA names: %s', strjoin(missing_names, ', '));
    end
    ST_mapped(calib_idx) = max(ST_vec(gsa_idx));
end
end
