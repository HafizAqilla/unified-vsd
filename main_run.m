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
DO_PLOTS            = true;
DO_OVERLAY          = true;
DO_GSA              = false;   % Run pre- and post-calibration GSA for paper
DO_FAST_CALIBRATION = false;  % Full run: Stage 1 + Stage 2 restarts

% DO_CALIBRATION: used by the post-surgery pipeline only.
%   true  — run fmincon calibration to fit post-op clinical targets.
%   false — run as a pure model prediction (baseline only, no optimisation).
%           Useful to see what the pre-surgery calibrated params predict
%           for the post-surgery state before any tuning.
DO_CALIBRATION = false;

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
%  STEP 1b — Warm-start: inject calibrated pre-surgery parameters
%
%  PURPOSE:
%   When running a post-surgery simulation, the patient's intrinsic
%   cardiovascular physiology (ventricular elastances, vascular R/C,
%   unstressed volumes) has already been personalised by the pre-surgery
%   calibration.  Injecting those calibrated values here ensures that the
%   post-surgery baseline simulation and subsequent calibration start from
%   a physiologically informed, patient-specific operating point — rather
%   than from the generic allometric-scaled values.
%
%  TRIGGER:
%   Activated by clinical.pre_surgery.CalibParams being a non-empty struct
%   (set in run_post_surgery.m by loading params_calibrated_pre_surgery_*.mat).
%
%  FIELDS COPIED (cardiac/vascular only — R.vsd intentionally excluded):
%   Ventricular elastances  : E.LV.EA, E.LV.EB, E.RV.EA, E.RV.EB
%   Atrial elastances       : E.LA.EA, E.LA.EB, E.RA.EA, E.RA.EB
%   Unstressed volumes      : V0.LV, V0.RV, V0.LA, V0.RA
%   Systemic resistances    : R.SAR, R.SC, R.SVEN
%   Pulmonary resistances   : R.PAR, R.PCOX, R.PCNO, R.PVEN
%   Compliances             : C.SAR, C.SVEN, C.PAR, C.PVEN
%
%  R.vsd is EXCLUDED: params_from_clinical (Step 2) sets it to 1e6 for
%  post_surgery — copying the pre-surgery finite value would re-open the shunt.
%% =====================================================================
has_calib_params = strcmp(scenario, 'post_surgery') && ...
                   isfield(clinical, 'pre_surgery')  && ...
                   isfield(clinical.pre_surgery, 'CalibParams') && ...
                   isstruct(clinical.pre_surgery.CalibParams)   && ...
                   ~isempty(fieldnames(clinical.pre_surgery.CalibParams));

if has_calib_params
    fprintf('\n=== [Step 1b] Warm-start: injecting pre-surgery calibrated params ===\n');
    pc = clinical.pre_surgery.CalibParams;   % shorthand

    % ---- Ventricular elastances -----------------------------------------
    warm_fields_E_LV = {'EA', 'EB'};
    for wf = warm_fields_E_LV
        if isfield(pc.E.LV, wf{1}); params0.E.LV.(wf{1}) = pc.E.LV.(wf{1}); end
    end
    warm_fields_E_RV = {'EA', 'EB'};
    for wf = warm_fields_E_RV
        if isfield(pc.E.RV, wf{1}); params0.E.RV.(wf{1}) = pc.E.RV.(wf{1}); end
    end

    % ---- Atrial elastances ----------------------------------------------
    warm_fields_E_LA = {'EA', 'EB'};
    for wf = warm_fields_E_LA
        if isfield(pc.E.LA, wf{1}); params0.E.LA.(wf{1}) = pc.E.LA.(wf{1}); end
    end
    warm_fields_E_RA = {'EA', 'EB'};
    for wf = warm_fields_E_RA
        if isfield(pc.E.RA, wf{1}); params0.E.RA.(wf{1}) = pc.E.RA.(wf{1}); end
    end

    % ---- Unstressed volumes ---------------------------------------------
    warm_V0 = {'LV', 'RV', 'LA', 'RA'};
    for wf = warm_V0
        if isfield(pc.V0, wf{1}); params0.V0.(wf{1}) = pc.V0.(wf{1}); end
    end

    % ---- Systemic resistances -------------------------------------------
    warm_R_sys = {'SAR', 'SC', 'SVEN'};
    for wf = warm_R_sys
        if isfield(pc.R, wf{1}); params0.R.(wf{1}) = pc.R.(wf{1}); end
    end

    % ---- Pulmonary resistances ------------------------------------------
    % R.vsd intentionally EXCLUDED — set by params_from_clinical in Step 2.
    warm_R_pul = {'PAR', 'PCOX', 'PCNO', 'PVEN'};
    for wf = warm_R_pul
        if isfield(pc.R, wf{1}); params0.R.(wf{1}) = pc.R.(wf{1}); end
    end

    % ---- Compliances ----------------------------------------------------
    warm_C = {'SAR', 'SVEN', 'PAR', 'PVEN'};
    for wf = warm_C
        if isfield(pc.C, wf{1}); params0.C.(wf{1}) = pc.C.(wf{1}); end
    end

    fprintf('[main_run] Warm-start applied. Key injected values:\n');
    fprintf('  E.LV: EA=%.4f  EB=%.4f  |  E.RV: EA=%.4f  EB=%.4f\n', ...
            params0.E.LV.EA, params0.E.LV.EB, params0.E.RV.EA, params0.E.RV.EB);
    fprintf('  V0.LV=%.2f mL  V0.RV=%.2f mL\n', params0.V0.LV, params0.V0.RV);
    fprintf('  R.SAR=%.4f  R.PAR=%.4f  C.SAR=%.4f  C.PAR=%.4f\n', ...
            params0.R.SAR, params0.R.PAR, params0.C.SAR, params0.C.PAR);
else
    if strcmp(scenario, 'post_surgery')
        fprintf('\n[main_run] No CalibParams found — post-surgery model will start from allometric scaling.\n');
    end
end

%% =====================================================================
%  STEP 2 — Map clinical measurements (HR, SVR, PVR, R_VSD)
%% =====================================================================
fprintf('\n=== [Step 2/10] Mapping clinical measurements (%.1fs elapsed) ===\n', toc(run_timer));
params0 = params_from_clinical(params0, clinical, scenario);

%% =====================================================================
%  POST-SURGERY PIPELINE (simple — returns early before GSA steps)
%
%  For post_surgery the workflow follows the lighter-weight pattern:
%    Step 3  — Baseline simulation     (always runs)
%    Step 4  — Calibration (optional)  (controlled by DO_CALIBRATION)
%    Step 5  — Validation report
%    Step 6  — Plots
%    Step 7  — GSA (optional)          (controlled by DO_GSA)
%    Step 8  — Timestamped save
%
%  Pre-surgery continues below (complex GSA-masked calibration pipeline).
%% =====================================================================
if strcmp(scenario, 'post_surgery')

    %% -----------------------------------------------------------------
    %  POST-SURGERY STEP 3 — Baseline simulation
    %% -----------------------------------------------------------------
    fprintf('\n=== [Post-Surgery Step 3] Baseline simulation (%.1fs elapsed) ===\n', toc(run_timer));
    sim_base     = integrate_system(params0);
    metrics_base = compute_clinical_indices(sim_base, params0);
    fprintf('[main_run] Baseline complete.\n');

    %% -----------------------------------------------------------------
    %  POST-SURGERY STEP 4 — Calibration (optional)
    %
    %  Set DO_CALIBRATION = false to use this as a pure model prediction:
    %  the pre-surgery calibrated parameters are injected in Step 1b and
    %  the baseline simulation above IS the model prediction.
    %  Set DO_CALIBRATION = true to further tune to post-op targets.
    %% -----------------------------------------------------------------
    params_cal  = params0;     % default: calibrated = baseline (no optimisation)
    metrics_cal = [];
    calib_out   = [];
    sim_cal     = sim_base;    % default: same simulation

    if DO_CALIBRATION
        fprintf('\n=== [Post-Surgery Step 4] Calibration — fmincon (%.1fs elapsed) ===\n', toc(run_timer));
        [params_cal, calib_out] = run_calibration(params0, clinical, scenario, [], DO_FAST_CALIBRATION, []);

        sim_cal     = integrate_system(params_cal);
        metrics_cal = compute_clinical_indices(sim_cal, params_cal);
        fprintf('[main_run] Calibration complete. Best J = %.6f\n', calib_out.fbest);

        fprintf('\n--- Calibrated parameter changes ---\n');
        param_tbl = table(calib_out.names(:), calib_out.x0(:), calib_out.xbest(:), ...
            (calib_out.xbest(:) - calib_out.x0(:)) ./ abs(calib_out.x0(:)) * 100, ...
            'VariableNames', {'Parameter', 'Initial', 'Calibrated', 'Change_pct'});
        disp(param_tbl);
    else
        fprintf('\n[main_run] DO_CALIBRATION=false — running as model prediction (no optimisation).\n');
        fprintf('[main_run] Baseline simulation IS the post-surgery prediction.\n');
    end

    %% -----------------------------------------------------------------
    %  POST-SURGERY STEP 5 — Validation report
    %% -----------------------------------------------------------------
    fprintf('\n=== [Post-Surgery Step 5] Validation report (%.1fs elapsed) ===\n', toc(run_timer));
    results_dir_ps = fullfile(root, 'results', 'tables');
    if ~exist(results_dir_ps, 'dir'), mkdir(results_dir_ps); end
    report = validation_report(clinical, metrics_base, metrics_cal, scenario, ...
                               'ResultsDir', results_dir_ps);

    %% -----------------------------------------------------------------
    %  POST-SURGERY STEP 6 — Plots
    %% -----------------------------------------------------------------
    fprintf('\n=== [Post-Surgery Step 6] Plots (%.1fs elapsed) ===\n', toc(run_timer));
    if DO_PLOTS
        plotting_tools(sim_base, params0, 'Baseline (post-surgery prediction)', scenario);
        if DO_CALIBRATION && ~isempty(metrics_cal)
            plotting_tools(sim_cal, params_cal, 'Calibrated', scenario);
            if DO_OVERLAY
                plot_overlay_comparison(sim_base, params0, 'Prediction', ...
                    sim_cal, params_cal, 'Calibrated', scenario);
            end
        end
    end

    %% -----------------------------------------------------------------
    %  POST-SURGERY STEP 7 — GSA (optional)
    %% -----------------------------------------------------------------
    if DO_GSA
        fprintf('\n=== [Post-Surgery Step 7] PCE GSA (%.1fs elapsed) ===\n', toc(run_timer));
        params_for_gsa = params_cal;   % use calibrated (or baseline if DO_CALIBRATION=false)
        gsa_ps_cfg = gsa_pce_setup(params_for_gsa, scenario);
        gsa_ps_out = gsa_run_pce(gsa_ps_cfg, params_for_gsa);

        gsa_dir_ps = fullfile(root, 'results', 'gsa');
        if ~exist(gsa_dir_ps, 'dir'), mkdir(gsa_dir_ps); end
        ts_ps     = datestr(now, 'yyyymmdd_HHMMSS');
        gsa_ps_fname = fullfile(gsa_dir_ps, ...
                                sprintf('gsa_pce_post_surgery_%s.mat', ts_ps));
        gsa_ps_save.scenario    = scenario;
        gsa_ps_save.clinical    = clinical;
        gsa_ps_save.params_gsa  = params_for_gsa;
        gsa_ps_save.gsa_cfg     = gsa_ps_cfg;
        gsa_ps_save.gsa_out     = gsa_ps_out;
        save(gsa_ps_fname, '-struct', 'gsa_ps_save');
        fprintf('[main_run] GSA results saved to:\n          %s\n', gsa_ps_fname);

        gsa_summary_ps = make_gsa_summary_table(gsa_ps_out);
        disp(gsa_summary_ps);
        make_gsa_matrix_table(gsa_ps_out, 0.1, true);
    else
        fprintf('\n[main_run] Post-surgery GSA skipped (DO_GSA=false).\n');
    end

    %% -----------------------------------------------------------------
    %  POST-SURGERY STEP 8 — Save calibrated (or prediction) parameters
    %% -----------------------------------------------------------------
    fprintf('\n=== [Post-Surgery Step 8] Saving results (%.1fs elapsed) ===\n', toc(run_timer));
    ts_ps2 = datestr(now, 'yyyymmdd_HHMMSS');
    fname_ps = sprintf('params_calibrated_post_surgery_%s.mat', ts_ps2);
    if DO_CALIBRATION
        save(fullfile(results_dir_ps, fname_ps), 'params_cal', 'calib_out', 'report');
    else
        % Prediction run: save baseline params labelled as prediction
        params_cal = params0;  %#ok<NASGU> % alias for consistent struct field name
        save(fullfile(results_dir_ps, fname_ps), 'params_cal', 'report');
    end
    fprintf('[main_run] Post-surgery results saved to:\n          %s\n', ...
            fullfile(results_dir_ps, fname_ps));

    fprintf('\n[main_run] Done. Scenario: post_surgery  |  Total time: %.1f s\n', toc(run_timer));
    return;   % early return — skip the pre-surgery GSA pipeline below
end

%% =====================================================================
%  STEP 3 — Baseline simulation  (PRE-SURGERY PIPELINE CONTINUES)
%% =====================================================================
fprintf('\n=== [Step 3/10] Baseline simulation (%.1fs elapsed) ===\n', toc(run_timer));
sim_base     = integrate_system(params0);
metrics_base = compute_clinical_indices(sim_base, params0);
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

% ---- Early save: persist calibrated parameters before long post-processing ----
% Saved immediately so params_calibrated_pre_surgery_*.mat exists even if
% Steps 7-10 (GSA, validation, plotting) crash later in the pipeline.
results_dir_early = fullfile(root, 'results', 'tables');
if ~exist(results_dir_early, 'dir'), mkdir(results_dir_early); end
ts_early   = datestr(now, 'yyyymmdd_HHMMSS');
fname_early = sprintf('params_calibrated_pre_surgery_%s.mat', ts_early);
save(fullfile(results_dir_early, fname_early), ...
     'params_cal', 'calib_out', 'optMask', 'sobol_ST_threshold');
fprintf('[main_run] Early parameter save (post-calibration):\n          %s\n', ...
        fullfile(results_dir_early, fname_early));

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
            plot_simulation_comparison(sim_base, params0, sim_cal, params_cal, 'Calibrated');
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
calib_diag.rmse_baseline   = report.rmse_baseline;
calib_diag.rmse_cal        = report.rmse_cal;
calib_diag.table_baseline  = report.table_baseline;   % per-metric baseline errors
calib_diag.table_cal       = report.table_cal;        % per-metric calibrated errors
calib_diag.primary_gate    = report.primary_gate;

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



