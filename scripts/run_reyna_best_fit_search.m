function result = run_reyna_best_fit_search(start_package_path, options)
% RUN_REYNA_BEST_FIT_SEARCH
% -----------------------------------------------------------------------
% Searches for a lower-error Reyna pre-surgery candidate using direct ODE
% simulation, log-space multi-starts, and a 10 percent target-gate polish.
%
% This script is intentionally experimental. It does not overwrite the main
% calibration pipeline. It exports the candidate, metric residuals, and
% parameter movement so a physiologist can inspect whether the improved fit
% is numerically and clinically defensible before promoting the method.
%
% INPUTS:
%   start_package_path - optional saved params_best_candidate MAT path    [-]
%   options            - optional struct controlling search budget        [-]
%
% OUTPUTS:
%   result             - struct with candidate, metrics, and output paths [-]
%
% ASSUMPTIONS:
%   - Reyna's RVEDV target remains a consistency-check target because the
%     stroke-volume audit marks it inconsistent with the catheter flow block.
%   - The 10 percent gate is a polishing target, not a claim that all
%     measurements are accurate to 10 percent.
%
% REFERENCES:
%   [1] Schiavazzi et al. (2017). Patient-specific parameter estimation in
%       single-ventricle lumped circulation models under uncertainty.
%   [2] docs/clinical_data_dictionary.md.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-15
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1
    start_package_path = '';
end
if nargin < 2 || isempty(options)
    options = struct();
end
options = default_options(options);

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

scenario = 'pre_surgery';
clinical = patient_reyna();
case_profile = build_case_calibration_profile(clinical, scenario);

if isempty(start_package_path)
    start_package_path = find_best_reyna_candidate(root);
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_reyna_best_fit_search', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end

fprintf('[run_reyna_best_fit_search] Start candidate:\n  %s\n', start_package_path);
params_start = load_candidate_params(start_package_path);
[params_ref, params_scaled, params_seeded, registry_context] = ...
    build_reyna_reference_context(clinical, scenario, case_profile);

[primary_metrics, primary_selection_table] = ...
    select_primary_metrics(clinical, [], scenario, case_profile);
calib_seed = calibration_param_sets( ...
    scenario, params_seeded, [], primary_metrics, case_profile, registry_context);
targets = build_target_table(clinical, scenario, case_profile);

start_eval = evaluate_candidate(params_start, clinical, scenario, targets);
fprintf(['[run_reyna_best_fit_search] Start primary RMSE %.4f, full RMSE %.4f, ', ...
    'max primary error %.2f%%.\n'], start_eval.rmse_primary, ...
    start_eval.rmse_full, start_eval.max_primary_abs_error_pct);

search_sets = build_search_sets(calib_seed);
best_result = make_round_result('start_candidate', {}, [], Inf, ...
    params_start, start_eval, [], []);

for set_idx = 1:numel(search_sets)
    set_cfg = search_sets{set_idx};
    names = set_cfg.names(ismember(set_cfg.names, calib_seed.names_all));
    names = names(:)';
    if isempty(names)
        continue;
    end

    [lb, ub] = bounds_for_names(calib_seed, names);
    if isfield(set_cfg, 'boundRelaxation')
        [lb, ub] = relax_bounds(lb, ub, names, set_cfg.boundRelaxation);
    end

    round_result = optimize_search_set(set_cfg.label, names, lb, ub, ...
        params_start, params_scaled, case_profile, clinical, scenario, ...
        targets, options, start_eval);

    if is_better_candidate(round_result.eval, best_result.eval)
        best_result = round_result;
    end
end

params_best = best_result.params;
best_eval = evaluate_candidate(params_best, clinical, scenario, targets);
comparison_table = build_comparison_table(start_eval, best_eval, targets);
parameter_table = build_parameter_table(params_start, params_best, ...
    params_scaled, calib_seed, best_result.names, case_profile, ...
    best_result.lb, best_result.ub);

comparison_csv = fullfile(tables_dir, 'reyna_best_fit_search_metrics.csv');
parameter_csv = fullfile(tables_dir, 'reyna_best_fit_search_parameters.csv');
primary_csv = fullfile(tables_dir, 'reyna_best_fit_search_primary_selection.csv');
manifest_path = fullfile(run_dir, 'run_manifest.txt');
mat_path = fullfile(mat_dir, 'reyna_best_fit_search_result.mat');

writetable(comparison_table, comparison_csv);
writetable(parameter_table, parameter_csv);
writetable(primary_selection_table, primary_csv);

result = struct();
result.run_dir = run_dir;
result.start_package_path = start_package_path;
result.round_label = best_result.label;
result.params_start = params_start;
result.params_best = params_best;
result.start_eval = start_eval;
result.best_eval = best_eval;
result.comparison_table = comparison_table;
result.parameter_table = parameter_table;
result.primary_metrics = primary_metrics;
result.primary_selection_table = primary_selection_table;
result.case_profile = case_profile;
result.calib_seed = calib_seed;
result.options = options;
result.paths = struct('comparison_csv', comparison_csv, ...
    'parameter_csv', parameter_csv, 'primary_csv', primary_csv, ...
    'manifest', manifest_path, 'mat', mat_path);

save(mat_path, 'result');
write_manifest(manifest_path, result);

fprintf('\n[run_reyna_best_fit_search] Best round: %s\n', best_result.label);
fprintf('[run_reyna_best_fit_search] Primary RMSE %.4f -> %.4f | full RMSE %.4f -> %.4f\n', ...
    start_eval.rmse_primary, best_eval.rmse_primary, ...
    start_eval.rmse_full, best_eval.rmse_full);
fprintf('[run_reyna_best_fit_search] Exports:\n  %s\n  %s\n  %s\n', ...
    comparison_csv, parameter_csv, manifest_path);
disp(comparison_table);
end

function options = default_options(options)
defaults = struct();
defaults.gateTarget = 0.10;              % [-] desired residual gate
defaults.gateWidth = 0.020;              % [-] hinge sharpness
defaults.nRandomStarts = 14;             % [count] log-space starts
defaults.randomSeed = 1515;              % [-] reproducibility seed
defaults.patternMaxEvaluations = 160;    % [count]
defaults.patternMaxIterations = 55;      % [count]
defaults.fminconMaxEvaluations = 520;    % [count]
defaults.fminconMaxIterations = 130;     % [count]
defaults.fastSimulation = true;          % [-] slightly relaxed SS during search
defaults.verboseOptimizer = false;       % [-]
defaults.objectiveMode = 'gate';         % 'gate' | 'rmse'

fields = fieldnames(defaults);
for idx = 1:numel(fields)
    name = fields{idx};
    if ~isfield(options, name) || isempty(options.(name))
        options.(name) = defaults.(name);
    end
end
end

function search_sets = build_search_sets(calib_seed)
active_names = calib_seed.names_all(calib_seed.mask);
preferred_all = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.LV.EB','V0.LV', ...
    'E.RV.EA','E.RV.EB','V0.RV','E.LA.EA','E.RA.EA','vsd.Cd'};
preferred_core = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.LV.EB','V0.LV', ...
    'E.LA.EA','E.RA.EA','vsd.Cd'};
preferred_chamber = {'E.LV.EA','E.LV.EB','V0.LV', ...
    'E.RV.EA','E.RV.EB','V0.RV','E.LA.EA','E.RA.EA'};

search_sets = {
    struct('label', 'all_active_relaxed_10pct', ...
        'names', {intersect(preferred_all, active_names, 'stable')}, ...
        'boundRelaxation', make_relaxation('moderate'))
    struct('label', 'flow_pressure_chamber_10pct', ...
        'names', {intersect(preferred_core, active_names, 'stable')}, ...
        'boundRelaxation', make_relaxation('moderate'))
    struct('label', 'chamber_atrial_10pct', ...
        'names', {intersect(preferred_chamber, active_names, 'stable')}, ...
        'boundRelaxation', make_relaxation('chamber'))
    };
end

function relaxation = make_relaxation(mode)
switch mode
    case 'chamber'
        relaxation.names = {'E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB', ...
            'V0.LV','V0.RV','E.LA.EA','E.RA.EA'};
        relaxation.lowerFactor = [0.75, 0.80, 0.85, 0.80, ...
            0.85, 0.85, 0.80, 0.80];
        relaxation.upperFactor = [1.10, 1.15, 1.20, 1.20, ...
            1.10, 1.10, 1.20, 1.20];
    otherwise
        relaxation.names = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
            'C.SAR','C.PAR','E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB', ...
            'V0.LV','V0.RV','E.LA.EA','E.RA.EA','vsd.Cd'};
        relaxation.lowerFactor = [0.80, 0.80, 0.80, ...
            0.75, 0.75, 0.75, 0.80, 0.85, 0.80, ...
            0.85, 0.85, 0.80, 0.80, 0.85];
        relaxation.upperFactor = [1.15, 1.25, 1.15, ...
            1.45, 1.45, 1.10, 1.15, 1.25, 1.20, ...
            1.10, 1.10, 1.20, 1.20, 1.15];
end
end

function round_result = optimize_search_set(label, names, lb, ub, ...
    params_start, params_scaled, case_profile, clinical, scenario, ...
    targets, options, start_eval)
% OPTIMIZE_SEARCH_SET - log-space multi-start plus patternsearch/SQP.
fprintf('\n[run_reyna_best_fit_search] Round %s\n', label);
fprintf('  Parameters: %s\n', strjoin(names, ', '));

x0 = values_from_params(params_start, params_scaled, names, case_profile);
x0 = min(max(x0(:), lb(:)), ub(:));
starts = build_log_space_starts(x0, lb(:), ub(:), options);

best_x = x0;
best_f = Inf;
for start_idx = 1:size(starts, 2)
    x_init = starts(:, start_idx);
    obj = @(x) search_objective(x, params_start, params_scaled, names, ...
        case_profile, clinical, scenario, targets, start_eval, ...
        options, lb(:), ub(:));

    [x_candidate, f_candidate] = local_solve(obj, x_init, lb(:), ub(:), options);
    if f_candidate < best_f
        best_f = f_candidate;
        best_x = x_candidate(:);
        candidate_params = apply_values(params_start, params_scaled, ...
            names, best_x, case_profile);
        candidate_eval = evaluate_candidate(candidate_params, clinical, ...
            scenario, targets);
        fprintf('  New best from start %02d: J %.4g | primary RMSE %.4f | max %.2f%%\n', ...
            start_idx, best_f, candidate_eval.rmse_primary, ...
            candidate_eval.max_primary_abs_error_pct);
    end
end

params_round = apply_values(params_start, params_scaled, names, best_x, case_profile);
round_eval = evaluate_candidate(params_round, clinical, scenario, targets);
round_result = make_round_result(label, names, best_x, best_f, ...
    params_round, round_eval, lb(:), ub(:));
fprintf('  Finished %s: primary RMSE %.4f, full RMSE %.4f, gate failures %d.\n', ...
    label, round_eval.rmse_primary, round_eval.rmse_full, ...
    round_eval.n_primary_over_gate);
end

function [x_best, f_best] = local_solve(obj, x0, lb, ub, options)
x_best = x0(:);
f_best = obj(x_best);
display_mode = 'off';
if options.verboseOptimizer
    display_mode = 'iter';
end

if exist('patternsearch', 'file') == 2
    try
        ps_opts = optimoptions('patternsearch', ...
            'Display', display_mode, ...
            'MaxFunctionEvaluations', options.patternMaxEvaluations, ...
            'MaxIterations', options.patternMaxIterations, ...
            'InitialMeshSize', 0.08, ...
            'MeshTolerance', 1e-3, ...
            'UseCompletePoll', true, ...
            'UseCompleteSearch', false);
        [x_ps, f_ps] = patternsearch(obj, x_best, [], [], [], [], ...
            lb, ub, [], ps_opts);
        if isfinite(f_ps) && f_ps < f_best
            x_best = x_ps(:);
            f_best = f_ps;
        end
    catch ME
        fprintf('  patternsearch skipped: %s\n', ME.message);
    end
end

if exist('fmincon', 'file') == 2
    try
        span = max(ub - lb, 1e-12);
        fc_opts = optimoptions('fmincon', ...
            'Algorithm', 'sqp', ...
            'Display', display_mode, ...
            'FiniteDifferenceType', 'central', ...
            'FiniteDifferenceStepSize', max(1e-4, 0.0015 * median(span)), ...
            'MaxFunctionEvaluations', options.fminconMaxEvaluations, ...
            'MaxIterations', options.fminconMaxIterations, ...
            'OptimalityTolerance', 1e-5, ...
            'StepTolerance', 1e-5);
        [x_fc, f_fc] = fmincon(obj, x_best, [], [], [], [], ...
            lb, ub, [], fc_opts);
        if isfinite(f_fc) && f_fc < f_best
            x_best = x_fc(:);
            f_best = f_fc;
        end
    catch ME
        fprintf('  fmincon skipped: %s\n', ME.message);
    end
end
end

function starts = build_log_space_starts(x0, lb, ub, options)
rng(options.randomSeed, 'twister');
n_params = numel(x0);
n_starts = max(1, round(options.nRandomStarts));
starts = zeros(n_params, n_starts);
starts(:, 1) = x0(:);

log_lb = log(max(lb(:), 1e-12));
log_ub = log(max(ub(:), 1e-12));
log_x0 = log(max(x0(:), 1e-12));
for start_idx = 2:n_starts
    jitter = 0.35 * randn(n_params, 1);
    if mod(start_idx, 4) == 0
        raw = rand(n_params, 1);
        trial = log_lb + raw .* (log_ub - log_lb);
    else
        trial = log_x0 + jitter .* max(log_ub - log_lb, 1e-6);
    end
    starts(:, start_idx) = exp(min(max(trial, log_lb), log_ub));
end
end

function J = search_objective(x, params_base, params_scaled, names, ...
    case_profile, clinical, scenario, targets, start_eval, options, lb, ub)
params = apply_values(params_base, params_scaled, names, x(:), case_profile);
if options.fastSimulation && isfield(params, 'sim')
    params.sim.nCyclesSteady = min(max(params.sim.nCyclesSteady, 28), 55);
    params.sim.ss_tol_P = max(params.sim.ss_tol_P, 0.35);
    params.sim.ss_tol_V = max(params.sim.ss_tol_V, 0.35);
end

eval = evaluate_candidate(params, clinical, scenario, targets);
if ~eval.valid
    J = 1e8;
    return;
end

gate = options.gateTarget;
gate_width = options.gateWidth;
primary_rel = abs(eval.primary_error_pct(:)) / 100;
calibration_rel = abs(eval.calibration_error_pct(:)) / 100;
soft_rel = abs(eval.soft_error_pct(:)) / 100;
consistency_rel = abs(eval.consistency_error_pct(:)) / 100;

if strcmpi(char(options.objectiveMode), 'rmse')
    J_gate = 0.35 * sum(eval.primary_weights(:) .* ...
        (max(0, primary_rel - gate) / gate_width).^2);
    J_primary = 95.0 * sum(eval.primary_weights(:) .* primary_rel.^2);
    J_cal = 14.0 * sum(eval.calibration_weights(:) .* calibration_rel.^2);
    J_soft = 2.5 * sum(soft_rel.^2);
    J_consistency = 0.10 * sum(min(consistency_rel, 0.35).^2);
else
    J_gate = sum(eval.primary_weights(:) .* ...
        (max(0, primary_rel - gate) / gate_width).^2);
    J_primary = 12.0 * sum(eval.primary_weights(:) .* primary_rel.^2);
    J_cal = 3.5 * sum(eval.calibration_weights(:) .* calibration_rel.^2);
    J_soft = 0.9 * sum(soft_rel.^2);
    J_consistency = 0.18 * sum(min(consistency_rel, 0.35).^2);
end

% Keep already-good pressure-flow anchors from degrading while the search
% trades chamber-volume and function errors.
protected = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin','LVEDV'};
protected_now = abs_error_for_metrics(eval.error_table, protected) / 100;
protected_start = abs_error_for_metrics(start_eval.error_table, protected) / 100;
protected_ceiling = min(0.12, max(0.055, protected_start + 0.020));
if strcmpi(char(options.objectiveMode), 'rmse')
    J_protect = 2.5 * sum((max(0, protected_now - protected_ceiling) / 0.020).^2);
else
    J_protect = 18.0 * sum((max(0, protected_now - protected_ceiling) / 0.020).^2);
end

% Stronger walls for the current dominant residuals.
if strcmpi(char(options.objectiveMode), 'rmse')
    J_walls = metric_wall(eval.error_table, 'LVEF', gate, 4.0) + ...
        metric_wall(eval.error_table, 'LVESV', gate, 3.0) + ...
        metric_wall(eval.error_table, 'LVEDV', gate, 2.0) + ...
        metric_wall(eval.error_table, 'LAP_mean', gate, 2.0) + ...
        metric_wall(eval.error_table, 'RVESV', gate, 1.5);
else
    J_walls = metric_wall(eval.error_table, 'LVEF', gate, 22.0) + ...
        metric_wall(eval.error_table, 'LVESV', gate, 18.0) + ...
        metric_wall(eval.error_table, 'LVEDV', gate, 12.0) + ...
        metric_wall(eval.error_table, 'LAP_mean', gate, 10.0) + ...
        metric_wall(eval.error_table, 'RVESV', gate, 8.0);
end

% Log-space local regularization and boundary avoidance.
x = x(:);
lb = lb(:);
ub = ub(:);
x0 = values_from_params(params_base, params_scaled, names, case_profile);
x0 = min(max(x0(:), lb), ub);
J_reg = 0.10 * sum((log(max(x, 1e-12) ./ max(x0, 1e-12))).^2);
span = max(ub - lb, 1e-12);
dist_lower = (x - lb) ./ span;
dist_upper = (ub - x) ./ span;
J_bound = 8.0 * sum(max(0, 0.045 - dist_lower).^2 + ...
    max(0, 0.045 - dist_upper).^2);

J_validity = 0;
if ~eval.physiology_valid
    J_validity = 1e5 + 1e4 * numel(eval.validity.failed_flags);
end

J = J_gate + J_primary + J_cal + J_soft + J_consistency + ...
    J_protect + J_walls + J_reg + J_bound + J_validity;
end

function penalty = metric_wall(error_table, metric_name, gate, weight)
err_rel = abs(signed_error_for_metric(error_table, metric_name)) / 100;
if ~isfinite(err_rel)
    penalty = 0;
else
    penalty = weight * (max(0, err_rel - gate) / 0.020)^2;
end
end

function eval = evaluate_candidate(params, clinical, scenario, targets)
% EVALUATE_CANDIDATE - simulate, score target errors, and check validity.
eval = struct();
eval.valid = false;
eval.physiology_valid = false;
eval.metrics = struct();
eval.validity = struct('is_valid', false, 'failed_flags', {{}});
eval.rmse_primary = Inf;
eval.rmse_full = Inf;
eval.rmse_hard = Inf;
eval.rmse_soft = Inf;
eval.error_table = table();
eval.primary_error_pct = Inf(0, 1);
eval.calibration_error_pct = Inf(0, 1);
eval.soft_error_pct = Inf(0, 1);
eval.consistency_error_pct = Inf(0, 1);
eval.primary_weights = Inf(0, 1);
eval.calibration_weights = Inf(0, 1);
eval.max_primary_abs_error_pct = Inf;
eval.n_primary_over_gate = Inf;

try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
    validity = evaluate_simulation_validity(sim, params, metrics, scenario, clinical);
catch
    return;
end

model_col = nan(height(targets), 1);
for idx = 1:height(targets)
    metric_name = targets.Metric{idx};
    if isfield(metrics, metric_name) && isfinite(metrics.(metric_name))
        model_col(idx) = metrics.(metric_name);
    end
end

clinical_col = targets.Clinical;
valid = isfinite(clinical_col) & isfinite(model_col);
error_pct_col = nan(height(targets), 1);
error_pct_col(valid) = 100 * (model_col(valid) - clinical_col(valid)) ./ ...
    max(abs(clinical_col(valid)), 1e-9);

eval.valid = any(valid);
eval.physiology_valid = validity.is_valid;
eval.metrics = metrics;
eval.validity = validity;
eval.error_table = table(targets.Metric, targets.Unit, clinical_col, ...
    model_col, error_pct_col, targets.Tier, targets.IncludedInCalibration, ...
    targets.IncludedInPrimaryRMSE, targets.Weight, ...
    'VariableNames', {'Metric','Unit','Clinical','Model','Error_pct', ...
    'Tier','IncludedInCalibration','IncludedInPrimaryRMSE','Weight'});

primary_mask = valid & targets.IncludedInPrimaryRMSE;
cal_mask = valid & targets.IncludedInCalibration;
hard_mask = valid & strcmp(targets.Tier, 'hard');
soft_mask = valid & strcmp(targets.Tier, 'soft');
consistency_mask = valid & strcmp(targets.Tier, 'consistency_check_only');

eval.rmse_primary = rmse_for_mask(error_pct_col, primary_mask);
eval.rmse_full = rmse_for_mask(error_pct_col, valid);
eval.rmse_hard = rmse_for_mask(error_pct_col, hard_mask);
eval.rmse_soft = rmse_for_mask(error_pct_col, soft_mask);
eval.primary_error_pct = error_pct_col(primary_mask);
eval.calibration_error_pct = error_pct_col(cal_mask);
eval.soft_error_pct = error_pct_col(soft_mask);
eval.consistency_error_pct = error_pct_col(consistency_mask);
eval.primary_weights = targets.Weight(primary_mask);
eval.calibration_weights = targets.Weight(cal_mask);
eval.max_primary_abs_error_pct = max(abs(eval.primary_error_pct));
eval.n_primary_over_gate = sum(abs(eval.primary_error_pct) > 10.0);
end

function value = rmse_for_mask(error_pct_col, mask)
if any(mask)
    value = sqrt(mean((error_pct_col(mask) / 100).^2));
else
    value = NaN;
end
end

function targets = build_target_table(clinical, scenario, case_profile)
target_struct = get_calibration_targets(scenario, clinical);
target_tiers = case_profile.targetTiers.table;
n_targets = numel(target_struct);

metric_col = {target_struct.Metric}';
clinical_col = [target_struct.ClinicalValue]';
unit_col = {target_struct.Unit}';
tier_col = repmat({'unavailable'}, n_targets, 1);
included_cal = false(n_targets, 1);
included_primary = false(n_targets, 1);
weight_col = ones(n_targets, 1);

for idx = 1:n_targets
    hit = find(strcmp(target_tiers.Metric, metric_col{idx}), 1, 'first');
    if ~isempty(hit)
        tier_col{idx} = target_tiers.Tier{hit};
        included_cal(idx) = logical(target_tiers.IncludedInCalibration(hit));
        included_primary(idx) = logical(target_tiers.IncludedInPrimaryRMSE(hit));
    end
    weight_col(idx) = target_weight(metric_col{idx}, tier_col{idx}, case_profile);
end

targets = table(metric_col, clinical_col, unit_col, tier_col, ...
    included_cal, included_primary, weight_col, ...
    'VariableNames', {'Metric','Clinical','Unit','Tier', ...
    'IncludedInCalibration','IncludedInPrimaryRMSE','Weight'});
end

function weight = target_weight(metric_name, tier, case_profile)
switch tier
    case 'hard'
        weight = 1.00;
    case 'soft'
        weight = 0.45;
    case 'validation_only'
        weight = 0.35;
    otherwise
        weight = 0.15;
end
if isfield(case_profile, 'metricWeightOverrides') && ...
        isfield(case_profile.metricWeightOverrides, metric_name)
    weight = weight * case_profile.metricWeightOverrides.(metric_name);
end
if isfield(case_profile, 'targetTiers') && ...
        isfield(case_profile.targetTiers, 'weights') && ...
        isfield(case_profile.targetTiers.weights, metric_name)
    weight = weight * case_profile.targetTiers.weights.(metric_name);
end
end

function tf = is_better_candidate(candidate_eval, incumbent_eval)
if candidate_eval.rmse_primary < incumbent_eval.rmse_primary - 1e-5
    tf = true;
    return;
end
if candidate_eval.rmse_primary > incumbent_eval.rmse_primary + 1e-5
    tf = false;
    return;
end
if candidate_eval.n_primary_over_gate < incumbent_eval.n_primary_over_gate
    tf = true;
    return;
end
if candidate_eval.n_primary_over_gate > incumbent_eval.n_primary_over_gate
    tf = false;
    return;
end
tf = candidate_eval.rmse_full < incumbent_eval.rmse_full;
end

function tbl = build_comparison_table(start_eval, best_eval, targets)
start_tbl = start_eval.error_table;
best_tbl = best_eval.error_table;
[is_common, idx_best] = ismember(start_tbl.Metric, best_tbl.Metric);
idx_start = find(is_common);
idx_best = idx_best(is_common);

start_abs = abs(start_tbl.Error_pct(idx_start));
best_abs = abs(best_tbl.Error_pct(idx_best));
tbl = table(start_tbl.Metric(idx_start), start_tbl.Unit(idx_start), ...
    start_tbl.Clinical(idx_start), start_tbl.Model(idx_start), ...
    start_tbl.Error_pct(idx_start), best_tbl.Model(idx_best), ...
    best_tbl.Error_pct(idx_best), start_abs, best_abs, ...
    best_abs - start_abs, best_abs <= 10.0, ...
    start_tbl.Tier(idx_start), start_tbl.IncludedInCalibration(idx_start), ...
    start_tbl.IncludedInPrimaryRMSE(idx_start), ...
    'VariableNames', {'Metric','Unit','Clinical','StartModel', ...
    'StartError_pct','BestModel','BestError_pct','StartAbsError_pct', ...
    'BestAbsError_pct','DeltaAbsError_pct','Pass_10pct','Tier', ...
    'IncludedInCalibration','IncludedInPrimaryRMSE'});

% Preserve target order from get_calibration_targets.
if height(tbl) ~= height(targets)
    return;
end
end

function tbl = build_parameter_table(params_start, params_best, reference_params, ...
    calib_seed, names, case_profile, used_lb, used_ub)
if isempty(names)
    tbl = table();
    return;
end
registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
start_values = values_from_params(params_start, reference_params, names, case_profile);
best_values = values_from_params(params_best, reference_params, names, case_profile);
if nargin < 7 || isempty(used_lb)
    used_lb = registry.lb(loc);
end
if nargin < 8 || isempty(used_ub)
    used_ub = registry.ub(loc);
end
tbl = table(names(:), registry.baseline_scaled(loc), start_values(:), ...
    best_values(:), best_values(:) ./ max(abs(start_values(:)), 1e-12), ...
    registry.lb(loc), registry.ub(loc), used_lb(:), used_ub(:), ...
    'VariableNames', {'Parameter','BaselineScaled','StartValue', ...
    'BestValue','RatioToStart','RegistryLowerBound','RegistryUpperBound', ...
    'UsedLowerBound','UsedUpperBound'});
end

function [lb, ub] = bounds_for_names(calib_seed, names)
registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
if any(loc == 0)
    missing = names(loc == 0);
    error('run_reyna_best_fit_search:missingRegistryName', ...
        'Missing registry names: %s', strjoin(missing, ', '));
end
lb = registry.lb(loc);
ub = registry.ub(loc);
end

function [lb, ub] = relax_bounds(lb, ub, names, relaxation)
for idx = 1:numel(relaxation.names)
    hit = strcmp(names, relaxation.names{idx});
    if ~any(hit)
        continue;
    end
    lb(hit) = lb(hit) * relaxation.lowerFactor(idx);
    ub(hit) = ub(hit) * relaxation.upperFactor(idx);
end
lb = max(lb, 1e-9);
ub = max(ub, lb * 1.001);
end

function values = values_from_params(params, reference_params, names, case_profile)
values = nan(numel(names), 1);
for idx = 1:numel(names)
    values(idx) = get_calibration_param_value( ...
        params, reference_params, names{idx}, case_profile);
end
end

function params = apply_values(params_base, reference_params, names, values, case_profile)
params = params_base;
for idx = 1:numel(names)
    params = set_calibration_param_value( ...
        params, reference_params, names{idx}, values(idx), case_profile);
end
end

function abs_error_pct = abs_error_for_metrics(error_table, metric_names)
abs_error_pct = nan(numel(metric_names), 1);
for idx = 1:numel(metric_names)
    hit = strcmp(error_table.Metric, metric_names{idx});
    if any(hit)
        abs_error_pct(idx) = abs(error_table.Error_pct(find(hit, 1, 'first')));
    end
end
abs_error_pct(~isfinite(abs_error_pct)) = 0;
end

function error_pct = signed_error_for_metric(error_table, metric_name)
error_pct = NaN;
hit = strcmp(error_table.Metric, metric_name);
if any(hit)
    error_pct = error_table.Error_pct(find(hit, 1, 'first'));
end
end

function round_result = make_round_result(label, names, x, objective, ...
    params, eval, lb, ub)
round_result = struct('label', label, 'names', {names}, 'x', x, ...
    'objective', objective, 'params', params, 'eval', eval, ...
    'lb', lb, 'ub', ub);
end

function params = load_candidate_params(start_package_path)
data = load(start_package_path);
if isfield(data, 'best_candidate') && isfield(data.best_candidate, 'params')
    params = data.best_candidate.params;
elseif isfield(data, 'scientific_candidate') && isfield(data.scientific_candidate, 'params')
    params = data.scientific_candidate.params;
elseif isfield(data, 'accepted_candidate') && isfield(data.accepted_candidate, 'params')
    params = data.accepted_candidate.params;
elseif isfield(data, 'params_cal')
    params = data.params_cal;
elseif isfield(data, 'params_best')
    params = data.params_best;
else
    error('run_reyna_best_fit_search:missingParams', ...
        'No candidate parameter struct found in %s.', start_package_path);
end
end

function best_candidate_path = find_best_reyna_candidate(root)
listing = dir(fullfile(root, 'results', 'runs', '*reyna_pre_surgery', ...
    'mat', 'params_best_candidate_pre_surgery.mat'));
if isempty(listing)
    error('run_reyna_best_fit_search:noBestCandidate', ...
        'No Reyna best-candidate package found under results/runs.');
end

best_rmse = Inf;
best_idx = 1;
for idx = 1:numel(listing)
    run_dir = fileparts(fileparts(listing(idx).folder));
    manifest_path = fullfile(run_dir, 'run_manifest.txt');
    candidate_rmse = parse_manifest_number(manifest_path, 'BestCandidateRMSE');
    if ~isfinite(candidate_rmse)
        candidate_rmse = listing(idx).datenum;
    end
    if candidate_rmse < best_rmse
        best_rmse = candidate_rmse;
        best_idx = idx;
    end
end
best_candidate_path = fullfile(listing(best_idx).folder, listing(best_idx).name);
end

function value = parse_manifest_number(path, key)
value = NaN;
if exist(path, 'file') ~= 2
    return;
end
lines = readlines(path);
prefix = sprintf('%s:', key);
hit = startsWith(strtrim(lines), prefix);
if ~any(hit)
    return;
end
text = extractAfter(strtrim(lines(find(hit, 1, 'first'))), strlength(prefix));
value = str2double(strtrim(text));
end

function [params_ref, params_scaled, params_seeded, registry_context] = ...
    build_reyna_reference_context(clinical, scenario, case_profile)
params_ref = default_parameters();
patient = struct();
patient.age_years = clinical.common.age_years;
patient.age_days = clinical.common.age_years * 365.25;
patient.weight_kg = clinical.common.weight_kg;
patient.height_cm = clinical.common.height_cm;
patient.sex = clinical.common.sex;
patient.maturation_mode = 'normal';
patient.scaling_mode = case_profile.preferredScalingMode;
if isfield(clinical.common, 'BSA') && isfinite(clinical.common.BSA)
    patient.BSA = clinical.common.BSA;
end
params_scaled = apply_scaling(params_ref, patient);
params_seeded = params_from_clinical(params_scaled, clinical, ...
    scenario, params_scaled, case_profile);
registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
end

function write_manifest(path, result)
fid = fopen(path, 'w');
if fid < 0
    error('run_reyna_best_fit_search:manifestOpenFailed', ...
        'Unable to write manifest: %s', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Reyna Best Fit Search\n');
fprintf(fid, '=====================\n');
fprintf(fid, 'StartPackage: %s\n', result.start_package_path);
fprintf(fid, 'RoundLabel: %s\n', result.round_label);
fprintf(fid, 'GateTarget: %.4f\n', result.options.gateTarget);
fprintf(fid, 'ObjectiveMode: %s\n', char(result.options.objectiveMode));
fprintf(fid, 'RandomStarts: %d\n', result.options.nRandomStarts);
fprintf(fid, 'StartPrimaryRMSE: %.6f\n', result.start_eval.rmse_primary);
fprintf(fid, 'BestPrimaryRMSE: %.6f\n', result.best_eval.rmse_primary);
fprintf(fid, 'StartFullRMSE: %.6f\n', result.start_eval.rmse_full);
fprintf(fid, 'BestFullRMSE: %.6f\n', result.best_eval.rmse_full);
fprintf(fid, 'StartHardRMSE: %.6f\n', result.start_eval.rmse_hard);
fprintf(fid, 'BestHardRMSE: %.6f\n', result.best_eval.rmse_hard);
fprintf(fid, 'StartSoftRMSE: %.6f\n', result.start_eval.rmse_soft);
fprintf(fid, 'BestSoftRMSE: %.6f\n', result.best_eval.rmse_soft);
fprintf(fid, 'StartMaxPrimaryAbsErrorPct: %.6f\n', ...
    result.start_eval.max_primary_abs_error_pct);
fprintf(fid, 'BestMaxPrimaryAbsErrorPct: %.6f\n', ...
    result.best_eval.max_primary_abs_error_pct);
fprintf(fid, 'StartPrimaryFailuresOver10pct: %d\n', ...
    result.start_eval.n_primary_over_gate);
fprintf(fid, 'BestPrimaryFailuresOver10pct: %d\n', ...
    result.best_eval.n_primary_over_gate);
fprintf(fid, 'PhysiologyValid: %d\n', result.best_eval.physiology_valid);
if ~result.best_eval.physiology_valid
    fprintf(fid, 'FailedValidityFlags: %s\n', ...
        strjoin(result.best_eval.validity.failed_flags, ', '));
end
fprintf(fid, 'ComparisonCSV: %s\n', result.paths.comparison_csv);
fprintf(fid, 'ParameterCSV: %s\n', result.paths.parameter_csv);
fprintf(fid, 'MatFile: %s\n', result.paths.mat);
end

function project_path = build_clean_project_path(root)
project_paths = strsplit(genpath(root), pathsep);
project_paths = project_paths(~cellfun('isempty', project_paths));
is_shadow = contains(project_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
    contains(project_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
    contains(project_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, project_paths);
project_path = strjoin(project_paths(~is_shadow & is_existing), pathsep);
end
