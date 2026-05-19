function result = run_reyna_15pct_gate_polish(start_package_path)
% RUN_REYNA_15PCT_GATE_POLISH
% -----------------------------------------------------------------------
% Polishes the saved Reyna scientific candidate against a 15% per-metric
% gate. This is a post-hoc, direct-ODE experiment: it starts from the best
% saved candidate, penalizes only clinically measured targets that remain
% outside the 15% band, and protects already acceptable pressure-flow
% targets from degradation.
%
% INPUTS:
%   start_package_path - optional params_best_candidate MAT file          [-]
%
% OUTPUTS:
%   result             - struct with output paths, metrics, and params    [-]
%
% REFERENCES:
%   [1] docs/reyna_pak_dipo_systemic_flow_revision.md
%   [2] docs/clinical_data_dictionary.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-12
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

scenario = 'pre_surgery';
clinical = patient_reyna();
case_profile = build_case_calibration_profile(clinical, scenario);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_reyna_15pct_gate_polish', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
mkdir(tables_dir);
mkdir(mat_dir);

if nargin < 1 || isempty(start_package_path)
    start_package_path = find_latest_reyna_best_candidate(root);
end

fprintf('[run_reyna_15pct_gate_polish] Loading start candidate:\n  %s\n', start_package_path);
loaded = load(start_package_path, 'best_candidate');
if ~isfield(loaded, 'best_candidate')
    error('run_reyna_15pct_gate_polish:missingCandidate', ...
        'MAT file does not contain best_candidate.');
end

params_start = loaded.best_candidate.params;
[params_ref, params_scaled, params_seeded, registry_context] = ...
    build_reyna_reference_context(clinical, scenario, case_profile); %#ok<ASGLU>
[primary_metrics, ~] = select_primary_metrics(clinical, [], scenario, case_profile);
calib_seed = calibration_param_sets( ...
    scenario, params_seeded, [], primary_metrics, case_profile, registry_context);

important_metrics = {'CO_Lmin','SAP_max','SAP_min','RVEDV','LVEF','LVESV'};
protected_metrics = {'PAP_mean','QpQs','RAP_mean','SAP_mean','SVR', ...
    'LVEDV','RVESV','LAP_mean'};
targets = build_target_table(clinical, scenario);

start_eval = evaluate_candidate(params_start, clinical, scenario, targets, ...
    important_metrics, protected_metrics);
fprintf('[run_reyna_15pct_gate_polish] Start RMSE %.4f, max important error %.2f%%, gate failures %d.\n', ...
    start_eval.rmse, start_eval.max_important_abs_error_pct, start_eval.n_important_over_gate);

base_names = {'group.R_sys_scale','R.SVEN','C.SAR', ...
    'E.LV.EA','E.LV.EB','V0.LV', ...
    'E.RV.EA','E.RV.EB','V0.RV'};
candidate_sets = {
    struct('label', 'requested_9_parameter_set', 'names', {base_names})
    struct('label', 'requested_set_plus_vsd_cd', 'names', {[base_names, {'vsd.Cd'}]})
    struct('label', 'rv_bound_relaxation_probe', ...
        'names', {[base_names, {'vsd.Cd'}]}, ...
        'boundRelaxations', make_rv_bound_relaxation())
    };

best_result = make_empty_round_result();
best_result.params = params_start;
best_result.eval = start_eval;
best_result.label = 'start_candidate';
best_result.x = [];
best_result.names = {};
best_result.objective = Inf;

for set_idx = 1:numel(candidate_sets)
    if set_idx > 1 && best_result.eval.n_measured_over_gate == 0
        break;
    end

    set_cfg = candidate_sets{set_idx};
    names = set_cfg.names(ismember(set_cfg.names, calib_seed.names_all));
    if isempty(names)
        warning('run_reyna_15pct_gate_polish:noNames', ...
            'Candidate set %s has no active parameters after registry filtering.', set_cfg.label);
        continue;
    end

    bound_relaxations = struct();
    if isfield(set_cfg, 'boundRelaxations')
        bound_relaxations = set_cfg.boundRelaxations;
    end

    round_result = optimize_gate_round( ...
        set_cfg.label, names, params_start, params_scaled, case_profile, ...
        clinical, scenario, calib_seed, targets, important_metrics, ...
        protected_metrics, start_eval, bound_relaxations);

    if is_better_gate_result(round_result.eval, best_result.eval)
        best_result = round_result;
    end
end

params_polished = best_result.params;
polished_eval = evaluate_candidate(params_polished, clinical, scenario, targets, ...
    important_metrics, protected_metrics);

comparison_tbl = build_comparison_table(start_eval, polished_eval, targets);
parameter_tbl = build_parameter_table(params_start, params_polished, params_scaled, ...
    calib_seed, best_result.names, case_profile, best_result.lb, best_result.ub);

comparison_csv = fullfile(tables_dir, 'reyna_15pct_gate_polish_comparison.csv');
parameter_csv = fullfile(tables_dir, 'reyna_15pct_gate_polish_parameters.csv');
manifest_path = fullfile(run_dir, 'run_manifest.txt');
mat_path = fullfile(mat_dir, 'reyna_15pct_gate_polish_result.mat');

writetable(comparison_tbl, comparison_csv);
writetable(parameter_tbl, parameter_csv);

result = struct();
result.run_dir = run_dir;
result.start_package_path = start_package_path;
result.round_label = best_result.label;
result.params_start = params_start;
result.params_polished = params_polished;
result.start_eval = start_eval;
result.polished_eval = polished_eval;
result.comparison_table = comparison_tbl;
result.parameter_table = parameter_tbl;
result.important_metrics = important_metrics;
result.protected_metrics = protected_metrics;
result.primary_metrics = primary_metrics;
result.case_profile = case_profile;
result.calib_seed = calib_seed;
result.paths = struct('comparison_csv', comparison_csv, ...
    'parameter_csv', parameter_csv, 'manifest', manifest_path, ...
    'mat', mat_path);

save(mat_path, 'result');
write_manifest(manifest_path, result);

fprintf('\n[run_reyna_15pct_gate_polish] Best round: %s\n', best_result.label);
fprintf('[run_reyna_15pct_gate_polish] Polished RMSE %.4f, max measured error %.2f%%, measured gate failures %d.\n', ...
    polished_eval.rmse, polished_eval.max_measured_abs_error_pct, ...
    polished_eval.n_measured_over_gate);
fprintf('[run_reyna_15pct_gate_polish] Exported:\n  %s\n  %s\n  %s\n', ...
    comparison_csv, parameter_csv, manifest_path);
disp(comparison_tbl);
end

function round_result = optimize_gate_round(label, names, params_start, params_scaled, ...
    case_profile, clinical, scenario, calib_seed, targets, important_metrics, ...
    protected_metrics, start_eval, bound_relaxations)
% OPTIMIZE_GATE_ROUND - pattern-search then fmincon local polish.
registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
if any(loc == 0)
    error('run_reyna_15pct_gate_polish:missingRegistryName', ...
        'At least one polish parameter is missing from the registry.');
end

lb = registry.lb(loc);
ub = registry.ub(loc);
if nargin >= 12 && isstruct(bound_relaxations) && isfield(bound_relaxations, 'names')
    [lb, ub] = apply_bound_relaxations(lb, ub, names, bound_relaxations);
end
x0 = values_from_params(params_start, params_scaled, names, case_profile);
x0 = min(max(x0(:), lb(:)), ub(:));
span = max(ub(:) - lb(:), 1e-12);

fprintf('\n[run_reyna_15pct_gate_polish] Round %s\n', label);
fprintf('  Active names: %s\n', strjoin(names, ', '));

obj = @(x) gate_objective(x, params_start, params_scaled, names, case_profile, ...
    clinical, scenario, targets, important_metrics, protected_metrics, ...
    start_eval, x0, lb(:), ub(:));

best_x = x0;
best_f = obj(x0);

if exist('patternsearch', 'file') == 2
    ps_opts = optimoptions('patternsearch', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 220, ...
        'MaxIterations', 60, ...
        'InitialMeshSize', 0.10, ...
        'MeshTolerance', 1e-3, ...
        'UseCompletePoll', true, ...
        'UseCompleteSearch', false);
    try
        [x_ps, f_ps] = patternsearch(obj, x0, [], [], [], [], lb(:), ub(:), [], ps_opts);
        if f_ps < best_f
            best_x = x_ps(:);
            best_f = f_ps;
        end
    catch ME
        fprintf('[run_reyna_15pct_gate_polish] patternsearch skipped: %s\n', ME.message);
    end
end

if exist('fmincon', 'file') == 2
    fc_opts = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'FiniteDifferenceType', 'central', ...
        'FiniteDifferenceStepSize', max(1e-4, 0.002 * median(span)), ...
        'MaxFunctionEvaluations', 700, ...
        'MaxIterations', 120, ...
        'OptimalityTolerance', 1e-5, ...
        'StepTolerance', 1e-5);
    try
        [x_fc, f_fc] = fmincon(obj, best_x, [], [], [], [], lb(:), ub(:), [], fc_opts);
        if f_fc < best_f
            best_x = x_fc(:);
            best_f = f_fc;
        end
    catch ME
        fprintf('[run_reyna_15pct_gate_polish] fmincon skipped: %s\n', ME.message);
    end
end

params_round = apply_values(params_start, params_scaled, names, best_x, case_profile);
round_eval = evaluate_candidate(params_round, clinical, scenario, targets, ...
    important_metrics, protected_metrics);

round_result = make_empty_round_result();
round_result.label = label;
round_result.names = names;
round_result.x = best_x;
round_result.objective = best_f;
round_result.params = params_round;
round_result.eval = round_eval;
round_result.lb = lb(:);
round_result.ub = ub(:);
fprintf('[run_reyna_15pct_gate_polish] Round %s finished: RMSE %.4f, max measured %.2f%%, failures %d.\n', ...
    label, round_eval.rmse, round_eval.max_measured_abs_error_pct, ...
    round_eval.n_measured_over_gate);
end

function J = gate_objective(x, params_base, params_scaled, names, case_profile, ...
    clinical, scenario, targets, important_metrics, protected_metrics, ...
    start_eval, x0, lb, ub)
% GATE_OBJECTIVE - hinge objective focused on the 15% clinical gate.
params = apply_values(params_base, params_scaled, names, x(:), case_profile);
params.sim.nCyclesSteady = min(max(params.sim.nCyclesSteady, 30), 60);
params.sim.ss_tol_P = max(params.sim.ss_tol_P, 0.5);
params.sim.ss_tol_V = max(params.sim.ss_tol_V, 0.5);

eval = evaluate_candidate(params, clinical, scenario, targets, ...
    important_metrics, protected_metrics);
if ~eval.valid
    J = 1e6;
    return;
end

important_weight = metric_weight_vector(important_metrics, struct( ...
    'CO_Lmin', 18.0, 'SAP_max', 7.0, 'SAP_min', 8.0, ...
    'RVEDV', 10.0, 'LVEF', 5.5, 'LVESV', 4.0));
protected_weight = metric_weight_vector(protected_metrics, struct( ...
    'PAP_mean', 6.0, 'QpQs', 7.0, 'RAP_mean', 5.0, 'SAP_mean', 6.0, ...
    'SVR', 6.0, 'LVEDV', 3.5, 'RVESV', 3.0, 'LAP_mean', 2.0));

gate = 0.15;                    % [-] 15% target
gate_width = 0.025;             % [-] sharpness around target
important_rel = eval.important_abs_error_pct(:) / 100;
important_excess = max(0, important_rel - gate);
J = sum(important_weight(:) .* (important_excess / gate_width).^2);

% Tie-breaker: after crossing the gate, prefer smaller residuals, but gently.
J = J + 0.04 * sum(important_weight(:) .* (important_rel / gate).^2);

measured_abs_rel = abs(eval.error_table.Error_pct(isfinite(eval.error_table.Error_pct))) / 100;
measured_excess = max(0, measured_abs_rel - gate);
J = J + 2.0 * sum((measured_excess / gate_width).^2);

co_error_rel = signed_error_for_metric(eval.error_table, 'CO_Lmin') / 100;
if isfinite(co_error_rel) && co_error_rel < -gate
    J = J + 45.0 * ((-gate - co_error_rel) / 0.020)^2;
end

protected_rel = eval.protected_abs_error_pct(:) / 100;
start_protected = start_eval.protected_abs_error_pct(:) / 100;
protected_ceiling = min(0.15, max(0.06, start_protected + 0.015));
protected_excess = max(0, protected_rel - protected_ceiling);
J = J + sum(protected_weight(:) .* (protected_excess / 0.020).^2);
protected_degradation_pct = max(0, ...
    max(eval.protected_abs_error_pct(:) - start_eval.protected_abs_error_pct(:)));

rmse_excess = max(0, eval.rmse - start_eval.rmse - 0.003);
J = J + 35.0 * (rmse_excess / 0.010)^2;

% Keep the polish local unless the 15% gate justifies movement.
x = x(:);
x0 = x0(:);
lb = lb(:);
ub = ub(:);
J = J + 0.12 * sum((log(max(x, 1e-12) ./ max(x0, 1e-12))).^2);

span = max(ub - lb, 1e-12);
dist_lower = (x - lb) ./ span;
dist_upper = (ub - x) ./ span;
J = J + 25.0 * sum(max(0, 0.06 - dist_lower).^2 + max(0, 0.06 - dist_upper).^2);

if protected_degradation_pct > 2.0
    J = J + 20.0 * ((protected_degradation_pct - 2.0) / 2.0)^2;
end
end

function eval = evaluate_candidate(params, clinical, scenario, targets, ...
    important_metrics, protected_metrics)
% EVALUATE_CANDIDATE - direct simulation metrics and gate summaries.
eval = struct();
eval.valid = false;
eval.metrics = struct();
eval.rmse = Inf;
eval.error_table = table();
eval.important_abs_error_pct = Inf(numel(important_metrics), 1);
eval.protected_abs_error_pct = Inf(numel(protected_metrics), 1);
eval.max_important_abs_error_pct = Inf;
eval.n_important_over_gate = Inf;
eval.max_measured_abs_error_pct = Inf;
eval.n_measured_over_gate = Inf;
eval.max_protected_degradation_pct = 0;

try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
catch
    return;
end

metric_col = targets.Metric;
clinical_col = targets.Clinical;
model_col = nan(height(targets), 1);
error_pct_col = nan(height(targets), 1);
for idx = 1:height(targets)
    metric_name = metric_col{idx};
    if isfield(metrics, metric_name) && isfinite(metrics.(metric_name))
        model_col(idx) = metrics.(metric_name);
    end
end

valid = isfinite(clinical_col) & isfinite(model_col);
error_pct_col(valid) = 100 * (model_col(valid) - clinical_col(valid)) ./ ...
    max(abs(clinical_col(valid)), 1e-9);
if any(valid)
    eval.rmse = sqrt(mean((error_pct_col(valid) / 100).^2));
end

eval.valid = any(valid);
eval.metrics = metrics;
eval.error_table = table(metric_col, clinical_col, model_col, error_pct_col, ...
    'VariableNames', {'Metric','Clinical','Model','Error_pct'});
eval.important_abs_error_pct = abs_error_for_metrics(eval.error_table, important_metrics);
eval.protected_abs_error_pct = abs_error_for_metrics(eval.error_table, protected_metrics);
eval.max_important_abs_error_pct = max(eval.important_abs_error_pct);
eval.n_important_over_gate = sum(eval.important_abs_error_pct > 15.0);
measured_abs_error_pct = abs(error_pct_col(valid));
eval.max_measured_abs_error_pct = max(measured_abs_error_pct);
eval.n_measured_over_gate = sum(measured_abs_error_pct > 15.0);
end

function tf = is_better_gate_result(candidate_eval, incumbent_eval)
% IS_BETTER_GATE_RESULT - prioritize measured gate failures, then max error.
if candidate_eval.n_measured_over_gate < incumbent_eval.n_measured_over_gate
    tf = true;
    return;
end
if candidate_eval.n_measured_over_gate > incumbent_eval.n_measured_over_gate
    tf = false;
    return;
end
if candidate_eval.max_measured_abs_error_pct < ...
        incumbent_eval.max_measured_abs_error_pct - 0.10
    tf = true;
    return;
end
if candidate_eval.max_measured_abs_error_pct > ...
        incumbent_eval.max_measured_abs_error_pct + 0.10
    tf = false;
    return;
end
tf = candidate_eval.rmse < incumbent_eval.rmse;
end

function tbl = build_target_table(clinical, scenario)
targets = get_calibration_targets(scenario, clinical);
metric_col = {targets.Metric}';
clinical_col = [targets.ClinicalValue]';
unit_col = {targets.Unit}';
tbl = table(metric_col, clinical_col, unit_col, ...
    'VariableNames', {'Metric','Clinical','Unit'});
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

function weights = metric_weight_vector(metric_names, weight_struct)
weights = ones(numel(metric_names), 1);
for idx = 1:numel(metric_names)
    metric_name = metric_names{idx};
    if isfield(weight_struct, metric_name)
        weights(idx) = weight_struct.(metric_name);
    end
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
abs_error_pct(~isfinite(abs_error_pct)) = Inf;
end

function error_pct = signed_error_for_metric(error_table, metric_name)
error_pct = NaN;
hit = strcmp(error_table.Metric, metric_name);
if any(hit)
    error_pct = error_table.Error_pct(find(hit, 1, 'first'));
end
end

function tbl = build_comparison_table(start_eval, polished_eval, targets)
start_tbl = start_eval.error_table;
polished_tbl = polished_eval.error_table;
[is_common, idx_polished] = ismember(start_tbl.Metric, polished_tbl.Metric);
idx_start = find(is_common);
idx_polished = idx_polished(is_common);

metric_col = start_tbl.Metric(idx_start);
unit_col = strings(numel(idx_start), 1);
for idx = 1:numel(idx_start)
    target_hit = strcmp(targets.Metric, metric_col{idx});
    if any(target_hit)
        unit_col(idx) = string(targets.Unit{find(target_hit, 1, 'first')});
    end
end

start_abs = abs(start_tbl.Error_pct(idx_start));
polished_abs = abs(polished_tbl.Error_pct(idx_polished));
tbl = table(metric_col, cellstr(unit_col), start_tbl.Clinical(idx_start), ...
    start_tbl.Model(idx_start), start_tbl.Error_pct(idx_start), ...
    polished_tbl.Model(idx_polished), polished_tbl.Error_pct(idx_polished), ...
    start_abs, polished_abs, polished_abs - start_abs, polished_abs <= 15.0, ...
    'VariableNames', {'Metric','Unit','Clinical','StartModel','StartError_pct', ...
    'PolishedModel','PolishedError_pct','StartAbsError_pct','PolishedAbsError_pct', ...
    'DeltaAbsError_pct','Pass_15pct'});
end

function tbl = build_parameter_table(params_start, params_polished, reference_params, ...
    calib_seed, names, case_profile, used_lb, used_ub)
if isempty(names)
    tbl = table();
    return;
end

registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
if nargin < 7 || isempty(used_lb)
    used_lb = registry.lb(loc);
end
if nargin < 8 || isempty(used_ub)
    used_ub = registry.ub(loc);
end
start_values = values_from_params(params_start, reference_params, names, case_profile);
polished_values = values_from_params(params_polished, reference_params, names, case_profile);
ratio_to_start = polished_values ./ max(abs(start_values), 1e-12);
tbl = table(names(:), registry.baseline_scaled(loc), start_values(:), ...
    polished_values(:), ratio_to_start(:), registry.lb(loc), registry.ub(loc), ...
    used_lb(:), used_ub(:), ...
    'VariableNames', {'Parameter','BaselineScaled','StartValue','PolishedValue', ...
    'RatioToStart','RegistryLowerBound','RegistryUpperBound', ...
    'UsedLowerBound','UsedUpperBound'});
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
params_seeded = params_from_clinical(params_scaled, clinical, scenario, params_scaled, case_profile);
registry_context = struct('params_adult', params_ref, 'params_scaled', params_scaled);
end

function relaxation = make_rv_bound_relaxation()
% MAKE_RV_BOUND_RELAXATION - explicit probe for RV volume-bound contact.
relaxation = struct();
relaxation.names = {'E.RV.EA','V0.RV'};
relaxation.lowerFactor = [1.00, 0.90];
relaxation.upperFactor = [1.25, 1.00];
end

function [lb, ub] = apply_bound_relaxations(lb, ub, names, relaxation)
% APPLY_BOUND_RELAXATIONS - adjust selected experimental bounds.
for idx = 1:numel(relaxation.names)
    hit = strcmp(names, relaxation.names{idx});
    if ~any(hit)
        continue;
    end
    lb(hit) = lb(hit) * relaxation.lowerFactor(idx);
    ub(hit) = ub(hit) * relaxation.upperFactor(idx);
end
end

function empty = make_empty_round_result()
empty = struct('label', '', 'names', {{}}, 'x', [], 'objective', Inf, ...
    'params', struct(), 'eval', struct(), 'lb', [], 'ub', []);
end

function write_manifest(path, result)
fid = fopen(path, 'w');
if fid < 0
    error('run_reyna_15pct_gate_polish:manifestOpenFailed', ...
        'Unable to write manifest: %s', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Reyna 15%% Gate Polish\n');
fprintf(fid, '=====================\n');
fprintf(fid, 'StartPackage: %s\n', result.start_package_path);
fprintf(fid, 'RoundLabel: %s\n', result.round_label);
fprintf(fid, 'StartRMSE: %.6f\n', result.start_eval.rmse);
fprintf(fid, 'PolishedRMSE: %.6f\n', result.polished_eval.rmse);
fprintf(fid, 'StartMaxImportantAbsErrorPct: %.6f\n', ...
    result.start_eval.max_important_abs_error_pct);
fprintf(fid, 'PolishedMaxImportantAbsErrorPct: %.6f\n', ...
    result.polished_eval.max_important_abs_error_pct);
fprintf(fid, 'StartImportantFailuresOver15pct: %d\n', ...
    result.start_eval.n_important_over_gate);
fprintf(fid, 'PolishedImportantFailuresOver15pct: %d\n', ...
    result.polished_eval.n_important_over_gate);
fprintf(fid, 'StartMeasuredFailuresOver15pct: %d\n', ...
    result.start_eval.n_measured_over_gate);
fprintf(fid, 'PolishedMeasuredFailuresOver15pct: %d\n', ...
    result.polished_eval.n_measured_over_gate);
fprintf(fid, 'StartMaxMeasuredAbsErrorPct: %.6f\n', ...
    result.start_eval.max_measured_abs_error_pct);
fprintf(fid, 'PolishedMaxMeasuredAbsErrorPct: %.6f\n', ...
    result.polished_eval.max_measured_abs_error_pct);
fprintf(fid, 'ImportantMetrics: %s\n', strjoin(result.important_metrics, ', '));
fprintf(fid, 'ProtectedMetrics: %s\n', strjoin(result.protected_metrics, ', '));
fprintf(fid, 'ComparisonCSV: %s\n', result.paths.comparison_csv);
fprintf(fid, 'ParameterCSV: %s\n', result.paths.parameter_csv);
fprintf(fid, 'MatFile: %s\n', result.paths.mat);
end

function best_candidate_path = find_latest_reyna_best_candidate(root)
listing = dir(fullfile(root, 'results', 'runs', '*reyna_pre_surgery', ...
    'mat', 'params_best_candidate_pre_surgery.mat'));
if isempty(listing)
    error('run_reyna_15pct_gate_polish:noBestCandidate', ...
        'No Reyna best-candidate package found under results/runs.');
end
[~, order] = sort([listing.datenum], 'descend');
best_candidate_path = fullfile(listing(order(1)).folder, listing(order(1)).name);
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
