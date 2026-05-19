function result = run_reyna_pressure_flow_grid_polish(start_package_path, options)
% RUN_REYNA_PRESSURE_FLOW_GRID_POLISH
% -----------------------------------------------------------------------
% Runs a deterministic local grid around a Reyna candidate to improve the
% already-plausible catheter pressure-flow targets without hiding the
% known echo-volume conflict.
%
% INPUTS:
%   start_package_path - optional MAT path containing best_candidate      [-]
%   options            - optional struct with factor grids                [-]
%
% OUTPUTS:
%   result             - struct with best candidate and exported tables   [-]
%
% ASSUMPTIONS:
%   - The pressure-flow block is judged separately from echo volume targets
%     because the pre-calibration audit found a critical stroke-volume
%     inconsistency between catheter-derived flow and echo-derived volumes.
%   - RVEDV remains consistency-check-only and is reported transparently.
%
% REFERENCES:
%   [1] docs/reyna_targeted_grid_polish_20260515.md.
%   [2] docs/calibration_data_governance_notes.md.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-15
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1 || isempty(start_package_path)
    start_package_path = default_start_path();
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
[params_scaled, params_seeded, registry_context] = ...
    build_reyna_reference_context(clinical, scenario, case_profile);
[primary_metrics, primary_selection_table] = ...
    select_primary_metrics(clinical, [], scenario, case_profile);
calib_seed = calibration_param_sets( ...
    scenario, params_seeded, [], primary_metrics, case_profile, registry_context);
targets = build_target_table(clinical, scenario, case_profile);

loaded = load(start_package_path, 'best_candidate');
params_start = loaded.best_candidate.params;
start_eval = evaluate_candidate(params_start, clinical, scenario, targets, options);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_reyna_pressure_flow_grid_polish', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end

names = {'group.R_sys_scale','R.SVEN','vsd.Cd','group.R_pul_scale'};
base_values = values_from_params(params_start, params_scaled, names, case_profile);
[lb, ub] = bounds_for_names(calib_seed, names);

fprintf('[run_reyna_pressure_flow_grid_polish] Start pressure-flow RMSE %.4f, primary RMSE %.4f\n', ...
    start_eval.rmse_pressure_flow, start_eval.rmse_primary);
fprintf('[run_reyna_pressure_flow_grid_polish] Grid sizes: Rsys=%d, Rven=%d, Cd=%d, Rpul=%d\n', ...
    numel(options.R_sys_factors), numel(options.R_SVEN_factors), ...
    numel(options.vsd_Cd_factors), numel(options.R_pul_factors));

rows = {};
best_eval = start_eval;
params_best = params_start;
best_values = base_values;
best_score = candidate_score(start_eval, options);
eval_count = 0;

for i_sys = 1:numel(options.R_sys_factors)
    for i_sven = 1:numel(options.R_SVEN_factors)
        for i_cd = 1:numel(options.vsd_Cd_factors)
            for i_pul = 1:numel(options.R_pul_factors)
                factors = [options.R_sys_factors(i_sys), ...
                    options.R_SVEN_factors(i_sven), ...
                    options.vsd_Cd_factors(i_cd), ...
                    options.R_pul_factors(i_pul)];
                trial_values = min(max(base_values(:) .* factors(:), lb(:)), ub(:));
                params_trial = apply_values(params_start, params_scaled, ...
                    names, trial_values, case_profile);
                trial_eval = evaluate_candidate(params_trial, clinical, scenario, targets, options);
                eval_count = eval_count + 1;
                rows(end + 1, :) = trial_row(factors, trial_values, trial_eval); %#ok<AGROW>

                trial_score = candidate_score(trial_eval, options);
                if trial_score < best_score
                    best_score = trial_score;
                    best_eval = trial_eval;
                    params_best = params_trial;
                    best_values = trial_values;
                    fprintf('  New best %03d | pressure-flow %.4f | primary %.4f | over5 %d | over10 %d\n', ...
                        eval_count, best_eval.rmse_pressure_flow, best_eval.rmse_primary, ...
                        best_eval.n_pressure_flow_over5, best_eval.n_primary_over_gate);
                end
            end
        end
    end
end

grid_table = cell2table(rows, 'VariableNames', grid_variable_names());
comparison_table = build_comparison_table(start_eval, best_eval);
parameter_table = table(names(:), base_values(:), best_values(:), ...
    best_values(:) ./ max(abs(base_values(:)), 1e-12), lb(:), ub(:), ...
    'VariableNames', {'Parameter','StartValue','BestValue', ...
    'RatioToStart','LowerBound','UpperBound'});

grid_csv = fullfile(tables_dir, 'reyna_pressure_flow_grid_all_trials.csv');
comparison_csv = fullfile(tables_dir, 'reyna_pressure_flow_grid_metrics.csv');
parameter_csv = fullfile(tables_dir, 'reyna_pressure_flow_grid_parameters.csv');
primary_csv = fullfile(tables_dir, 'reyna_pressure_flow_primary_selection.csv');
manifest_path = fullfile(run_dir, 'run_manifest.txt');
mat_path = fullfile(mat_dir, 'reyna_pressure_flow_grid_result.mat');

writetable(grid_table, grid_csv);
writetable(comparison_table, comparison_csv);
writetable(parameter_table, parameter_csv);
writetable(primary_selection_table, primary_csv);

best_candidate = struct();
best_candidate.params = params_best;
best_candidate.metrics = best_eval.metrics;
best_candidate.validity = best_eval.validity;
best_candidate.report_table = comparison_table;
best_candidate.primary_rmse = best_eval.rmse_primary;
best_candidate.full_rmse = best_eval.rmse_full;
best_candidate.pressure_flow_rmse = best_eval.rmse_pressure_flow;
best_candidate.parameter_table = parameter_table;

result = struct();
result.run_dir = run_dir;
result.start_package_path = start_package_path;
result.params_start = params_start;
result.params_best = params_best;
result.start_eval = start_eval;
result.best_eval = best_eval;
result.grid_table = grid_table;
result.comparison_table = comparison_table;
result.parameter_table = parameter_table;
result.primary_metrics = primary_metrics;
result.primary_selection_table = primary_selection_table;
result.best_candidate = best_candidate;
result.options = options;
result.paths = struct('grid_csv', grid_csv, 'comparison_csv', comparison_csv, ...
    'parameter_csv', parameter_csv, 'primary_csv', primary_csv, ...
    'manifest', manifest_path, 'mat', mat_path);

save(mat_path, 'result', 'best_candidate');
write_manifest(manifest_path, result);

fprintf('[run_reyna_pressure_flow_grid_polish] Done. Pressure-flow RMSE %.4f -> %.4f\n', ...
    start_eval.rmse_pressure_flow, best_eval.rmse_pressure_flow);
fprintf('[run_reyna_pressure_flow_grid_polish] Exported:\n  %s\n  %s\n', ...
    comparison_csv, parameter_csv);
disp(comparison_table);
end

function options = default_options(options)
defaults = struct();
defaults.R_sys_factors = [0.78, 0.88, 1.00, 1.08];
defaults.R_SVEN_factors = [0.85, 1.00, 1.15];
defaults.vsd_Cd_factors = [1.00, 1.08, 1.16, 1.20];
defaults.R_pul_factors = [0.90, 1.00, 1.12];
defaults.pressure_flow_metrics = {'RAP_mean','PAP_mean','SAP_min', ...
    'SAP_max','SAP_mean','QpQs','SVR','CO_Lmin'};
defaults.primary_gate_pct = 10.0;
defaults.pressure_flow_goal_pct = 5.0;
defaults.max_primary_over10 = 3;
fields = fieldnames(defaults);
for idx = 1:numel(fields)
    field_name = fields{idx};
    if ~isfield(options, field_name) || isempty(options.(field_name))
        options.(field_name) = defaults.(field_name);
    end
end
end

function path = default_start_path()
root = fileparts(fileparts(mfilename('fullpath')));
path = fullfile(root, 'results', 'runs', ...
    '20260515_173742_reyna_targeted_grid_polish', 'mat', ...
    'reyna_targeted_grid_result.mat');
end

function score = candidate_score(eval, options)
if ~eval.valid || ~eval.physiology_valid
    score = Inf;
    return;
end
gate_excess = max(0, eval.max_primary_abs_error_pct - options.primary_gate_pct) / 100; % [-]
if eval.n_primary_over_gate > options.max_primary_over10
    primary_gate_penalty = 0.10 * (eval.n_primary_over_gate - options.max_primary_over10);
else
    primary_gate_penalty = 0;
end
score = eval.rmse_pressure_flow + 0.03 * eval.n_pressure_flow_over5 + ...
    0.20 * eval.pressure_flow_excess5_rmse + ...
    0.08 * eval.rmse_primary + 0.01 * eval.n_primary_over_gate + ...
    0.03 * gate_excess + primary_gate_penalty;
end

function names = grid_variable_names()
names = {'R_sys_factor','R_SVEN_factor','vsd_Cd_factor','R_pul_factor', ...
    'R_sys_value','R_SVEN_value','vsd_Cd_value','R_pul_value', ...
    'PressureFlowRMSE','PressureFlowOver5pct','PressureFlowExcess5RMSE', ...
    'PrimaryRMSE','FullRMSE','HardRMSE','SoftRMSE', ...
    'MaxPrimaryAbsErrorPct','PrimaryFailuresOver10pct','PhysiologyValid', ...
    'RAP_mean_ErrorPct','PAP_mean_ErrorPct','SAP_min_ErrorPct', ...
    'SAP_max_ErrorPct','SAP_mean_ErrorPct','QpQs_ErrorPct', ...
    'SVR_ErrorPct','CO_Lmin_ErrorPct','LAP_mean_ErrorPct', ...
    'LVEDV_ErrorPct','LVESV_ErrorPct','RVEDV_ErrorPct','RVESV_ErrorPct', ...
    'LVEF_ErrorPct'};
end

function row = trial_row(factors, values, eval)
metric_names = {'RAP_mean','PAP_mean','SAP_min','SAP_max','SAP_mean', ...
    'QpQs','SVR','CO_Lmin','LAP_mean','LVEDV','LVESV', ...
    'RVEDV','RVESV','LVEF'};
row = num2cell([factors(:)', values(:)', eval.rmse_pressure_flow, ...
    eval.n_pressure_flow_over5, eval.pressure_flow_excess5_rmse, ...
    eval.rmse_primary, eval.rmse_full, eval.rmse_hard, eval.rmse_soft, ...
    eval.max_primary_abs_error_pct, eval.n_primary_over_gate, ...
    double(eval.physiology_valid)]);
for idx = 1:numel(metric_names)
    row{end + 1} = signed_error_for_metric(eval.error_table, metric_names{idx}); %#ok<AGROW>
end
end

function eval = evaluate_candidate(params, clinical, scenario, targets, options)
eval = struct();
eval.valid = false;
eval.physiology_valid = false;
eval.metrics = struct();
eval.validity = struct('is_valid', false, 'failed_flags', {{}});
eval.rmse_primary = Inf;
eval.rmse_full = Inf;
eval.rmse_hard = Inf;
eval.rmse_soft = Inf;
eval.rmse_pressure_flow = Inf;
eval.pressure_flow_excess5_rmse = Inf;
eval.error_table = table();
eval.max_primary_abs_error_pct = Inf;
eval.n_primary_over_gate = Inf;
eval.n_pressure_flow_over5 = Inf;

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

primary_mask = valid & targets.IncludedInPrimaryRMSE;
hard_mask = valid & strcmp(targets.Tier, 'hard');
soft_mask = valid & strcmp(targets.Tier, 'soft');
pressure_flow_mask = valid & ismember(targets.Metric, options.pressure_flow_metrics(:));
pressure_flow_abs = abs(error_pct_col(pressure_flow_mask));
pressure_flow_excess = max(0, pressure_flow_abs - options.pressure_flow_goal_pct) / 100;

eval.valid = any(valid);
eval.physiology_valid = validity.is_valid;
eval.metrics = metrics;
eval.validity = validity;
eval.error_table = table(targets.Metric, targets.Unit, clinical_col, ...
    model_col, error_pct_col, targets.Tier, targets.IncludedInCalibration, ...
    targets.IncludedInPrimaryRMSE, ...
    'VariableNames', {'Metric','Unit','Clinical','Model','Error_pct', ...
    'Tier','IncludedInCalibration','IncludedInPrimaryRMSE'});
eval.rmse_primary = rmse_for_mask(error_pct_col, primary_mask);
eval.rmse_full = rmse_for_mask(error_pct_col, valid);
eval.rmse_hard = rmse_for_mask(error_pct_col, hard_mask);
eval.rmse_soft = rmse_for_mask(error_pct_col, soft_mask);
eval.rmse_pressure_flow = rmse_for_mask(error_pct_col, pressure_flow_mask);
eval.pressure_flow_excess5_rmse = sqrt(mean(pressure_flow_excess.^2));
eval.max_primary_abs_error_pct = max(abs(error_pct_col(primary_mask)));
eval.n_primary_over_gate = sum(abs(error_pct_col(primary_mask)) > options.primary_gate_pct);
eval.n_pressure_flow_over5 = sum(pressure_flow_abs > options.pressure_flow_goal_pct);
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
tier_tbl = case_profile.targetTiers.table;
n_targets = numel(target_struct);
metric_col = {target_struct.Metric}';
clinical_col = [target_struct.ClinicalValue]';
unit_col = {target_struct.Unit}';
tier_col = repmat({'unavailable'}, n_targets, 1);
included_cal = false(n_targets, 1);
included_primary = false(n_targets, 1);
for idx = 1:n_targets
    hit = find(strcmp(tier_tbl.Metric, metric_col{idx}), 1, 'first');
    if ~isempty(hit)
        tier_col{idx} = tier_tbl.Tier{hit};
        included_cal(idx) = logical(tier_tbl.IncludedInCalibration(hit));
        included_primary(idx) = logical(tier_tbl.IncludedInPrimaryRMSE(hit));
    end
end
targets = table(metric_col, clinical_col, unit_col, tier_col, ...
    included_cal, included_primary, ...
    'VariableNames', {'Metric','Clinical','Unit','Tier', ...
    'IncludedInCalibration','IncludedInPrimaryRMSE'});
end

function tbl = build_comparison_table(start_eval, best_eval)
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
    best_abs - start_abs, best_abs <= 5.0, best_abs <= 10.0, ...
    start_tbl.Tier(idx_start), start_tbl.IncludedInCalibration(idx_start), ...
    start_tbl.IncludedInPrimaryRMSE(idx_start), ...
    'VariableNames', {'Metric','Unit','Clinical','StartModel', ...
    'StartError_pct','BestModel','BestError_pct','StartAbsError_pct', ...
    'BestAbsError_pct','DeltaAbsError_pct','Pass_5pct','Pass_10pct', ...
    'Tier','IncludedInCalibration','IncludedInPrimaryRMSE'});
end

function [lb, ub] = bounds_for_names(calib_seed, names)
registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
if any(loc == 0)
    missing = names(loc == 0);
    error('run_reyna_pressure_flow_grid_polish:missingRegistryName', ...
        'Missing registry names: %s', strjoin(missing, ', '));
end
lb = registry.lb(loc);
ub = registry.ub(loc);
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

function error_pct = signed_error_for_metric(error_table, metric_name)
error_pct = NaN;
hit = strcmp(error_table.Metric, metric_name);
if any(hit)
    error_pct = error_table.Error_pct(find(hit, 1, 'first'));
end
end

function [params_scaled, params_seeded, registry_context] = ...
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
params_seeded = params_from_clinical(params_scaled, clinical, scenario, ...
    params_scaled, case_profile);
registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
end

function write_manifest(path, result)
fid = fopen(path, 'w');
if fid < 0
    error('run_reyna_pressure_flow_grid_polish:manifestOpenFailed', ...
        'Unable to write manifest: %s', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Reyna Pressure-Flow Grid Polish\n');
fprintf(fid, '===============================\n');
fprintf(fid, 'StartPackage: %s\n', result.start_package_path);
fprintf(fid, 'StartPressureFlowRMSE: %.6f\n', result.start_eval.rmse_pressure_flow);
fprintf(fid, 'BestPressureFlowRMSE: %.6f\n', result.best_eval.rmse_pressure_flow);
fprintf(fid, 'StartPressureFlowOver5pct: %d\n', result.start_eval.n_pressure_flow_over5);
fprintf(fid, 'BestPressureFlowOver5pct: %d\n', result.best_eval.n_pressure_flow_over5);
fprintf(fid, 'StartPrimaryRMSE: %.6f\n', result.start_eval.rmse_primary);
fprintf(fid, 'BestPrimaryRMSE: %.6f\n', result.best_eval.rmse_primary);
fprintf(fid, 'StartFullRMSE: %.6f\n', result.start_eval.rmse_full);
fprintf(fid, 'BestFullRMSE: %.6f\n', result.best_eval.rmse_full);
fprintf(fid, 'BestMaxPrimaryAbsErrorPct: %.6f\n', ...
    result.best_eval.max_primary_abs_error_pct);
fprintf(fid, 'BestPrimaryFailuresOver10pct: %d\n', ...
    result.best_eval.n_primary_over_gate);
fprintf(fid, 'PhysiologyValid: %d\n', result.best_eval.physiology_valid);
fprintf(fid, 'ComparisonCSV: %s\n', result.paths.comparison_csv);
fprintf(fid, 'ParameterCSV: %s\n', result.paths.parameter_csv);
fprintf(fid, 'GridCSV: %s\n', result.paths.grid_csv);
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
