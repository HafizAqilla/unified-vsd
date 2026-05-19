function result = run_reyna_targeted_grid_polish(start_package_path, options)
% RUN_REYNA_TARGETED_GRID_POLISH
% -----------------------------------------------------------------------
% Runs a deterministic local grid around the strongest Reyna best-candidate
% residual directions identified by scan_reyna_targeted_sensitivity.m.
%
% INPUTS:
%   start_package_path - optional params_best_candidate MAT path          [-]
%   options            - optional struct with factor grids                [-]
%
% OUTPUTS:
%   result             - struct with best candidate and exported tables   [-]
%
% ASSUMPTIONS:
%   - C.SAR is allowed to move inside registry bounds because catheter pulse
%     pressure and Fick flow constrain systemic arterial compliance.
%   - RVEDV remains a consistency-check target and is reported, not fitted.
%
% REFERENCES:
%   [1] Stergiopulos et al. (1994). Arterial compliance from pressure-flow.
%   [2] docs/clinical_data_dictionary.md.
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
start_eval = evaluate_candidate(params_start, clinical, scenario, targets);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_reyna_targeted_grid_polish', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end

names = {'C.SAR','E.LV.EB','V0.LV','E.LV.EA'};
base_values = values_from_params(params_start, params_scaled, names, case_profile);
[lb, ub] = bounds_for_names(calib_seed, names);
lb = max(1e-9, lb(:) .* options.bound_lower_factors(:));
ub = max(lb * 1.001, ub(:) .* options.bound_upper_factors(:));

fprintf('[run_reyna_targeted_grid_polish] Start primary RMSE %.4f, full RMSE %.4f\n', ...
    start_eval.rmse_primary, start_eval.rmse_full);
fprintf('[run_reyna_targeted_grid_polish] Grid sizes: C.SAR=%d, E.LV.EB=%d, V0.LV=%d, E.LV.EA=%d\n', ...
    numel(options.C_SAR_factors), numel(options.E_LV_EB_factors), ...
    numel(options.V0_LV_factors), numel(options.E_LV_EA_factors));

rows = {};
best_eval = start_eval;
params_best = params_start;
best_values = base_values;
best_score = candidate_score(start_eval);
eval_count = 0;

for i_c = 1:numel(options.C_SAR_factors)
    for i_eb = 1:numel(options.E_LV_EB_factors)
        for i_v0 = 1:numel(options.V0_LV_factors)
            for i_ea = 1:numel(options.E_LV_EA_factors)
                factors = [options.C_SAR_factors(i_c), ...
                    options.E_LV_EB_factors(i_eb), ...
                    options.V0_LV_factors(i_v0), ...
                    options.E_LV_EA_factors(i_ea)];
                trial_values = min(max(base_values(:) .* factors(:), lb(:)), ub(:));
                params_trial = apply_values(params_start, params_scaled, ...
                    names, trial_values, case_profile);
                trial_eval = evaluate_candidate(params_trial, clinical, scenario, targets);
                eval_count = eval_count + 1;
                rows(end + 1, :) = trial_row(factors, trial_values, trial_eval); %#ok<AGROW>

                trial_score = candidate_score(trial_eval);
                if trial_score < best_score
                    best_score = trial_score;
                    best_eval = trial_eval;
                    params_best = params_trial;
                    best_values = trial_values;
                    fprintf('  New best %03d | primary %.4f | full %.4f | max %.2f%%\n', ...
                        eval_count, best_eval.rmse_primary, best_eval.rmse_full, ...
                        best_eval.max_primary_abs_error_pct);
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

grid_csv = fullfile(tables_dir, 'reyna_targeted_grid_all_trials.csv');
comparison_csv = fullfile(tables_dir, 'reyna_targeted_grid_metrics.csv');
parameter_csv = fullfile(tables_dir, 'reyna_targeted_grid_parameters.csv');
primary_csv = fullfile(tables_dir, 'reyna_targeted_grid_primary_selection.csv');
manifest_path = fullfile(run_dir, 'run_manifest.txt');
mat_path = fullfile(mat_dir, 'reyna_targeted_grid_result.mat');

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

fprintf('[run_reyna_targeted_grid_polish] Done. Primary RMSE %.4f -> %.4f\n', ...
    start_eval.rmse_primary, best_eval.rmse_primary);
fprintf('[run_reyna_targeted_grid_polish] Exported:\n  %s\n  %s\n', ...
    comparison_csv, parameter_csv);
disp(comparison_table);
end

function options = default_options(options)
defaults = struct();
defaults.C_SAR_factors = [0.60, 0.65, 0.70, 0.75, 0.80];
defaults.E_LV_EB_factors = [1.00, 1.08, 1.15, 1.22];
defaults.V0_LV_factors = [1.00, 1.05, 1.10];
defaults.E_LV_EA_factors = [0.94, 1.00];
defaults.bound_lower_factors = [1.00, 1.00, 1.00, 1.00];
defaults.bound_upper_factors = [1.00, 1.00, 1.00, 1.00];
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
    '20260514_155253_reyna_pre_surgery', 'mat', ...
    'params_best_candidate_pre_surgery.mat');
end

function score = candidate_score(eval)
if ~eval.valid || ~eval.physiology_valid
    score = Inf;
    return;
end
gate_excess = max(0, eval.max_primary_abs_error_pct - 10.0) / 100; % [-]
score = eval.rmse_primary + 0.12 * eval.rmse_full + ...
    0.010 * eval.n_primary_over_gate + 0.020 * gate_excess;
end

function names = grid_variable_names()
names = {'C_SAR_factor','E_LV_EB_factor','V0_LV_factor','E_LV_EA_factor', ...
    'C_SAR_value','E_LV_EB_value','V0_LV_value','E_LV_EA_value', ...
    'PrimaryRMSE','FullRMSE','HardRMSE','SoftRMSE','MaxPrimaryAbsErrorPct', ...
    'PrimaryFailuresOver10pct','PhysiologyValid','LAP_mean_ErrorPct', ...
    'LVEF_ErrorPct','LVESV_ErrorPct','LVEDV_ErrorPct','CO_Lmin_ErrorPct', ...
    'SAP_min_ErrorPct','SAP_mean_ErrorPct','SAP_max_ErrorPct','QpQs_ErrorPct','PAP_mean_ErrorPct', ...
    'RAP_mean_ErrorPct','SVR_ErrorPct'};
end

function row = trial_row(factors, values, eval)
metric_names = {'LAP_mean','LVEF','LVESV','LVEDV','CO_Lmin', ...
    'SAP_min','SAP_mean','SAP_max','QpQs','PAP_mean','RAP_mean','SVR'};
row = num2cell([factors(:)', values(:)', eval.rmse_primary, ...
    eval.rmse_full, eval.rmse_hard, eval.rmse_soft, ...
    eval.max_primary_abs_error_pct, eval.n_primary_over_gate, ...
    double(eval.physiology_valid)]);
for idx = 1:numel(metric_names)
    row{end + 1} = signed_error_for_metric(eval.error_table, metric_names{idx}); %#ok<AGROW>
end
end

function eval = evaluate_candidate(params, clinical, scenario, targets)
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

primary_mask = valid & targets.IncludedInPrimaryRMSE;
hard_mask = valid & strcmp(targets.Tier, 'hard');
soft_mask = valid & strcmp(targets.Tier, 'soft');

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
eval.max_primary_abs_error_pct = max(abs(error_pct_col(primary_mask)));
eval.n_primary_over_gate = sum(abs(error_pct_col(primary_mask)) > 10.0);
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
    best_abs - start_abs, best_abs <= 10.0, ...
    start_tbl.Tier(idx_start), start_tbl.IncludedInCalibration(idx_start), ...
    start_tbl.IncludedInPrimaryRMSE(idx_start), ...
    'VariableNames', {'Metric','Unit','Clinical','StartModel', ...
    'StartError_pct','BestModel','BestError_pct','StartAbsError_pct', ...
    'BestAbsError_pct','DeltaAbsError_pct','Pass_10pct','Tier', ...
    'IncludedInCalibration','IncludedInPrimaryRMSE'});
end

function [lb, ub] = bounds_for_names(calib_seed, names)
registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
if any(loc == 0)
    missing = names(loc == 0);
    error('run_reyna_targeted_grid_polish:missingRegistryName', ...
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
    error('run_reyna_targeted_grid_polish:manifestOpenFailed', ...
        'Unable to write manifest: %s', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Reyna Targeted Grid Polish\n');
fprintf(fid, '==========================\n');
fprintf(fid, 'StartPackage: %s\n', result.start_package_path);
fprintf(fid, 'StartPrimaryRMSE: %.6f\n', result.start_eval.rmse_primary);
fprintf(fid, 'BestPrimaryRMSE: %.6f\n', result.best_eval.rmse_primary);
fprintf(fid, 'StartFullRMSE: %.6f\n', result.start_eval.rmse_full);
fprintf(fid, 'BestFullRMSE: %.6f\n', result.best_eval.rmse_full);
fprintf(fid, 'StartHardRMSE: %.6f\n', result.start_eval.rmse_hard);
fprintf(fid, 'BestHardRMSE: %.6f\n', result.best_eval.rmse_hard);
fprintf(fid, 'StartSoftRMSE: %.6f\n', result.start_eval.rmse_soft);
fprintf(fid, 'BestSoftRMSE: %.6f\n', result.best_eval.rmse_soft);
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
