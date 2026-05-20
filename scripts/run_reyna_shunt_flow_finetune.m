function result = run_reyna_shunt_flow_finetune(start_package_path, options)
% RUN_REYNA_SHUNT_FLOW_FINETUNE
% -----------------------------------------------------------------------
% Locally polishes Reyna pre-surgery shunt-flow consistency while preserving
% the low-RMSE pressure-flow fit.
%
% INPUTS:
%   start_package_path - MAT file containing best_candidate.params        [-]
%   options            - struct with local grid factors                   [-]
%
% OUTPUTS:
%   result             - struct with start/best metrics and export paths  [-]
%
% ASSUMPTIONS:
%   - Q_shunt_Lmin is a derived consistency target from Qp - Qs, so it is
%     weighted behind the direct catheter anchors RAP/PAP/SAP/QpQs/Qs.
%   - VSD geometry is refreshed from the active clinical record before
%     polishing so the candidate matches the current data dictionary.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md, Reyna CO comparator note.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-20
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
params_start = align_start_candidate(loaded.best_candidate.params, clinical);
start_eval = evaluate_candidate(params_start, clinical, scenario, targets, options);

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_reyna_shunt_flow_finetune', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end

names = {'group.R_sys_scale','R.SVEN','vsd.Cd','group.R_pul_scale'};
base_values = values_from_params(params_start, params_scaled, names, case_profile);
[lb, ub] = bounds_for_names(calib_seed, names);

fprintf('[run_reyna_shunt_flow_finetune] Start direct RMSE %.4f | Qshunt error %.2f%%\n', ...
    start_eval.rmse_direct, signed_error_for_metric(start_eval.error_table, 'Q_shunt_Lmin'));

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
                    fprintf('  New best %03d | direct %.4f | Qshunt %.2f%% | QpQs %.2f%% | CO %.2f%%\n', ...
                        eval_count, best_eval.rmse_direct, ...
                        signed_error_for_metric(best_eval.error_table, 'Q_shunt_Lmin'), ...
                        signed_error_for_metric(best_eval.error_table, 'QpQs'), ...
                        signed_error_for_metric(best_eval.error_table, 'CO_Lmin'));
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

grid_csv = fullfile(tables_dir, 'reyna_shunt_flow_finetune_all_trials.csv');
comparison_csv = fullfile(tables_dir, 'reyna_shunt_flow_finetune_metrics.csv');
parameter_csv = fullfile(tables_dir, 'reyna_shunt_flow_finetune_parameters.csv');
primary_csv = fullfile(tables_dir, 'reyna_shunt_flow_primary_selection.csv');
manifest_path = fullfile(run_dir, 'run_manifest.txt');
mat_path = fullfile(mat_dir, 'reyna_shunt_flow_finetune_result.mat');

writetable(grid_table, grid_csv);
writetable(comparison_table, comparison_csv);
writetable(parameter_table, parameter_csv);
writetable(primary_selection_table, primary_csv);

best_candidate = struct();
best_candidate.label = 'shunt_flow_finetune';
best_candidate.params = params_best;
best_candidate.metrics = best_eval.metrics;
best_candidate.validity = best_eval.validity;
best_candidate.report_table = comparison_table;
best_candidate.primary_rmse = best_eval.rmse_primary;
best_candidate.direct_rmse = best_eval.rmse_direct;
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

fprintf('[run_reyna_shunt_flow_finetune] Done. Direct RMSE %.4f -> %.4f | Qshunt %.2f%% -> %.2f%%\n', ...
    start_eval.rmse_direct, best_eval.rmse_direct, ...
    signed_error_for_metric(start_eval.error_table, 'Q_shunt_Lmin'), ...
    signed_error_for_metric(best_eval.error_table, 'Q_shunt_Lmin'));
disp(comparison_table);
end

function options = default_options(options)
defaults = struct();
defaults.R_sys_factors = [0.90, 0.96, 1.00, 1.04];
defaults.R_SVEN_factors = [0.90, 1.00, 1.10];
defaults.vsd_Cd_factors = [1.00, 1.10, 1.20, 1.30];
defaults.R_pul_factors = [0.92, 1.00, 1.08];
defaults.direct_metrics = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
defaults.primary_gate_pct = 12.0;
defaults.shunt_weight = 0.55;
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
preferred_path = fullfile(root, 'results', 'runs', ...
    '20260514_155253_reyna_pre_surgery', 'mat', ...
    'params_best_candidate_pre_surgery.mat');
if exist(preferred_path, 'file')
    path = preferred_path;
    return;
end
path = latest_existing_result_file(root, '*_reyna_pre_surgery', ...
    fullfile('mat', 'params_best_candidate_pre_surgery.mat'), ...
    'Reyna pre-surgery best-candidate package');
end

function score = candidate_score(eval, options)
if ~eval.valid || ~eval.physiology_valid
    score = Inf;
    return;
end
qshunt_err_rel = abs(signed_error_for_metric(eval.error_table, 'Q_shunt_Lmin')) / 100;
direct_gate_excess = max(0, eval.max_direct_abs_error_pct - options.primary_gate_pct) / 100;
score = eval.rmse_direct + options.shunt_weight * qshunt_err_rel + ...
    0.35 * direct_gate_excess^2 + 0.02 * eval.n_direct_over_gate;
end

function names = grid_variable_names()
names = {'R_sys_factor','R_SVEN_factor','vsd_Cd_factor','R_pul_factor', ...
    'R_sys_value','R_SVEN_value','vsd_Cd_value','R_pul_value', ...
    'DirectRMSE','PrimaryRMSE','FullRMSE','HardRMSE','SoftRMSE', ...
    'MaxDirectAbsErrorPct','DirectFailuresOverGate','PhysiologyValid', ...
    'RAP_mean_ErrorPct','PAP_mean_ErrorPct','SAP_mean_ErrorPct', ...
    'QpQs_ErrorPct','CO_Lmin_ErrorPct','Q_shunt_Lmin_ErrorPct', ...
    'PAP_min_ErrorPct','PAP_max_ErrorPct','SAP_min_ErrorPct','SAP_max_ErrorPct'};
end

function row = trial_row(factors, values, eval)
metric_names = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin', ...
    'Q_shunt_Lmin','PAP_min','PAP_max','SAP_min','SAP_max'};
row = num2cell([factors(:)', values(:)', eval.rmse_direct, ...
    eval.rmse_primary, eval.rmse_full, eval.rmse_hard, eval.rmse_soft, ...
    eval.max_direct_abs_error_pct, eval.n_direct_over_gate, ...
    double(eval.physiology_valid)]);
for idx = 1:numel(metric_names)
    row{end + 1} = signed_error_for_metric(eval.error_table, metric_names{idx}); %#ok<AGROW>
end
end

function eval = evaluate_candidate(params, clinical, scenario, targets, options)
eval = empty_eval();
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
direct_mask = valid & ismember(targets.Metric, options.direct_metrics(:));
direct_abs = abs(error_pct_col(direct_mask));

eval.valid = any(valid);
eval.physiology_valid = validity.is_valid;
eval.metrics = metrics;
eval.validity = validity;
eval.error_table = table(targets.Metric, targets.Unit, clinical_col, ...
    model_col, error_pct_col, targets.Tier, targets.IncludedInCalibration, ...
    targets.IncludedInPrimaryRMSE, ...
    'VariableNames', {'Metric','Unit','Clinical','Model','Error_pct', ...
    'Tier','IncludedInCalibration','IncludedInPrimaryRMSE'});
eval.rmse_direct = rmse_for_mask(error_pct_col, direct_mask);
eval.rmse_primary = rmse_for_mask(error_pct_col, primary_mask);
eval.rmse_full = rmse_for_mask(error_pct_col, valid);
eval.rmse_hard = rmse_for_mask(error_pct_col, hard_mask);
eval.rmse_soft = rmse_for_mask(error_pct_col, soft_mask);
eval.max_direct_abs_error_pct = max(direct_abs);
eval.n_direct_over_gate = sum(direct_abs > options.primary_gate_pct);
end

function eval = empty_eval()
eval = struct();
eval.valid = false;
eval.physiology_valid = false;
eval.metrics = struct();
eval.validity = struct('is_valid', false, 'failed_flags', {{}}); 
eval.rmse_direct = Inf;
eval.rmse_primary = Inf;
eval.rmse_full = Inf;
eval.rmse_hard = Inf;
eval.rmse_soft = Inf;
eval.error_table = table();
eval.max_direct_abs_error_pct = Inf;
eval.n_direct_over_gate = Inf;
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
    error('run_reyna_shunt_flow_finetune:missingRegistryName', ...
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
if isempty(error_table)
    return;
end
hit = strcmp(error_table.Metric, metric_name);
if any(hit)
    error_pct = error_table.Error_pct(find(hit, 1, 'first'));
end
end

function params = align_start_candidate(params, clinical)
params.HR = clinical.common.HR;
src = clinical.pre_surgery;
if isfield(src, 'VSD_mode') && ~isempty(src.VSD_mode)
    params.vsd.mode = char(src.VSD_mode);
end
if isfield(src, 'VSD_diameter_mm') && isfinite(src.VSD_diameter_mm)
    params.vsd.diameter_mm = src.VSD_diameter_mm;
    params.vsd.area_mm2 = pi * (src.VSD_diameter_mm / 2)^2;
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
    error('run_reyna_shunt_flow_finetune:manifestOpenFailed', ...
        'Unable to write manifest: %s', path);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, 'Reyna Shunt-Flow Finetune\n');
fprintf(fid, '=========================\n');
fprintf(fid, 'StartPackage: %s\n', result.start_package_path);
fprintf(fid, 'StartDirectRMSE: %.6f\n', result.start_eval.rmse_direct);
fprintf(fid, 'BestDirectRMSE: %.6f\n', result.best_eval.rmse_direct);
fprintf(fid, 'StartPrimaryRMSE: %.6f\n', result.start_eval.rmse_primary);
fprintf(fid, 'BestPrimaryRMSE: %.6f\n', result.best_eval.rmse_primary);
fprintf(fid, 'StartFullRMSE: %.6f\n', result.start_eval.rmse_full);
fprintf(fid, 'BestFullRMSE: %.6f\n', result.best_eval.rmse_full);
fprintf(fid, 'StartQshuntErrorPct: %.6f\n', ...
    signed_error_for_metric(result.start_eval.error_table, 'Q_shunt_Lmin'));
fprintf(fid, 'BestQshuntErrorPct: %.6f\n', ...
    signed_error_for_metric(result.best_eval.error_table, 'Q_shunt_Lmin'));
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

function path = latest_existing_result_file(root, run_pattern, relative_file, label)
% LATEST_EXISTING_RESULT_FILE - newest local run artifact matching pattern.
runs_dir = fullfile(root, 'results', 'runs');
matches = dir(fullfile(runs_dir, run_pattern));
path = '';
best_datenum = -Inf;                           % [datenum]

for idx = 1:numel(matches)
    if ~matches(idx).isdir
        continue;
    end
    candidate_path = fullfile(matches(idx).folder, matches(idx).name, relative_file);
    info = dir(candidate_path);
    if isempty(info)
        continue;
    end
    if info.datenum > best_datenum
        best_datenum = info.datenum;           % [datenum]
        path = candidate_path;                 % [char]
    end
end

if isempty(path)
    error('run_reyna_shunt_flow_finetune:missingDefaultStartPath', ...
        'No %s found under %s. Pass start_package_path explicitly.', ...
        label, runs_dir);
end
end
