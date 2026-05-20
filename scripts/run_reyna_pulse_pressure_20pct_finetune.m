function result = run_reyna_pulse_pressure_20pct_finetune(start_result_path, options)
% RUN_REYNA_PULSE_PRESSURE_20PCT_FINETUNE
% -----------------------------------------------------------------------
% Runs a local +/-20 percent polish around the Reyna shunt-flow candidate,
% focusing on systemic and pulmonary pulse pressure while protecting the
% direct pressure-flow and shunt-flow fit.
%
% INPUTS:
%   start_result_path - MAT file containing best_candidate.params        [-]
%   options           - optional factor grids and score weights          [-]
%
% OUTPUTS:
%   result            - struct with best candidate and exported tables   [-]
%
% ASSUMPTIONS:
%   - C.SAR is the main lever for systemic arterial pulse pressure.
%   - C.PAR is the main lever for pulmonary arterial pulse pressure.
%   - Q_shunt_Lmin is a derived target and must remain protected, but it
%     should not overpower direct catheter mean-pressure anchors.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-20
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1 || isempty(start_result_path)
    start_result_path = default_start_path();
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

loaded = load(start_result_path, 'best_candidate');
params_start = loaded.best_candidate.params;
start_eval = evaluate_candidate(params_start, clinical, scenario);

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_reyna_pulse_pressure_20pct_finetune', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end

names = {'C.SAR','C.PAR','group.R_sys_scale','R.SVEN','vsd.Cd','group.R_pul_scale'};
base_values = values_from_params(params_start, params_scaled, names, case_profile);
[lb, ub] = bounds_for_names(calib_seed, names);

fprintf('[run_reyna_pulse_pressure_20pct_finetune] Start pulse RMSE %.4f | direct RMSE %.4f | Qshunt %.2f%%\n', ...
    start_eval.rmse_pulse, start_eval.rmse_direct, start_eval.Q_shunt_error_pct);

rows = {};
best_eval = start_eval;
params_best = params_start;
best_values = base_values;
best_score = candidate_score(start_eval, options);
eval_count = 0;

for i_csar = 1:numel(options.C_SAR_factors)
    for i_cpar = 1:numel(options.C_PAR_factors)
        for i_rsys = 1:numel(options.R_sys_factors)
            for i_sven = 1:numel(options.R_SVEN_factors)
                for i_cd = 1:numel(options.vsd_Cd_factors)
                    for i_rpul = 1:numel(options.R_pul_factors)
                        factors = [options.C_SAR_factors(i_csar), ...
                            options.C_PAR_factors(i_cpar), ...
                            options.R_sys_factors(i_rsys), ...
                            options.R_SVEN_factors(i_sven), ...
                            options.vsd_Cd_factors(i_cd), ...
                            options.R_pul_factors(i_rpul)];
                        trial_values = min(max(base_values(:) .* factors(:), lb(:)), ub(:));
                        params_trial = apply_values(params_start, params_scaled, ...
                            names, trial_values, case_profile);
                        trial_eval = evaluate_candidate(params_trial, clinical, scenario);
                        eval_count = eval_count + 1;
                        rows(end + 1, :) = trial_row(factors, trial_values, trial_eval); %#ok<AGROW>

                        trial_score = candidate_score(trial_eval, options);
                        if trial_score < best_score
                            best_score = trial_score;
                            best_eval = trial_eval;
                            params_best = params_trial;
                            best_values = trial_values;
                            fprintf('  New best %03d | pulse %.4f | direct %.4f | Qshunt %.2f%% | SAPpp %.2f%% | PAPpp %.2f%%\n', ...
                                eval_count, best_eval.rmse_pulse, best_eval.rmse_direct, ...
                                best_eval.Q_shunt_error_pct, best_eval.SAP_pulse_error_pct, ...
                                best_eval.PAP_pulse_error_pct);
                        end
                    end
                end
            end
        end
    end
end

grid_table = cell2table(rows, 'VariableNames', grid_variable_names());
comparison_table = comparison_table_from_evals(start_eval, best_eval);
parameter_table = table(names(:), base_values(:), best_values(:), ...
    best_values(:) ./ max(abs(base_values(:)), 1e-12), lb(:), ub(:), ...
    'VariableNames', {'Parameter','StartValue','BestValue', ...
    'RatioToStart','LowerBound','UpperBound'});

grid_csv = fullfile(tables_dir, 'reyna_pulse_pressure_20pct_all_trials.csv');
comparison_csv = fullfile(tables_dir, 'reyna_pulse_pressure_20pct_metrics.csv');
parameter_csv = fullfile(tables_dir, 'reyna_pulse_pressure_20pct_parameters.csv');
primary_csv = fullfile(tables_dir, 'reyna_pulse_pressure_20pct_primary_selection.csv');
manifest_path = fullfile(run_dir, 'run_manifest.txt');
mat_path = fullfile(mat_dir, 'reyna_pulse_pressure_20pct_result.mat');

writetable(grid_table, grid_csv);
writetable(comparison_table, comparison_csv);
writetable(parameter_table, parameter_csv);
writetable(primary_selection_table, primary_csv);

best_candidate = struct();
best_candidate.label = 'pulse_pressure_20pct_finetune';
best_candidate.params = params_best;
best_candidate.metrics = best_eval.metrics;
best_candidate.validity = best_eval.validity;
best_candidate.report_table = comparison_table;
best_candidate.parameter_table = parameter_table;
best_candidate.rmse_direct = best_eval.rmse_direct;
best_candidate.rmse_pulse = best_eval.rmse_pulse;

result = struct();
result.run_dir = run_dir;
result.start_result_path = start_result_path;
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

fprintf('[run_reyna_pulse_pressure_20pct_finetune] Done. Pulse RMSE %.4f -> %.4f | direct RMSE %.4f -> %.4f\n', ...
    start_eval.rmse_pulse, best_eval.rmse_pulse, ...
    start_eval.rmse_direct, best_eval.rmse_direct);
disp(comparison_table);
end

function options = default_options(options)
defaults = struct();
defaults.C_SAR_factors = [0.80, 0.85, 0.90, 1.00];
defaults.C_PAR_factors = [1.00, 1.10, 1.20];
defaults.R_sys_factors = [1.00, 1.05];
defaults.R_SVEN_factors = [0.95, 1.00, 1.05];
defaults.vsd_Cd_factors = [0.95, 1.00, 1.05];
defaults.R_pul_factors = [0.95, 1.00, 1.05];
defaults.direct_weight = 1.00;
defaults.pulse_weight = 0.75;
defaults.shunt_weight = 0.35;
defaults.direct_gate_pct = 12.0;
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
path = latest_existing_result_file(root, '*_reyna_shunt_flow_finetune', ...
    fullfile('mat', 'reyna_shunt_flow_finetune_result.mat'), ...
    'Reyna shunt-flow finetune result');
end

function score = candidate_score(eval, options)
if ~eval.valid || ~eval.physiology_valid
    score = Inf;
    return;
end
direct_excess = max(0, eval.max_direct_abs_error_pct - options.direct_gate_pct) / 100;
score = options.direct_weight * eval.rmse_direct + ...
    options.pulse_weight * eval.rmse_pulse + ...
    options.shunt_weight * abs(eval.Q_shunt_error_pct) / 100 + ...
    0.45 * direct_excess^2;
end

function names = grid_variable_names()
names = {'C_SAR_factor','C_PAR_factor','R_sys_factor','R_SVEN_factor', ...
    'vsd_Cd_factor','R_pul_factor','C_SAR_value','C_PAR_value', ...
    'R_sys_value','R_SVEN_value','vsd_Cd_value','R_pul_value', ...
    'DirectRMSE','PulseRMSE','MaxDirectAbsErrorPct','PhysiologyValid', ...
    'RAP_mean_ErrorPct','PAP_mean_ErrorPct','SAP_mean_ErrorPct', ...
    'QpQs_ErrorPct','CO_Lmin_ErrorPct','Q_shunt_Lmin_ErrorPct', ...
    'SAP_pulse_ErrorPct','PAP_pulse_ErrorPct','SAP_min_ErrorPct', ...
    'SAP_max_ErrorPct','PAP_min_ErrorPct','PAP_max_ErrorPct'};
end

function row = trial_row(factors, values, eval)
row = num2cell([factors(:)', values(:)', eval.rmse_direct, eval.rmse_pulse, ...
    eval.max_direct_abs_error_pct, double(eval.physiology_valid), ...
    eval.RAP_mean_error_pct, eval.PAP_mean_error_pct, eval.SAP_mean_error_pct, ...
    eval.QpQs_error_pct, eval.CO_Lmin_error_pct, eval.Q_shunt_error_pct, ...
    eval.SAP_pulse_error_pct, eval.PAP_pulse_error_pct, ...
    eval.SAP_min_error_pct, eval.SAP_max_error_pct, ...
    eval.PAP_min_error_pct, eval.PAP_max_error_pct]);
end

function eval = evaluate_candidate(params, clinical, scenario)
eval = empty_eval();
try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
    validity = evaluate_simulation_validity(sim, params, metrics, scenario, clinical);
catch
    return;
end

src = clinical.(scenario);
targets = struct();
targets.RAP_mean = src.RAP_mean_mmHg;
targets.PAP_mean = src.PAP_mean_mmHg;
targets.SAP_mean = src.SAP_mean_mmHg;
targets.QpQs = src.QpQs;
targets.CO_Lmin = src.CO_Lmin;
targets.Q_shunt_Lmin = src.Q_shunt_Lmin;
targets.SAP_min = src.SAP_dia_mmHg;
targets.SAP_max = src.SAP_sys_mmHg;
targets.PAP_min = src.PAP_dia_mmHg;
targets.PAP_max = src.PAP_sys_mmHg;
targets.SAP_pulse = src.SAP_sys_mmHg - src.SAP_dia_mmHg;
targets.PAP_pulse = src.PAP_sys_mmHg - src.PAP_dia_mmHg;

direct_names = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
direct_errors = nan(numel(direct_names), 1);
for idx = 1:numel(direct_names)
    direct_errors(idx) = error_pct(metrics.(direct_names{idx}), targets.(direct_names{idx}));
end

eval.valid = true;
eval.physiology_valid = validity.is_valid;
eval.metrics = metrics;
eval.validity = validity;
eval.rmse_direct = sqrt(mean((direct_errors / 100).^2));
eval.RAP_mean_error_pct = direct_errors(1);
eval.PAP_mean_error_pct = direct_errors(2);
eval.SAP_mean_error_pct = direct_errors(3);
eval.QpQs_error_pct = direct_errors(4);
eval.CO_Lmin_error_pct = direct_errors(5);
eval.Q_shunt_error_pct = error_pct(metrics.Q_shunt_Lmin, targets.Q_shunt_Lmin);
eval.SAP_min_error_pct = error_pct(metrics.SAP_min, targets.SAP_min);
eval.SAP_max_error_pct = error_pct(metrics.SAP_max, targets.SAP_max);
eval.PAP_min_error_pct = error_pct(metrics.PAP_min, targets.PAP_min);
eval.PAP_max_error_pct = error_pct(metrics.PAP_max, targets.PAP_max);
eval.SAP_pulse_error_pct = error_pct(metrics.SAP_pulse, targets.SAP_pulse);
eval.PAP_pulse_error_pct = error_pct(metrics.PAP_pulse, targets.PAP_pulse);
eval.rmse_pulse = sqrt(mean(([eval.SAP_pulse_error_pct; ...
    eval.PAP_pulse_error_pct] / 100).^2));
eval.max_direct_abs_error_pct = max(abs(direct_errors));
end

function eval = empty_eval()
eval = struct();
eval.valid = false;
eval.physiology_valid = false;
eval.metrics = struct();
eval.validity = struct('is_valid', false, 'failed_flags', {{}}); 
eval.rmse_direct = Inf;
eval.rmse_pulse = Inf;
eval.max_direct_abs_error_pct = Inf;
eval.RAP_mean_error_pct = Inf;
eval.PAP_mean_error_pct = Inf;
eval.SAP_mean_error_pct = Inf;
eval.QpQs_error_pct = Inf;
eval.CO_Lmin_error_pct = Inf;
eval.Q_shunt_error_pct = Inf;
eval.SAP_min_error_pct = Inf;
eval.SAP_max_error_pct = Inf;
eval.PAP_min_error_pct = Inf;
eval.PAP_max_error_pct = Inf;
eval.SAP_pulse_error_pct = Inf;
eval.PAP_pulse_error_pct = Inf;
end

function value = error_pct(model_value, clinical_value)
value = 100 * (model_value - clinical_value) / max(abs(clinical_value), 1e-9);
end

function tbl = comparison_table_from_evals(start_eval, best_eval)
metric = {'RAP_mean'; 'PAP_mean'; 'SAP_mean'; 'QpQs'; 'CO_Lmin'; ...
    'Q_shunt_Lmin'; 'SAP_pulse'; 'PAP_pulse'; 'SAP_min'; 'SAP_max'; ...
    'PAP_min'; 'PAP_max'};
start_error = [start_eval.RAP_mean_error_pct; start_eval.PAP_mean_error_pct; ...
    start_eval.SAP_mean_error_pct; start_eval.QpQs_error_pct; ...
    start_eval.CO_Lmin_error_pct; start_eval.Q_shunt_error_pct; ...
    start_eval.SAP_pulse_error_pct; start_eval.PAP_pulse_error_pct; ...
    start_eval.SAP_min_error_pct; start_eval.SAP_max_error_pct; ...
    start_eval.PAP_min_error_pct; start_eval.PAP_max_error_pct];
best_error = [best_eval.RAP_mean_error_pct; best_eval.PAP_mean_error_pct; ...
    best_eval.SAP_mean_error_pct; best_eval.QpQs_error_pct; ...
    best_eval.CO_Lmin_error_pct; best_eval.Q_shunt_error_pct; ...
    best_eval.SAP_pulse_error_pct; best_eval.PAP_pulse_error_pct; ...
    best_eval.SAP_min_error_pct; best_eval.SAP_max_error_pct; ...
    best_eval.PAP_min_error_pct; best_eval.PAP_max_error_pct];
tbl = table(metric, start_error, best_error, abs(start_error), abs(best_error), ...
    abs(best_error) - abs(start_error), ...
    'VariableNames', {'Metric','StartError_pct','BestError_pct', ...
    'StartAbsError_pct','BestAbsError_pct','DeltaAbsError_pct'});
end

function [lb, ub] = bounds_for_names(calib_seed, names)
registry = calib_seed.parameterRegistry;
[~, loc] = ismember(names, registry.name);
if any(loc == 0)
    missing = names(loc == 0);
    error('run_reyna_pulse_pressure_20pct_finetune:missingRegistryName', ...
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
    error('run_reyna_pulse_pressure_20pct_finetune:manifestOpenFailed', ...
        'Unable to write manifest: %s', path);
end
cleaner = onCleanup(@() fclose(fid));
fprintf(fid, 'Reyna Pulse-Pressure 20pct Finetune\n');
fprintf(fid, '===================================\n');
fprintf(fid, 'StartResult: %s\n', result.start_result_path);
fprintf(fid, 'StartDirectRMSE: %.6f\n', result.start_eval.rmse_direct);
fprintf(fid, 'BestDirectRMSE: %.6f\n', result.best_eval.rmse_direct);
fprintf(fid, 'StartPulseRMSE: %.6f\n', result.start_eval.rmse_pulse);
fprintf(fid, 'BestPulseRMSE: %.6f\n', result.best_eval.rmse_pulse);
fprintf(fid, 'StartQshuntErrorPct: %.6f\n', result.start_eval.Q_shunt_error_pct);
fprintf(fid, 'BestQshuntErrorPct: %.6f\n', result.best_eval.Q_shunt_error_pct);
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
    error('run_reyna_pulse_pressure_20pct_finetune:missingDefaultStartPath', ...
        'No %s found under %s. Pass start_result_path explicitly.', ...
        label, runs_dir);
end
end
