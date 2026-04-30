%% diagnose_rmse_plateau.m
% DIAGNOSE_RMSE_PLATEAU
% -----------------------------------------------------------------------
% Lightweight diagnostic for the current calibration plateau. Reads the
% latest saved calibration package and summarizes:
%   - grouped RMSE contributions
%   - worst calibrated residuals
%   - failed validity gates
%   - active parameters parked on bounds
%
% USAGE:
%   >> run('scripts/diagnose_rmse_plateau.m')
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-30
% VERSION:  1.0
% -----------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

scenario = 'pre_surgery';
results_file = fullfile(root, 'results', 'tables', sprintf('params_calibrated_%s.mat', scenario));
if ~isfile(results_file)
    error('diagnose_rmse_plateau:missingResultsFile', ...
        'Could not find %s. Run main_run first.', results_file);
end

data = load(results_file, 'report', 'calib_out', 'validity_base', 'validity_cal');
report = data.report;
calib_out = data.calib_out;
validity_base = data.validity_base;
validity_cal = data.validity_cal;

fprintf('============================================================\n');
fprintf('RMSE PLATEAU DIAGNOSTIC  |  scenario=%s\n', scenario);
fprintf('============================================================\n');
fprintf('Overall RMSE baseline:   %.4f\n', report.rmse_baseline);
fprintf('Overall RMSE calibrated: %.4f\n', report.rmse_cal);
fprintf('\n');

group_tbl = grouped_rmse_table(report.table_cal);
fprintf('--- Grouped RMSE contribution ---\n');
disp(group_tbl);

fprintf('--- Worst calibrated residuals ---\n');
disp(report.sorted_errors(1:min(8, height(report.sorted_errors)), :));

fprintf('--- Validity gates (baseline) ---\n');
print_validity_summary(validity_base);
fprintf('--- Validity gates (calibrated) ---\n');
print_validity_summary(validity_cal);

fprintf('--- Active parameters near bounds ---\n');
bound_tbl = active_bound_hits(calib_out);
if isempty(bound_tbl)
    fprintf('No active parameters are within 2%% of their lower/upper bounds.\n');
else
    disp(bound_tbl);
end

fprintf('--- Objective breakdown (top weighted contributors) ---\n');
if isfield(calib_out, 'objective_breakdown') && ~isempty(calib_out.objective_breakdown)
    breakdown = calib_out.objective_breakdown;
    breakdown = sortrows(breakdown, 'WeightedContribution', 'descend');
    disp(breakdown(1:min(8, height(breakdown)), :));
else
    fprintf('No objective breakdown stored in results package.\n');
end

function tbl = grouped_rmse_table(table_cal)
metric_group = strings(height(table_cal), 1);
for i = 1:height(table_cal)
    metric_group(i) = classify_metric(table_cal.Metric{i});
end

groups = unique(metric_group, 'stable');
n_groups = numel(groups);
tbl = table('Size', [n_groups, 4], ...
    'VariableTypes', {'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Group', 'NMetrics', 'RMSE', 'MeanAbsError_pct'});

for i = 1:n_groups
    mask = metric_group == groups(i) & ~isnan(table_cal.Clinical) & ~isnan(table_cal.Calibrated);
    clinical = table_cal.Clinical(mask);
    model = table_cal.Calibrated(mask);
    err_pct = abs(table_cal.Error_pct(mask));
    tbl.Group(i) = groups(i);
    tbl.NMetrics(i) = sum(mask);
    if any(mask)
        rel_err = (model - clinical) ./ max(abs(clinical), 1e-9);
        tbl.RMSE(i) = sqrt(mean(rel_err.^2));
        tbl.MeanAbsError_pct(i) = mean(err_pct);
    else
        tbl.RMSE(i) = NaN;
        tbl.MeanAbsError_pct(i) = NaN;
    end
end
end

function group_name = classify_metric(metric_name)
if ismember(metric_name, {'QpQs', 'PVR', 'SVR', 'CO_Lmin', 'VSD_frac_pct'})
    group_name = "flow_resistance";
elseif contains(metric_name, 'EF')
    group_name = "function";
elseif contains(metric_name, 'EDV') || contains(metric_name, 'ESV')
    group_name = "volume";
else
    group_name = "pressure";
end
end

function print_validity_summary(validity)
if validity.is_valid
    fprintf('PASS | no failed flags\n');
else
    fprintf('FAIL | penalty=%.1f | flags=%s\n', ...
        validity.penalty, strjoin(validity.failed_flags, ', '));
end
if isfield(validity, 'bound_summary')
    disp(validity.bound_summary);
end
fprintf('\n');
end

function tbl = active_bound_hits(calib_out)
if ~all(isfield(calib_out, {'names', 'xbest', 'lb', 'ub'}))
    tbl = [];
    return;
end

span = max(calib_out.ub - calib_out.lb, 1e-9);
dist_to_lb = (calib_out.xbest - calib_out.lb) ./ span;
dist_to_ub = (calib_out.ub - calib_out.xbest) ./ span;
mask = dist_to_lb <= 0.02 | dist_to_ub <= 0.02;
if ~any(mask)
    tbl = [];
    return;
end

side = repmat("mid", numel(calib_out.names), 1);
side(dist_to_lb <= 0.02) = "lower";
side(dist_to_ub <= 0.02) = "upper";
tbl = table(string(calib_out.names(mask)), calib_out.xbest(mask), ...
    calib_out.lb(mask), calib_out.ub(mask), side(mask), ...
    'VariableNames', {'Parameter', 'Value', 'LowerBound', 'UpperBound', 'NearestSide'});
end

function project_path = build_clean_project_path(project_root)
project_paths = strsplit(genpath(project_root), pathsep);
project_paths = project_paths(~cellfun('isempty', project_paths));
is_shadow = contains(project_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, project_paths);
project_path = strjoin(project_paths(~is_shadow & is_existing), pathsep);
end
