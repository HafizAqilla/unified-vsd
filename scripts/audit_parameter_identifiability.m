function audit = audit_parameter_identifiability(varargin)
% AUDIT_PARAMETER_IDENTIFIABILITY
% -----------------------------------------------------------------------
% Local finite-difference identifiability audit for case-profile calibration.
%
% For each benchmark case, this script perturbs the active calibration
% parameters around the clinical-seeded baseline and ranks parameters by
% their normalized sensitivity to the metrics allowed by the case profile.
% It also flags highly correlated sensitivity columns, which indicate
% parameter pairs that are difficult to estimate independently.
%
% INPUTS:
%   Name/value options:
%     Cases             - cell array: reyna | razka | patient_profile_A  [-]
%     Scenario          - 'pre_surgery' | 'post_surgery'                 [-]
%     RelativeStep      - fractional parameter perturbation              [-]
%     Fidelity          - 'moderate' or 'strict' solver settings         [-]
%     WriteTables       - write CSV summary tables                       [-]
%     MaxParameters     - optional cap for smoke tests                   [-]
%
% OUTPUTS:
%   audit             - struct array with tables and sensitivity matrix  [-]
%
% ASSUMPTIONS:
%   - The audit is local to the parameter vector being tested.
%   - Moderate fidelity uses the same relaxed steady-state settings as the
%     calibration search; use Fidelity='strict' before final claims.
%
% REFERENCES:
%   [1] Brun et al. (2001). Practical identifiability of dynamic models.
%   [2] Saltelli et al. (2010). Global sensitivity analysis.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.1
% -----------------------------------------------------------------------

project_root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(project_root));

opts = parse_options(project_root, varargin{:});
case_names = cellstr(opts.Cases);
n_cases = numel(case_names);
audit = repmat(empty_case_audit(), n_cases, 1);

for case_idx = 1:n_cases
    clinical = load_case_profile(case_names{case_idx});
    case_profile = build_case_calibration_profile(clinical, opts.Scenario);
    params0 = build_clinical_seed_params(clinical, opts.Scenario, case_profile);
    params0 = apply_audit_fidelity(params0, opts.Fidelity);

    calib = calibration_param_sets(opts.Scenario, params0, [], {}, case_profile);
    [metric_names, clinical_values] = available_metric_targets( ...
        clinical, opts.Scenario, calib.metricFields);
    parameter_names = calib.names(:)';
    if isfinite(opts.MaxParameters)
        n_keep = min(numel(parameter_names), opts.MaxParameters);
        parameter_names = parameter_names(1:n_keep);
    end

    fprintf('\n[audit_parameter_identifiability] %s | mode=%s | %d params | %d metrics\n', ...
        case_profile.case_id, case_profile.mode, numel(parameter_names), numel(metric_names));

    case_audit = evaluate_case_identifiability(params0, calib, case_profile, ...
        metric_names, clinical_values, parameter_names, opts);
    audit(case_idx) = case_audit;

    if opts.WriteTables
        write_case_tables(case_audit, opts.OutputDir);
    end
end
end

function case_audit = evaluate_case_identifiability(params0, calib, case_profile, ...
    metric_names, clinical_values, parameter_names, opts)
% EVALUATE_CASE_IDENTIFIABILITY - compute local sensitivities for one case.
base_metrics = safe_metric_vector(params0, metric_names);
n_metrics = numel(metric_names);
n_parameters = numel(parameter_names);
sensitivity_matrix = nan(n_metrics, n_parameters);

for parameter_idx = 1:n_parameters
    parameter_name = parameter_names{parameter_idx};
    x0 = get_calibration_param_value(params0, calib.referenceParams, parameter_name, case_profile);
    step_abs = opts.RelativeStep * max(abs(x0), opts.ParameterScaleFloor);
    [x_minus, x_plus] = bounded_step(x0, step_abs, calib, parameter_name);

    params_minus = set_calibration_param_value( ...
        params0, calib.referenceParams, parameter_name, x_minus, case_profile);
    params_plus = set_calibration_param_value( ...
        params0, calib.referenceParams, parameter_name, x_plus, case_profile);
    y_minus = safe_metric_vector(params_minus, metric_names);
    y_plus = safe_metric_vector(params_plus, metric_names);

    delta_x_rel = (x_plus - x_minus) / max(abs(x0), opts.ParameterScaleFloor);
    metric_scale = max([abs(base_metrics(:)), abs(clinical_values(:)), ...
        opts.MetricScaleFloor * ones(n_metrics, 1)], [], 2);
    sensitivity_matrix(:, parameter_idx) = ((y_plus(:) - y_minus(:)) ./ metric_scale) ./ ...
        max(abs(delta_x_rel), eps);
end

case_audit = empty_case_audit();
case_audit.case_id = case_profile.case_id;
case_audit.case_mode = case_profile.mode;
case_audit.scenario = case_profile.scenario;
case_audit.metric_names = metric_names(:);
case_audit.parameter_names = parameter_names(:);
case_audit.clinical_values = clinical_values(:);
case_audit.base_values = base_metrics(:);
case_audit.sensitivity_matrix = sensitivity_matrix;
case_audit.parameter_table = build_parameter_table(case_audit, opts.WeakScoreThreshold);
case_audit.metric_table = build_metric_table(case_audit);
case_audit.correlation_table = build_correlation_table(case_audit, opts.CorrelationThreshold);
case_audit.recommendation_table = build_recommendation_table(case_audit);
print_case_summary(case_audit);
end

function parameter_table = build_parameter_table(case_audit, weak_score_threshold)
% BUILD_PARAMETER_TABLE - rank parameters by normalized sensitivity score.
sensitivity_matrix = case_audit.sensitivity_matrix;
parameter_names = case_audit.parameter_names;
metric_names = case_audit.metric_names;
n_parameters = numel(parameter_names);

score_col = nan(n_parameters, 1);
max_abs_col = nan(n_parameters, 1);
dominant_metric_col = cell(n_parameters, 1);
weak_col = false(n_parameters, 1);

for parameter_idx = 1:n_parameters
    column_values = sensitivity_matrix(:, parameter_idx);
    finite_values = column_values(isfinite(column_values));
    if isempty(finite_values)
        dominant_metric_col{parameter_idx} = '';
        weak_col(parameter_idx) = true;
        continue;
    end
    score_col(parameter_idx) = norm(finite_values, 2) / sqrt(numel(finite_values));
    [max_abs_col(parameter_idx), metric_idx] = max(abs(column_values));
    dominant_metric_col{parameter_idx} = metric_names{metric_idx};
    weak_col(parameter_idx) = score_col(parameter_idx) < weak_score_threshold;
end

parameter_table = table(parameter_names(:), repmat({case_audit.case_mode}, n_parameters, 1), ...
    score_col, max_abs_col, dominant_metric_col, weak_col, ...
    'VariableNames', {'Parameter', 'CaseMode', 'IdentifiabilityScore', ...
    'MaxAbsSensitivity', 'DominantMetric', 'WeaklyIdentifiable'});
parameter_table = sortrows(parameter_table, 'IdentifiabilityScore', 'descend');
end

function metric_table = build_metric_table(case_audit)
% BUILD_METRIC_TABLE - rank metric dependence on active parameters.
sensitivity_matrix = case_audit.sensitivity_matrix;
n_metrics = numel(case_audit.metric_names);
max_abs_col = nan(n_metrics, 1);
dominant_parameter_col = cell(n_metrics, 1);

for metric_idx = 1:n_metrics
    row_values = sensitivity_matrix(metric_idx, :);
    [max_abs_col(metric_idx), parameter_idx] = max(abs(row_values));
    if isempty(parameter_idx) || isnan(max_abs_col(metric_idx))
        dominant_parameter_col{metric_idx} = '';
    else
        dominant_parameter_col{metric_idx} = case_audit.parameter_names{parameter_idx};
    end
end

metric_table = table(case_audit.metric_names(:), case_audit.clinical_values(:), ...
    case_audit.base_values(:), max_abs_col, dominant_parameter_col, ...
    'VariableNames', {'Metric', 'ClinicalValue', 'BaselineValue', ...
    'MaxAbsSensitivity', 'DominantParameter'});
metric_table = sortrows(metric_table, 'MaxAbsSensitivity', 'descend');
end

function correlation_table = build_correlation_table(case_audit, correlation_threshold)
% BUILD_CORRELATION_TABLE - flag highly similar parameter sensitivity columns.
parameter_names = case_audit.parameter_names;
sensitivity_matrix = case_audit.sensitivity_matrix;
rows = {};

for i = 1:numel(parameter_names)
    for j = i + 1:numel(parameter_names)
        corr_ij = column_correlation(sensitivity_matrix(:, i), sensitivity_matrix(:, j));
        if isfinite(corr_ij) && abs(corr_ij) >= correlation_threshold
            rows(end + 1, :) = {parameter_names{i}, parameter_names{j}, corr_ij, abs(corr_ij)}; %#ok<AGROW>
        end
    end
end

if isempty(rows)
    correlation_table = table(cell(0, 1), cell(0, 1), zeros(0, 1), zeros(0, 1), ...
        'VariableNames', {'ParameterA', 'ParameterB', 'Correlation', 'AbsCorrelation'});
else
    correlation_table = cell2table(rows, ...
        'VariableNames', {'ParameterA', 'ParameterB', 'Correlation', 'AbsCorrelation'});
    correlation_table = sortrows(correlation_table, 'AbsCorrelation', 'descend');
end
end

function recommendation_table = build_recommendation_table(case_audit)
% BUILD_RECOMMENDATION_TABLE - convert audit findings into calibration actions.
rows = {};
weak_rows = case_audit.parameter_table(case_audit.parameter_table.WeaklyIdentifiable, :);
for weak_idx = 1:height(weak_rows)
    rows(end + 1, :) = {case_audit.case_id, case_audit.case_mode, ...
        'weak_parameter', weak_rows.Parameter{weak_idx}, '', ...
        'freeze_or_strengthen_prior', ...
        sprintf('Score %.4f; dominant metric %s.', ...
        weak_rows.IdentifiabilityScore(weak_idx), weak_rows.DominantMetric{weak_idx})}; %#ok<AGROW>
end

for corr_idx = 1:height(case_audit.correlation_table)
    corr_row = case_audit.correlation_table(corr_idx, :);
    rows(end + 1, :) = {case_audit.case_id, case_audit.case_mode, ...
        'correlated_pair', corr_row.ParameterA{1}, corr_row.ParameterB{1}, ...
        'group_or_choose_representative', ...
        sprintf('Abs sensitivity correlation %.4f.', corr_row.AbsCorrelation)}; %#ok<AGROW>
end

if isempty(rows)
    recommendation_table = table(cell(0, 1), cell(0, 1), cell(0, 1), ...
        cell(0, 1), cell(0, 1), cell(0, 1), cell(0, 1), ...
        'VariableNames', {'CaseID', 'CaseMode', 'FindingType', 'ParameterA', ...
        'ParameterB', 'RecommendedAction', 'Rationale'});
else
    recommendation_table = cell2table(rows, ...
        'VariableNames', {'CaseID', 'CaseMode', 'FindingType', 'ParameterA', ...
        'ParameterB', 'RecommendedAction', 'Rationale'});
end
end

function [metric_names, clinical_values] = available_metric_targets(clinical, scenario, allowed_metrics)
% AVAILABLE_METRIC_TARGETS - return case-profile metrics with finite targets.
targets = get_calibration_targets(scenario, clinical);
target_names = {targets.Metric};
clinical_all = [targets.ClinicalValue];
is_available = isfinite(clinical_all) & ismember(target_names, allowed_metrics);
metric_names = target_names(is_available);
clinical_values = clinical_all(is_available);
end

function metric_values = safe_metric_vector(params, metric_names)
% SAFE_METRIC_VECTOR - run model and extract requested metrics.
try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
catch
    metric_values = nan(numel(metric_names), 1);
    return;
end

metric_values = nan(numel(metric_names), 1);
for metric_idx = 1:numel(metric_names)
    metric_name = metric_names{metric_idx};
    if isfield(metrics, metric_name) && isfinite(metrics.(metric_name))
        metric_values(metric_idx) = metrics.(metric_name);
    end
end
end

function params = build_clinical_seed_params(clinical, scenario, case_profile)
% BUILD_CLINICAL_SEED_PARAMS - reproduce main_run scaling and clinical mapping.
params_ref = default_parameters();
patient.age_years = clinical.common.age_years;       % [years]
patient.age_days = clinical.common.age_years * 365.25; % [days]
patient.weight_kg = clinical.common.weight_kg;       % [kg]
patient.height_cm = clinical.common.height_cm;       % [cm]
patient.sex = clinical.common.sex;                   % [-]
patient.maturation_mode = 'normal';                  % [-]
if isfield(clinical.common, 'maturation_mode') && ~isempty(clinical.common.maturation_mode)
    patient.maturation_mode = clinical.common.maturation_mode;
end
if isfield(clinical.common, 'BSA') && isfinite(clinical.common.BSA)
    patient.BSA = clinical.common.BSA;               % [m^2]
end
params_scaled = apply_scaling(params_ref, patient);
params = params_from_clinical(params_scaled, clinical, scenario, params_scaled, case_profile);
end

function params = apply_audit_fidelity(params, fidelity)
% APPLY_AUDIT_FIDELITY - choose documented solver fidelity for audit runs.
switch lower(char(fidelity))
    case 'strict'
        return;
    case 'moderate'
        params.sim.nCyclesSteady = min(max(params.sim.nCyclesSteady, 30), 60);
        params.sim.ss_tol_P = max(params.sim.ss_tol_P, 0.5);
        params.sim.ss_tol_V = max(params.sim.ss_tol_V, 0.5);
    otherwise
        error('audit_parameter_identifiability:unknownFidelity', ...
            'Fidelity must be ''moderate'' or ''strict''.');
end
end

function [x_minus, x_plus] = bounded_step(x0, step_abs, calib, parameter_name)
% BOUNDED_STEP - perturb parameter while respecting calibration bounds.
name_idx = find(strcmp(calib.names, parameter_name), 1, 'first');
if isempty(name_idx)
    x_minus = x0 - step_abs;
    x_plus = x0 + step_abs;
    return;
end
x_minus = max(calib.lb(name_idx), x0 - step_abs);
x_plus = min(calib.ub(name_idx), x0 + step_abs);
if x_minus == x_plus
    x_minus = x0;
    x_plus = x0;
end
end

function corr_value = column_correlation(a, b)
% COLUMN_CORRELATION - Pearson correlation for finite sensitivity entries.
finite_mask = isfinite(a) & isfinite(b);
if nnz(finite_mask) < 3
    corr_value = NaN;
    return;
end
a = a(finite_mask) - mean(a(finite_mask));
b = b(finite_mask) - mean(b(finite_mask));
denom = norm(a, 2) * norm(b, 2);
if denom <= eps
    corr_value = NaN;
else
    corr_value = (a' * b) / denom;
end
end

function write_case_tables(case_audit, output_dir)
% WRITE_CASE_TABLES - persist CSV tables for review and reproducibility.
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
safe_case_id = regexprep(lower(case_audit.case_id), '[^a-z0-9_]+', '_');
prefix = fullfile(output_dir, sprintf('%s_%s_%s', ...
    datestr(now, 'yyyymmdd_HHMMSS'), safe_case_id, case_audit.case_mode));
writetable(case_audit.parameter_table, [prefix '_parameters.csv']);
writetable(case_audit.metric_table, [prefix '_metrics.csv']);
writetable(case_audit.correlation_table, [prefix '_correlations.csv']);
writetable(case_audit.recommendation_table, [prefix '_recommendations.csv']);
fprintf('[audit_parameter_identifiability] Wrote tables: %s_*.csv\n', prefix);
end

function print_case_summary(case_audit)
% PRINT_CASE_SUMMARY - concise console summary for one audited case.
fprintf('--- Parameter identifiability ranking: %s (%s) ---\n', ...
    case_audit.case_id, case_audit.case_mode);
disp(case_audit.parameter_table(1:min(8, height(case_audit.parameter_table)), :));

weak_count = nnz(case_audit.parameter_table.WeaklyIdentifiable);
fprintf('Weakly identifiable parameters: %d/%d\n', weak_count, height(case_audit.parameter_table));
if ~isempty(case_audit.correlation_table)
    fprintf('Highly correlated parameter pairs: %d\n', height(case_audit.correlation_table));
    disp(case_audit.correlation_table(1:min(8, height(case_audit.correlation_table)), :));
else
    fprintf('Highly correlated parameter pairs: 0\n');
end
if ~isempty(case_audit.recommendation_table)
    fprintf('Actionable audit recommendations: %d\n', height(case_audit.recommendation_table));
end
end

function clinical = load_case_profile(case_name)
% LOAD_CASE_PROFILE - map short case names to clinical profile functions.
switch lower(char(case_name))
    case {'reyna', 'patient_reyna'}
        clinical = patient_reyna();
    case {'razka', 'patient_profile_razka'}
        clinical = patient_profile_Razka();
    case {'patient_profile_a', 'profile_a', 'a'}
        clinical = patient_profile_A();
    otherwise
        error('audit_parameter_identifiability:unknownCase', ...
            'Unknown case name: %s', char(case_name));
end
end

function opts = parse_options(project_root, varargin)
parser = inputParser();
addParameter(parser, 'Cases', {'reyna', 'razka', 'patient_profile_A'}, ...
    @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(parser, 'Scenario', 'pre_surgery', @(x) ischar(x) || isstring(x));
addParameter(parser, 'RelativeStep', 0.03, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'Fidelity', 'moderate', @(x) ischar(x) || isstring(x));
addParameter(parser, 'WriteTables', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'MaxParameters', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'WeakScoreThreshold', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(parser, 'CorrelationThreshold', 0.95, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
addParameter(parser, 'MetricScaleFloor', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'ParameterScaleFloor', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(parser, 'OutputDir', fullfile(project_root, 'results', 'identifiability'), ...
    @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
if ischar(opts.Cases) || isstring(opts.Cases)
    opts.Cases = cellstr(opts.Cases);
end
opts.Scenario = char(opts.Scenario);
opts.Fidelity = char(opts.Fidelity);
opts.OutputDir = char(opts.OutputDir);
opts.WriteTables = logical(opts.WriteTables);
end

function case_audit = empty_case_audit()
case_audit = struct( ...
    'case_id', '', ...
    'case_mode', '', ...
    'scenario', '', ...
    'metric_names', {{}}, ...
    'parameter_names', {{}}, ...
    'clinical_values', [], ...
    'base_values', [], ...
    'sensitivity_matrix', [], ...
    'parameter_table', table(), ...
    'metric_table', table(), ...
    'correlation_table', table(), ...
    'recommendation_table', table());
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
