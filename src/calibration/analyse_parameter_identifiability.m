function report = analyse_parameter_identifiability(params_best, scenario, calib_out, gate_table, results_dir)
% ANALYSE_PARAMETER_IDENTIFIABILITY
% -----------------------------------------------------------------------
% Reports whether the retained active parameter set is JOINTLY identifiable,
% not merely individually sensitive.
%
% Motivation
% ----------
% GSA/Sobol screening ranks parameters one at a time by their marginal
% influence on each output. It cannot detect that two retained parameters
% are collinear -- for example C.SAR and group.R_sys_scale trading off along
% the systemic RC time constant that governs pulse-pressure shape (identified
% in docs/reyna_zhang_fullmetric_results_20260828.md section 3.5). Two
% collinear parameters can each look individually sensitive while carrying
% almost no independent information once the other is known.
%
% Method
% ------
% At the calibrated operating point, build the scaled sensitivity matrix
%   S(i,j) = (d y_i / d theta_j) * (theta_j / sigma_i)
% via central finite differences (1% relative step) over every governed
% metric i and active parameter j, using params_best as its own calibration
% reference so the perturbation is exactly a 1% multiplicative step around
% the current operating point regardless of the parameter's underlying
% representation (physical value or group scale).
%
% Then report:
%   - cond(S)                    condition number; >1e3 indicates near-dependence
%   - pairwise column correlation of S; |rho| > 0.9 flags a collinear pair
%   - per-parameter column norm; a near-zero column is a parameter GSA
%     should have already dropped
%
% This function REPORTS ONLY. It does not remove or fix any parameter --
% that judgement (which member of a collinear pair to hold fixed, and at
% what value) requires physiological reasoning that belongs to the user.
%
% INPUTS:
%   params_best  - calibrated parameter struct (physical, post-calibration) [-]
%   scenario     - 'pre_surgery' | 'post_surgery'                          [-]
%   calib_out    - run_calibration output struct; must carry .names        [-]
%                  (active parameter names) and .caseProfile
%   gate_table   - .table field from export_full_metric_gate, must carry
%                  Metric, Sigma, InPrimaryRMSE                            [-]
%   results_dir  - output directory for exported CSVs [char/string]        [-]
%
% OUTPUTS:
%   report - struct with:
%       .parameter_table    per-parameter table: Parameter, ColumnNorm,
%                            Inactive, MaxAbsCorrelation, MostCorrelatedWith,
%                            Flagged
%       .pair_table          every parameter pair with |rho| and a
%                             Collinear flag
%       .condition_number    cond(S)
%       .n_metrics, .n_parameters
%       .csv_path             per-parameter CSV path
%       .pair_csv_path        pair-table CSV path
%       .sensitivity_matrix   the raw S matrix [n_metrics x n_parameters]
%
% ASSUMPTIONS:
%   - Only rows with InPrimaryRMSE = true in gate_table are used as the
%     metric set i, consistent with the governed acceptance mask elsewhere.
%   - A metric or parameter that fails to evaluate (simulation error, or a
%     non-finite forward value) is dropped from S with a console warning
%     rather than aborting the whole analysis.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 3)
%   [2] docs/reyna_zhang_fullmetric_results_20260828.md section 3.5
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 5
    results_dir = '';
end

report = empty_report();

if ~isstruct(calib_out) || ~isfield(calib_out, 'names') || isempty(calib_out.names)
    warning('analyse_parameter_identifiability:noActiveParameters', ...
        'calib_out.names is empty; skipping identifiability analysis.');
    return;
end
if ~istable(gate_table) || ~all(ismember({'Metric','Sigma','InPrimaryRMSE'}, ...
        gate_table.Properties.VariableNames))
    warning('analyse_parameter_identifiability:missingGateColumns', ...
        'gate_table is missing required columns; skipping identifiability analysis.');
    return;
end

case_profile = struct();
if isfield(calib_out, 'caseProfile')
    case_profile = calib_out.caseProfile;
end

governed = gate_table(logical(gate_table.InPrimaryRMSE), :);
metric_names = governed.Metric(:)';
metric_sigma = governed.Sigma(:)';
n_metrics = numel(metric_names);
param_names = calib_out.names(:)';
n_parameters = numel(param_names);

if n_metrics == 0 || n_parameters == 0
    return;
end

relative_step = 0.01;
S = nan(n_metrics, n_parameters);
skipped_params = {};

for j = 1:n_parameters
    name = param_names{j};
    x0 = get_calibration_param_value(params_best, params_best, name, case_profile);
    if ~isfinite(x0) || x0 == 0
        skipped_params{end+1} = name; %#ok<AGROW>
        continue;
    end
    delta = abs(x0) * relative_step;

    params_plus = set_calibration_param_value(params_best, params_best, ...
        name, x0 + delta, case_profile);
    params_minus = set_calibration_param_value(params_best, params_best, ...
        name, x0 - delta, case_profile);

    y_plus = evaluate_metric_set(params_plus, metric_names);
    y_minus = evaluate_metric_set(params_minus, metric_names);

    for i = 1:n_metrics
        if isfinite(y_plus(i)) && isfinite(y_minus(i))
            dydtheta = (y_plus(i) - y_minus(i)) / (2 * delta);
            S(i, j) = dydtheta * (x0 / max(metric_sigma(i), 1e-9));
        end
    end
end

if ~isempty(skipped_params)
    warning('analyse_parameter_identifiability:degenerateParameters', ...
        ['The following active parameters had a non-finite or zero ', ...
         'operating-point value and were skipped: %s'], ...
        strjoin(skipped_params, ', '));
end

report.sensitivity_matrix = S;
report.n_metrics = n_metrics;
report.n_parameters = n_parameters;

tables_out = identifiability_tables_from_matrix(S, param_names);
report.parameter_table = tables_out.parameter_table;
report.pair_table = tables_out.pair_table;
report.condition_number = tables_out.condition_number;

if isempty(report.parameter_table)
    return;
end

if ~isempty(results_dir)
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    report.csv_path = fullfile(results_dir, ...
        sprintf('parameter_identifiability_%s.csv', scenario));
    writetable(report.parameter_table, report.csv_path);
    report.pair_csv_path = fullfile(results_dir, ...
        sprintf('parameter_identifiability_pairs_%s.csv', scenario));
    writetable(report.pair_table, report.pair_csv_path);
end

print_identifiability_report(report, scenario);
end

% =========================================================================
function report = empty_report()
report = struct( ...
    'parameter_table', table(), ...
    'pair_table', table(), ...
    'condition_number', NaN, ...
    'n_metrics', 0, ...
    'n_parameters', 0, ...
    'csv_path', '', ...
    'pair_csv_path', '', ...
    'sensitivity_matrix', []);
end

% =========================================================================
function y = evaluate_metric_set(params, metric_names)
% EVALUATE_METRIC_SET - forward-simulate params and read the requested
% metrics, mirroring the sim overrides applied in objective_calibration.m
% for consistency (evaluate_metrics local function).
y = nan(1, numel(metric_names));
try
    params.sim.nCyclesSteady = min(max(params.sim.nCyclesSteady, 30), 60);
    params.sim.ss_tol_P = max(params.sim.ss_tol_P, 0.5);
    params.sim.ss_tol_V = max(params.sim.ss_tol_V, 0.5);
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
catch
    return;
end
for i = 1:numel(metric_names)
    mf = metric_names{i};
    if isfield(metrics, mf) && isfinite(metrics.(mf))
        y(i) = metrics.(mf);
    end
end
end

% =========================================================================
function print_identifiability_report(report, scenario)
fprintf('\n--- PARAMETER IDENTIFIABILITY — %s ---\n', ...
    upper(strrep(scenario, '_', ' ')));
fprintf(['  Scaled sensitivity matrix S(i,j) = dy_i/dtheta_j * theta_j/sigma_i, ', ...
    '%d metrics x %d parameters.\n'], report.n_metrics, ...
    height(report.parameter_table));
if isfinite(report.condition_number)
    fprintf('  Condition number: %.3g', report.condition_number);
    if report.condition_number > 1e3
        fprintf('  [near-dependence: cond > 1e3]');
    end
    fprintf('\n');
else
    fprintf('  Condition number: unavailable\n');
end
disp(report.parameter_table);

if isempty(report.pair_table) || height(report.pair_table) == 0
    return;
end
collinear = report.pair_table(logical(report.pair_table.Collinear), :);
if height(collinear) == 0
    fprintf('  [PASS] No parameter pair exceeds |rho| = 0.9.\n');
else
    for idx = 1:height(collinear)
        fprintf(2, '  [COLLINEAR] %s <-> %s  (rho = %+.3f)\n', ...
            collinear.ParameterA{idx}, collinear.ParameterB{idx}, ...
            collinear.Correlation(idx));
    end
end
if ~isempty(report.csv_path)
    fprintf('  Exported: %s\n', report.csv_path);
end
if ~isempty(report.pair_csv_path)
    fprintf('  Exported: %s\n', report.pair_csv_path);
end
end
