function status = classify_calibration_run(report, plausibility_metrics, thresholds)
% CLASSIFY_CALIBRATION_RUN
% -----------------------------------------------------------------------
% Classifies a calibration run using both output fit and parameter
% plausibility rather than RMSE alone.
%
% INPUTS:
%   report                - validation_report output struct              [-]
%   plausibility_metrics  - output of evaluate_parameter_plausibility    [-]
%   thresholds            - optional thresholds struct with fields:
%       .excellent_error_pct   default 5                                [%]
%       .primary_error_pct     default 10                               [%]
%       .secondary_error_pct   default 15                               [%]
%       .max_warning_fraction  default 0.20 (reported, not rejecting)   [-]
%       .reject_on_any_fail    default true                             [-]
%
% OUTPUTS:
%   status               - struct with fields:
%       .label                ACCEPT | PROMISING_NEAR_MISS | ...
%       .fit_ok               logical
%       .plausibility_ok      logical
%       .primary_fit_ok       logical
%       .secondary_fit_ok     logical
%       .n_primary_fail       count
%       .n_secondary_fail     count
%       .n_warning            count
%       .n_fail               count
%       .warning_fraction     fraction
%       .summary              human-readable sentence
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3 || isempty(thresholds)
    thresholds = default_thresholds();
else
    thresholds = merge_thresholds(default_thresholds(), thresholds);
end
if nargin < 2 || isempty(plausibility_metrics)
    plausibility_metrics = struct();
end

status = struct();
status.thresholds = thresholds;

[primary_excellent_ok, n_primary_excellent_fail] = evaluate_primary_fit(report, thresholds.excellent_error_pct);
[primary_fit_ok, n_primary_fail] = evaluate_primary_fit(report, thresholds.primary_error_pct);
[secondary_fit_ok, n_secondary_fail] = evaluate_secondary_fit(report, thresholds.secondary_error_pct);

status.primary_excellent_ok = primary_excellent_ok;
status.primary_fit_ok = primary_fit_ok;
status.secondary_fit_ok = secondary_fit_ok;
status.fit_ok = primary_fit_ok && secondary_fit_ok;
status.n_primary_excellent_fail = n_primary_excellent_fail;
status.n_primary_fail = n_primary_fail;
status.n_secondary_fail = n_secondary_fail;
status.max_primary_abs_error = compute_max_abs_error(report, 'primary');
status.max_secondary_abs_error = compute_max_abs_error(report, 'secondary');
status.rmse_improvement_frac = compute_rmse_improvement(report);
status.rmse_not_worse = compute_rmse_not_worse(report);
status.strong_rmse_improvement = status.rmse_improvement_frac >= ...
    thresholds.accept_rmse_improvement_frac;

status.n_warning = 0;
status.n_fail = 0;
status.warning_fraction = 0;
if isfield(plausibility_metrics, 'n_warning')
    status.n_warning = plausibility_metrics.n_warning;
end
if isfield(plausibility_metrics, 'n_fail')
    status.n_fail = plausibility_metrics.n_fail;
end
if isfield(plausibility_metrics, 'warning_fraction')
    status.warning_fraction = plausibility_metrics.warning_fraction;
end

if thresholds.reject_on_any_fail
    no_plausibility_fail = (status.n_fail == 0);
else
    no_plausibility_fail = true;
end
status.plausibility_ok = no_plausibility_fail;
status.plausibility_warning_ok = status.warning_fraction <= thresholds.max_warning_fraction;

% Acceptance must cover every metric the governed RMSE is computed over, not
% only the five selected primaries. A run that passes 5 of 5 while failing a
% graded metric at 24% is not an accepted candidate.
[status.governed_gate_ok, status.n_governed_fail, status.governed_gate_total, ...
    status.governed_gate_failures] = evaluate_governed_gate(report);

% Chi-squared goodness-of-fit, recorded for visibility only (PRD
% reyna_statistical_calibration_v1 Phase 2). NOT used to gate ACCEPT in this
% phase: changing the objective (Phase 1) and the acceptance rule in the same
% change would make neither attributable. Observe chi2/N across runs first.
[status.chi2_per_obs, status.chi2_interpretation] = extract_chi2_summary(report);

if status.fit_ok && status.governed_gate_ok && status.plausibility_ok && ...
        status.rmse_not_worse
    status.label = 'ACCEPT';
elseif is_promising_near_miss(status, thresholds)
    status.label = 'PROMISING_NEAR_MISS';
elseif status.fit_ok && ~status.plausibility_ok
    status.label = 'OUTPUT_FIT_ONLY';
elseif ~status.fit_ok && status.plausibility_ok
    status.label = 'PHYSIOLOGICAL_BUT_POOR_FIT';
else
    status.label = 'REJECT';
end

status.summary = sprintf(['%s | excellent5_fail=%d | primary10_fail=%d | ', ...
    'governed_gate=%d/%d | secondary15_fail=%d | plausibility_warning=%d | ', ...
    'plausibility_fail=%d | RMSE_improvement=%.1f%%'], ...
    status.label, status.n_primary_excellent_fail, status.n_primary_fail, ...
    status.governed_gate_total - status.n_governed_fail, status.governed_gate_total, ...
    status.n_secondary_fail, status.n_warning, status.n_fail, ...
    100 * status.rmse_improvement_frac);
if status.n_governed_fail > 0
    status.summary = sprintf('%s | governed_gate_failures=%s', ...
        status.summary, strjoin(status.governed_gate_failures, ','));
end
if isfinite(status.chi2_per_obs)
    status.summary = sprintf('%s | chi2_per_obs=%.2f (%s)', ...
        status.summary, status.chi2_per_obs, status.chi2_interpretation);
end
end

% =========================================================================
function [chi2_per_obs, interpretation] = extract_chi2_summary(report)
chi2_per_obs = NaN;
interpretation = 'unavailable';
if ~isstruct(report) || ~isfield(report, 'chi_squared') || ...
        ~isstruct(report.chi_squared)
    return;
end
chi2_per_obs = report.chi_squared.chi2_per_obs;
interpretation = report.chi_squared.interpretation;
end

% =========================================================================
function [gate_ok, n_fail, n_total, failures] = evaluate_governed_gate(report)
% EVALUATE_GOVERNED_GATE - acceptance across the whole governed RMSE mask.
%
% Backward compatible: a report without a full-metric gate (older callers,
% or a bare validation_report invocation) neither passes nor blocks, so the
% historical five-metric behaviour is preserved.
gate_ok = true;
n_fail = 0;
n_total = 0;
failures = {};

if ~isstruct(report) || ~isfield(report, 'full_metric_gate') || ...
        ~isstruct(report.full_metric_gate) || ...
        ~isfield(report.full_metric_gate, 'table') || ...
        isempty(report.full_metric_gate.table)
    return;
end

tbl = report.full_metric_gate.table;
governed = logical(tbl.InPrimaryRMSE);
n_total = nnz(governed);
if n_total == 0
    return;
end

failing = governed & ~logical(tbl.WithinGate);
n_fail = nnz(failing);
failures = tbl.Metric(failing)';
gate_ok = (n_fail == 0);
end

function [fit_ok, n_fail] = evaluate_primary_fit(report, threshold_pct)
fit_ok = false;
n_fail = 0;
if ~isstruct(report) || ~isfield(report, 'primary_gate') || isempty(report.primary_gate)
    return;
end

gate_tbl = report.primary_gate;
valid_rows = ~isnan(gate_tbl.AbsError_pct);
if ~any(valid_rows)
    return;
end
n_fail = sum(abs(gate_tbl.AbsError_pct(valid_rows)) > threshold_pct);
fit_ok = (n_fail == 0);
end

function [fit_ok, n_fail] = evaluate_secondary_fit(report, threshold_pct)
fit_ok = true;
n_fail = 0;
if ~isstruct(report) || ~isfield(report, 'table_cal') || isempty(report.table_cal)
    return;
end
if ~isfield(report, 'primary_metrics')
    primary_metrics = {};
else
    primary_metrics = report.primary_metrics(:)';
end
if ~isfield(report, 'validation_holdout_metrics')
    holdout_metrics = {};
else
    holdout_metrics = report.validation_holdout_metrics(:)';
end

tbl = report.table_cal;
is_secondary = ~ismember(tbl.Metric, primary_metrics) & ...
    ~ismember(tbl.Metric, holdout_metrics);
if ismember('IncludedInCalibration', tbl.Properties.VariableNames)
    is_secondary = is_secondary & tbl.IncludedInCalibration;
end
valid_rows = is_secondary & ~isnan(tbl.Error_pct);
if ~any(valid_rows)
    return;
end

n_fail = sum(abs(tbl.Error_pct(valid_rows)) > threshold_pct);
fit_ok = (n_fail == 0);
end

function thresholds = default_thresholds()
thresholds = struct( ...
    'excellent_error_pct', 5, ...
    'primary_error_pct', 10, ...
    'secondary_error_pct', 15, ...
    'max_warning_fraction', 0.20, ...
    'reject_on_any_fail', true, ...
    'accept_rmse_improvement_frac', 0.20, ...
    'promising_rmse_improvement_frac', 0.20, ...
    'promising_max_primary_fail', 3, ...
    'promising_max_secondary_fail', 6, ...
    'promising_max_primary_error_pct', 15, ...
    'promising_max_secondary_error_pct', 25);
end

function merged = merge_thresholds(defaults, overrides)
merged = defaults;
fields = fieldnames(overrides);
for idx = 1:numel(fields)
    merged.(fields{idx}) = overrides.(fields{idx});
end
end

function tf = is_promising_near_miss(status, thresholds)
tf = false;
if status.fit_ok
    return;
end

no_plausibility_fail = (status.n_fail == 0);

tf = no_plausibility_fail && ...
    status.rmse_improvement_frac >= thresholds.promising_rmse_improvement_frac && ...
    status.n_primary_fail <= thresholds.promising_max_primary_fail && ...
    status.n_secondary_fail <= thresholds.promising_max_secondary_fail && ...
    status.max_primary_abs_error <= thresholds.promising_max_primary_error_pct && ...
    status.max_secondary_abs_error <= thresholds.promising_max_secondary_error_pct;
end

function max_abs_error = compute_max_abs_error(report, gate_type)
max_abs_error = Inf;
if ~isstruct(report)
    return;
end

switch gate_type
    case 'primary'
        if ~isfield(report, 'primary_gate') || isempty(report.primary_gate)
            return;
        end
        tbl = report.primary_gate;
        if ~ismember('AbsError_pct', tbl.Properties.VariableNames)
            return;
        end
        values = abs(tbl.AbsError_pct(~isnan(tbl.AbsError_pct)));
    otherwise
        if ~isfield(report, 'table_cal') || isempty(report.table_cal)
            return;
        end
        tbl = report.table_cal;
        if ~ismember('Error_pct', tbl.Properties.VariableNames)
            return;
        end
        if ~isfield(report, 'primary_metrics')
            primary_metrics = {};
        else
            primary_metrics = report.primary_metrics(:)';
        end
        is_secondary = ~ismember(tbl.Metric, primary_metrics);
        if ismember('IncludedInCalibration', tbl.Properties.VariableNames)
            is_secondary = is_secondary & tbl.IncludedInCalibration;
        end
        values = abs(tbl.Error_pct(is_secondary & ~isnan(tbl.Error_pct)));
end

if isempty(values)
    return;
end
max_abs_error = max(values);
end

function improvement_frac = compute_rmse_improvement(report)
improvement_frac = -Inf;
if ~isstruct(report) || ~isfield(report, 'rmse_baseline') || ~isfield(report, 'rmse_cal')
    return;
end
if ~isfinite(report.rmse_baseline) || ~isfinite(report.rmse_cal) || report.rmse_baseline <= 0
    return;
end
improvement_frac = (report.rmse_baseline - report.rmse_cal) / report.rmse_baseline;
end

function tf = compute_rmse_not_worse(report)
tf = false;
if ~isstruct(report) || ~isfield(report, 'rmse_baseline') || ~isfield(report, 'rmse_cal')
    return;
end
if ~isfinite(report.rmse_baseline) || ~isfinite(report.rmse_cal)
    return;
end
tf = report.rmse_cal <= report.rmse_baseline + 1e-9;
end
