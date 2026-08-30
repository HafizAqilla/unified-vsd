function gate = export_full_metric_gate(report, scenario, results_dir, gate_pct)
% EXPORT_FULL_METRIC_GATE
% -----------------------------------------------------------------------
% Builds and exports the per-metric acceptance table across every clinical
% validation target, not only the selected primary metrics.
%
% The governed primary RMSE and the five-metric primary gate both describe
% subsets of the available clinical evidence. This export makes the whole
% distribution visible so an acceptance claim can be stated as "n of N"
% against a disclosed denominator.
%
% INPUTS:
%   report      - struct from validation_report                          [-]
%   scenario    - 'pre_surgery' | 'post_surgery'                         [-]
%   results_dir - output directory for the exported CSV [char/string]    [-]
%   gate_pct    - optional patient-acceptance band, default 10           [%]
%
% OUTPUTS:
%   gate        - struct with:
%       .table            per-metric table (all finite-comparator rows)
%       .n_total          rows with a finite clinical comparator      [count]
%       .n_within_gate    rows within the acceptance band             [count]
%       .n_primary_total  rows inside the governed primary RMSE mask  [count]
%       .n_primary_within rows inside the mask and within the band    [count]
%       .worst_metric     metric name with the largest absolute error [char]
%       .worst_abs_pct    that absolute error                            [%]
%       .csv_path         written file path, empty when not exported  [char]
%
% ASSUMPTIONS:
%   - Rows without a finite clinical comparator are excluded; they cannot
%     be scored against a patient measurement.
%   - Objective membership is read from the target-tier governance table,
%     so a metric that is graded but never fitted is visible as such.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 0)
%   [2] docs/calibration_data_governance_notes.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3
    results_dir = '';
end
if nargin < 4 || isempty(gate_pct) || ~isfinite(gate_pct)
    gate_pct = 10;
end

gate = empty_gate(gate_pct);

src_tbl = select_source_table(report);
if isempty(src_tbl)
    return;
end

value_col = 'Calibrated';
if ~ismember(value_col, src_tbl.Properties.VariableNames)
    value_col = 'Baseline';
end

mask = isfinite(src_tbl.Clinical);
rows = src_tbl(mask, :);
if isempty(rows)
    return;
end

abs_err = abs(rows.Error_pct);
within_gate = ~isnan(abs_err) & abs_err <= gate_pct;
within_excellent = ~isnan(abs_err) & abs_err <= 0.5 * gate_pct;

in_objective = logical(rows.IncludedInCalibration);
in_primary = logical(rows.IncludedInPrimaryRMSE);

sigma = resolve_row_sigma(rows);
z_score = (rows.(value_col) - rows.Clinical) ./ sigma;
z_score_sq = z_score .^ 2;

gate_tbl = table( ...
    rows.Metric, rows.Unit, rows.Tier, rows.Clinical, rows.(value_col), ...
    rows.Error_pct, abs_err, sigma, z_score, z_score_sq, in_objective, in_primary, ...
    within_gate, within_excellent, rows.Flag, ...
    'VariableNames', {'Metric','Unit','Tier','Clinical','Model', ...
    'Error_pct','AbsError_pct','Sigma','ZScore','ZScoreSquared', ...
    'InObjective','InPrimaryRMSE','WithinGate','WithinExcellent','Flag'});

% Sorted worst-first by |z|, the statistically meaningful ranking (PRD
% reyna_statistical_calibration_v1 Phase 2): a metric known to +-0.5 units
% that misses by 1 unit is a worse fit than one known to +-10 units that
% misses by 5, even though both may show the same or a smaller percent error.
gate_tbl = sortrows(gate_tbl, 'ZScoreSquared', 'descend', ...
    'MissingPlacement', 'last');

gate.table = gate_tbl;
gate.gate_pct = gate_pct;
gate.n_total = height(gate_tbl);
gate.n_within_gate = nnz(gate_tbl.WithinGate);
gate.n_within_excellent = nnz(gate_tbl.WithinExcellent);
gate.n_primary_total = nnz(gate_tbl.InPrimaryRMSE);
gate.n_primary_within = nnz(gate_tbl.InPrimaryRMSE & gate_tbl.WithinGate);
gate.n_graded_not_fitted = nnz(gate_tbl.InPrimaryRMSE & ~gate_tbl.InObjective);
gate.outside_gate_metrics = gate_tbl.Metric(~gate_tbl.WithinGate)';

finite_err = gate_tbl.AbsError_pct(~isnan(gate_tbl.AbsError_pct));
if ~isempty(finite_err)
    [gate.worst_abs_pct, worst_ix] = max(gate_tbl.AbsError_pct);
    gate.worst_metric = gate_tbl.Metric{worst_ix};
end

if ~isempty(results_dir)
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    gate.csv_path = fullfile(results_dir, ...
        sprintf('full_metric_gate_%s.csv', scenario));
    writetable(gate_tbl, gate.csv_path);
end

print_full_metric_gate(gate, scenario);

end

% =========================================================================
function sigma = resolve_row_sigma(rows)
% RESOLVE_ROW_SIGMA - measurement uncertainty per row, same resolution order
% as build_target_sigma_map in build_case_calibration_profile.m:
% UncertaintyAbs, then UncertaintyFraction*|Clinical|, then a 10% fallback.
% Kept independent of that function (rather than calling it) because this
% export works directly off the validation_report table, which already
% carries UncertaintyAbs/UncertaintyFraction per row.
n = height(rows);
sigma = nan(n, 1);

has_abs_col = ismember('UncertaintyAbs', rows.Properties.VariableNames);
has_frac_col = ismember('UncertaintyFraction', rows.Properties.VariableNames);

for idx = 1:n
    abs_sigma = NaN;
    if has_abs_col
        abs_sigma = rows.UncertaintyAbs(idx);
    end
    frac_sigma = NaN;
    if has_frac_col && isfinite(rows.UncertaintyFraction(idx)) && ...
            rows.UncertaintyFraction(idx) > 0
        frac_sigma = abs(rows.Clinical(idx)) * rows.UncertaintyFraction(idx);
    end

    if isfinite(abs_sigma) && abs_sigma > 0
        sigma(idx) = abs_sigma;
    elseif isfinite(frac_sigma) && frac_sigma > 0
        sigma(idx) = frac_sigma;
    else
        sigma(idx) = 0.10 * abs(rows.Clinical(idx));
    end
end
sigma = max(sigma, 1e-9);
end

% =========================================================================
function gate = empty_gate(gate_pct)
gate = struct( ...
    'table', table(), ...
    'gate_pct', gate_pct, ...
    'n_total', 0, ...
    'n_within_gate', 0, ...
    'n_within_excellent', 0, ...
    'n_primary_total', 0, ...
    'n_primary_within', 0, ...
    'n_graded_not_fitted', 0, ...
    'outside_gate_metrics', {{}}, ...
    'worst_metric', '', ...
    'worst_abs_pct', NaN, ...
    'csv_path', '');
end

% =========================================================================
function src_tbl = select_source_table(report)
% SELECT_SOURCE_TABLE - prefer the calibrated comparison table.
src_tbl = table();
if ~isstruct(report)
    return;
end
if isfield(report, 'table_cal') && ~isempty(report.table_cal)
    src_tbl = report.table_cal;
elseif isfield(report, 'table_baseline') && ~isempty(report.table_baseline)
    src_tbl = report.table_baseline;
end
end

% =========================================================================
function print_full_metric_gate(gate, scenario)
% PRINT_FULL_METRIC_GATE - console summary with a disclosed denominator.
fprintf('\n--- FULL METRIC %g%% GATE — %s ---\n', ...
    gate.gate_pct, upper(strrep(scenario, '_', ' ')));
if gate.n_total == 0
    fprintf('  No clinical validation targets available.\n');
    return;
end

fprintf('  All clinical targets:      %d of %d within %g%%\n', ...
    gate.n_within_gate, gate.n_total, gate.gate_pct);
fprintf('  Governed primary RMSE set: %d of %d within %g%%\n', ...
    gate.n_primary_within, gate.n_primary_total, gate.gate_pct);
fprintf('  Within %g%% excellent band: %d of %d\n', ...
    0.5 * gate.gate_pct, gate.n_within_excellent, gate.n_total);

if gate.n_graded_not_fitted > 0
    fprintf(2, '  [GOVERNANCE] %d metric(s) counted in primary RMSE but absent from the objective.\n', ...
        gate.n_graded_not_fitted);
end

if isempty(gate.outside_gate_metrics)
    fprintf('  [PASS] Every clinical validation target is within %g%%.\n', gate.gate_pct);
else
    for idx = 1:numel(gate.outside_gate_metrics)
        name = gate.outside_gate_metrics{idx};
        row = find(strcmp(gate.table.Metric, name), 1, 'first');
        fprintf(2, '  [OUTSIDE] %-14s |error| = %6.2f%%  (tier %s)\n', ...
            name, gate.table.AbsError_pct(row), gate.table.Tier{row});
    end
end

if ~isempty(gate.worst_metric)
    fprintf('  Worst metric: %s (|error| = %.2f%%)\n', ...
        gate.worst_metric, gate.worst_abs_pct);
end
if ~isempty(gate.csv_path)
    fprintf('  Exported: %s\n', gate.csv_path);
end
end
