function summary_table = compare_run_metrics(reference_csv, candidate_csv)
% COMPARE_RUN_METRICS
% -----------------------------------------------------------------------
% Compare key calibrated metrics between two validation CSV exports.
%
% INPUTS:
%   reference_csv - path to a baseline/reference validation CSV file
%   candidate_csv - path to a candidate validation CSV file
%
% OUTPUTS:
%   summary_table - table of selected metrics and error deltas
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-02
% VERSION:  1.0
% -----------------------------------------------------------------------

reference_table = readtable(reference_csv);
candidate_table = readtable(candidate_csv);

key_metrics = { ...
    'RAP_mean'
    'PAP_min'
    'PAP_mean'
    'SAP_mean'
    'QpQs'
    'SVR'
    'CO_Lmin'
    'LVEDV'
    'LVESV'
    'LVEF'
    };

n_metrics = numel(key_metrics);
metric_col = cell(n_metrics, 1);
reference_err_col = nan(n_metrics, 1);
candidate_err_col = nan(n_metrics, 1);
delta_err_col = nan(n_metrics, 1);

for i = 1:n_metrics
    metric_name = key_metrics{i};
    metric_col{i} = metric_name;

    ref_idx = find(strcmp(reference_table.Metric, metric_name), 1, 'first');
    cand_idx = find(strcmp(candidate_table.Metric, metric_name), 1, 'first');

    if ~isempty(ref_idx)
        reference_err_col(i) = reference_table.Error_pct(ref_idx);
    end
    if ~isempty(cand_idx)
        candidate_err_col(i) = candidate_table.Error_pct(cand_idx);
    end

    if isfinite(reference_err_col(i)) && isfinite(candidate_err_col(i))
        delta_err_col(i) = abs(candidate_err_col(i)) - abs(reference_err_col(i));
    end
end

summary_table = table(metric_col, reference_err_col, candidate_err_col, delta_err_col, ...
    'VariableNames', {'Metric', 'ReferenceError_pct', 'CandidateError_pct', 'AbsoluteDelta_pct'});

disp(summary_table);
end
