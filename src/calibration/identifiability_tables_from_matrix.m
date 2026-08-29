function out = identifiability_tables_from_matrix(S, param_names)
% IDENTIFIABILITY_TABLES_FROM_MATRIX
% -----------------------------------------------------------------------
% Pure math: condition number, column correlation, and per-parameter/
% per-pair tables from an already-built scaled sensitivity matrix S.
%
% Kept separate from analyse_parameter_identifiability.m (which drives the
% forward ODE model to build S) so this reduction step can be tested
% directly against synthetic matrices -- duplicated columns, orthogonal
% columns, a zero column -- without paying for a simulation.
%
% INPUTS:
%   S           - scaled sensitivity matrix, columns = parameters,
%                 rows = metrics. May contain NaN for a metric/parameter
%                 pair that failed to evaluate.               [n_m x n_p]
%   param_names - cellstr of parameter names, one per column of S  [1xn_p]
%
% OUTPUTS:
%   out - struct with:
%       .parameter_table    Parameter, ColumnNorm, Inactive,
%                            MaxAbsCorrelation, MostCorrelatedWith, Flagged
%       .pair_table          ParameterA, ParameterB, Correlation, Collinear
%       .condition_number    cond(S) over the columns with any finite entry
%
% A column that is entirely NaN (every perturbation for that parameter
% failed to evaluate) is dropped entirely, not zero-filled, so it cannot be
% falsely reported as "inactive": inactive means "evaluated but has no
% effect on any governed metric", which is a different finding from
% "could not be evaluated".
%
% "Inactive" (near-zero column norm) is a flag, not automatic removal: PRD
% reyna_statistical_calibration_v1 Phase 3 is report-only. A near-zero
% column here should already have been dropped by GSA screening; seeing one
% survive to this stage is itself a finding worth investigating.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 3)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------
out = struct('parameter_table', table(), 'pair_table', table(), ...
    'condition_number', NaN);

if isempty(S) || isempty(param_names)
    return;
end

valid_cols = ~all(isnan(S), 1);
S_valid = S(:, valid_cols);
S_valid(isnan(S_valid)) = 0;
valid_names = param_names(valid_cols);

if isempty(S_valid) || size(S_valid, 2) < 1
    return;
end

out.condition_number = safe_cond(S_valid);

column_norm = sqrt(sum(S_valid.^2, 1));
positive_norms = column_norm(column_norm > 0);
if isempty(positive_norms)
    inactive_threshold = 1e-6;
else
    inactive_threshold = 1e-6 * max(positive_norms);
end
is_inactive = column_norm <= inactive_threshold;

corr_mat = safe_column_correlation(S_valid);

n_valid = numel(valid_names);
max_abs_corr = nan(1, n_valid);
most_corr_with = repmat({''}, 1, n_valid);
flagged = false(1, n_valid);
pair_rows = {};
for a = 1:n_valid
    others = setdiff(1:n_valid, a);
    if isempty(others)
        continue;
    end
    [best_val, best_ix] = max(abs(corr_mat(a, others)));
    max_abs_corr(a) = best_val;
    most_corr_with{a} = valid_names{others(best_ix)};
    flagged(a) = best_val > 0.9;
end
for a = 1:n_valid
    for b = (a+1):n_valid
        pair_rows(end+1, :) = {valid_names{a}, valid_names{b}, ...
            corr_mat(a, b), abs(corr_mat(a, b)) > 0.9}; %#ok<AGROW>
    end
end

out.parameter_table = table(valid_names(:), column_norm(:), is_inactive(:), ...
    max_abs_corr(:), most_corr_with(:), flagged(:), ...
    'VariableNames', {'Parameter','ColumnNorm','Inactive', ...
    'MaxAbsCorrelation','MostCorrelatedWith','Flagged'});

if isempty(pair_rows)
    out.pair_table = table(cell(0,1), cell(0,1), zeros(0,1), false(0,1), ...
        'VariableNames', {'ParameterA','ParameterB','Correlation','Collinear'});
else
    out.pair_table = cell2table(pair_rows, ...
        'VariableNames', {'ParameterA','ParameterB','Correlation','Collinear'});
    out.pair_table = sortrows(out.pair_table, 'Correlation', 'descend', ...
        'ComparisonMethod', 'abs');
end
end

% =========================================================================
function c = safe_cond(S)
try
    c = cond(S);
catch
    c = NaN;
end
end

function corr_mat = safe_column_correlation(S)
% SAFE_COLUMN_CORRELATION - pairwise Pearson correlation between columns.
% A constant column (zero variance) would make corrcoef return NaN; guard
% it by returning zero correlation for such columns instead of propagating
% NaN into the flagging logic.
n = size(S, 2);
try
    corr_mat = corrcoef(S);
    corr_mat(isnan(corr_mat)) = 0;
catch
    corr_mat = eye(n);
end
end
