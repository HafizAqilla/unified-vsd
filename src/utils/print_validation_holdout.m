function holdout_report = print_validation_holdout(gate_table)
% PRINT_VALIDATION_HOLDOUT
% -----------------------------------------------------------------------
% Console block and summary struct for genuinely independent validation
% holdout metrics: excluded from both fitting and the governed RMSE/
% chi-squared, so a prediction against them is an out-of-sample test.
%
% Deliberately distinct from a derived/algebraic exclusion tier
% ('derived_validation'). A metric computed from other already-fitted
% targets (e.g. SVR = (SAP_mean - RAP_mean) / CO_Lmin) carries no
% information independent of what the objective already sees, so it
% cannot serve as a holdout no matter how well it lands. See
% docs/reyna_statistical_calibration_prd.md Phase 5 for why SVR is not
% designated a holdout despite the PRD's original suggestion.
%
% INPUTS:
%   gate_table - the .table field from export_full_metric_gate, must
%                carry a 'Tier' column (and, when present, 'WithinGate')  [-]
%
% OUTPUTS:
%   holdout_report - struct with:
%       .table          rows tiered 'validation_holdout'
%       .n_total        count of such rows
%       .n_within_gate  count of those within the acceptance band
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 5)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------
holdout_report = struct('table', table(), 'n_total', 0, 'n_within_gate', 0);

if ~istable(gate_table) || ~ismember('Tier', gate_table.Properties.VariableNames)
    return;
end

rows = gate_table(strcmp(gate_table.Tier, 'validation_holdout'), :);
holdout_report.table = rows;
holdout_report.n_total = height(rows);

fprintf('\n--- VALIDATION HOLDOUT ---\n');
if holdout_report.n_total == 0
    fprintf(['  No metric is currently designated validation_holdout. ', ...
        'Every finite target is either fitted or reported as derived/', ...
        'consistency-only evidence, not held out as an independent ', ...
        'prediction test.\n']);
    return;
end

if ismember('WithinGate', rows.Properties.VariableNames)
    holdout_report.n_within_gate = nnz(logical(rows.WithinGate));
end
disp(rows);
fprintf('  Held out of both the objective and the governed RMSE: %d of %d within the gate band.\n', ...
    holdout_report.n_within_gate, holdout_report.n_total);
end
