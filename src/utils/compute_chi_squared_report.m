function report = compute_chi_squared_report(gate_table, n_active_parameters)
% COMPUTE_CHI_SQUARED_REPORT
% -----------------------------------------------------------------------
% Reduces the governed-metric gate table to a chi-squared goodness-of-fit
% statistic that cannot be improved by redefining the reported denominator.
%
% Motivation
% ----------
% A percentage acceptance gate treats a pressure recorded with three
% identical repeats the same as one known only to +-20%. It is not a
% statistical criterion. Chi-squared, computed over residuals normalised by
% each metric's own declared measurement uncertainty, is: it answers whether
% the fit is as good as the data allows (chi2/N ~ 1), worse than the data
% supports (chi2/N >> 1), or fitting below the noise floor (chi2/N << 1).
%
% This is the standard discrepancy-principle criterion from ill-posed
% inverse-problem theory (Morozov): stop fitting once residuals match the
% measurement noise; fitting tighter than that is fitting noise, not
% physics.
%
% INPUTS:
%   gate_table          - the .table field from export_full_metric_gate,
%                          must carry Sigma, ZScore, ZScoreSquared,
%                          InPrimaryRMSE                                 [-]
%   n_active_parameters - count of fitted (active) parameters          [int]
%
% OUTPUTS:
%   report - struct with:
%       .chi2            sum of z^2 over governed rows                   [-]
%       .n_obs           governed row count                          [count]
%       .n_skipped       governed rows excluded for non-finite z    [count]
%       .n_parameters    echoed input                                [count]
%       .dof             max(n_obs - n_parameters, 0)                [count]
%       .chi2_per_obs    chi2 / n_obs                                    [-]
%       .chi2_reduced    chi2 / dof, NaN when dof <= 0                   [-]
%       .interpretation  'underfit' | 'consistent' | 'overfit'        [char]
%       .dof_note        'insufficient_dof' when dof <= 2, else
%                        'sufficient_dof' (never '', to keep the struct a
%                        valid single-row table for struct2table)       [char]
%       .worst_metric    metric name with largest |z|, or 'none' when
%                        unavailable (never '')                         [char]
%       .worst_z         that z value                                    [-]
%
% ASSUMPTIONS:
%   - Only rows with InPrimaryRMSE = true are counted: this is the governed
%     acceptance set, consistent with export_full_metric_gate.
%   - A reduced chi-squared computed at dof <= 2 is not statistically
%     stable; report it, but flag it via dof_note rather than presenting it
%     alone.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 2)
%   [2] Morozov's discrepancy principle (ill-posed inverse problems)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------

report = empty_report(n_active_parameters);

if nargin < 1 || isempty(gate_table) || ~istable(gate_table)
    return;
end
required_cols = {'ZScore', 'ZScoreSquared', 'InPrimaryRMSE', 'Metric'};
if ~all(ismember(required_cols, gate_table.Properties.VariableNames))
    return;
end

governed = logical(gate_table.InPrimaryRMSE);
governed_tbl = gate_table(governed, :);

finite_mask = isfinite(governed_tbl.ZScoreSquared);
report.n_skipped = nnz(~finite_mask);
scored_tbl = governed_tbl(finite_mask, :);

report.n_obs = height(scored_tbl);
if report.n_obs == 0
    return;
end

report.chi2 = sum(scored_tbl.ZScoreSquared);
report.dof = max(report.n_obs - report.n_parameters, 0);
report.chi2_per_obs = report.chi2 / report.n_obs;
if report.dof > 0
    report.chi2_reduced = report.chi2 / report.dof;
else
    report.chi2_reduced = NaN;
end
if report.dof <= 2
    report.dof_note = 'insufficient_dof';
end

report.interpretation = classify_chi2_per_obs(report.chi2_per_obs);

[max_z2, worst_ix] = max(scored_tbl.ZScoreSquared);
report.worst_metric = scored_tbl.Metric{worst_ix};
report.worst_z = sign(scored_tbl.ZScore(worst_ix)) * sqrt(max_z2);

print_chi_squared_report(report);
end

% =========================================================================
function report = empty_report(n_active_parameters)
if nargin < 1 || isempty(n_active_parameters) || ~isfinite(n_active_parameters)
    n_active_parameters = 0;
end
report = struct( ...
    'chi2', NaN, ...
    'n_obs', 0, ...
    'n_skipped', 0, ...
    'n_parameters', round(n_active_parameters), ...
    'dof', 0, ...
    'chi2_per_obs', NaN, ...
    'chi2_reduced', NaN, ...
    'interpretation', 'unavailable', ...
    'dof_note', 'sufficient_dof', ...
    'worst_metric', 'none', ...
    'worst_z', NaN);
% dof_note and worst_metric are never left as '' (0x0 char): struct2table
% (used by export_chi_squared_summary in validation_report.m) treats an
% empty-char field as having 0 rows against every other field's 1, and
% throws "different numbers of rows" -- this crashed a completed 6-start
% calibration run before its results could be saved (2026-08-29). Every
% char field must carry a non-empty placeholder so the struct is always a
% valid single-row table.
end

% =========================================================================
function label = classify_chi2_per_obs(chi2_per_obs)
% CLASSIFY_CHI2_PER_OBS - discrepancy-principle bands.
%   > 2.0        underfit   (model or data inconsistent)
%   0.5 - 2.0    consistent (residuals match measurement noise)
%   < 0.5        overfit    (fitting below the noise floor)
if ~isfinite(chi2_per_obs)
    label = 'unavailable';
elseif chi2_per_obs > 2.0
    label = 'underfit';
elseif chi2_per_obs < 0.5
    label = 'overfit';
else
    label = 'consistent';
end
end

% =========================================================================
function print_chi_squared_report(report)
fprintf('\n--- CHI-SQUARED ---\n');
fprintf('  N (governed observations) : %d', report.n_obs);
if report.n_skipped > 0
    fprintf('  (%d skipped: non-finite z)', report.n_skipped);
end
fprintf('\n');
fprintf('  p (active parameters)     : %d\n', report.n_parameters);
fprintf('  dof = N - p               : %d', report.dof);
if strcmp(report.dof_note, 'insufficient_dof')
    fprintf('  [dof <= 2: reduced chi2 is not statistically stable]');
end
fprintf('\n');
fprintf('  chi2                      : %.3f\n', report.chi2);
fprintf('  chi2 / N                  : %.3f\n', report.chi2_per_obs);
if isfinite(report.chi2_reduced)
    fprintf('  chi2 / dof (reduced)      : %.3f\n', report.chi2_reduced);
else
    fprintf('  chi2 / dof (reduced)      : n/a (dof <= 0)\n');
end
fprintf('  interpretation             : %s\n', report.interpretation);
switch report.interpretation
    case 'underfit'
        fprintf(['  [chi2/N > 2.0] Residuals exceed measurement noise: ', ...
            'model or data inconsistent.\n']);
    case 'overfit'
        fprintf(['  [chi2/N < 0.5] Fitting below the noise floor: ', ...
            'overfitting, or declared uncertainties are too generous.\n']);
    case 'consistent'
        fprintf(['  [0.5 <= chi2/N <= 2.0] Residuals are consistent with ', ...
            'declared measurement noise.\n']);
end
if ~isempty(report.worst_metric)
    fprintf('  Worst metric (by |z|)     : %s (z = %+.2f)\n', ...
        report.worst_metric, report.worst_z);
end
end
