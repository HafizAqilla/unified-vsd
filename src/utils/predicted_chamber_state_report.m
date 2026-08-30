function report = predicted_chamber_state_report(metrics, clinical, scenario, excluded_evidence, results_dir)
% PREDICTED_CHAMBER_STATE_REPORT
% -----------------------------------------------------------------------
% Reports model-predicted chamber volumes and ejection fractions as
% PREDICTIONS rather than fitted targets, screens them against broad
% paediatric plausibility ranges, and cross-checks their direction against
% post-operative echo evidence where that evidence exists.
%
% Motivation
% ----------
% Reyna's chamber block is H+1 post-operative echo, so it is not a valid
% pre-operative comparator and has been removed as a calibration target.
% That makes the pre-operative chamber state a genuine model output with no
% comparator, which raises the reviewer question directly: what are the
% pre-procedure volumes, and would they be higher or lower than the measured
% post-closure ones?
%
% Physiology used for the direction check (pre-operative, restrictive
% left-to-right VSD, versus the same patient after closure):
%   - LVEDV: expected HIGHER pre-operatively. The shunt raises pulmonary
%     blood flow, so pulmonary venous return volume-loads the left ventricle;
%     closure removes that load.
%   - LVEF: expected HIGHER OR EQUAL pre-operatively. Part of LV ejection
%     goes into the lower-resistance right ventricle and pulmonary bed, so
%     effective afterload is reduced; closure restores full systemic
%     afterload.
%   - RVEDV / RVESV: NOT directionally specified. In a VSD the shunt is
%     ejected in systole largely into the pulmonary artery, so the right
%     ventricle is typically pressure-loaded rather than volume-loaded, and
%     the direction depends on defect size and pulmonary vascular state.
%     Reported without a directional expectation.
%
% INPUTS:
%   metrics           - struct from compute_clinical_indices (calibrated)  [-]
%   clinical          - unified clinical struct                            [-]
%   scenario          - 'pre_surgery' | 'post_surgery'                     [-]
%   excluded_evidence - recipe.excluded_evidence struct, or []             [-]
%   results_dir       - output directory for the exported CSV        [char]
%
% OUTPUTS:
%   report - struct with:
%       .table          per-chamber prediction, range screen, direction check
%       .n_within_range count inside the paediatric soft range      [count]
%       .n_total        chambers reported                           [count]
%       .direction_ok   all specified directions consistent       [logical]
%       .csv_path       written file path, empty when not exported   [char]
%
% ASSUMPTIONS:
%   - Reference ranges are broad screening bands, not normative paediatric
%     percentiles; passing them is a plausibility check, not validation.
%   - The post-operative comparison is a DIRECTION check only. It is not an
%     accuracy claim, because the two values describe different loading
%     states of the same patient.
%
% REFERENCES:
%   [1] src/utils/clinical_reference_ranges.m
%   [2] docs/reyna_zhang_fullmetric_results_20260828.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 4
    excluded_evidence = [];
end
if nargin < 5
    results_dir = '';
end

report = struct('table', table(), 'n_within_range', 0, 'n_total', 0, ...
    'direction_ok', true, 'csv_path', '');

chamber_metrics = {'LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVEF', 'RVEF'};
evidence_fields = {'LVEDV_mL', 'LVESV_mL', 'RVEDV_mL', 'RVESV_mL', 'LVEF', 'RVEF'};
expected_dir = {'higher', 'unspecified', 'unspecified', 'unspecified', ...
    'higher_or_equal', 'unspecified'};

ranges = clinical_reference_ranges(scenario, clinical);

metric_col = {};
unit_col = {};
model_col = [];
soft_low_col = [];
soft_high_col = [];
within_col = logical([]);
ref_col = [];
delta_col = [];
expected_col = {};
consistent_col = {};

for idx = 1:numel(chamber_metrics)
    name = chamber_metrics{idx};
    if ~isfield(metrics, name) || ~isfinite(metrics.(name))
        continue;
    end

    [soft_low, soft_high, unit] = lookup_range(ranges, name);
    model_value = metrics.(name);

    ref_value = NaN;
    if isstruct(excluded_evidence) && isfield(excluded_evidence, evidence_fields{idx})
        ref_value = excluded_evidence.(evidence_fields{idx});
    end

    delta_pct = NaN;
    if isfinite(ref_value) && abs(ref_value) > 1e-9
        delta_pct = 100 * (model_value - ref_value) / abs(ref_value);
    end

    [consistency, is_bad] = direction_consistency(expected_dir{idx}, delta_pct);
    if is_bad
        report.direction_ok = false;
    end

    metric_col{end+1, 1} = name;                 %#ok<AGROW>
    unit_col{end+1, 1} = unit;                   %#ok<AGROW>
    model_col(end+1, 1) = model_value;           %#ok<AGROW>
    soft_low_col(end+1, 1) = soft_low;           %#ok<AGROW>
    soft_high_col(end+1, 1) = soft_high;         %#ok<AGROW>
    within_col(end+1, 1) = isfinite(soft_low) && isfinite(soft_high) && ...
        model_value >= soft_low && model_value <= soft_high;  %#ok<AGROW>
    ref_col(end+1, 1) = ref_value;               %#ok<AGROW>
    delta_col(end+1, 1) = delta_pct;             %#ok<AGROW>
    expected_col{end+1, 1} = expected_dir{idx};  %#ok<AGROW>
    consistent_col{end+1, 1} = consistency;      %#ok<AGROW>
end

if isempty(metric_col)
    return;
end

report.table = table(metric_col, unit_col, model_col, soft_low_col, ...
    soft_high_col, within_col, ref_col, delta_col, expected_col, consistent_col, ...
    'VariableNames', {'Metric','Unit','ModelPredicted','PaedSoftLow', ...
    'PaedSoftHigh','WithinPaedRange','PostOpReference','Delta_vs_PostOp_pct', ...
    'ExpectedDirection','DirectionCheck'});
report.n_total = height(report.table);
report.n_within_range = nnz(within_col);

if ~isempty(results_dir)
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    report.csv_path = fullfile(results_dir, ...
        sprintf('predicted_chamber_state_%s.csv', scenario));
    writetable(report.table, report.csv_path);
end

print_report(report, scenario, excluded_evidence);
end

% =========================================================================
function [soft_low, soft_high, unit] = lookup_range(ranges, metric_name)
soft_low = NaN;
soft_high = NaN;
unit = '';
if isempty(ranges)
    return;
end
idx = find(strcmp(ranges.Metric, metric_name), 1, 'first');
if isempty(idx)
    return;
end
soft_low = ranges.SoftLow(idx);
soft_high = ranges.SoftHigh(idx);
unit = ranges.Unit{idx};
end

% =========================================================================
function [label, is_bad] = direction_consistency(expected, delta_pct)
% DIRECTION_CONSISTENCY - compare predicted direction with physiology.
% A small band around zero counts as "equal" so noise is not read as a
% direction violation.
tol_pct = 2.0;
is_bad = false;

if strcmp(expected, 'unspecified')
    label = 'not_specified';
    return;
end
if ~isfinite(delta_pct)
    label = 'no_reference';
    return;
end

switch expected
    case 'higher'
        if delta_pct > tol_pct
            label = 'consistent';
        else
            label = 'INCONSISTENT';
            is_bad = true;
        end
    case 'higher_or_equal'
        if delta_pct > -tol_pct
            label = 'consistent';
        else
            label = 'INCONSISTENT';
            is_bad = true;
        end
    otherwise
        label = 'not_specified';
end
end

% =========================================================================
function print_report(report, scenario, excluded_evidence)
fprintf('\n--- PREDICTED CHAMBER STATE — %s ---\n', ...
    upper(strrep(scenario, '_', ' ')));
fprintf(['  These are model PREDICTIONS with no valid comparator in this ', ...
    'scenario.\n  They are screened for plausibility, not validated.\n\n']);
disp(report.table);

fprintf('  Within paediatric screening range: %d of %d\n', ...
    report.n_within_range, report.n_total);

if isstruct(excluded_evidence) && isfield(excluded_evidence, 'timing')
    fprintf(['  Reference column is %s evidence, retained for DIRECTION ', ...
        'only.\n  It is not an accuracy target: the two values describe ', ...
        'different loading states.\n'], excluded_evidence.timing);
end

if report.direction_ok
    fprintf('  [PASS] Every specified direction matches shunt physiology.\n');
else
    bad = report.table.Metric(strcmp(report.table.DirectionCheck, 'INCONSISTENT'));
    for idx = 1:numel(bad)
        fprintf(2, ['  [DIRECTION] %s moves against the expected ', ...
            'pre-operative shunt loading.\n'], bad{idx});
    end
end
end
