function findings = assert_evidence_timing_governance(clinical, scenario, target_tiers, recipe)
% ASSERT_EVIDENCE_TIMING_GOVERNANCE
% -----------------------------------------------------------------------
% Blocks evidence recorded at one surgical timepoint from silently fitting a
% scenario at another.
%
% The Reyna LV/RV volume block is H+1 post-operative echo. patient_reyna()
% marks it unavailable for pre-surgery fitting, while the calibration recipe
% re-injects it. Four of those rows are demoted to consistency-check-only,
% but any row left at hard/soft tier is actively fitted, and override_IC
% additionally seeds chamber V0/E from the same block. This check makes both
% channels explicit and refuses the cross-timing case unless the recipe opts
% in deliberately.
%
% INPUTS:
%   clinical     - unified clinical struct                               [-]
%   scenario     - 'pre_surgery' | 'post_surgery'                        [-]
%   target_tiers - build_target_tiers output                             [-]
%   recipe       - calibration recipe struct (may be empty)              [-]
%
% OUTPUTS:
%   findings     - struct with:
%       .table          per-metric evidence timing and tier
%       .violations     metric names fitted across timepoints    [cellstr]
%       .allowed        true when the recipe opted in            [logical]
%
% ERRORS:
%   assert_evidence_timing_governance:crossTimingEvidence when a metric whose
%   evidence timing disagrees with the scenario carries a fitted tier and the
%   recipe has not set allow_cross_timing_evidence.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 5)
%   [2] config/patient_reyna.m lines 81-91
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------

findings = struct('table', table(), 'violations', {{}}, 'allowed', false);

if nargin < 4 || ~isstruct(recipe) || ~isfield(recipe, 'evidence_timing') || ...
        ~isstruct(recipe.evidence_timing)
    return;
end

findings.allowed = isfield(recipe, 'allow_cross_timing_evidence') && ...
    isequal(recipe.allow_cross_timing_evidence, true);

targets = get_calibration_targets(scenario, clinical);
scenario_timing = scenario_expected_timing(scenario);

metric_col = {};
field_col = {};
timing_col = {};
tier_col = {};
fitted_col = logical([]);

for idx = 1:numel(targets)
    field_name = targets(idx).ClinicalField;
    if ~isfield(recipe.evidence_timing, field_name)
        continue;
    end
    timing = char(recipe.evidence_timing.(field_name));
    tier = tier_for_metric(target_tiers, targets(idx).Metric);
    fitted = ismember(tier, {'hard', 'soft'});

    metric_col{end+1, 1} = targets(idx).Metric; %#ok<AGROW>
    field_col{end+1, 1} = field_name;           %#ok<AGROW>
    timing_col{end+1, 1} = timing;              %#ok<AGROW>
    tier_col{end+1, 1} = tier;                  %#ok<AGROW>
    fitted_col(end+1, 1) = fitted;              %#ok<AGROW>
end

if isempty(metric_col)
    return;
end

cross_timing = ~strcmp(timing_col, scenario_timing) & ...
    ~strcmp(timing_col, 'derived');
findings.table = table(metric_col, field_col, timing_col, tier_col, ...
    fitted_col, cross_timing, ...
    'VariableNames', {'Metric','ClinicalField','EvidenceTiming','Tier', ...
    'Fitted','CrossTiming'});

violation_mask = cross_timing & fitted_col;
findings.violations = metric_col(violation_mask)';

if ~any(violation_mask) || findings.allowed
    return;
end

error('assert_evidence_timing_governance:crossTimingEvidence', ...
    ['Scenario "%s" fits targets whose evidence was recorded at a ', ...
     'different timepoint: %s. Demote them to consistency_only/', ...
     'validation_holdout, or set recipe.allow_cross_timing_evidence = true ', ...
     'and state the limitation in the results report.'], ...
    scenario, strjoin(findings.violations, ', '));
end

% =========================================================================
function timing = scenario_expected_timing(scenario)
switch char(scenario)
    case 'pre_surgery'
        timing = 'pre_operative';
    case 'post_surgery'
        timing = 'post_operative_H1';
    otherwise
        timing = 'unknown';
end
end

function tier = tier_for_metric(target_tiers, metric_name)
tier = 'unavailable';
if ~isstruct(target_tiers) || ~isfield(target_tiers, 'table') || ...
        isempty(target_tiers.table)
    return;
end
idx = find(strcmp(target_tiers.table.Metric, metric_name), 1, 'first');
if ~isempty(idx)
    tier = target_tiers.table.Tier{idx};
end
end
