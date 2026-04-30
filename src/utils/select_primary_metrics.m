function [primary_metrics, selection_table] = select_primary_metrics(clinical, gsa_out, scenario)
% SELECT_PRIMARY_METRICS
% -----------------------------------------------------------------------
% Selects five patient-specific primary calibration metrics.
%
% The selection keeps clinically mandatory paediatric VSD targets fixed when
% available, then fills the remaining slots using measurement reliability and
% front-end GSA support.
%
% INPUTS:
%   clinical         - unified clinical struct
%   gsa_out          - initial GSA output struct (may be empty)
%   scenario         - 'pre_surgery' or 'post_surgery'
%
% OUTPUTS:
%   primary_metrics  - cell array of selected model metric names
%   selection_table  - table documenting scores and reasons
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

targets = get_calibration_targets(scenario, clinical);
n_targets = numel(targets);

metric_col = cell(n_targets, 1);
available_col = false(n_targets, 1);
mandatory_col = false(n_targets, 1);
candidate_col = false(n_targets, 1);
reliability_col = cell(n_targets, 1);
gsa_score_col = zeros(n_targets, 1);
reliability_score_col = zeros(n_targets, 1);
total_score_col = zeros(n_targets, 1);
selected_col = false(n_targets, 1);

for i = 1:n_targets
    metric_col{i} = targets(i).Metric;
    available_col(i) = ~isnan(targets(i).ClinicalValue);
    mandatory_col(i) = targets(i).MandatoryPrimary;
    candidate_col(i) = targets(i).CandidatePrimary;
    reliability_col{i} = targets(i).Reliability;
    reliability_score_col(i) = reliability_score(targets(i).Reliability);
    gsa_score_col(i) = metric_gsa_score(gsa_out, targets(i).Metric);
end

if max(gsa_score_col) > 0
    gsa_score_col = gsa_score_col ./ max(gsa_score_col);
end

total_score_col = reliability_score_col + gsa_score_col;

eligible_mandatory = available_col & mandatory_col;
selected_col(eligible_mandatory) = true;

slots_remaining = max(0, 5 - nnz(selected_col));
eligible_candidates = available_col & candidate_col & ~selected_col;
candidate_idx = find(eligible_candidates);
[~, order] = sort(total_score_col(candidate_idx), 'descend');
candidate_idx = candidate_idx(order);
candidate_idx = candidate_idx(1:min(slots_remaining, numel(candidate_idx)));
selected_col(candidate_idx) = true;

% If sparse clinical data leave fewer than five metrics, fill with any
% calibrated available metric so each run still has a reproducible target set.
slots_remaining = max(0, 5 - nnz(selected_col));
if slots_remaining > 0
    fallback_idx = find(available_col & [targets.UseForCalibration]' & ~selected_col);
    [~, order] = sort(total_score_col(fallback_idx), 'descend');
    fallback_idx = fallback_idx(order);
    fallback_idx = fallback_idx(1:min(slots_remaining, numel(fallback_idx)));
    selected_col(fallback_idx) = true;
end

primary_metrics = metric_col(selected_col);
selection_table = table(metric_col, available_col, mandatory_col, candidate_col, ...
    reliability_col, reliability_score_col, gsa_score_col, total_score_col, selected_col, ...
    'VariableNames', {'Metric', 'Available', 'MandatoryPrimary', 'CandidatePrimary', ...
                      'Reliability', 'ReliabilityScore', 'GsaScore', 'TotalScore', 'Selected'});

end


function score = reliability_score(reliability)
switch lower(strtrim(reliability))
    case 'high'
        score = 3;
    case 'moderate'
        score = 2;
    case 'derived'
        score = 1;
    otherwise
        score = 0;
end
end


function score = metric_gsa_score(gsa_out, metric_name)
score = 0;
if isempty(gsa_out) || ~isstruct(gsa_out) || ~isfield(gsa_out, metric_name)
    return;
end
metric_out = gsa_out.(metric_name);
if isstruct(metric_out) && isfield(metric_out, 'ST')
    ST = metric_out.ST(:);
    ST = ST(isfinite(ST));
    if ~isempty(ST)
        score = max(abs(ST));
    end
end
end
