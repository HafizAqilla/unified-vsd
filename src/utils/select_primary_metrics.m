function [primary_metrics, selection_table] = select_primary_metrics(clinical, gsa_out, scenario, caseProfile)
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
%   caseProfile      - optional evidence-aware calibration profile struct
%
% OUTPUTS:
%   primary_metrics  - cell array of selected model metric names
%   selection_table  - table documenting scores and reasons
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 4 || isempty(caseProfile)
    caseProfile = struct();
end

targets = get_calibration_targets(scenario, clinical);
n_targets = numel(targets);
allowed_metrics = allowed_metric_names(caseProfile);
max_primary_metrics = 5;   % [count] historical default for full-data cases
if isfield(caseProfile, 'maxPrimaryMetrics') && isfinite(caseProfile.maxPrimaryMetrics)
    max_primary_metrics = max(1, round(caseProfile.maxPrimaryMetrics));
end

metric_col = cell(n_targets, 1);
available_col = false(n_targets, 1);
allowed_col = true(n_targets, 1);
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
    allowed_col(i) = isempty(allowed_metrics) || ismember(targets(i).Metric, allowed_metrics);
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

eligible_mandatory = available_col & allowed_col & mandatory_col;
selected_col(eligible_mandatory) = true;

slots_remaining = max(0, max_primary_metrics - nnz(selected_col));
preferred_metrics = preferred_metric_names(caseProfile);
if slots_remaining > 0 && ~isempty(preferred_metrics)
    preferred_idx = preferred_metric_indices(metric_col, available_col, allowed_col, ...
        selected_col, preferred_metrics);
    preferred_idx = preferred_idx(1:min(slots_remaining, numel(preferred_idx)));
    selected_col(preferred_idx) = true;
end

slots_remaining = max(0, max_primary_metrics - nnz(selected_col));
eligible_candidates = available_col & allowed_col & candidate_col & ~selected_col;
candidate_idx = find(eligible_candidates);
[~, order] = sort(total_score_col(candidate_idx), 'descend');
candidate_idx = candidate_idx(order);
candidate_idx = candidate_idx(1:min(slots_remaining, numel(candidate_idx)));
selected_col(candidate_idx) = true;

% If sparse clinical data leave fewer than five metrics, fill with any
% calibrated available metric so each run still has a reproducible target set.
slots_remaining = max(0, max_primary_metrics - nnz(selected_col));
if slots_remaining > 0
    fallback_idx = find(available_col & allowed_col & [targets.UseForCalibration]' & ~selected_col);
    [~, order] = sort(total_score_col(fallback_idx), 'descend');
    fallback_idx = fallback_idx(order);
    fallback_idx = fallback_idx(1:min(slots_remaining, numel(fallback_idx)));
    selected_col(fallback_idx) = true;
end

primary_metrics = metric_col(selected_col);
selection_table = table(metric_col, available_col, allowed_col, mandatory_col, candidate_col, ...
    reliability_col, reliability_score_col, gsa_score_col, total_score_col, selected_col, ...
    'VariableNames', {'Metric', 'Available', 'AllowedByCaseProfile', 'MandatoryPrimary', 'CandidatePrimary', ...
                      'Reliability', 'ReliabilityScore', 'GsaScore', 'TotalScore', 'Selected'});

end

function allowed_metrics = allowed_metric_names(caseProfile)
allowed_metrics = {};
if isfield(caseProfile, 'allowedMetricFields') && ~isempty(caseProfile.allowedMetricFields)
    allowed_metrics = caseProfile.allowedMetricFields(:)';
end
if isfield(caseProfile, 'validationHoldoutMetrics') && ~isempty(caseProfile.validationHoldoutMetrics)
    allowed_metrics = setdiff(allowed_metrics, caseProfile.validationHoldoutMetrics(:)', 'stable');
end
end

function preferred_metrics = preferred_metric_names(caseProfile)
preferred_metrics = {};
if isfield(caseProfile, 'preferredPrimaryMetrics') && ~isempty(caseProfile.preferredPrimaryMetrics)
    preferred_metrics = caseProfile.preferredPrimaryMetrics(:)';
end
end

function preferred_idx = preferred_metric_indices(metric_col, available_col, allowed_col, selected_col, preferred_metrics)
preferred_idx = [];
for i = 1:numel(preferred_metrics)
    idx = find(strcmp(metric_col, preferred_metrics{i}) & available_col & allowed_col & ~selected_col, 1);
    if ~isempty(idx)
        preferred_idx(end + 1) = idx; %#ok<AGROW>
    end
end
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
