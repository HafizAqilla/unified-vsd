function target_config = build_target_tiers(clinical, scenario, audit, config)
% BUILD_TARGET_TIERS
% -----------------------------------------------------------------------
% Builds explicit calibration target tiers from clinical availability and
% pre-calibration consistency audit results.
%
% INPUTS:
%   clinical  - unified clinical struct from config/ patient profiles    [-]
%   scenario  - scenario string: 'pre_surgery' | 'post_surgery'          [-]
%   audit     - audit_clinical_consistency output struct                 [-]
%   config    - optional target-tier policy struct                       [-]
%
% OUTPUTS:
%   target_config - struct with hard/soft/report-only target lists,
%                   weights, inclusion flags, and an audit table         [-]
%
% ASSUMPTIONS:
%   - Targets excluded from calibration remain visible in validation.
%   - Consistency-only, derived-validation, and validation-holdout targets
%     are excluded from primary RMSE but included in full RMSE for
%     transparent reporting.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%   [2] docs/calibration_data_governance_notes.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-14
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(scenario)
    scenario = 'pre_surgery';
end
if nargin < 3 || isempty(audit)
    audit = audit_clinical_consistency(clinical, scenario);
end
if nargin < 4 || isempty(config)
    config = default_target_tier_config();
else
    config = merge_struct(default_target_tier_config(), config);
end
config = apply_scenario_target_tier_config(config, scenario);

targets = get_calibration_targets(scenario, clinical);      % [-]
metric_names = {targets.Metric};                            % [cellstr]
clinical_values = [targets.ClinicalValue];                  % [mixed units]
available = isfinite(clinical_values);                      % [-]
available_metrics = metric_names(available);                % [cellstr]

hard = intersect(config.hard, available_metrics, 'stable');
soft = intersect(config.soft, available_metrics, 'stable');
consistency_only = intersect(config.consistency_only, available_metrics, 'stable');
derived_validation = intersect(config.derived_validation, available_metrics, 'stable');
validation_holdout = intersect(config.validation_holdout, available_metrics, 'stable');
primary_rmse_holdout = intersect(config.primary_rmse_holdout, available_metrics, 'stable');
consistency_reasons = struct();
holdout_reasons = struct();

if isfield(audit, 'recommended_target_tier_changes')
    changes = audit.recommended_target_tier_changes;
    for idx = 1:numel(changes)
        metric_name = changes(idx).metric;
        if ~ismember(metric_name, available_metrics)
            continue;
        end
        if strcmpi(changes(idx).tier, 'consistency_check_only')
            consistency_only = unique([consistency_only, {metric_name}], 'stable');
            consistency_reasons.(metric_name) = changes(idx).reason;
        end
    end
end

% EF is mathematically derived from LVEDV and LVESV. When all three are
% available, fitting EF in addition to both volumes double-counts the same
% echo measurement block. Keep it visible for validation, but do not let it
% act as an independent calibration/RMSE anchor.
lv_volume_fit_metrics = intersect({'LVEDV','LVESV'}, ...
    unique([hard, soft], 'stable'), 'stable');
if all(ismember({'LVEDV','LVESV','LVEF'}, available_metrics)) && ...
        numel(lv_volume_fit_metrics) == 2
    consistency_only = unique([consistency_only, {'LVEF'}], 'stable');
    consistency_reasons.LVEF = ['EF is directly derived from LVEDV and ', ...
        'LVESV; excluded from fitting to avoid double-counting echo volumes.'];
end

% In sparse catheterisation records without echo volume anchors, direct
% ventricular EDPs are clinically visible but weakly identifiable in this
% lumped model. Report them as holdout rows instead of letting them dominate
% a pressure-flow calibration.
if ismember('LVEDP', available_metrics) && ...
        ~any(ismember({'LVEDV','LVESV'}, available_metrics))
    validation_holdout = unique([validation_holdout, {'LVEDP'}], 'stable');
    holdout_reasons.LVEDP = ['LVEDP is reported as validation holdout because ', ...
        'no LV volume anchor is available to identify diastolic stiffness/preload.'];
end
if ismember('RVEDP', available_metrics) && ...
        ~any(ismember({'RVEDV','RVESV'}, available_metrics))
    validation_holdout = unique([validation_holdout, {'RVEDP'}], 'stable');
    holdout_reasons.RVEDP = ['RVEDP is reported as validation holdout because ', ...
        'no RV volume anchor is available to identify diastolic stiffness/preload.'];
end

hard = setdiff(hard, consistency_only, 'stable');
hard = setdiff(hard, derived_validation, 'stable');
hard = setdiff(hard, validation_holdout, 'stable');
soft = setdiff(soft, consistency_only, 'stable');
soft = setdiff(soft, derived_validation, 'stable');
soft = setdiff(soft, validation_holdout, 'stable');
soft = setdiff(soft, hard, 'stable');

included_in_calibration = unique([hard, soft], 'stable');
included_in_primary_rmse = setdiff(available_metrics, ...
    unique([consistency_only, derived_validation, validation_holdout, ...
    primary_rmse_holdout], 'stable'), 'stable');
excluded_from_primary_rmse = unique([consistency_only, derived_validation, ...
    validation_holdout, primary_rmse_holdout], 'stable');

weights = struct();
for idx = 1:numel(hard)
    weights.(hard{idx}) = config.hard_weight_multiplier;
end
for idx = 1:numel(soft)
    weights.(soft{idx}) = config.soft_weight_multiplier;
end
specific_names = fieldnames(config.metric_weight_multipliers);
for idx = 1:numel(specific_names)
    metric_name = specific_names{idx};
    if ismember(metric_name, included_in_calibration)
        weights.(metric_name) = config.metric_weight_multipliers.(metric_name);
    end
end

target_config = struct();
target_config.policy = config.policy_name;
target_config.hard = hard;
target_config.soft = soft;
target_config.consistency_only = consistency_only;
target_config.derived_validation = derived_validation;
target_config.validation_holdout = validation_holdout;
target_config.primary_rmse_holdout = primary_rmse_holdout;
target_config.included_in_calibration = included_in_calibration;
target_config.included_in_primary_rmse = included_in_primary_rmse;
target_config.excluded_from_primary_rmse = excluded_from_primary_rmse;
target_config.weights = weights;
target_config.audit_summary = audit.summary;
target_config.consistency_reasons = consistency_reasons;
target_config.holdout_reasons = holdout_reasons;
target_config.table = build_tier_table(targets, hard, soft, ...
    consistency_only, derived_validation, validation_holdout, ...
    primary_rmse_holdout, audit, consistency_reasons, holdout_reasons);
end

function config = default_target_tier_config()
config = struct();
config.policy_name = 'flow_volume_consistency_governance_v1';
config.hard = {'CO_Lmin','QpQs','PAP_mean','SAP_mean','RAP_mean', ...
    'LVEDV','LVESV','LVEF'};
config.soft = {'Q_shunt_Lmin','SAP_max','SAP_min','RVESV'};
config.consistency_only = {};
config.derived_validation = {'PVR','SVR'};
config.validation_holdout = {};
config.primary_rmse_holdout = {};
config.hard_weight_multiplier = 1.00;
config.soft_weight_multiplier = 0.45;
config.metric_weight_multipliers = struct( ...
    'CO_Lmin', 1.10, ...
    'QpQs', 1.00, ...
    'Q_shunt_Lmin', 1.20, ...
    'PAP_mean', 0.90, ...
    'SAP_mean', 1.05, ...
    'RAP_mean', 0.90, ...
    'LVEDV', 0.80, ...
    'LVESV', 0.85, ...
    'LVEF', 0.85, ...
    'SAP_max', 0.45, ...
    'SAP_min', 0.40, ...
    'RVESV', 0.45);
end

function config = apply_scenario_target_tier_config(config, scenario)
% APPLY_SCENARIO_TARGET_TIER_CONFIG - scenario-specific target governance.
if ~strcmp(char(scenario), 'post_surgery')
    return;
end

config.hard = setdiff(config.hard, {'LVESV'}, 'stable');
config.soft = unique([config.soft, {'RVEDV','RVESV'}], 'stable');
config.validation_holdout = unique([config.validation_holdout, {'LVESV'}], 'stable');
config.metric_weight_multipliers.RVEDV = 0.45;
config.metric_weight_multipliers.RVESV = 0.45;
end

function tier_table = build_tier_table(targets, hard, soft, consistency_only, ...
    derived_validation, validation_holdout, primary_rmse_holdout, audit, ...
    consistency_reasons, holdout_reasons)
n_targets = numel(targets);
metric_col = cell(n_targets, 1);
tier_col = cell(n_targets, 1);
included_cal_col = false(n_targets, 1);
included_primary_rmse_col = false(n_targets, 1);
flag_col = cell(n_targets, 1);
reason_col = cell(n_targets, 1);

for idx = 1:n_targets
    metric_name = targets(idx).Metric;
    metric_col{idx} = metric_name;
    flag_col{idx} = 'none';
    reason_col{idx} = 'none';

    if ismember(metric_name, consistency_only)
        tier_col{idx} = 'consistency_check_only';
        included_cal_col(idx) = false;
        included_primary_rmse_col(idx) = false;
        [flag_col{idx}, reason_col{idx}] = consistency_flag(metric_name, audit, consistency_reasons);
    elseif ismember(metric_name, derived_validation)
        tier_col{idx} = 'derived_validation';
        included_cal_col(idx) = false;
        included_primary_rmse_col(idx) = false;
        flag_col{idx} = 'derived_validation';
        reason_col{idx} = 'Derived from source pressures/flows; retained for validation only.';
    elseif ismember(metric_name, validation_holdout)
        tier_col{idx} = 'validation_holdout';
        included_cal_col(idx) = false;
        included_primary_rmse_col(idx) = false;
        flag_col{idx} = 'validation_holdout';
        reason_col{idx} = holdout_reason(metric_name, holdout_reasons);
    elseif ismember(metric_name, hard)
        tier_col{idx} = 'hard';
        included_cal_col(idx) = true;
        included_primary_rmse_col(idx) = ~ismember(metric_name, primary_rmse_holdout);
        if ~included_primary_rmse_col(idx)
            flag_col{idx} = 'primary_rmse_holdout';
            reason_col{idx} = 'Used during calibration/reporting but excluded from governed primary RMSE.';
        end
    elseif ismember(metric_name, soft)
        tier_col{idx} = 'soft';
        included_cal_col(idx) = true;
        included_primary_rmse_col(idx) = ~ismember(metric_name, primary_rmse_holdout);
        if ~included_primary_rmse_col(idx)
            flag_col{idx} = 'primary_rmse_holdout';
            reason_col{idx} = 'Used during calibration/reporting but excluded from governed primary RMSE.';
        end
    elseif isfinite(targets(idx).ClinicalValue)
        tier_col{idx} = 'validation_only';
        included_cal_col(idx) = false;
        included_primary_rmse_col(idx) = true;
    else
        tier_col{idx} = 'unavailable';
        included_cal_col(idx) = false;
        included_primary_rmse_col(idx) = false;
    end
end

tier_table = table(metric_col, tier_col, included_cal_col, ...
    included_primary_rmse_col, flag_col, reason_col, ...
    'VariableNames', {'Metric','Tier','IncludedInCalibration', ...
    'IncludedInPrimaryRMSE','Flag','Reason'});
end

function reason = holdout_reason(metric_name, holdout_reasons)
reason = 'Target is retained as transparent validation holdout.';
if isstruct(holdout_reasons) && isfield(holdout_reasons, metric_name)
    reason = holdout_reasons.(metric_name);
end
end

function [flag, reason] = consistency_flag(metric_name, audit, consistency_reasons)
flag = 'consistency_check_only';
reason = 'Target is retained for validation but excluded from fitting.';
if isstruct(consistency_reasons) && isfield(consistency_reasons, metric_name)
    reason = consistency_reasons.(metric_name);
end
if ~isstruct(audit) || ~isfield(audit, 'recommended_target_tier_changes')
    return;
end
changes = audit.recommended_target_tier_changes;
for idx = 1:numel(changes)
    if strcmp(changes(idx).metric, metric_name)
        flag = 'inconsistent_or_unverified';
        reason = changes(idx).reason;
        return;
    end
end
end

function merged = merge_struct(defaults, overrides)
merged = defaults;
fields = fieldnames(overrides);
for idx = 1:numel(fields)
    merged.(fields{idx}) = overrides.(fields{idx});
end
end
