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
derived_reasons = struct();

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

% Post-operative BP echo records may report LVESV and EF directly, while
% LVEDV, stroke volume, and CO are calculated from those same source fields.
% Fit the independent echo pair and keep algebraic consequences as derived
% validation rows, avoiding repeated leverage from one measurement block.
if strcmp(char(scenario), 'post_surgery')
    [derived_validation, hard, soft, validation_holdout, derived_reasons] = ...
        apply_post_op_echo_derived_tiers(clinical, targets, ...
        derived_validation, hard, soft, validation_holdout, derived_reasons);
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
    primary_rmse_holdout, audit, consistency_reasons, holdout_reasons, ...
    derived_reasons);
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

config.soft = unique([config.soft, {'RVEDV','RVESV'}], 'stable');
config.metric_weight_multipliers.RVEDV = 0.45;
config.metric_weight_multipliers.RVESV = 0.45;
end

function [derived_validation, hard, soft, validation_holdout, derived_reasons] = ...
    apply_post_op_echo_derived_tiers(clinical, targets, derived_validation, ...
    hard, soft, validation_holdout, derived_reasons)
% APPLY_POST_OP_ECHO_DERIVED_TIERS - classify algebraic BP echo derivatives.
src = clinical.post_surgery;                                % [-]
metric_names = {targets.Metric};                            % [cellstr]

has_lvesv = has_finite_field(src, 'LVESV_mL');              % [-]
has_ef = has_finite_field(src, 'EF') || has_finite_field(src, 'LVEF'); % [-]
if ~(has_lvesv && has_ef)
    return;
end

hard = unique([hard, {'LVESV','LVEF'}], 'stable');
validation_holdout = setdiff(validation_holdout, {'LVESV','LVEF'}, 'stable');

LVESV_mL = finite_field(src, 'LVESV_mL');                   % [mL]
EF = finite_field(src, 'EF');                               % [fraction]
if ~isfinite(EF)
    EF = finite_field(src, 'LVEF');                         % [fraction]
end
LVEDV_derived_mL = LVESV_mL / max(1 - EF, 1e-9);            % [mL]

if has_finite_field(src, 'LVEDV_mL') && ...
        values_match(src.LVEDV_mL, LVEDV_derived_mL)
    derived_validation = unique([derived_validation, {'LVEDV'}], 'stable');
    hard = setdiff(hard, {'LVEDV'}, 'stable');
    soft = setdiff(soft, {'LVEDV'}, 'stable');
    derived_reasons.LVEDV = ['LVEDV is algebraically derived from ', ...
        'post-operative BP echo LVESV and EF.'];
end

if has_finite_field(src, 'CO_Lmin')
    SV_lv_mL = LVEDV_derived_mL - LVESV_mL;                 % [mL/beat]
    HR_bpm = scenario_HR_bpm(clinical);                     % [bpm]
    CO_derived_Lmin = SV_lv_mL * HR_bpm / 1000;             % [L/min]
    if values_match(src.CO_Lmin, CO_derived_Lmin) || ...
            ~has_independent_CO_comparator(src)
        derived_validation = unique([derived_validation, {'CO_Lmin'}], 'stable');
        hard = setdiff(hard, {'CO_Lmin'}, 'stable');
        soft = setdiff(soft, {'CO_Lmin'}, 'stable');
        derived_reasons.CO_Lmin = ['CO_Lmin is echo-derived from ', ...
            'LV stroke volume and heart rate, not an independent flow measurement.'];
    end
end

if ~ismember('LVEF', metric_names)
    hard = setdiff(hard, {'LVEF'}, 'stable');
end
end

function tier_table = build_tier_table(targets, hard, soft, consistency_only, ...
    derived_validation, validation_holdout, primary_rmse_holdout, audit, ...
    consistency_reasons, holdout_reasons, derived_reasons)
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
        reason_col{idx} = derived_validation_reason(metric_name, derived_reasons);
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

function reason = derived_validation_reason(metric_name, derived_reasons)
reason = 'Derived from source measurements; retained for validation only.';
if ismember(metric_name, {'PVR','SVR'})
    reason = 'Derived from source pressures/flows; retained for validation only.';
end
if isstruct(derived_reasons) && isfield(derived_reasons, metric_name)
    reason = derived_reasons.(metric_name);
end
end

function reason = holdout_reason(metric_name, holdout_reasons)
reason = 'Target is retained as transparent validation holdout.';
if isstruct(holdout_reasons) && isfield(holdout_reasons, metric_name)
    reason = holdout_reasons.(metric_name);
end
end

function value = finite_field(src, field_name)
value = NaN;
if isstruct(src) && isfield(src, field_name) && isnumeric(src.(field_name)) && ...
        isscalar(src.(field_name)) && isfinite(src.(field_name))
    value = src.(field_name);
end
end

function tf = has_finite_field(src, field_name)
tf = isfinite(finite_field(src, field_name));
end

function tf = values_match(observed, expected)
abs_tol = 1e-6;                                              % [same unit]
rel_tol = 1e-4;                                              % [-]
tf = isfinite(observed) && isfinite(expected) && ...
    abs(observed - expected) <= max(abs_tol, rel_tol * max(abs(expected), 1));
end

function HR_bpm = scenario_HR_bpm(clinical)
HR_bpm = NaN;                                                % [bpm]
if isfield(clinical, 'post_surgery') && ...
        has_finite_field(clinical.post_surgery, 'HR')
    HR_bpm = clinical.post_surgery.HR;                       % [bpm]
elseif isfield(clinical, 'common') && has_finite_field(clinical.common, 'HR')
    HR_bpm = clinical.common.HR;                             % [bpm]
end
end

function tf = has_independent_CO_comparator(src)
tf = false;
if ~isstruct(src) || ~isfield(src, 'CO_comparator') || isempty(src.CO_comparator)
    return;
end
comparator = lower(strtrim(char(src.CO_comparator)));        % [-]
tf = any(strcmp(comparator, {'qs_lmin','fick','thermodilution','cath'}));
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
