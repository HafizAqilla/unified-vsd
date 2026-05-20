function J = objective_calibration(x, params0, clinical, calib, scenario, pce_surrogate)
% OBJECTIVE_CALIBRATION
% -----------------------------------------------------------------------
% Tiered calibration objective:
%   J = J_primary + lambda_secondary * J_secondary + J_reg + J_invalid
%
% Primary metrics are normalised to the 5% target.
% Secondary metrics are normalised to the 10% target.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  3.1
% -----------------------------------------------------------------------

if nargin < 6
    pce_surrogate = [];
end

params = params0;
for i = 1:numel(calib.names)
    params = set_calibration_param_value( ...
        params, calib.referenceParams, calib.names{i}, x(i), calib.caseProfile);
end

[metrics, sim, validity_penalty] = evaluate_metrics(params, x, clinical, calib, scenario, pce_surrogate);
if isempty(metrics)
    J = validity_penalty;
    return;
end

targets = get_calibration_targets(scenario, clinical);
target_names = {targets.Metric};
systemic_bundle = build_systemic_bundle(clinical, scenario, targets);
pressure_bundle = build_pressure_waveform_bundle(clinical, scenario);
shunt_bundle = build_shunt_fraction_bundle(clinical, scenario);

J_primary = 0;
J_secondary = 0;
J_clinical_guard = 0;
for k = 1:numel(calib.metricFields)
    mf = calib.metricFields{k};
    if is_consistency_only_metric(calib, mf)
        continue;
    end
    idx = find(strcmp(target_names, mf), 1, 'first');
    if isempty(idx)
        continue;
    end
    y_clin = targets(idx).ClinicalValue;
    if isnan(y_clin) || ~isfield(metrics, mf) || ~isfinite(metrics.(mf))
        continue;
    end

    if should_route_metric_to_systemic_bundle(mf, systemic_bundle)
        % MAP, RAP, Qs, and derived SVR describe one coupled systemic load.
        % Route them through one uncertainty-aware bundle to prevent
        % double-counting the same clinical inconsistency.
        continue;
    end

    err_rel = abs(metrics.(mf) - y_clin) / max(abs(y_clin), 1e-6);
    weight = calib.weights.(mf);
    tier = calibration_metric_tier(calib, mf);
    if strcmp(tier, 'hard') || (~strcmp(tier, 'soft') && ismember(mf, calib.primaryMetrics))
        J_primary = J_primary + weight * (err_rel / calib.primaryTarget)^2;
    else
        J_secondary = J_secondary + weight * (err_rel / calib.secondaryTarget)^2;
    end

    J_clinical_guard = J_clinical_guard + clinical_guard_penalty(mf, err_rel, metrics.(mf), y_clin, calib);
end

J_reg = 0;
if isfield(calib, 'regLambda') && calib.regLambda > 0
    J_reg = calib.regLambda * sum(((x(:) - calib.x0(:)) ./ max(abs(calib.x0(:)), 1e-6)).^2);
end

J_systemic = systemic_bundle_penalty(metrics, systemic_bundle, calib);
J_pressure = pressure_waveform_penalty(metrics, pressure_bundle, calib);
J_shunt = shunt_fraction_penalty(metrics, shunt_bundle, calib);
J_param_plausibility = parameter_plausibility_penalty(x, calib);
J_boundary = boundary_plausibility_penalty(x, calib);

J = J_primary + calib.secondaryLambda * J_secondary + J_systemic + ...
    J_pressure + J_shunt + J_clinical_guard + J_reg + ...
    J_param_plausibility + J_boundary + validity_penalty;

if ~isempty(sim) && ~sim.ss_reached
    J = J + calib.invalidPenaltyScale;
end

function tf = is_consistency_only_metric(calib, metric_name)
tf = strcmp(calibration_metric_tier(calib, metric_name), 'consistency_check_only');
end

function tier = calibration_metric_tier(calib, metric_name)
tier = '';
if ~isfield(calib, 'targetTiers') || ~isstruct(calib.targetTiers) || ...
        ~isfield(calib.targetTiers, 'table') || isempty(calib.targetTiers.table)
    return;
end
tier_tbl = calib.targetTiers.table;
idx = find(strcmp(tier_tbl.Metric, metric_name), 1, 'first');
if ~isempty(idx)
    tier = tier_tbl.Tier{idx};
end
end

function penalty = parameter_plausibility_penalty(x, calib)
% PARAMETER_PLAUSIBILITY_PENALTY - penalize drift from scaled baseline anchors.
penalty = 0;
if ~isfield(calib, 'parameterRegistryActive') || isempty(calib.parameterRegistryActive)
    return;
end
lambda = 0.5;
if isfield(calib, 'paramPlausibilityLambda') && isfinite(calib.paramPlausibilityLambda)
    lambda = calib.paramPlausibilityLambda;
end
if lambda <= 0
    return;
end

registry = calib.parameterRegistryActive;
baseline = registry.baseline_scaled(:);
x = x(:);
valid = isfinite(x) & isfinite(baseline) & x > 0 & baseline > 0;
if ~any(valid)
    return;
end

weights = parameter_confidence_weights(registry.confidence_level);
ratio = x(valid) ./ baseline(valid);
penalty = lambda * sum(weights(valid) .* (log(ratio)).^2);
end

function penalty = boundary_plausibility_penalty(x, calib)
% BOUNDARY_PLAUSIBILITY_PENALTY - discourage optimizer from parking on bounds.
penalty = 0;
if ~isfield(calib, 'parameterRegistryActive') || isempty(calib.parameterRegistryActive)
    return;
end
lambda = 20.0;
if isfield(calib, 'boundaryPlausibilityLambda') && isfinite(calib.boundaryPlausibilityLambda)
    lambda = calib.boundaryPlausibilityLambda;
end
if lambda <= 0
    return;
end

lb = calib.parameterRegistryActive.lb(:);
ub = calib.parameterRegistryActive.ub(:);
x = x(:);
span = max(ub - lb, 1e-12);
dist_lower = (x - lb) ./ span;
dist_upper = (ub - x) ./ span;
near_lower = max(0, 0.10 - dist_lower);
near_upper = max(0, 0.10 - dist_upper);
penalty = lambda * sum(near_lower.^2 + near_upper.^2);
end

function weights = parameter_confidence_weights(confidence_levels)
weights = ones(numel(confidence_levels), 1);
for confidence_idx = 1:numel(confidence_levels)
    switch lower(strtrim(confidence_levels{confidence_idx}))
        case 'high'
            weights(confidence_idx) = 2.0;
        case 'medium'
            weights(confidence_idx) = 1.0;
        otherwise
            weights(confidence_idx) = 0.4;
    end
end
end

function tf = should_route_metric_to_systemic_bundle(metric_name, bundle)
% SHOULD_ROUTE_METRIC_TO_SYSTEMIC_BUNDLE - avoid duplicate systemic penalties.
systemic_metrics = {'SAP_mean','RAP_mean','CO_Lmin','SVR'};
tf = bundle.has_full_systemic_targets && ismember(metric_name, systemic_metrics);
end

function bundle = build_systemic_bundle(clinical, scenario, targets)
% BUILD_SYSTEMIC_BUNDLE — collect coupled systemic targets for consistency checks.
bundle = struct();
bundle.has_full_systemic_targets = false;
bundle.Qs_target_Lmin = NaN;
bundle.SAP_mean_target_mmHg = NaN;
bundle.RAP_mean_target_mmHg = NaN;
bundle.SVR_target_WU = NaN;
bundle.Qs_sigma_Lmin = NaN;
bundle.SAP_mean_sigma_mmHg = NaN;
bundle.RAP_mean_sigma_mmHg = NaN;
bundle.SVR_sigma_WU = NaN;
bundle.Qs_floor_Lmin = NaN;
bundle.Qs_floor_sigma_Lmin = NaN;
bundle.SVR_ceiling_WU = NaN;
bundle.SVR_ceiling_sigma_WU = NaN;

if ~isstruct(clinical) || ~isfield(clinical, scenario)
    return;
end

src = clinical.(scenario);
required_fields = {'CO_Lmin', 'SAP_mean_mmHg', 'RAP_mean_mmHg'};
if ~all(isfield(src, required_fields))
    return;
end

target_vec = [src.CO_Lmin, src.SAP_mean_mmHg, src.RAP_mean_mmHg];
if any(isnan(target_vec))
    return;
end

bundle.has_full_systemic_targets = true;
bundle.Qs_target_Lmin = src.CO_Lmin;
bundle.SAP_mean_target_mmHg = src.SAP_mean_mmHg;
bundle.RAP_mean_target_mmHg = src.RAP_mean_mmHg;
bundle.SVR_target_WU = (src.SAP_mean_mmHg - src.RAP_mean_mmHg) / max(src.CO_Lmin, 1e-6);
bundle.Qs_sigma_Lmin = target_sigma(targets, 'CO_Lmin', src.CO_Lmin, 0.15);
bundle.SAP_mean_sigma_mmHg = target_sigma(targets, 'SAP_mean', src.SAP_mean_mmHg, 0.05);
bundle.RAP_mean_sigma_mmHg = target_sigma(targets, 'RAP_mean', src.RAP_mean_mmHg, 0.05);
bundle.SVR_sigma_WU = target_sigma(targets, 'SVR', bundle.SVR_target_WU, 0.15);
bundle.Qs_floor_Lmin = max(0, bundle.Qs_target_Lmin - bundle.Qs_sigma_Lmin);
bundle.Qs_floor_sigma_Lmin = max(0.50 * bundle.Qs_sigma_Lmin, 0.05 * bundle.Qs_target_Lmin);
bundle.SVR_ceiling_WU = bundle.SVR_target_WU + bundle.SVR_sigma_WU;
bundle.SVR_ceiling_sigma_WU = max(bundle.SVR_sigma_WU, 0.05 * bundle.SVR_target_WU);
if isfield(src, 'SVR_WU') && ~isnan(src.SVR_WU)
    bundle.SVR_documented_WU = src.SVR_WU;
else
    bundle.SVR_documented_WU = NaN;
end
end

function penalty = systemic_bundle_penalty(metrics, bundle, calib)
% SYSTEMIC_BUNDLE_PENALTY — couple Qs, SVR, and systemic-flow balance.
penalty = 0;
if ~bundle.has_full_systemic_targets
    return;
end

required_metrics = {'CO_Lmin','SAP_mean','RAP_mean','SVR'};
for required_idx = 1:numel(required_metrics)
    if ~isfield(metrics, required_metrics{required_idx}) || ...
            ~isfinite(metrics.(required_metrics{required_idx}))
        return;
    end
end

qs_z = normalized_target_residual(metrics.CO_Lmin, bundle.Qs_target_Lmin, bundle.Qs_sigma_Lmin);
sap_z = normalized_target_residual(metrics.SAP_mean, bundle.SAP_mean_target_mmHg, bundle.SAP_mean_sigma_mmHg);
rap_z = normalized_target_residual(metrics.RAP_mean, bundle.RAP_mean_target_mmHg, bundle.RAP_mean_sigma_mmHg);
svr_z = normalized_target_residual(metrics.SVR, bundle.SVR_target_WU, bundle.SVR_sigma_WU);

flow_balance_z = 0;
if isfield(metrics, 'Qao_Lmin') && isfield(metrics, 'Qs_Lmin')
    flow_balance_z = abs(metrics.Qao_Lmin - metrics.Qs_Lmin) / ...
        max(bundle.Qs_sigma_Lmin, 1e-6);
end

weights = systemic_load_weights(calib);
qs_deficit_z = max(0, (bundle.Qs_target_Lmin - metrics.CO_Lmin) / ...
    max(bundle.Qs_sigma_Lmin, 1e-6));
svr_high_z = max(0, (metrics.SVR - bundle.SVR_target_WU) / ...
    max(bundle.SVR_sigma_WU, 1e-6));
co_floor_z = max(0, (bundle.Qs_floor_Lmin - metrics.CO_Lmin) / ...
    max(bundle.Qs_floor_sigma_Lmin, 1e-6));
svr_ceiling_z = max(0, (metrics.SVR - bundle.SVR_ceiling_WU) / ...
    max(bundle.SVR_ceiling_sigma_WU, 1e-6));

penalty = penalty + weights.qs * qs_z^2;
penalty = penalty + weights.sap * sap_z^2;
penalty = penalty + weights.rap * rap_z^2;
penalty = penalty + weights.svr * svr_z^2;
penalty = penalty + weights.high_svr * svr_high_z^2;
penalty = penalty + weights.co_floor * co_floor_z^2;
penalty = penalty + weights.svr_ceiling * svr_ceiling_z^2;
penalty = penalty + weights.low_qs_high_svr * ...
    qs_deficit_z^2 * (1 + svr_high_z);
penalty = penalty + weights.low_flow_svr_wall * ...
    co_floor_z^2 * (1 + svr_ceiling_z)^2;
penalty = penalty + weights.flow_balance * flow_balance_z^2;
end

function weights = systemic_load_weights(calib)
weights = struct( ...
    'qs', 0.80, ...
    'sap', 0.75, ...
    'rap', 0.40, ...
    'svr', 0.55, ...
    'high_svr', 0.00, ...
    'low_qs_high_svr', 0.00, ...
    'co_floor', 0.00, ...
    'svr_ceiling', 0.00, ...
    'low_flow_svr_wall', 0.00, ...
    'flow_balance', 0.20);
if isfield(calib, 'caseProfile') && isstruct(calib.caseProfile) && ...
        isfield(calib.caseProfile, 'systemicLoadWeights') && ...
        isstruct(calib.caseProfile.systemicLoadWeights)
    override = calib.caseProfile.systemicLoadWeights;
    fields = fieldnames(override);
    for field_idx = 1:numel(fields)
        weights.(fields{field_idx}) = override.(fields{field_idx});
    end
end
end

function sigma = target_sigma(targets, metric_name, target_value, fallback_fraction)
% TARGET_SIGMA - resolve metric sigma from target metadata.
sigma = abs(target_value) * fallback_fraction;
if ~isfinite(sigma) || sigma <= 0
    sigma = max(abs(target_value), 1) * 0.10;
end

idx = find(strcmp({targets.Metric}, metric_name), 1, 'first');
if isempty(idx)
    sigma = max(sigma, 1e-6);
    return;
end

if isfield(targets, 'UncertaintyAbs') && isfinite(targets(idx).UncertaintyAbs) && ...
        targets(idx).UncertaintyAbs > 0
    sigma = targets(idx).UncertaintyAbs;
elseif isfield(targets, 'UncertaintyFraction') && ...
        isfinite(targets(idx).UncertaintyFraction) && targets(idx).UncertaintyFraction > 0
    sigma = abs(target_value) * targets(idx).UncertaintyFraction;
end
sigma = max(sigma, 1e-6);
end

function residual = normalized_target_residual(model_value, target_value, sigma)
% NORMALIZED_TARGET_RESIDUAL - absolute residual in uncertainty units.
residual = abs(model_value - target_value) / max(sigma, 1e-6);
end

function bundle = build_pressure_waveform_bundle(clinical, scenario)
% BUILD_PRESSURE_WAVEFORM_BUNDLE — collect systolic/diastolic pressure shape targets.
bundle = struct();
bundle.has_systemic_shape = false;
bundle.has_pulmonary_shape = false;

if ~isstruct(clinical) || ~isfield(clinical, scenario)
    return;
end

src = clinical.(scenario);
if all(isfield(src, {'SAP_sys_mmHg', 'SAP_dia_mmHg'})) && ...
        all(~isnan([src.SAP_sys_mmHg, src.SAP_dia_mmHg]))
    bundle.has_systemic_shape = true;
    bundle.SAP_min_target_mmHg = src.SAP_dia_mmHg;
    bundle.SAP_max_target_mmHg = src.SAP_sys_mmHg;
    bundle.SAP_pulse_target_mmHg = src.SAP_sys_mmHg - src.SAP_dia_mmHg;
end

if all(isfield(src, {'PAP_sys_mmHg', 'PAP_dia_mmHg'})) && ...
        all(~isnan([src.PAP_sys_mmHg, src.PAP_dia_mmHg]))
    bundle.has_pulmonary_shape = true;
    bundle.PAP_min_target_mmHg = src.PAP_dia_mmHg;
    bundle.PAP_max_target_mmHg = src.PAP_sys_mmHg;
    bundle.PAP_pulse_target_mmHg = src.PAP_sys_mmHg - src.PAP_dia_mmHg;
end
end

function penalty = pressure_waveform_penalty(metrics, bundle, calib)
% PRESSURE_WAVEFORM_PENALTY — pressure-shape penalties beyond mean pressures.
penalty = 0;

if bundle.has_systemic_shape
    sap_max_err_rel = abs(metrics.SAP_max - bundle.SAP_max_target_mmHg) / ...
        max(abs(bundle.SAP_max_target_mmHg), 1e-6);
    sap_min_err_rel = abs(metrics.SAP_min - bundle.SAP_min_target_mmHg) / ...
        max(abs(bundle.SAP_min_target_mmHg), 1e-6);
    sap_pulse_err_rel = abs(metrics.SAP_pulse - bundle.SAP_pulse_target_mmHg) / ...
        max(abs(bundle.SAP_pulse_target_mmHg), 1e-6);

    penalty = penalty + 0.8 * (sap_max_err_rel / calib.secondaryTarget)^2;
    penalty = penalty + 0.6 * (sap_min_err_rel / calib.secondaryTarget)^2;
    penalty = penalty + 0.5 * (sap_pulse_err_rel / calib.secondaryTarget)^2;
end

if bundle.has_pulmonary_shape
    pap_pulse_err_rel = abs(metrics.PAP_pulse - bundle.PAP_pulse_target_mmHg) / ...
        max(abs(bundle.PAP_pulse_target_mmHg), 1e-6);
    penalty = penalty + 0.2 * (pap_pulse_err_rel / calib.secondaryTarget)^2;
end
end

function bundle = build_shunt_fraction_bundle(clinical, scenario)
% BUILD_SHUNT_FRACTION_BUNDLE — derive a soft shunt-severity target from Qp/Qs.
bundle = struct();
bundle.has_target = false;
bundle.has_flow_consistency_target = false;
bundle.VSD_frac_target_pct = NaN;
bundle.Q_shunt_target_Lmin = NaN;
bundle.collapse_floor_scale = NaN;
bundle.penalty_weight = 0;
bundle.flow_consistency_weight = 0.15;

if ~strcmp(scenario, 'pre_surgery') || ~isstruct(clinical) || ~isfield(clinical, scenario)
    return;
end

src = clinical.(scenario);
if ~isfield(src, 'QpQs') || isnan(src.QpQs) || src.QpQs <= 1
    return;
end

if isfield(src, 'CO_Lmin') && isfinite(src.CO_Lmin)
    bundle.has_flow_consistency_target = true;
    bundle.Q_shunt_target_Lmin = src.CO_Lmin * (src.QpQs - 1);
    if isfield(src, 'Q_shunt_Lmin') && isfinite(src.Q_shunt_Lmin) && src.Q_shunt_Lmin > 0
        bundle.Q_shunt_target_Lmin = src.Q_shunt_Lmin;
        bundle.flow_consistency_weight = 0.65;
    end
end

% Apply this guard only when clinical data clearly indicates a sizeable
% left-to-right shunt. Mild shunts can otherwise lose their best
% pressure-flow basin to an overly aggressive shunt-fraction target.
has_documented_shunt_flow = isfield(src, 'Q_shunt_Lmin') && ...
    ~isnan(src.Q_shunt_Lmin) && src.Q_shunt_Lmin > 0;

if src.QpQs >= 1.5 && has_documented_shunt_flow
    bundle.has_target = true;
    bundle.collapse_floor_scale = 0.68;
    bundle.penalty_weight = 0.70;
elseif src.QpQs >= 1.10 && (~isfield(src, 'CO_Lmin') || isnan(src.CO_Lmin))
    % When Qp/Qs is measured but systemic flow is not, keep a gentle floor
    % under shunt fraction so calibration does not mimic a near-closed VSD.
    bundle.has_target = true;
    bundle.collapse_floor_scale = 0.45;
    bundle.penalty_weight = 0.25;
elseif ~bundle.has_flow_consistency_target
    return;
end

bundle.VSD_frac_target_pct = 100 * (1 - 1 / src.QpQs);
end

function penalty = shunt_fraction_penalty(metrics, bundle, calib)
% SHUNT_FRACTION_PENALTY — discourage shunt collapse inconsistent with Qp/Qs.
penalty = 0;
if ~bundle.has_target && ~bundle.has_flow_consistency_target
    return;
end

if bundle.has_target && isfield(metrics, 'VSD_frac_pct') && isfinite(metrics.VSD_frac_pct)
    collapse_floor_pct = bundle.collapse_floor_scale * bundle.VSD_frac_target_pct;
    if metrics.VSD_frac_pct < collapse_floor_pct
        collapse_rel = (collapse_floor_pct - metrics.VSD_frac_pct) / ...
            max(abs(bundle.VSD_frac_target_pct), 1e-6);
        penalty = penalty + bundle.penalty_weight * (collapse_rel / calib.secondaryTarget)^2;
    end
end

if bundle.has_flow_consistency_target && isfield(metrics, 'Qp_minus_Qs_Lmin') && ...
        isfield(metrics, 'Qvsd_Lmin') && isfinite(metrics.Qp_minus_Qs_Lmin) && ...
        isfinite(metrics.Qvsd_Lmin)
    accounting_gap_rel = abs(metrics.Qvsd_Lmin - metrics.Qp_minus_Qs_Lmin) / ...
        max(abs(bundle.Q_shunt_target_Lmin), 1e-6);
    target_gap_rel = abs(metrics.Qp_minus_Qs_Lmin - bundle.Q_shunt_target_Lmin) / ...
        max(abs(bundle.Q_shunt_target_Lmin), 1e-6);
    penalty = penalty + bundle.flow_consistency_weight * ...
        ((accounting_gap_rel / calib.secondaryTarget)^2 + ...
        0.25 * (target_gap_rel / calib.secondaryTarget)^2);
end
end

function penalty = clinical_guard_penalty(metric_name, err_rel, y_model, y_clin, calib)
penalty = 0;

pressure_metrics = {'RAP_mean','LAP_mean','PAP_min','PAP_max','PAP_mean', ...
                    'SAP_min','SAP_max','SAP_mean'};
volume_metrics = {'LVEDV','LVESV','RVEDV','RVESV'};
function_metrics = {'LVEF','RVEF'};
flow_metrics = {'CO_Lmin','SVR','PVR'};
guard = profile_volume_flow_guard(calib);

pressure_rel_thresh = 0.12;
pressure_scale = 35;
mean_pressure_rel_thresh = 0.08;
mean_pressure_scale = 55;

volume_rel_thresh = 0.15;
volume_scale = 30;
volume_upper_ratio = 1.25;
volume_upper_scale = 60;
preload_lower_ratio = 0.90;
preload_lower_scale = 45;

flow_rel_thresh = 0.20;
flow_scale = 15;
function_rel_thresh = 0.15;
function_scale = 25;

if ismember(metric_name, pressure_metrics)
    penalty = penalty + soft_excess_penalty(err_rel, pressure_rel_thresh, pressure_scale);
end
if ismember(metric_name, {'RAP_mean','LAP_mean','PAP_mean'})
    penalty = penalty + soft_excess_penalty(err_rel, mean_pressure_rel_thresh, mean_pressure_scale);
end
if ismember(metric_name, volume_metrics)
    penalty = penalty + soft_excess_penalty(err_rel, volume_rel_thresh, volume_scale);
    if y_model > volume_upper_ratio * y_clin
        ratio_excess = (y_model / max(y_clin, 1e-6)) - volume_upper_ratio;
        penalty = penalty + volume_upper_scale * ratio_excess^2;
    end
    if ismember(metric_name, {'LVEDV','RVEDV'}) && y_model < preload_lower_ratio * y_clin
        ratio_deficit = preload_lower_ratio - (y_model / max(y_clin, 1e-6));
        penalty = penalty + preload_lower_scale * ratio_deficit^2;
    end
end
if ismember(metric_name, function_metrics)
    penalty = penalty + soft_excess_penalty(err_rel, function_rel_thresh, function_scale);
end
if ismember(metric_name, flow_metrics)
    penalty = penalty + soft_excess_penalty(err_rel, flow_rel_thresh, flow_scale);
end
if guard.enabled
    ratio = y_model / max(y_clin, 1e-6);
    if strcmp(metric_name, 'RVEDV') && ratio > guard.RVEDV_upper_ratio
        penalty = penalty + guard.RVEDV_upper_scale * (ratio - guard.RVEDV_upper_ratio)^2;
    end
    if strcmp(metric_name, 'LVESV') && ratio < guard.LVESV_lower_ratio
        penalty = penalty + guard.LVESV_lower_scale * (guard.LVESV_lower_ratio - ratio)^2;
    end
    if strcmp(metric_name, 'LVEF') && ratio > guard.LVEF_upper_ratio
        penalty = penalty + guard.LVEF_upper_scale * (ratio - guard.LVEF_upper_ratio)^2;
    end
    if strcmp(metric_name, 'CO_Lmin') && y_model < y_clin && err_rel > guard.CO_low_rel_thresh
        penalty = penalty + guard.CO_low_scale * ...
            ((err_rel - guard.CO_low_rel_thresh) / guard.CO_low_rel_thresh)^2;
    end
end
end

function guard = profile_volume_flow_guard(calib)
guard = struct( ...
    'enabled', false, ...
    'RVEDV_upper_ratio', Inf, ...
    'RVEDV_upper_scale', 0, ...
    'LVESV_lower_ratio', -Inf, ...
    'LVESV_lower_scale', 0, ...
    'LVEF_upper_ratio', Inf, ...
    'LVEF_upper_scale', 0, ...
    'CO_low_rel_thresh', Inf, ...
    'CO_low_scale', 0);
if ~isfield(calib, 'caseProfile') || ~isfield(calib.caseProfile, 'volumeFlowGuard') || ...
        ~isstruct(calib.caseProfile.volumeFlowGuard)
    return;
end
profile_guard = calib.caseProfile.volumeFlowGuard;
names = fieldnames(guard);
for i = 1:numel(names)
    if isfield(profile_guard, names{i}) && ~isempty(profile_guard.(names{i}))
        guard.(names{i}) = profile_guard.(names{i});
    end
end
end

function penalty = soft_excess_penalty(err_rel, threshold, scale)
excess = max(0, err_rel - threshold);
penalty = scale * (excess / threshold)^2;
end

end

function [metrics, sim, penalty] = evaluate_metrics(params, ~, clinical, calib, scenario, pce_surrogate)
metrics = struct();
sim = [];
penalty = 0;

use_surrogate = ~isempty(pce_surrogate) && can_use_surrogate(calib.metricFields, pce_surrogate);
if use_surrogate
    params_surrogate = params;
    x_full = pce_surrogate.cfg.x0;
    for idx_full = 1:numel(pce_surrogate.cfg.names)
        x_full(idx_full) = get_calibration_param_value( ...
            params_surrogate, params_surrogate, pce_surrogate.cfg.names{idx_full}, struct());
    end
    x_full = min(max(x_full(:), pce_surrogate.cfg.lb(:)), pce_surrogate.cfg.ub(:));
    for k = 1:numel(calib.metricFields)
        mf = calib.metricFields{k};
        metrics.(mf) = uq_evalModel(pce_surrogate.(mf).surrogate, x_full(:)');
    end
    return;
end

params.sim.nCyclesSteady = min(max(params.sim.nCyclesSteady, 30), 60);
params.sim.ss_tol_P = max(params.sim.ss_tol_P, 0.5);
params.sim.ss_tol_V = max(params.sim.ss_tol_V, 0.5);

try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
    validity = evaluate_simulation_validity(sim, params, metrics, scenario, clinical);
    penalty = validity.penalty;
catch
    metrics = [];
    sim = [];
    penalty = calib.invalidPenaltyScale * 10;
end
end

function tf = can_use_surrogate(metric_fields, pce_surrogate)
tf = true;
for k = 1:numel(metric_fields)
    mf = metric_fields{k};
    if ~isfield(pce_surrogate, mf) || ~isfield(pce_surrogate.(mf), 'surrogate')
        tf = false;
        return;
    end
end
end
