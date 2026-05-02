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
% DATE:     2026-04-28
% VERSION:  3.0
% -----------------------------------------------------------------------

if nargin < 6
    pce_surrogate = [];
end

params = params0;
for i = 1:numel(calib.names)
    params = set_param_by_name(params, calib.names{i}, x(i));
end

[metrics, sim, validity_penalty] = evaluate_metrics(params, x, clinical, calib, scenario, pce_surrogate);
if isempty(metrics)
    J = validity_penalty;
    return;
end

targets = get_calibration_targets(scenario, clinical);
target_names = {targets.Metric};
systemic_bundle = build_systemic_bundle(clinical, scenario);
pressure_bundle = build_pressure_waveform_bundle(clinical, scenario);

J_primary = 0;
J_secondary = 0;
J_clinical_guard = 0;
for k = 1:numel(calib.metricFields)
    mf = calib.metricFields{k};
    idx = find(strcmp(target_names, mf), 1, 'first');
    if isempty(idx)
        continue;
    end
    y_clin = targets(idx).ClinicalValue;
    if isnan(y_clin) || ~isfield(metrics, mf) || ~isfinite(metrics.(mf))
        continue;
    end

    if strcmp(mf, 'SVR') && systemic_bundle.has_full_systemic_targets
        % When MAP, RAP, and Qs are all explicit clinical targets, SVR is
        % a derived systemic consistency check rather than an independent
        % free target. Skipping the direct SVR term avoids double-counting
        % the same mismatch while keeping a dedicated bundle penalty below.
        continue;
    end

    err_rel = abs(metrics.(mf) - y_clin) / max(abs(y_clin), 1e-6);
    weight = calib.weights.(mf);
    if ismember(mf, calib.primaryMetrics)
        J_primary = J_primary + weight * (err_rel / calib.primaryTarget)^2;
    else
        J_secondary = J_secondary + weight * (err_rel / calib.secondaryTarget)^2;
    end

    J_clinical_guard = J_clinical_guard + clinical_guard_penalty(mf, err_rel, metrics.(mf), y_clin);
end

J_reg = 0;
if isfield(calib, 'regLambda') && calib.regLambda > 0
    J_reg = calib.regLambda * sum(((x(:) - calib.x0(:)) ./ max(abs(calib.x0(:)), 1e-6)).^2);
end

J_systemic = systemic_bundle_penalty(metrics, systemic_bundle, calib);
J_pressure = pressure_waveform_penalty(metrics, pressure_bundle, calib);

J = J_primary + calib.secondaryLambda * J_secondary + J_systemic + ...
    J_pressure + J_clinical_guard + J_reg + validity_penalty;

if ~isempty(sim) && ~sim.ss_reached
    J = J + calib.invalidPenaltyScale;
end

function bundle = build_systemic_bundle(clinical, scenario)
% BUILD_SYSTEMIC_BUNDLE — collect coupled systemic targets for consistency checks.
bundle = struct();
bundle.has_full_systemic_targets = false;
bundle.Qs_target_Lmin = NaN;
bundle.SAP_mean_target_mmHg = NaN;
bundle.RAP_mean_target_mmHg = NaN;
bundle.SVR_target_WU = NaN;

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

qs_err_rel = abs(metrics.CO_Lmin - bundle.Qs_target_Lmin) / max(abs(bundle.Qs_target_Lmin), 1e-6);
svr_err_rel = abs(metrics.SVR - bundle.SVR_target_WU) / max(abs(bundle.SVR_target_WU), 1e-6);

flow_balance_rel = 0;
if isfield(metrics, 'Qao_Lmin') && isfield(metrics, 'Qs_Lmin')
    flow_balance_rel = abs(metrics.Qao_Lmin - metrics.Qs_Lmin) / ...
        max(abs(bundle.Qs_target_Lmin), 1e-6);
end

penalty = penalty + 1.4 * (qs_err_rel / calib.primaryTarget)^2;
penalty = penalty + 0.7 * (svr_err_rel / calib.secondaryTarget)^2;
penalty = penalty + 0.4 * (flow_balance_rel / calib.secondaryTarget)^2;
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

function penalty = clinical_guard_penalty(metric_name, err_rel, y_model, y_clin)
penalty = 0;

pressure_metrics = {'RAP_mean','LAP_mean','PAP_min','PAP_max','PAP_mean', ...
                    'SAP_min','SAP_max','SAP_mean'};
volume_metrics = {'LVEDV','LVESV','RVEDV','RVESV'};
flow_metrics = {'CO_Lmin','SVR','PVR'};

pressure_rel_thresh = 0.12;
pressure_scale = 35;
mean_pressure_rel_thresh = 0.08;
mean_pressure_scale = 55;

volume_rel_thresh = 0.15;
volume_scale = 30;
volume_upper_ratio = 1.25;
volume_upper_scale = 60;

flow_rel_thresh = 0.20;
flow_scale = 15;

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
end
if ismember(metric_name, flow_metrics)
    penalty = penalty + soft_excess_penalty(err_rel, flow_rel_thresh, flow_scale);
end
end

function penalty = soft_excess_penalty(err_rel, threshold, scale)
excess = max(0, err_rel - threshold);
penalty = scale * (excess / threshold)^2;
end

end

function [metrics, sim, penalty] = evaluate_metrics(params, x, clinical, calib, scenario, pce_surrogate)
metrics = struct();
sim = [];
penalty = 0;

use_surrogate = ~isempty(pce_surrogate) && can_use_surrogate(calib.metricFields, pce_surrogate);
if use_surrogate
    x_full = pce_surrogate.cfg.x0;
    for idx_full = 1:numel(pce_surrogate.cfg.names)
        idx = find(strcmp(calib.names, pce_surrogate.cfg.names{idx_full}), 1, 'first');
        if ~isempty(idx)
            x_full(idx_full) = x(idx);
        end
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

function params = set_param_by_name(params, name, value)
parts = strsplit(name, '.');
switch numel(parts)
    case 2
        params.(parts{1}).(parts{2}) = value;
    case 3
        params.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_param_by_name:unsupportedDepth', ...
            'Only 2- or 3-level notation supported. Got: %s', name);
end
end
