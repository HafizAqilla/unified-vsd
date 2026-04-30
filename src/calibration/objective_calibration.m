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

J = J_primary + calib.secondaryLambda * J_secondary + J_clinical_guard + J_reg + validity_penalty;

if ~isempty(sim) && ~sim.ss_reached
    J = J + calib.invalidPenaltyScale;
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
