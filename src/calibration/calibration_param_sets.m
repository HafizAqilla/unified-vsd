function calib = calibration_param_sets(scenario, params0, optMask, primaryMetrics)
% CALIBRATION_PARAM_SETS
% -----------------------------------------------------------------------
% Returns staged calibration configuration for pediatric VSD fitting.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  2.0
% -----------------------------------------------------------------------

if nargin < 3 || isempty(optMask)
    optMask = [];
end
if nargin < 4 || isempty(primaryMetrics)
    primaryMetrics = {};
end
if ischar(primaryMetrics)
    primaryMetrics = {primaryMetrics};
elseif isstring(primaryMetrics)
    primaryMetrics = cellstr(primaryMetrics);
end

calib = struct();
calib.primaryMetrics = primaryMetrics(:)';
calib.primaryTarget = params0.calibration.primary_target_pct / 100;
calib.secondaryTarget = params0.calibration.secondary_target_pct / 100;
calib.secondaryLambda = params0.calibration.secondary_lambda;
calib.invalidPenaltyScale = params0.calibration.invalid_penalty_scale;
calib.regLambda = 0;

switch scenario
    case 'pre_surgery'
        calib.names_all = {
            'R.SAR'
            'R.SC'
            'R.SVEN'
            'R.PAR'
            'R.PCOX'
            'R.PVEN'
            'C.SAR'
            'C.SVEN'
            'C.PAR'
            'C.PVEN'
            'E.LV.EA'
            'E.LV.EB'
            'E.RV.EA'
            'E.RV.EB'
            'E.LA.EA'
            'E.RA.EA'
            'V0.LV'
            'V0.RV'
            'R.vsd'
            };
        calib.metricFields = {
            'RAP_mean'
            'LAP_mean'
            'PAP_min'
            'PAP_max'
            'PAP_mean'
            'SAP_min'
            'SAP_max'
            'SAP_mean'
            'QpQs'
            'SVR'
            'PVR'
            'CO_Lmin'
            'LVEDV'
            'LVESV'
            'RVEDV'
            'RVESV'
            'LVEF'
            'RVEF'
            };
        calib.stageA.names = {'R.SAR','R.SC','R.SVEN','R.PAR','R.PCOX','R.PVEN', ...
                      'C.SAR','C.SVEN','C.PAR','C.PVEN','R.vsd'};
        calib.stageA.metricFields = {'RAP_mean','LAP_mean','SAP_mean','PAP_mean','SVR','PVR','QpQs','CO_Lmin'};
        calib.stageB.names = {'E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB','E.LA.EA','E.RA.EA','V0.LV','V0.RV'};
        calib.stageB.metricFields = {'RAP_mean','LAP_mean','LVEDV','LVESV','LVEF','RVEDV','RVESV','RVEF'};

        if isfield(params0, 'vsd') && isfield(params0.vsd, 'mode')
            if strcmpi(params0.vsd.mode, 'orifice_bidirectional')
                calib.names_all = calib.names_all(~strcmp(calib.names_all, 'R.vsd'));
                calib.stageA.names = calib.stageA.names(~strcmp(calib.stageA.names, 'R.vsd'));
            end
        end

    case 'post_surgery'
        calib.names_all = {
            'R.SAR'
            'R.SC'
            'R.SVEN'
            'R.PAR'
            'R.PCOX'
            'R.PVEN'
            'C.SAR'
            'C.SVEN'
            'C.PAR'
            'C.PVEN'
            'E.LV.EA'
            'E.LV.EB'
            'E.RV.EA'
            'E.RV.EB'
            'V0.LV'
            'V0.RV'
            };
        calib.metricFields = {
            'RAP_mean'
            'PAP_mean'
            'SAP_min'
            'SAP_max'
            'SAP_mean'
            'QpQs'
            'SVR'
            'PVR'
            'CO_Lmin'
            'LVEDV'
            'RVEDV'
            'LVEF'
            'RVEF'
            };
        calib.stageA.names = {'R.SAR','R.SC','R.SVEN','R.PAR','R.PCOX','R.PVEN', ...
                      'C.SAR','C.SVEN','C.PAR','C.PVEN'};
        calib.stageA.metricFields = {'RAP_mean','SAP_mean','PAP_mean','SVR','PVR','QpQs','CO_Lmin'};
        calib.stageB.names = {'E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB','V0.LV','V0.RV'};
        calib.stageB.metricFields = {'RAP_mean','LVEDV','LVEF','RVEDV','RVEF'};

    otherwise
        error('calibration_param_sets:unknownScenario', ...
            'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

calib.weights = make_metric_weights(calib.metricFields, calib.primaryMetrics);

calib.x0_all = pack_x(params0, calib.names_all);
[calib.lb_all, calib.ub_all] = bounds_from_names(params0, calib.names_all, scenario);

if isempty(optMask)
    optMask = true(numel(calib.names_all), 1);
end
optMask = logical(optMask(:));
if numel(optMask) ~= numel(calib.names_all)
    error('calibration_param_sets:maskSizeMismatch', ...
        'Mask length (%d) must match parameter count (%d).', ...
        numel(optMask), numel(calib.names_all));
end
if ~any(optMask)
    error('calibration_param_sets:emptyActiveSet', ...
        'Mask deactivates all parameters.');
end

calib.mask = optMask;
calib.names = calib.names_all(calib.mask);
calib.x0 = calib.x0_all(calib.mask);
calib.lb = calib.lb_all(calib.mask);
calib.ub = calib.ub_all(calib.mask);

end

function weights = make_metric_weights(metric_fields, primary_metrics)
weights = struct();
for k = 1:numel(metric_fields)
    weights.(metric_fields{k}) = 1.0;
end
for k = 1:numel(primary_metrics)
    mf = primary_metrics{k};
    if ismember(mf, metric_fields)
        weights.(mf) = 1.5;
    end
end
if isfield(weights, 'QpQs'), weights.QpQs = max(weights.QpQs, 1.5); end
if isfield(weights, 'PAP_mean'), weights.PAP_mean = max(weights.PAP_mean, 1.4); end
if isfield(weights, 'PVR'), weights.PVR = max(weights.PVR, 1.4); end
if isfield(weights, 'SAP_mean'), weights.SAP_mean = max(weights.SAP_mean, 1.3); end
if isfield(weights, 'RAP_mean'), weights.RAP_mean = max(weights.RAP_mean, 1.5); end
if isfield(weights, 'LAP_mean'), weights.LAP_mean = max(weights.LAP_mean, 1.5); end
if isfield(weights, 'LVEDV'), weights.LVEDV = max(weights.LVEDV, 1.4); end
if isfield(weights, 'RVEDV'), weights.RVEDV = max(weights.RVEDV, 1.4); end
end

function x = pack_x(params, names)
x = zeros(numel(names), 1);
for i = 1:numel(names)
    x(i) = get_param_by_name(params, names{i});
end
end

function [lb, ub] = bounds_from_names(params, names, scenario)
lb = zeros(numel(names), 1);
ub = zeros(numel(names), 1);

for i = 1:numel(names)
    x0 = get_param_by_name(params, names{i});
    nm = names{i};

    if startsWith(nm, 'R.')
        if strcmp(nm, 'R.vsd') && strcmp(scenario, 'pre_surgery')
            lb(i) = max(0.001, 0.05 * x0);
            ub(i) = min(500, 20 * max(x0, 0.01));
        else
            lb(i) = 0.25 * x0;
            ub(i) = 4.00 * x0;
        end
    elseif startsWith(nm, 'C.')
        lb(i) = 0.30 * x0;
        ub(i) = 6.00 * x0;
    elseif startsWith(nm, 'E.')
        lb(i) = 0.30 * x0;
        ub(i) = 4.00 * x0;
    elseif startsWith(nm, 'V0.')
        lb(i) = 0.30 * x0;
        ub(i) = 3.00 * x0;
    else
        lb(i) = 0.30 * x0;
        ub(i) = 4.00 * x0;
    end
end
end

function v = get_param_by_name(params, name)
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts)
    v = v.(parts{k});
end
end
