function calib = calibration_param_sets(scenario, params0, optMask, primaryMetrics, caseProfile, registryContext)
% CALIBRATION_PARAM_SETS
% -----------------------------------------------------------------------
% Returns staged calibration configuration for pediatric VSD fitting.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  3.0
% -----------------------------------------------------------------------

if nargin < 3 || isempty(optMask)
    optMask = [];
end
if nargin < 4 || isempty(primaryMetrics)
    primaryMetrics = {};
end
if nargin < 5 || isempty(caseProfile)
    caseProfile = struct();
end
if nargin < 6 || isempty(registryContext)
    registryContext = struct();
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
calib.paramPlausibilityLambda = 0.5;
calib.boundaryPlausibilityLambda = 20.0;
calib.caseProfile = caseProfile;
if isfield(caseProfile, 'targetTiers')
    calib.targetTiers = caseProfile.targetTiers;
else
    calib.targetTiers = struct();
end

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
            'C.PAR'
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
        calib.stageA.names = {'R.SAR','R.SC','R.SVEN','R.PAR','R.PCOX','R.PVEN','C.SAR','C.PAR','R.vsd'};
        calib.stageA.metricFields = {'RAP_mean','LAP_mean','SAP_min','SAP_max','SAP_mean', ...
            'PAP_mean','SVR','PVR','QpQs','CO_Lmin'};
        calib.stageB.names = {'E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB','E.LA.EA','E.RA.EA','V0.LV','V0.RV'};
        calib.stageB.metricFields = {'RAP_mean','LAP_mean','LVEDV','LVESV','LVEF','RVEDV','RVESV','RVEF'};

        if isfield(params0, 'vsd') && isfield(params0.vsd, 'mode')
            if strcmpi(params0.vsd.mode, 'orifice_bidirectional')
                calib.names_all(strcmp(calib.names_all, 'R.vsd')) = {'vsd.Cd'};
                calib.stageA.names(strcmp(calib.stageA.names, 'R.vsd')) = {'vsd.Cd'};
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
            'C.PAR'
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
        calib.stageA.names = {'R.SAR','R.SC','R.SVEN','R.PAR','R.PCOX','R.PVEN','C.SAR','C.PAR'};
        calib.stageA.metricFields = {'RAP_mean','SAP_min','SAP_max','SAP_mean', ...
            'PAP_mean','SVR','PVR','QpQs','CO_Lmin'};
        calib.stageB.names = {'E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB','V0.LV','V0.RV'};
        calib.stageB.metricFields = {'RAP_mean','LVEDV','LVEF','RVEDV','RVEF'};

    otherwise
        error('calibration_param_sets:unknownScenario', ...
            'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

calib = apply_case_profile_coupling(calib, caseProfile);
calib = apply_case_profile_metrics(calib, caseProfile);
if isfield(caseProfile, 'regLambda') && isfinite(caseProfile.regLambda)
    calib.regLambda = caseProfile.regLambda;
end
if isfield(caseProfile, 'paramPlausibilityLambdaScale') && ...
        isfinite(caseProfile.paramPlausibilityLambdaScale)
    calib.paramPlausibilityLambda = calib.paramPlausibilityLambda * ...
        caseProfile.paramPlausibilityLambdaScale;
end
if isfield(caseProfile, 'boundaryPlausibilityLambdaScale') && ...
        isfinite(caseProfile.boundaryPlausibilityLambdaScale)
    calib.boundaryPlausibilityLambda = calib.boundaryPlausibilityLambda * ...
        caseProfile.boundaryPlausibilityLambdaScale;
end
calib.primaryMetrics = intersect(calib.primaryMetrics, calib.metricFields, 'stable');
calib.weights = make_metric_weights(calib.metricFields, calib.primaryMetrics, caseProfile);
calib.referenceParams = resolve_group_reference_params(params0, registryContext);
calib.registryContext = registryContext;
calib.parameterRegistry = build_registry(params0, scenario, caseProfile, calib.names_all, registryContext);
validate_bounds(calib.parameterRegistry, scenario);
[calib.x0_all, calib.lb_all, calib.ub_all, calib.names_all, calib.parameterRegistry] = ...
    build_calibration_vector(calib.parameterRegistry, calib.names_all);

if isempty(optMask)
    optMask = true(numel(calib.names_all), 1);
end
optMask = logical(optMask(:));
if numel(optMask) ~= numel(calib.names_all)
    error('calibration_param_sets:maskSizeMismatch', ...
        'Mask length (%d) must match parameter count (%d).', ...
        numel(optMask), numel(calib.names_all));
end
optMask = apply_case_profile_parameter_mask(optMask, calib.names_all, caseProfile);
optMask = apply_case_profile_required_mask(optMask, calib.names_all, caseProfile);
if ~any(optMask)
    error('calibration_param_sets:emptyActiveSet', ...
        'Mask deactivates all parameters.');
end

calib.mask = optMask;
calib.names = calib.names_all(calib.mask);
calib.x0 = calib.x0_all(calib.mask);
calib.lb = calib.lb_all(calib.mask);
calib.ub = calib.ub_all(calib.mask);
calib.parameterRegistryActive = calib.parameterRegistry(calib.mask, :);

end

function calib = apply_case_profile_coupling(calib, case_profile)
if ~isfield(case_profile, 'coupledParameterGroups') || isempty(case_profile.coupledParameterGroups)
    return;
end
calib.names_all = collapse_coupled_parameter_names(calib.names_all, case_profile);
calib.stageA.names = collapse_coupled_parameter_names(calib.stageA.names, case_profile);
calib.stageB.names = collapse_coupled_parameter_names(calib.stageB.names, case_profile);
end

function calib = apply_case_profile_metrics(calib, case_profile)
if ~isfield(case_profile, 'allowedMetricFields') || isempty(case_profile.allowedMetricFields)
    return;
end
allowed = case_profile.allowedMetricFields(:)';
calib.metricFields = intersect(calib.metricFields, allowed, 'stable');
calib.stageA.metricFields = intersect(calib.stageA.metricFields, allowed, 'stable');
calib.stageB.metricFields = intersect(calib.stageB.metricFields, allowed, 'stable');
end

function optMask = apply_case_profile_parameter_mask(optMask, names_all, case_profile)
if ~isfield(case_profile, 'allowedFreeParameters') || isempty(case_profile.allowedFreeParameters)
    return;
end
allowed_mask = ismember(names_all(:), case_profile.allowedFreeParameters(:));
optMask = optMask & allowed_mask;
if ~any(optMask) && any(allowed_mask)
    first_allowed = find(allowed_mask, 1, 'first');
    optMask(first_allowed) = true;
end
end

function optMask = apply_case_profile_required_mask(optMask, names_all, case_profile)
if isfield(case_profile, 'systemicPolishEnabled') && case_profile.systemicPolishEnabled && ...
        isfield(case_profile, 'systemicPolishNames') && ~isempty(case_profile.systemicPolishNames)
    required_mask = ismember(names_all(:), case_profile.systemicPolishNames(:));
    optMask = optMask | required_mask;
    return;
end
end

function weights = make_metric_weights(metric_fields, primary_metrics, case_profile)
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
if isfield(weights, 'SAP_min'), weights.SAP_min = max(weights.SAP_min, 1.3); end
if isfield(weights, 'SAP_max'), weights.SAP_max = max(weights.SAP_max, 1.4); end
if isfield(weights, 'CO_Lmin'), weights.CO_Lmin = max(weights.CO_Lmin, 1.8); end
if isfield(weights, 'SVR'), weights.SVR = min(weights.SVR, 1.0); end
if isfield(weights, 'RAP_mean'), weights.RAP_mean = max(weights.RAP_mean, 1.5); end
if isfield(weights, 'LAP_mean'), weights.LAP_mean = max(weights.LAP_mean, 1.5); end
if isfield(weights, 'LVEDV'), weights.LVEDV = max(weights.LVEDV, 1.4); end
if isfield(weights, 'RVEDV'), weights.RVEDV = max(weights.RVEDV, 1.4); end

if isfield(case_profile, 'metricWeightOverrides') && isstruct(case_profile.metricWeightOverrides)
    override_names = fieldnames(case_profile.metricWeightOverrides);
    for k = 1:numel(override_names)
        mf = override_names{k};
        if isfield(weights, mf)
            weights.(mf) = weights.(mf) * case_profile.metricWeightOverrides.(mf);
        end
    end
end

if isfield(case_profile, 'targetTiers') && isstruct(case_profile.targetTiers)
    weights = apply_target_tier_weight_multipliers(weights, case_profile.targetTiers);
end
end

function weights = apply_target_tier_weight_multipliers(weights, target_tiers)
if ~isfield(target_tiers, 'table') || isempty(target_tiers.table)
    return;
end
tier_tbl = target_tiers.table;
for idx = 1:height(tier_tbl)
    metric_name = tier_tbl.Metric{idx};
    if ~isfield(weights, metric_name)
        continue;
    end
    switch tier_tbl.Tier{idx}
        case 'consistency_check_only'
            weights.(metric_name) = 0;
    end
end

if isfield(target_tiers, 'weights') && isstruct(target_tiers.weights)
    names = fieldnames(target_tiers.weights);
    for idx = 1:numel(names)
        metric_name = names{idx};
        if isfield(weights, metric_name)
            weights.(metric_name) = weights.(metric_name) * target_tiers.weights.(metric_name);
        end
    end
end
end

function registry = build_registry(params0, scenario, case_profile, names_all, registry_context)
params_adult = params0;
params_scaled = params0;
if isstruct(registry_context)
    if isfield(registry_context, 'params_adult') && ~isempty(registry_context.params_adult)
        params_adult = registry_context.params_adult;
    end
    if isfield(registry_context, 'params_scaled') && ~isempty(registry_context.params_scaled)
        params_scaled = registry_context.params_scaled;
    end
end

registry = build_parameter_registry( ...
    params_adult, params_scaled, params0, scenario, case_profile, names_all);
end

function reference_params = resolve_group_reference_params(params0, registry_context)
reference_params = params0;
if isstruct(registry_context) && isfield(registry_context, 'params_scaled') && ...
        ~isempty(registry_context.params_scaled)
    reference_params = registry_context.params_scaled;
end
end
