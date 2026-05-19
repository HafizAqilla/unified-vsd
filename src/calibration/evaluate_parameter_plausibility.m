function results = evaluate_parameter_plausibility(fitted_params, parameter_registry)
% EVALUATE_PARAMETER_PLAUSIBILITY
% -----------------------------------------------------------------------
% Evaluates calibrated parameter values against the centralized registry.
%
% INPUTS:
%   fitted_params       - fitted calibration vector in registry order     [-]
%   parameter_registry  - registry table from build_parameter_registry    [-]
%
% OUTPUTS:
%   results             - struct with fields:
%       .table               per-parameter plausibility table             [-]
%       .n_ok                number of OK parameters                      [-]
%       .n_warning           number of WARNING parameters                 [-]
%       .n_fail              number of FAIL parameters                    [-]
%       .warning_fraction    warning fraction among all parameters        [-]
%       .all_within_bounds   true if all fitted values lie in bounds      [-]
%
% FLAG LOGIC:
%   OK:
%       inside bounds and not near a bound
%   WARNING:
%       inside bounds but within 10% of lower or upper bound
%   FAIL:
%       outside bounds, non-finite, or wrong sign for positive-definite
%       parameter classes
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || ~istable(parameter_registry)
    error('evaluate_parameter_plausibility:invalidRegistry', ...
        'parameter_registry must be supplied as a table.');
end
if nargin < 1
    error('evaluate_parameter_plausibility:missingFittedParams', ...
        'fitted_params must be supplied.');
end

fitted_vec = fitted_params(:);
if numel(fitted_vec) ~= height(parameter_registry)
    error('evaluate_parameter_plausibility:sizeMismatch', ...
        'fitted_params length (%d) must match registry height (%d).', ...
        numel(fitted_vec), height(parameter_registry));
end

n_params = height(parameter_registry);
baseline_scaled = parameter_registry.baseline_scaled;
lb = parameter_registry.lb;
ub = parameter_registry.ub;
span = max(ub - lb, 1e-12);

absolute_change = fitted_vec - baseline_scaled;
ratio_to_baseline = fitted_vec ./ max(abs(baseline_scaled), 1e-12);
ratio_to_baseline(~isfinite(baseline_scaled) | abs(baseline_scaled) < 1e-12) = NaN;
log_ratio_to_baseline = nan(n_params, 1);
positive_ratio = ratio_to_baseline > 0 & isfinite(ratio_to_baseline);
log_ratio_to_baseline(positive_ratio) = log(ratio_to_baseline(positive_ratio));

within_bounds = fitted_vec >= lb & fitted_vec <= ub;
distance_to_lower_bound = (fitted_vec - lb) ./ span;
distance_to_upper_bound = (ub - fitted_vec) ./ span;
near_lower_bound = within_bounds & distance_to_lower_bound <= 0.10;
near_upper_bound = within_bounds & distance_to_upper_bound <= 0.10;

plausibility_flag = strings(n_params, 1);
for idx = 1:n_params
    name = parameter_registry.name{idx};
    fitted_value = fitted_vec(idx);
    positive_definite = startsWith(name, 'R.') || startsWith(name, 'C.') || ...
        startsWith(name, 'L.') || startsWith(name, 'E.') || ...
        startsWith(name, 'V0.') || strcmp(name, 'vsd.Cd') || ...
        startsWith(name, 'group.');

    if ~isfinite(fitted_value) || (positive_definite && fitted_value <= 0) || ~within_bounds(idx)
        plausibility_flag(idx) = "FAIL";
    elseif near_lower_bound(idx) || near_upper_bound(idx)
        plausibility_flag(idx) = "WARNING";
    else
        plausibility_flag(idx) = "OK";
    end
end

results = struct();
results.table = table( ...
    parameter_registry.name, parameter_registry.group, parameter_registry.unit, ...
    baseline_scaled, parameter_registry.seeded_value, fitted_vec, absolute_change, ...
    ratio_to_baseline, log_ratio_to_baseline, within_bounds, ...
    distance_to_lower_bound, distance_to_upper_bound, near_lower_bound, ...
    near_upper_bound, plausibility_flag, ...
    'VariableNames', {'Parameter','Group','Unit','BaselineScaled','SeededValue', ...
    'FittedValue','AbsoluteChange','RatioToBaseline','LogRatioToBaseline', ...
    'WithinBounds','DistanceToLowerBound','DistanceToUpperBound', ...
    'NearLowerBound','NearUpperBound','PlausibilityFlag'});

results.n_ok = sum(results.table.PlausibilityFlag == "OK");
results.n_warning = sum(results.table.PlausibilityFlag == "WARNING");
results.n_fail = sum(results.table.PlausibilityFlag == "FAIL");
results.warning_fraction = results.n_warning / max(n_params, 1);
results.all_within_bounds = all(results.table.WithinBounds);
end
