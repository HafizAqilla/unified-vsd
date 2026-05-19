function [x0, lb, ub, param_names, registry_subset] = build_calibration_vector(parameter_registry, active_params)
% BUILD_CALIBRATION_VECTOR
% -----------------------------------------------------------------------
% Builds an optimisation vector and bounds from the centralized parameter
% registry while preserving the requested parameter order.
%
% INPUTS:
%   parameter_registry - registry table from build_parameter_registry.m   [-]
%   active_params      - ordered calibration variable names               [cellstr]
%
% OUTPUTS:
%   x0                - starting vector from registry.seeded_value        [-]
%   lb                - lower bounds from registry                        [-]
%   ub                - upper bounds from registry                        [-]
%   param_names       - parameter names in optimisation order            [cellstr]
%   registry_subset   - registry rows matching param_names               [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1 || ~istable(parameter_registry)
    error('build_calibration_vector:invalidRegistry', ...
        'parameter_registry must be a table.');
end
if nargin < 2 || isempty(active_params)
    active_params = parameter_registry.name;
end
if ischar(active_params)
    active_params = {active_params};
elseif isstring(active_params)
    active_params = cellstr(active_params);
end

param_names = active_params(:);
[is_present, loc] = ismember(param_names, parameter_registry.name);
if ~all(is_present)
    missing_names = strjoin(param_names(~is_present), ', ');
    error('build_calibration_vector:missingParameters', ...
        'Registry is missing requested parameters: %s', missing_names);
end

registry_subset = parameter_registry(loc, :);
x0 = registry_subset.seeded_value;
lb = registry_subset.lb;
ub = registry_subset.ub;
param_names = registry_subset.name;
end
