function value = get_calibration_param_value(params, reference_params, name, case_profile)
% GET_CALIBRATION_PARAM_VALUE
% -----------------------------------------------------------------------
% Reads a calibration variable value from the parameter struct.
%
% Coupled calibration variables are represented as dimensionless scales
% relative to the supplied reference parameter set.
%
% INPUTS:
%   params            - parameter struct to read                        [-]
%   reference_params  - reference parameter struct for grouped scales   [-]
%   name              - calibration variable name                       [-]
%   case_profile      - calibration governance profile                  [-]
%
% OUTPUTS:
%   value             - scalar calibration value                        [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 4
    case_profile = struct();
end
if nargin < 2 || isempty(reference_params)
    reference_params = params;
end

group = get_coupled_parameter_group(name, case_profile);
if group.found
    scale_values = nan(numel(group.memberNames), 1);
    for member_idx = 1:numel(group.memberNames)
        member_name = group.memberNames{member_idx};
        current_value = get_nested_param(params, member_name);
        reference_value = get_nested_param(reference_params, member_name);
        scale_values(member_idx) = current_value / max(reference_value, 1e-12);
    end
    finite_scales = scale_values(isfinite(scale_values));
    if isempty(finite_scales)
        value = 1.0;
    else
        value = mean(finite_scales);
    end
    return;
end

value = get_nested_param(params, name);
end

function value = get_nested_param(params, name)
parts = strsplit(name, '.');
value = params;
for part_idx = 1:numel(parts)
    value = value.(parts{part_idx});
end
end
