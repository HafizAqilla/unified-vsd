function params = set_calibration_param_value(params, reference_params, name, value, case_profile)
% SET_CALIBRATION_PARAM_VALUE
% -----------------------------------------------------------------------
% Writes a calibration variable into the parameter struct.
%
% Coupled calibration variables apply the same dimensionless scale to all
% physical members relative to the supplied reference parameter set.
%
% INPUTS:
%   params            - parameter struct to modify                      [-]
%   reference_params  - reference parameter struct for grouped scales   [-]
%   name              - calibration variable name                       [-]
%   value             - scalar calibration value                        [-]
%   case_profile      - calibration governance profile                  [-]
%
% OUTPUTS:
%   params            - updated parameter struct                        [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 5
    case_profile = struct();
end
if nargin < 2 || isempty(reference_params)
    reference_params = params;
end

group = get_coupled_parameter_group(name, case_profile);
if group.found
    for member_idx = 1:numel(group.memberNames)
        member_name = group.memberNames{member_idx};
        reference_value = get_nested_param(reference_params, member_name);
        params = set_nested_param(params, member_name, reference_value * value);
    end
    params = enforce_vascular_rc_coupling(params, reference_params, case_profile);
    return;
end

params = set_nested_param(params, name, value);
params = enforce_vascular_rc_coupling(params, reference_params, case_profile);
end

function value = get_nested_param(params, name)
parts = strsplit(name, '.');
value = params;
for part_idx = 1:numel(parts)
    value = value.(parts{part_idx});
end
end

function params = set_nested_param(params, name, value)
parts = strsplit(name, '.');
switch numel(parts)
    case 2
        params.(parts{1}).(parts{2}) = value;
    case 3
        params.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_calibration_param_value:unsupportedDepth', ...
            'Only 2- or 3-level dot notation supported. Got: %s', name);
end
end
