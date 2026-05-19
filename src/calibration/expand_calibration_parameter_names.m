function member_names = expand_calibration_parameter_names(parameter_name, case_profile)
% EXPAND_CALIBRATION_PARAMETER_NAMES
% -----------------------------------------------------------------------
% Expands a calibration variable into the physical parameter names it
% controls. Ungrouped variables expand to themselves.
%
% INPUTS:
%   parameter_name - calibration variable name                          [-]
%   case_profile   - calibration governance profile                     [-]
%
% OUTPUTS:
%   member_names   - physical parameter names                           [cell]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.0
% -----------------------------------------------------------------------

group = get_coupled_parameter_group(parameter_name, case_profile);
if group.found
    member_names = group.memberNames(:)';
else
    member_names = {parameter_name};
end
end
