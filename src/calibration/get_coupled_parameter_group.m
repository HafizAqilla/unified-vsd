function group = get_coupled_parameter_group(parameter_name, case_profile)
% GET_COUPLED_PARAMETER_GROUP
% -----------------------------------------------------------------------
% Returns a coupled-parameter group definition for a calibration variable.
%
% INPUTS:
%   parameter_name - calibration variable name                          [-]
%   case_profile   - calibration governance profile                     [-]
%
% OUTPUTS:
%   group          - struct with fields:
%                    .found        logical                              [-]
%                    .name         group calibration name               [-]
%                    .memberNames  physical parameter names             [cell]
%                    .groupType    'R' | 'C' | ''                       [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.0
% -----------------------------------------------------------------------

group = struct( ...
    'found', false, ...
    'name', '', ...
    'memberNames', {{}}, ...
    'groupType', '');

if nargin < 2 || ~isstruct(case_profile) || ...
        ~isfield(case_profile, 'coupledParameterGroups') || ...
        isempty(case_profile.coupledParameterGroups)
    return;
end

group_defs = case_profile.coupledParameterGroups;
for group_idx = 1:numel(group_defs)
    if strcmp(group_defs(group_idx).name, parameter_name)
        group = group_defs(group_idx);
        group.found = true;
        return;
    end
end
end
