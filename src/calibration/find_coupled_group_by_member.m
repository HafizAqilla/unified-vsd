function group = find_coupled_group_by_member(member_name, case_profile)
% FIND_COUPLED_GROUP_BY_MEMBER
% -----------------------------------------------------------------------
% Returns the coupled-parameter group that contains a physical parameter.
%
% INPUTS:
%   member_name    - physical parameter name                            [-]
%   case_profile   - calibration governance profile                     [-]
%
% OUTPUTS:
%   group          - struct (see get_coupled_parameter_group)           [-]
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
    if ismember(member_name, group_defs(group_idx).memberNames)
        group = group_defs(group_idx);
        group.found = true;
        return;
    end
end
end
