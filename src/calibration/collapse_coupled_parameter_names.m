function names_out = collapse_coupled_parameter_names(names_in, case_profile)
% COLLAPSE_COUPLED_PARAMETER_NAMES
% -----------------------------------------------------------------------
% Rewrites physical parameter lists into case-profile coupled variables.
%
% INPUTS:
%   names_in       - physical calibration parameter names               [cell]
%   case_profile   - calibration governance profile                     [-]
%
% OUTPUTS:
%   names_out      - calibration variable names after grouping          [cell]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.0
% -----------------------------------------------------------------------

if isempty(names_in)
    names_out = names_in;
    return;
end

names_in = names_in(:)';
names_out = {};

for name_idx = 1:numel(names_in)
    current_name = names_in{name_idx};
    group = find_coupled_group_by_member(current_name, case_profile);
    if group.found
        collapsed_name = group.name;
    else
        collapsed_name = current_name;
    end
    if ~ismember(collapsed_name, names_out)
        names_out{end + 1} = collapsed_name; %#ok<AGROW>
    end
end
end
