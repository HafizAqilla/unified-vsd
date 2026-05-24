function params = apply_calibration_recipe_to_params(params, reference_params, recipe, case_profile)
% APPLY_CALIBRATION_RECIPE_TO_PARAMS
% -----------------------------------------------------------------------
% Apply explicit recipe-controlled parameter and initial-condition seeds.
%
% INPUTS:
%   params           - parameter struct after clinical seeding           [-]
%   reference_params - scaled reference parameters for grouped updates   [-]
%   recipe           - calibration recipe struct                         [-]
%   case_profile     - calibration governance profile                    [-]
%
% OUTPUTS:
%   params           - parameter struct with recipe seeds applied        [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-24
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 4
    case_profile = struct();
end
if nargin < 3 || ~isstruct(recipe) || isempty(fieldnames(recipe))
    return;
end
if nargin < 2 || isempty(reference_params)
    reference_params = params;
end

if isfield(recipe, 'fixed_parameter_values') && ...
        isstruct(recipe.fixed_parameter_values)
    params = apply_fixed_parameter_values( ...
        params, reference_params, recipe.fixed_parameter_values, case_profile);
end

if isfield(recipe, 'initial_conditions') && ...
        isstruct(recipe.initial_conditions) && ...
        isfield(recipe.initial_conditions, 'V')
    if ~isfield(params, 'ic') || ~isstruct(params.ic)
        params.ic = struct();
    end
    params.ic.V = recipe.initial_conditions.V(:);
end
end

function params = apply_fixed_parameter_values(params, reference_params, fixed_values, case_profile)
if ~isfield(fixed_values, 'names') || ~isfield(fixed_values, 'values')
    return;
end
names = fixed_values.names(:);
values = fixed_values.values(:);
if numel(names) ~= numel(values)
    error('apply_calibration_recipe_to_params:invalidFixedValues', ...
        'Fixed-parameter recipe names and values must have the same length.');
end
for idx = 1:numel(names)
    params = set_calibration_param_value( ...
        params, reference_params, char(names{idx}), values(idx), case_profile);
end
end
