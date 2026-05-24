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

scaling_mode = resolve_recipe_scaling_mode(params, reference_params);

if isfield(recipe, 'fixed_parameter_values') && ...
        isstruct(recipe.fixed_parameter_values) && ...
        applies_to_scaling_mode(recipe.fixed_parameter_values, scaling_mode)
    params = apply_fixed_parameter_values( ...
        params, reference_params, recipe.fixed_parameter_values, case_profile);
end

if isfield(recipe, 'initial_conditions') && ...
        isstruct(recipe.initial_conditions) && ...
        isfield(recipe.initial_conditions, 'V') && ...
        applies_to_scaling_mode(recipe.initial_conditions, scaling_mode)
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

function scaling_mode = resolve_recipe_scaling_mode(params, reference_params)
scaling_mode = '';
if isstruct(params) && isfield(params, 'scaling') && isstruct(params.scaling)
    scaling_mode = first_nonempty_scaling_mode(params.scaling);
end
if isempty(scaling_mode) && isstruct(reference_params) && ...
        isfield(reference_params, 'scaling') && isstruct(reference_params.scaling)
    scaling_mode = first_nonempty_scaling_mode(reference_params.scaling);
end
scaling_mode = normalize_scaling_mode(scaling_mode);
end

function scaling_mode = first_nonempty_scaling_mode(scaling)
scaling_mode = '';
if isfield(scaling, 'mode') && ~isempty(scaling.mode)
    scaling_mode = char(scaling.mode);
elseif isfield(scaling, 'requested_mode') && ~isempty(scaling.requested_mode)
    scaling_mode = char(scaling.requested_mode);
end
end

function tf = applies_to_scaling_mode(recipe_block, scaling_mode)
if ~isfield(recipe_block, 'scaling_modes') || isempty(recipe_block.scaling_modes)
    tf = true;
    return;
end
allowed_modes = normalize_scaling_mode_cell(recipe_block.scaling_modes);
tf = ismember(normalize_scaling_mode(scaling_mode), allowed_modes);
end

function modes = normalize_scaling_mode_cell(values)
if ischar(values)
    values = {values};
elseif isstring(values)
    values = cellstr(values);
end
modes = cell(size(values));
for idx = 1:numel(values)
    modes{idx} = normalize_scaling_mode(values{idx});
end
end

function mode = normalize_scaling_mode(mode)
mode = lower(strtrim(char(string(mode))));
switch mode
    case {'lundquist', 'lundqvist', 'bsa'}
        mode = 'lundquist_bsa';
end
end
