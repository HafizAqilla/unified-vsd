function [recipe, found] = load_calibration_recipe(clinical, scenario)
% LOAD_CALIBRATION_RECIPE
% -----------------------------------------------------------------------
% Load an explicit patient-scenario calibration recipe when one exists.
%
% INPUTS:
%   clinical - unified clinical struct                                  [-]
%   scenario - scenario string                                          [-]
%
% OUTPUTS:
%   recipe   - calibration recipe struct, or empty struct               [-]
%   found    - true when a matching recipe exists                       [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-24
% VERSION:  1.0
% -----------------------------------------------------------------------

recipe = struct();
found = false;

if nargin < 2 || isempty(scenario) || nargin < 1 || ~isstruct(clinical)
    return;
end

patient_label = resolve_patient_label(clinical);
if isempty(patient_label)
    return;
end

recipe_name = sprintf('%s_%s', patient_label, char(scenario));
recipe_name = matlab.lang.makeValidName(lower(recipe_name));
if exist(recipe_name, 'file') ~= 2
    return;
end

recipe = feval(recipe_name);
found = isstruct(recipe) && isfield(recipe, 'id') && ...
    isfield(recipe, 'scenario') && strcmpi(char(recipe.scenario), char(scenario));
end

function label = resolve_patient_label(clinical)
label = '';
if isfield(clinical, 'common')
    if isfield(clinical.common, 'patient_name') && ~isempty(clinical.common.patient_name)
        label = char(clinical.common.patient_name);
    elseif isfield(clinical.common, 'patient_id') && ~isempty(clinical.common.patient_id)
        label = char(clinical.common.patient_id);
    end
end
label = regexprep(lower(label), '[^a-z0-9]+', '_');
label = regexprep(label, '^_+|_+$', '');
end
