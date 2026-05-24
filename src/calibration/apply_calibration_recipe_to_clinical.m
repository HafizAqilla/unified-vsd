function clinical = apply_calibration_recipe_to_clinical(clinical, scenario, recipe)
% APPLY_CALIBRATION_RECIPE_TO_CLINICAL
% -----------------------------------------------------------------------
% Apply explicit recipe-controlled patient/scenario evidence overrides.
%
% INPUTS:
%   clinical - unified clinical struct                                  [-]
%   scenario - scenario string                                          [-]
%   recipe   - calibration recipe struct                                [-]
%
% OUTPUTS:
%   clinical - clinical struct with recipe evidence applied             [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-24
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3 || ~isstruct(recipe) || isempty(fieldnames(recipe))
    return;
end

if isfield(recipe, 'demographics') && isstruct(recipe.demographics)
    clinical.common = merge_struct_fields(clinical.common, recipe.demographics);
end

override_field = sprintf('%s_overrides', char(scenario));
if isfield(recipe, override_field) && isstruct(recipe.(override_field))
    if ~isfield(clinical, scenario) || ~isstruct(clinical.(scenario))
        clinical.(scenario) = struct();
    end
    clinical.(scenario) = merge_struct_fields( ...
        clinical.(scenario), recipe.(override_field));
end

if ~isfield(clinical, 'calibration_recipes') || ~isstruct(clinical.calibration_recipes)
    clinical.calibration_recipes = struct();
end
clinical.calibration_recipes.(char(scenario)) = recipe;
end

function dst = merge_struct_fields(dst, src)
names = fieldnames(src);
for idx = 1:numel(names)
    dst.(names{idx}) = src.(names{idx});
end
end
