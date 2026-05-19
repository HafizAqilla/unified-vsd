function [params, warm_start] = apply_post_surgery_warm_start(params, clinical, scenario)
% APPLY_POST_SURGERY_WARM_START
% -----------------------------------------------------------------------
% Copies calibrated pre-surgery cardiovascular parameters into a post-op run.
%
% Longer description: surgical VSD closure changes shunt geometry, but it
% does not instantly reset the patient's calibrated ventricular, vascular,
% and compliance parameters. This utility transfers those patient-specific
% parameters from clinical.pre_surgery.CalibParams into the post-surgery
% parameter seed before params_from_clinical closes the VSD.
%
% INPUTS:
%   params    - parameter struct after demographic scaling              [-]
%   clinical  - clinical struct with optional pre_surgery.CalibParams   [-]
%   scenario  - simulation scenario string                              [-]
%
% OUTPUTS:
%   params      - parameter struct with warm-start fields copied        [-]
%   warm_start  - audit struct describing copied fields                 [-]
%
% ASSUMPTIONS:
%   - VSD closure is handled later by params_from_clinical.
%   - R.vsd is deliberately excluded to avoid reopening the shunt.
%
% SIGN CONVENTIONS:
%   - No flows or pressure gradients are computed here.
%
% REFERENCES:
%   [1] docs/theory_notes.md. VSD scenario transition assumptions.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-18
% VERSION:  1.0
% -----------------------------------------------------------------------

warm_start = struct('applied', false, 'source', ...
    'clinical.pre_surgery.CalibParams', 'copiedFields', {{}});

if ~strcmp(scenario, 'post_surgery') || ~has_calib_params(clinical)
    return;
end

params_pre = clinical.pre_surgery.CalibParams;
copy_paths = {
    'E.LV.EA', 'E.LV.EB', 'E.RV.EA', 'E.RV.EB', ...
    'E.LA.EA', 'E.LA.EB', 'E.RA.EA', 'E.RA.EB', ...
    'V0.LV', 'V0.RV', 'V0.LA', 'V0.RA', ...
    'R.SAR', 'R.SC', 'R.SVEN', 'R.PAR', 'R.PCOX', 'R.PCNO', 'R.PVEN', ...
    'C.SAR', 'C.SVEN', 'C.PAR', 'C.PVEN'};

for path_idx = 1:numel(copy_paths)
    [params, did_copy] = copy_numeric_scalar(params, params_pre, copy_paths{path_idx});
    if did_copy
        warm_start.copiedFields{end + 1} = copy_paths{path_idx};
    end
end

warm_start.applied = ~isempty(warm_start.copiedFields);

end

function tf = has_calib_params(clinical)
% HAS_CALIB_PARAMS - true when a non-empty pre-op parameter seed exists.
tf = isstruct(clinical) && isfield(clinical, 'pre_surgery') && ...
    isfield(clinical.pre_surgery, 'CalibParams') && ...
    isstruct(clinical.pre_surgery.CalibParams) && ...
    ~isempty(fieldnames(clinical.pre_surgery.CalibParams));
end

function [target, did_copy] = copy_numeric_scalar(target, source, path_text)
% COPY_NUMERIC_SCALAR - copy a finite scalar at dot path if present in both structs.
did_copy = false;
path_parts = strsplit(path_text, '.');
if ~has_path(source, path_parts) || ~has_path(target, path_parts)
    return;
end
value = get_path_value(source, path_parts);
if ~(isnumeric(value) && isscalar(value) && isfinite(value))
    return;
end
target = set_path_value(target, path_parts, value);
did_copy = true;
end

function tf = has_path(s, path_parts)
% HAS_PATH - verify that all nested fields exist.
tf = isstruct(s);
for part_idx = 1:numel(path_parts)
    if ~tf || ~isfield(s, path_parts{part_idx})
        tf = false;
        return;
    end
    s = s.(path_parts{part_idx});
end
end

function value = get_path_value(s, path_parts)
% GET_PATH_VALUE - read a nested field value from a validated dot path.
value = s;
for part_idx = 1:numel(path_parts)
    value = value.(path_parts{part_idx});
end
end

function s = set_path_value(s, path_parts, value)
% SET_PATH_VALUE - set a two- or three-level nested field value.
switch numel(path_parts)
    case 2
        s.(path_parts{1}).(path_parts{2}) = value;
    case 3
        s.(path_parts{1}).(path_parts{2}).(path_parts{3}) = value;
    otherwise
        error('apply_post_surgery_warm_start:unsupportedPath', ...
            'Only two- or three-level parameter paths are supported: %s', ...
            strjoin(path_parts, '.'));
end
end
