function params = enforce_vascular_rc_coupling(params, reference_params, case_profile)
% ENFORCE_VASCULAR_RC_COUPLING
% -----------------------------------------------------------------------
% Recomputes vascular compliances from current vascular resistances using
% the BSA-scaled baseline as the anchoring reference.
%
% INPUTS:
%   params            - parameter struct to update                     [-]
%   reference_params  - BSA-scaled reference parameter struct          [-]
%   case_profile      - calibration governance profile                 [-]
%
% OUTPUTS:
%   params            - parameter struct with coupled vascular C terms  [-]
%
% REFERENCES:
%   [1] Kung et al. (2013). Journal of Biomechanics 46:423-429.
%   [2] Pennati and Fumero (2000). Annals of Biomedical Engineering.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3
    case_profile = struct();
end
if nargin < 2 || isempty(reference_params)
    reference_params = params;
end

rules = get_vascular_rc_coupling_rules(case_profile);
if isempty(rules)
    return;
end

for rule_idx = 1:numel(rules)
    rule = rules(rule_idx);
    R_ref = get_nested_param(reference_params, rule.resistanceName);
    C_ref = get_nested_param(reference_params, rule.complianceName);
    R_cur = get_nested_param(params, rule.resistanceName);

    if ~isfinite(R_ref) || ~isfinite(C_ref) || ~isfinite(R_cur) || R_ref <= 0 || C_ref <= 0 || R_cur <= 0
        continue;
    end

    resistance_scale = R_cur / max(R_ref, 1e-12);
    compliance_scale = resistance_scale ^ rule.exponent;
    compliance_scale = min(max(compliance_scale, rule.minScale), rule.maxScale);
    params = set_nested_param(params, rule.complianceName, C_ref * compliance_scale);
end

params.calibration.vascular_rc_coupling_source = rules(1).source;
params.calibration.vascular_rc_coupling_enabled = true;
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
        error('enforce_vascular_rc_coupling:unsupportedDepth', ...
            'Only 2- or 3-level dot notation supported. Got: %s', name);
end
end
