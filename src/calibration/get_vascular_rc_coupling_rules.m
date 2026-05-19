function rules = get_vascular_rc_coupling_rules(case_profile)
% GET_VASCULAR_RC_COUPLING_RULES
% -----------------------------------------------------------------------
% Returns vascular resistance-compliance coupling rules used to enforce
% literature-backed priors after BSA scaling.
%
% The coupling follows the resting-condition relation used by Kung et al.
% for BSA-adjusted lumped vascular beds:
%   C_i = C_i,BSA * (R_i / R_i,BSA)^(-4/3)
%
% INPUTS:
%   case_profile - calibration governance profile                      [-]
%
% OUTPUTS:
%   rules        - struct array describing R-C coupling               [-]
%
% REFERENCES:
%   [1] Kung et al. (2013). Journal of Biomechanics 46:423-429.
%   [2] Pennati and Fumero (2000). Annals of Biomedical Engineering.
%   [3] Baretta et al. (2011). Patient-specific BSA-adjusted LPN tuning.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1
    case_profile = struct();
end

enabled = true;
if isfield(case_profile, 'useVascularRcCoupling') && ~isempty(case_profile.useVascularRcCoupling)
    enabled = logical(case_profile.useVascularRcCoupling);
end
if ~enabled
    rules = struct('resistanceName', {}, 'complianceName', {}, ...
        'exponent', {}, 'minScale', {}, 'maxScale', {}, 'source', {});
    return;
end

rules = repmat(struct( ...
    'resistanceName', '', ...
    'complianceName', '', ...
    'exponent', -4/3, ...
    'minScale', 0.05, ...
    'maxScale', 8.0, ...
    'source', 'Kung2013_PennatiFumero2000_resting_RC_relation'), 4, 1);

% Large-artery compliances (C.SAR, C.PAR) are anchored separately from
% direct pulse-pressure formulas, so exact R-C coupling is applied only to
% the distributed capillary/venous beds.
rules(1).resistanceName = 'R.SC';
rules(1).complianceName = 'C.SC';

rules(2).resistanceName = 'R.SVEN';
rules(2).complianceName = 'C.SVEN';

rules(3).resistanceName = 'R.PCOX';
rules(3).complianceName = 'C.PCOX';

rules(4).resistanceName = 'R.PVEN';
rules(4).complianceName = 'C.PVEN';
end
