function params = apply_maturation(params, age_days, maturation_mode)
% APPLY_MATURATION
% -----------------------------------------------------------------------
% Adds pediatric age-dependent vascular adaptation on top of size scaling.
%
% INPUTS:
%   params            - parameter struct after physiological scaling [-]
%   age_days          - patient age                                 [days]
%   maturation_mode   - 'normal' | 'pvr_fixed_day3' |
%                       'pvr_fixed_day30' | 'none'
%
% OUTPUTS:
%   params            - parameter struct with vascular maturation    [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3 || isempty(maturation_mode)
    maturation_mode = 'normal';
end

maturation_mode = lower(char(maturation_mode));
age_days = max(0, age_days);

switch maturation_mode
    case 'none'
        params.maturation.mode = maturation_mode;
        params.maturation.age_days = age_days;
        params.maturation.pvr_factor = 1.0;
        params.maturation.svr_factor = 1.0;
        return;
    case 'normal'
        pvr_age_days = age_days;
    case 'pvr_fixed_day3'
        pvr_age_days = 3;
    case 'pvr_fixed_day30'
        pvr_age_days = 30;
    otherwise
        error('apply_maturation:unknownMode', ...
            'Unsupported maturation mode: %s', maturation_mode);
end

% Pulmonary resistance is high at birth and falls rapidly over the first
% weeks, while systemic resistance rises more gradually through infancy.
pvr_factor = 1.0 + 2.5 * exp(-pvr_age_days / 21.0);
svr_factor = 0.70 + 0.30 * (1.0 - exp(-age_days / 90.0));

params.R.PAR  = params.R.PAR  * pvr_factor;
params.R.PCOX = params.R.PCOX * pvr_factor;
params.R.PCNO = params.R.PCNO * pvr_factor;
params.R.PVEN = params.R.PVEN * pvr_factor;

params.C.PAR  = params.C.PAR  / sqrt(pvr_factor);
params.C.PCOX = params.C.PCOX / sqrt(pvr_factor);
params.C.PCNO = params.C.PCNO / sqrt(pvr_factor);
params.C.PVEN = params.C.PVEN / sqrt(pvr_factor);

params.R.SAR  = params.R.SAR  * svr_factor;
params.R.SC   = params.R.SC   * svr_factor;
params.R.SVEN = params.R.SVEN * svr_factor;

params.C.SAR  = params.C.SAR  / max(svr_factor, 1e-6);
params.C.SC   = params.C.SC   / max(svr_factor, 1e-6);
params.C.SVEN = params.C.SVEN / max(svr_factor, 1e-6);

params.maturation.mode = maturation_mode;
params.maturation.age_days = age_days;
params.maturation.pvr_factor = pvr_factor;
params.maturation.svr_factor = svr_factor;
params.maturation.pvr_reference_age_days = pvr_age_days;

end
