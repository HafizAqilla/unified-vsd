function mask = create_optimization_mask(ST, threshold, parameterNames, caseProfile)
% CREATE_OPTIMIZATION_MASK
% -----------------------------------------------------------------------
% Builds a boolean optimisation mask from Sobol total-order indices.
%
% This function performs parameter screening prior to optimisation:
% parameters with low total sensitivity are frozen, while influential
% parameters are marked active.
%
% INPUTS:
%   ST         - Sobol total-order indices [nParam x nMetric] or [nParam x 1]
%                Values are dimensionless and expected in [0, 1].
%   threshold  - activation threshold for ST [dimensionless]
%                Typical value: 0.10
%   parameterNames - optional calibration parameter names [nParam x 1]
%   caseProfile    - optional evidence-aware calibration profile struct
%
% OUTPUTS:
%   mask       - logical column vector [nParam x 1]
%                true  = active parameter (optimise)
%                false = frozen parameter
%
% ASSUMPTIONS:
%   - If ST is multi-metric, conservative aggregation is used:
%     ST_agg(i) = max_j ST(i,j).
%   - Any NaN/Inf/negative ST entries are treated as 0 sensitivity.
%
% REFERENCES:
%   [1] Saltelli et al. (2010). Computer Physics Communications 181:259-270.
%   [2] Jansen (1999). Computer Physics Communications 117:35-43.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-16
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(threshold)
    threshold = 0.10;   % [dimensionless] default screening threshold
end
if nargin < 3
    parameterNames = {};
end
if nargin < 4
    caseProfile = struct();
end

validateattributes(ST, {'double', 'single'}, {'2d', 'nonempty'}, ...
    mfilename, 'ST', 1);
validateattributes(threshold, {'double', 'single'}, ...
    {'scalar', 'real', 'finite', '>=', 0, '<=', 1}, ...
    mfilename, 'threshold', 2);

ST = double(ST);

% Sanitize non-physical numerical artefacts from Monte-Carlo estimation.
ST(~isfinite(ST)) = 0;
ST(ST < 0)        = 0;

% Aggregate across metrics conservatively: keep parameter active if
% it is influential for any metric.
if size(ST, 2) > 1
    ST_agg = max(ST, [], 2);   % [dimensionless]
else
    ST_agg = ST(:);            % [dimensionless]
end

mask = ST_agg >= threshold;    % logical [nParam x 1]
mask = apply_case_profile_mask(mask, ST_agg, parameterNames, caseProfile);

% Safety guard: never return an all-false mask.
if ~any(mask)
    [~, idx_max] = max(ST_agg);
    mask(idx_max) = true;
end

mask = logical(mask(:));

end

function mask = apply_case_profile_mask(mask, ST_agg, parameterNames, caseProfile)
% APPLY_CASE_PROFILE_MASK — restrict active parameters by case evidence.
if isempty(parameterNames) || ~isfield(caseProfile, 'allowedFreeParameters') || ...
        isempty(caseProfile.allowedFreeParameters)
    return;
end
parameterNames = parameterNames(:);
allowed_mask = ismember(parameterNames, caseProfile.allowedFreeParameters(:));
if numel(allowed_mask) ~= numel(mask) || ~any(allowed_mask)
    return;
end

mask = mask & allowed_mask;
if ~any(mask)
    ST_allowed = ST_agg;
    ST_allowed(~allowed_mask) = -Inf;
    [~, idx_max] = max(ST_allowed);
    mask(idx_max) = true;
end
end
