function J = objective_calibration(x, params0, clinical, calib, scenario)
% OBJECTIVE_CALIBRATION
% -----------------------------------------------------------------------
% Weighted normalised least-squares objective for fmincon.
%
%   J = Σ_k  w_k · ((y_k(x) − y_k^{clin}) / y_k^{clin})²
%     + λ · Σ_i  ((x_i − x0_i) / x0_i)²
%
% INPUTS:
%   x         - parameter vector (free variables, ordered per calib.names)
%   params0   - baseline parameter struct (NOT modified in-place)
%   clinical  - unified clinical struct
%   calib     - configuration struct from calibration_param_sets.m
%   scenario  - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   J         - scalar objective value  (dimensionless)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

%% 1. Update parameter struct with trial x
params = params0;
for i = 1:numel(calib.names)
    params = set_param_by_name(params, calib.names{i}, x(i));
end

%% 2. Simulate to steady state
try
    sim = integrate_system(params);
catch
    J = 1e9;   % large penalty for solver failure
    return;
end

%% 3. Compute metrics
metrics = compute_clinical_indices(sim, params);

%% 4. Select clinical target sub-section
switch scenario
    case 'pre_surgery',  src = clinical.pre_surgery;
    case 'post_surgery', src = clinical.post_surgery;
end

%% 5. Build clinical lookup table (field mapping same as validation_report.m)
% Map model metric name → clinical struct field name
field_map = build_field_map(scenario);

%% 6. Weighted squared errors
J = 0;
for k = 1:numel(calib.metricFields)
    mf  = calib.metricFields{k};
    w   = calib.weights.(mf);

    % Retrieve clinical value via field map
    if isfield(field_map, mf) && isfield(src, field_map.(mf))
        y_clin = src.(field_map.(mf));
    elseif isfield(src, mf)
        y_clin = src.(mf);
    else
        continue;
    end

    if isnan(y_clin), continue; end

    if isfield(metrics, mf)
        y_model = metrics.(mf);
        denom   = max(abs(y_clin), 1e-6);
        J = J + w * ((y_model - y_clin) / denom)^2;
    end
end

%% 7. Regularisation: pull parameters toward initial guess
if isfield(calib, 'regLambda') && calib.regLambda > 0
    x0  = calib.x0(:);
    J   = J + calib.regLambda * sum(((x(:) - x0) ./ max(abs(x0), 1e-6)).^2);
end

end  % objective_calibration

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function params = set_param_by_name(params, name, value)
% SET_PARAM_BY_NAME — assign value to nested struct field (dot notation)
parts = strsplit(name, '.');
switch numel(parts)
    case 2, params.(parts{1}).(parts{2}) = value;
    case 3, params.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_param_by_name:unsupportedDepth', ...
              'Only 2- or 3-level notation supported. Got: %s', name);
end
end

function fmap = build_field_map(scenario)
% BUILD_FIELD_MAP — maps model metric name → clinical sub-struct field name
%   This must stay in sync with validation_report.m metric_defs.
fmap = struct();
switch scenario
    case 'pre_surgery'
        fmap.RAP_mean  = 'RAP_mean_mmHg';
        fmap.PAP_min   = 'PAP_sys_mmHg';
        fmap.PAP_max   = 'PAP_sys_mmHg';
        fmap.PAP_mean  = 'PAP_mean_mmHg';
        fmap.QpQs      = 'QpQs';
        fmap.PVR       = 'PVR_WU';
        fmap.SVR       = 'SVR_WU';
        fmap.LVEDV     = 'LVEDV_mL';
        fmap.LVESV     = 'LVESV_mL';
        fmap.RVEDV     = 'RVEDV_mL';
        fmap.RVESV     = 'RVESV_mL';
        fmap.LVEF      = 'LVEF';
    case 'post_surgery'
        fmap.SAP_min   = 'SAP_sys_mmHg';
        fmap.SAP_max   = 'SAP_sys_mmHg';
        fmap.SAP_mean  = 'MAP_mmHg';
        fmap.SVR       = 'SVR_WU';
        fmap.PVR       = 'PVR_WU';
        fmap.LVEF      = 'LVEF';
        fmap.RVEF      = 'RVEF';
        fmap.QpQs      = 'QpQs';
        fmap.LVEDV     = 'LVEDV_mL';
        fmap.RVEDV     = 'RVEDV_mL';
        fmap.RAP_mean  = 'RAP_mean_mmHg';
        fmap.PAP_mean  = 'PAP_mean_mmHg';
end
end
