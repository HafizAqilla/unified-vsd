function calib = calibration_param_sets(scenario, params0)
% CALIBRATION_PARAM_SETS
% -----------------------------------------------------------------------
% Returns the scenario-specific calibration configuration struct.
%
% The free-parameter list depends on the clinical scenario:
%   'pre_surgery'  — R_VSD, pulmonary resistances, ventricular elastances.
%                    These capture the dominant pre-surgery pathology
%                    (shunt magnitude and pulmonary hypertension).
%   'post_surgery' — systemic/pulmonary resistances, compliances, and
%                    ventricular elastances.  R_VSD is EXCLUDED (fixed at
%                    large value; the septal defect is closed).
%
% The calibration objective function, fmincon options, and convergence
% criteria are IDENTICAL for both scenarios — only the free-parameter list
% and metric targets change.
%
% INPUTS:
%   scenario  - 'pre_surgery' | 'post_surgery'
%   params0   - current (scaled) parameter struct (used to set x0 & bounds)
%
% OUTPUTS:
%   calib     - struct with fields:
%     .names           cell array of parameter names (dot-notation)
%     .metricFields    cell array of metric field names to match
%     .weights         struct of per-metric weights
%     .regLambda       regularisation strength
%     .x0              initial guess vector
%     .lb, .ub         lower / upper bounds vectors
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

calib = struct();

switch scenario

    %% ==================================================================
    case 'pre_surgery'
    % Free parameters: VSD resistance (primary pathology), pulmonary
    % circuit resistances, ventricular elastances.
    % ================================================================
    % Resistances are pre-conditioned analytically by params_from_clinical
    % (Ohm's law from SVR/PVR and VSD gradient/flow).  R.vsd is included
    % here because the mean-gradient conversion factor adds uncertainty.
    % The 4 elastance parameters are genuinely underdetermined.
    calib.names = {
        'R.SAR'       % Systemic arterial resistance
        'R.SVEN'      % Systemic venous resistance
        'R.PAR'       % Pulmonary arterial resistance
        'R.PVEN'      % Pulmonary venous resistance
        'C.SVEN'      % Systemic venous compliance
        'E.LV.EA'     % LV active elastance
        'E.LV.EB'     % LV passive elastance
        'E.RV.EA'     % RV active elastance
        'V0.LV'       % LV unstressed volume
        'V0.RV'       % RV unstressed volume
        'R.vsd'       % VSD shunt resistance
        };

    calib.metricFields = {
        'RAP_mean'
        'PAP_mean'
        'SAP_mean'
        'QpQs'
        'PVR'
        'SVR'
        'LVEDV'
        'LVESV'
        'RVEDV'
        'RVESV'
        'LVEF'
        };

    % Weights — emphasise volume/EF metrics (elastances are the free params)
    calib.weights = struct();
    for k = 1:numel(calib.metricFields)
        calib.weights.(calib.metricFields{k}) = 1.0;
    end
    calib.weights.LVEF     = 4.0;   % primary elastance target
    calib.weights.LVEDV    = 3.0;   % LV volume overload
    calib.weights.LVESV    = 3.0;
    calib.weights.RVEDV    = 2.5;
    calib.weights.RVESV    = 2.5;
    calib.weights.SAP_mean = 2.0;   % MAP — secondary check
    calib.weights.PAP_mean = 1.5;   % should already be near-target from Ohm's law
    calib.weights.PVR      = 1.5;
    calib.weights.QpQs     = 1.5;

    %% ==================================================================
    case 'post_surgery'
    % Free parameters: systemic/pulmonary circuit resistances and
    % compliances, ventricular elastances.  R_VSD is NOT free.
    % ================================================================
    calib.names = {
        'R.SAR'       % Systemic arterial resistance
        'R.SVEN'      % Systemic venous resistance
        'R.PAR'       % Pulmonary arterial resistance
        'R.PCOX'      % Pulmonary capillary resistance
        'R.PVEN'      % Pulmonary venous resistance
        'C.SAR'       % Systemic arterial compliance
        'C.PAR'       % Pulmonary arterial compliance
        'E.LV.EA'     % LV active elastance
        'E.LV.EB'     % LV passive elastance
        'E.RV.EA'     % RV active elastance
        'E.RV.EB'     % RV passive elastance
        };

    calib.metricFields = {
        'SAP_min'
        'SAP_max'
        'SAP_mean'
        'SVR'
        'PVR'
        'LVEF'
        'RVEF'
        'QpQs'
        'LVEDV'
        'RVEDV'
        'RAP_mean'
        'PAP_mean'
        };

    calib.weights = struct();
    for k = 1:numel(calib.metricFields)
        calib.weights.(calib.metricFields{k}) = 1.0;
    end
    calib.weights.SAP_mean = 3.0;
    calib.weights.SVR      = 2.5;
    calib.weights.LVEF     = 3.0;
    calib.weights.RVEF     = 2.5;
    calib.weights.QpQs     = 4.0;   % should recover to ~1.0 post-surgery

    otherwise
        error('calibration_param_sets:unknownScenario', ...
              'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

%% Regularisation
%  Set to 0: the analytical pre-conditioning in params_from_clinical already
%  provides a good starting point.  Non-zero lambda was preventing the
%  optimiser from moving parameters far enough to match clinical targets.
calib.regLambda = 0;

%% Initial guess and bounds from current params0
calib.x0 = pack_x(params0, calib.names);
[calib.lb, calib.ub] = bounds_from_names(params0, calib.names, scenario);

end  % calibration_param_sets

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function x = pack_x(params, names)
% PACK_X — extract parameter values into a column vector by dot-notation name
x = zeros(numel(names), 1);
for i = 1:numel(names)
    x(i) = get_param_by_name(params, names{i});
end
end

function [lb, ub] = bounds_from_names(params, names, scenario)
% BOUNDS_FROM_NAMES — physiological bounds per parameter type
lb = zeros(numel(names), 1);
ub = zeros(numel(names), 1);

for i = 1:numel(names)
    x0  = get_param_by_name(params, names{i});
    nm  = names{i};

    if startsWith(nm, 'R.')
        if strcmp(nm, 'R.vsd') && strcmp(scenario, 'pre_surgery')
            lb(i) = max(0.001, 0.05 * x0);
            ub(i) = min(500,   20.0 * x0);
        else
            lb(i) = 0.3 * x0;
            ub(i) = 3.0 * x0;
        end
    elseif startsWith(nm, 'C.')
        lb(i) = 0.5 * x0;
        ub(i) = 2.0 * x0;
    elseif startsWith(nm, 'E.')
        lb(i) = 0.3 * x0;
        ub(i) = 3.0 * x0;
    else
        lb(i) = 0.5 * x0;
        ub(i) = 2.0 * x0;
    end
end
end

function v = get_param_by_name(params, name)
% GET_PARAM_BY_NAME — resolve 'A.B.C' dot notation into nested struct access
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts)
    v = v.(parts{k});
end
end
