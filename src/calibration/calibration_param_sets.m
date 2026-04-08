function calib = calibration_param_sets(scenario, params0, optMask)
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
%   optMask   - optional logical mask for active parameters [nParam x 1]
%               true  = active in optimisation
%               false = frozen (kept at baseline)
%
% OUTPUTS:
%   calib     - struct with fields:
%     .names_all       full cell array of parameter names (dot-notation)
%     .names           active-only parameter names
%     .metricFields    cell array of metric field names to match
%     .weights         struct of per-metric weights
%     .regLambda       regularisation strength
%     .mask            logical active-parameter mask
%     .x0_all          full initial guess vector
%     .lb_all, .ub_all full lower / upper bounds vectors
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
    calib.names_all = calib.names;

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

    % Weights — drive optimizer toward the largest post-calibration residuals.
    % Post-calibration error profile:
    %   Very bad (>30%):  PAP_mean (-38%), PVR (-57%), LVEDV (-44%), RVEDV (-36%)
    %   Significant:      RAP_mean (+35%), LVESV (-38%)
    %   Good (<10%):      SAP_mean (1.6%), QpQs (-1.1%), LVEF (-5.7%)
    % Primary metrics (QpQs, SAP_mean, LVEF) have smooth-scale amplification
    % in objective_calibration.m so their base weight can be moderate here.
    calib.weights = struct();
    for k = 1:numel(calib.metricFields)
        calib.weights.(calib.metricFields{k}) = 1.0;
    end
    calib.weights.LVEDV    = 5.0;   % major residual: -44%; top priority
    calib.weights.RVEDV    = 4.5;   % major residual: -36%
    calib.weights.PAP_mean = 3.5;   % major residual: -38%; was 1.5
    calib.weights.PVR      = 3.0;   % major residual: -57%; was 1.5
    calib.weights.LVESV    = 3.0;   % significant residual: -38%
    calib.weights.LVEF     = 3.0;   % primary; smooth scaling handles constraint; was 4.0
    calib.weights.RAP_mean = 2.0;   % significant residual: +35%; was 1.0
    calib.weights.QpQs     = 2.0;   % primary; was 1.5
    calib.weights.SAP_mean = 2.0;   % primary; MAP check
    calib.weights.SVR      = 1.5;   % borderline: -5.7%; was 1.0
    calib.weights.RVESV    = 1.5;   % minor residual: -11%; was 2.5

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
    calib.names_all = calib.names;

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

%% Initial guess and bounds from current params0 (full parameter set)
calib.x0_all = pack_x(params0, calib.names_all);
[calib.lb_all, calib.ub_all] = bounds_from_names(params0, calib.names_all, scenario);

%% Optional optimisation mask (Batch 2 GSA bridge)
if nargin < 3 || isempty(optMask)
    optMask = true(numel(calib.names_all), 1);
end

optMask = logical(optMask(:));
if numel(optMask) ~= numel(calib.names_all)
    error('calibration_param_sets:maskSizeMismatch', ...
          'Mask length (%d) must match parameter count (%d).', ...
          numel(optMask), numel(calib.names_all));
end
if ~any(optMask)
    error('calibration_param_sets:emptyActiveSet', ...
          'Mask deactivates all parameters. At least one active parameter is required.');
end

calib.mask = optMask;

% Active subset passed to the optimizer.
calib.names = calib.names_all(calib.mask);
calib.x0    = calib.x0_all(calib.mask);
calib.lb    = calib.lb_all(calib.mask);
calib.ub    = calib.ub_all(calib.mask);

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
        elseif ismember(nm, {'R.SAR', 'R.SVEN'})
            % Systemic resistances are analytically pre-conditioned; tighten
            % to prevent large excursions away from the Ohm's-law estimate.
            lb(i) = 0.4 * x0;
            ub(i) = 2.5 * x0;
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
    elseif startsWith(nm, 'V0.')
        % Widened upper bound: LVEDV and RVEDV are under-predicted by ~40%.
        % The optimizer needs room to substantially increase unstressed volumes.
        lb(i) = 0.2 * x0;
        ub(i) = 4.0 * x0;
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
