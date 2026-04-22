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
%     .names           cell array of ACTIVE parameter names (dot-notation)
%     .names_all       cell array of ALL candidate names (superset of .names)
%     .mask            logical column vector: true where param is active
%     .metricFields    cell array of metric field names to match
%     .weights         struct of per-metric weights
%     .regLambda       regularisation strength
%     .x0              initial guess vector for ACTIVE params only
%     .x0_all          initial guess vector for ALL candidate params
%     .lb, .ub         lower / upper bounds vectors (ACTIVE params only)
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
        'R.SVEN'      % Systemic venous resistance
        'R.PAR'       % Pulmonary arterial resistance (must be included)
        'R.PVEN'      % Pulmonary venous resistance (must be included)
        'C.SVEN'      % Systemic venous compliance (must be included)
        'E.LV.EA'     % LV active elastance
        'E.LV.EB'     % LV passive elastance
        'E.RV.EA'     % RV active elastance
        'E.RV.EB'     % RV passive elastance
        'V0.RV'       % RV unstressed volume
        'V0.LV'       % LV unstressed volume — large LVESV mismatch observed
        'V0.LA'       % LA unstressed volume — controls LA filling → LV preload → LVEDV/LVESV
        'R.vsd'       % VSD shunt resistance (must be included)
        'R.SC'        % Systemic arteriolar resistance   — primary SVR contributor
        'R.PCOX'      % Pulmonary arteriolar resistance
        'C.PVEN'      % Pulmonary venous compliance
        'C.PAR'       % Pulmonary arterial compliance    — controls PAP_min
        'R.SAR'       % Aortic resistance               — controls MAP
        'C.SAR'       % Aortic (Windkessel) compliance  — controls pulse pressure
        };

    calib.metricFields = {
        'RAP_mean'
        'RAP_max'
        'RAP_min'
        'PAP_min'
        'PAP_max'
        'PAP_mean'
        'SAP_min'
        'SAP_max'
        'SAP_mean'
        'QpQs'
        'LVEDV'
        'LVESV'
        'RVEDV'
        'RVESV'
        };

    % Weights — prioritise mean pressures and QpQs as primary targets.
    % Systolic peaks (SAP_max, PAP_max) are the dominant clinical readouts
    % for VSD severity; they receive the highest weight.
    % Diastolic values (SAP_min, PAP_min) constrain vascular compliance
    % and receive standard weight.
    % RAP_mean is a secondary check (usually near-target from Ohm's law).
    calib.weights = struct();
    for k = 1:numel(calib.metricFields)
        calib.weights.(calib.metricFields{k}) = 1.0;   % default all to 1.0
    end
    % --- Primary haemodynamic anchors (pressures + flow ratio) ----------
    % SAP_min/SAP_max together constrain both contractility (peak) and
    % vascular compliance (trough).  Both must be fitted simultaneously;
    % previous run showed ±10–12% error on both — boost equally.
    calib.weights.SAP_max  = 5.0;   % systolic BP — primary VSD target (+10.6% error)
    calib.weights.SAP_min  = 5.0;   % diastolic BP — pulse-pressure pair (-12.3% error)
    calib.weights.SAP_mean = 3.0;   % MAP — key haemodynamic anchor (+5.4% error)
    calib.weights.PAP_max  = 4.0;   % PA systolic — primary indicator of PH (-10.5% error)
    calib.weights.PAP_mean = 3.0;   % mean PAP — pulmonary hypertension marker
    calib.weights.PAP_min  = 3.0;   % PA diastolic — compliance constraint
    calib.weights.QpQs     = 6.0;   % shunt ratio — #1 VSD severity metric (-9.7% error)
    calib.weights.RAP_mean = 3.0;   % venous pressure — no longer corrupted by NaN RAP_min
    calib.weights.RAP_max  = 2.0;   % venous peak pressure
    calib.weights.RAP_min  = 0.0;   % excluded: NaN clinical target — weight must be 0
    % --- Volume metrics (large errors observed, boost significantly) -----
    calib.weights.RVESV    = 7.0;   % RV ESV — residual +21% error
    calib.weights.RVEDV    = 4.0;   % RV EDV — near target (-2%)
    calib.weights.LVESV    = 7.0;   % LV ESV — -35% error, V0.LA added to address this
    calib.weights.LVEDV    = 3.0;   % LV EDV — near target (-1.8%)



    %% ==================================================================
    case 'post_surgery'
    % Free parameters: systemic/pulmonary circuit resistances and
    % compliances, ventricular elastances.  R_VSD is NOT free.
    % ================================================================
    calib.names = {
        'R.SVEN'      % Systemic venous resistance
        'R.PAR'       % Pulmonary arterial resistance (must be included)
        'R.PVEN'      % Pulmonary venous resistance (must be included)
        'C.SVEN'      % Systemic venous compliance (must be included)
        'E.LV.EA'     % LV active elastance
        'E.LV.EB'     % LV passive elastance
        'E.RV.EA'     % RV active elastance (second )
        'E.RV.EB'     % RV passive elastance
        'V0.RV'       % RV unstressed volume
        'R.SC'        % Systemic arteriolar resistance   — primary SVR contributor
        'R.PCOX'      % Pulmonary arteriolar resistance
        'C.PVEN'      % Pulmonary venous compliance
        'R.SAR'       % Aortic resistance               — controls MAP
        'C.SAR'       % Aortic (Windkessel) compliance  — controls pulse pressure
        };

    calib.metricFields = {
        'PAP_max'
        'PAP_min'
        'PAP_mean'
        'SAP_max'
        'SAP_min'
        'SAP_mean'
        'RAP_mean'
        'RAP_max'
        'RAP_min'
        'QpQs'
        'LVEDV'
        'LVESV'
        'RVEDV'
        'RVESV'
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

%% Traceability fields — required by run_calibration.m for calib_out
%
%  names_all : the complete candidate set from which .names was drawn.
%              Currently identical to .names (all candidates are active).
%              Reserved for future partial-freeze workflows where some
%              candidates may be frozen at their params0 values.
%
%  mask      : logical column vector (numel = numel(names_all)).
%              mask(i) == true  → param i is free (included in x vector).
%              mask(i) == false → param i is frozen at params0 value.
%              Currently all-true because every candidate is free.
%
%  x0_all    : initial-guess values for every candidate param (same order
%              as names_all).  Enables full audit of what was optimised
%              vs what was held constant.
calib.names_all = calib.names;               % [-]  all candidate names
calib.mask      = true(numel(calib.names_all), 1);  % [-]  all active
calib.x0_all    = calib.x0;                 % [-]  full initial vector

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

    if startsWith(nm, 'V0.')
        % Unstressed volumes: moderate range to preserve blood conservation
        lb(i) = 0.30 * x0;
        ub(i) = 3.00 * x0;
    elseif startsWith(nm, 'R.')
        if strcmp(nm, 'R.vsd') && strcmp(scenario, 'pre_surgery')
            % VSD resistance: very wide range — the Gorlin estimate has high
            % uncertainty; the true R.vsd is the primary free variable.
            lb(i) = max(0.001, 0.02 * x0);   % down to 2% of x0
            ub(i) = min(500,   30.0 * x0);   % up to 30× x0
        else
            % General resistances: wider than before (0.2×–5×) to allow
            % paediatric scaling correction beyond allometric estimates.
            lb(i) = 0.20 * x0;   % was 0.3×
            ub(i) = 5.00 * x0;   % was 3.0×
        end
    elseif startsWith(nm, 'C.')
        if strcmp(nm, 'C.SAR')
            % C.SAR (aortic Windkessel compliance): allometric scaling gives
            % C.SAR ≈ 0.038 mL/mmHg for a 13.4 kg child — far too stiff.
            % Clinical aortic compliance for a 3-year-old is typically
            % 0.3–1.5 mL/mmHg.  Previous run still showed ±12% SAP error,
            % meaning the optimizer did not reach the physiological range.
            % Upper bound widened to 40× to guarantee the full clinical
            % range is accessible from the allometric starting point.
            lb(i) = 0.10 * x0;    % allow slightly smaller if needed
            ub(i) = 40.0 * x0;   % widened: must cover 0.038 → 1.5 mL/mmHg
        else
            % All other compliances: wider (0.3×–4×)
            lb(i) = 0.30 * x0;
            ub(i) = 4.00 * x0;
        end
    elseif startsWith(nm, 'E.')
        lb(i) = 0.20 * x0;   % was 0.3×
        ub(i) = 5.00 * x0;   % was 3.0×
    else
        lb(i) = 0.30 * x0;
        ub(i) = 4.00 * x0;
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