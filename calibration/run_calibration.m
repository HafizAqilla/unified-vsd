function [params_best, calib_out] = run_calibration(params0, clinical, scenario)
% RUN_CALIBRATION
% -----------------------------------------------------------------------
% Calibrates a scenario-specific subset of model parameters to clinical
% haemodynamic targets using fmincon (gradient-based, constrained).
%
% The free-parameter list and metric targets are determined by 'scenario';
% the optimisation engine itself is identical for both scenarios.
%
% INPUTS:
%   params0   - baseline (scaled) parameter struct
%   clinical  - unified clinical struct (from patient_template.m)
%   scenario  - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   params_best - parameter struct with calibrated values
%   calib_out   - struct with convergence summary:
%                   .names    free parameter names
%                   .x0       initial values
%                   .xbest    optimal values
%                   .fbest    optimal objective value
%                   .lb, .ub  bounds used
%
% CALIBRATION OBJECTIVE  (weighted normalised least squares):
%   J(x) = Σ_k  w_k · ((y_k(x) − y_k^{clin}) / y_k^{clin})²
%         + λ · Σ_i  ((x_i − x0_i) / x0_i)²          (regularisation)
%
% SOLVER:
%   fmincon with 'interior-point' algorithm.
%   MultiStart (10 random starts) is used if the Global Optimization
%   Toolbox is available; otherwise a single fmincon run is performed.
%
% REFERENCES:
%   [1] objective_calibration.m — objective function.
%   [2] calibration_param_sets.m — scenario-specific free-parameter lists.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

%% Retrieve scenario-specific calibration configuration
calib = calibration_param_sets(scenario, params0);

%% Objective function handle
obj = @(x) objective_calibration(x, params0, clinical, calib, scenario);

%% Solver options
opts = optimoptions('fmincon', ...
    'Algorithm',           'interior-point', ...
    'Display',             'iter-detailed', ...
    'MaxFunctionEvaluations', 600, ...
    'MaxIterations',       200, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance',       1e-8);

problem = createOptimProblem('fmincon', ...
    'x0',       calib.x0, ...
    'objective', obj, ...
    'lb',       calib.lb, ...
    'ub',       calib.ub, ...
    'options',  opts);

%% Run optimisation
try
    ms    = MultiStart('Display', 'off');
    [xbest, fbest] = run(ms, problem, 10);
    fprintf('[run_calibration] MultiStart completed. Best J = %.6f\n', fbest);
catch
    warning('run_calibration:noGOT', ...
            'Global Optimization Toolbox not available; using single fmincon run.');
    [xbest, fbest] = fmincon(problem);
    fprintf('[run_calibration] fmincon completed. Best J = %.6f\n', fbest);
end

%% Unpack best parameters
params_best = params0;
for i = 1:numel(calib.names)
    params_best = set_param_by_name(params_best, calib.names{i}, xbest(i));
end

%% Collect output summary
calib_out         = struct();
calib_out.names   = calib.names;
calib_out.x0      = calib.x0;
calib_out.xbest   = xbest;
calib_out.fbest   = fbest;
calib_out.lb      = calib.lb;
calib_out.ub      = calib.ub;
calib_out.scenario = scenario;

end  % run_calibration

% =========================================================================
%  LOCAL HELPER
% =========================================================================

function params = set_param_by_name(params, name, value)
% SET_PARAM_BY_NAME — assign value to nested struct field via dot notation
parts = strsplit(name, '.');
switch numel(parts)
    case 2
        params.(parts{1}).(parts{2}) = value;
    case 3
        params.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_param_by_name:unsupportedDepth', ...
              'Only 2- or 3-level dot notation supported. Got: %s', name);
end
end
