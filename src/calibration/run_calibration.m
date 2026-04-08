function [params_best, calib_out] = run_calibration(params0, clinical, scenario, optMask)
% RUN_CALIBRATION
% -----------------------------------------------------------------------
% Calibrates a scenario-specific subset of model parameters to clinical
% haemodynamic targets using a 2-stage bounded optimisation strategy.
%
% STRATEGY:
%   Stage 1 — Single start from the analytical baseline (x0).
%     Uses a soft L2 regularisation term (regLambda = 0.005) to prevent
%     the optimizer from drifting far from the physics-based starting point
%     during the initial search.  Interior-point + L-BFGS.
%
%   Stage 2 — Diverse multi-start restarts (n_restart points).
%     Regularisation is disabled so the optimizer can explore freely.
%     Starting points cover lower-quartile, upper-quartile, and one uniform
%     random sample of [lb, ub].  The best solution across all starts
%     (evaluated on the unregularised objective) is kept.
%
%   Selection — lowest unregularised J among all Stage 1 + Stage 2 runs,
%     restricted to runs with exitflag >= -1 (convergence or iteration limit).
%
% RATIONALE:
%   Single-start local optimization is vulnerable to poor basins when the
%   parameter landscape has multiple local minima (confirmed by run-to-run
%   variability in prior calibration results).  Two-stage strategy retains
%   the benefits of the analytical starting point while adding escape
%   capability via diverse restarts.
%
% INPUTS:
%   params0   - baseline (scaled, pre-conditioned by params_from_clinical)
%   clinical  - unified clinical struct
%   scenario  - 'pre_surgery' | 'post_surgery'
%   optMask   - optional logical mask of active parameters [nParam x 1]
%
% OUTPUTS:
%   params_best - parameter struct with optimized parameters
%   calib_out   - struct with convergence summary
%
% REFERENCES:
%   [1] Byrd, Lu, Nocedal (1995). A Limited Memory Algorithm for Bound
%       Constrained Optimization. SIAM J. Scientific Computing 16(5).
%   [2] objective_calibration.m — weighted least-squares + guard terms
%   [3] calibration_param_sets.m — free-parameter lists and bounds
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-08
% VERSION:  4.0
% -----------------------------------------------------------------------

if nargin < 4
    optMask = [];
end

%% Retrieve scenario-specific calibration configuration
calib = calibration_param_sets(scenario, params0, optMask);
d     = numel(calib.names);

if exist('fmincon', 'file') ~= 2
    error('run_calibration:missingFmincon', ...
          ['fmincon is not available. Install Optimization Toolbox ' ...
           'to run calibration.']);
end

% -----------------------------------------------------------------------
%  Objective handles
%    obj1 — Stage 1: soft regularisation (lambda = 0.005) stabilises the
%           search near the physics-based starting point.
%    obj2 — Stage 2 + evaluation: no regularisation; used for fair
%           comparison across all starts and for final reporting.
% -----------------------------------------------------------------------
calib1           = calib;
calib1.regLambda = 0.005;   % soft pull toward x0 in Stage 1
obj1 = @(x) objective_calibration(x, params0, clinical, calib1, scenario);

calib2           = calib;
calib2.regLambda = 0;       % exploration in Stage 2; unregularised evaluation
obj2 = @(x) objective_calibration(x, params0, clinical, calib2, scenario);

% Unregularised objective at baseline — used for improvement reporting.
J0 = obj2(calib.x0);
fprintf('[run_calibration] J0 (unregularised, baseline) = %.6f\n', J0);

% -----------------------------------------------------------------------
%  Stage 1 — interior-point L-BFGS from analytical baseline
% -----------------------------------------------------------------------
fprintf('[run_calibration] Stage 1: fmincon interior-point + L-BFGS (%d active params)...\n', d);

opts_s1 = optimoptions('fmincon', ...
    'Algorithm',                'interior-point', ...
    'HessianApproximation',     'lbfgs',          ...
    'FiniteDifferenceStepSize', 1e-5,             ...
    'FiniteDifferenceType',     'forward',        ...
    'Display',                  'iter-detailed',  ...
    'MaxFunctionEvaluations',   4000,             ...
    'MaxIterations',            300,              ...
    'OptimalityTolerance',      1e-6,             ...
    'StepTolerance',            1e-8);

[x_s1, ~, exitflag_s1, output_s1] = fmincon( ...
    obj1, calib.x0, [], [], [], [], calib.lb, calib.ub, [], opts_s1);

% Re-evaluate Stage 1 result on the unregularised objective for fair comparison.
f_s1 = obj2(x_s1);
fprintf('[run_calibration] Stage 1 done.  J (unreg) = %.6f | exitflag = %d\n', ...
        f_s1, exitflag_s1);

% Initialise best-so-far with Stage 1 result.
xbest  = x_s1;
fbest  = f_s1;
best_stage = 1;
best_restart = 0;

% -----------------------------------------------------------------------
%  Stage 2 — diverse multi-start restarts (no regularisation)
% -----------------------------------------------------------------------
n_restart = 3;
rng(2026, 'combRecursive');   % fixed seed for reproducibility

% Starting points: lower-quartile, upper-quartile, one uniform random.
starts = [calib.lb + 0.25*(calib.ub - calib.lb), ...
          calib.lb + 0.75*(calib.ub - calib.lb), ...
          calib.lb + rand(d, 1).*(calib.ub - calib.lb)];

opts_s2 = optimoptions('fmincon', ...
    'Algorithm',                'interior-point', ...
    'HessianApproximation',     'lbfgs',          ...
    'FiniteDifferenceStepSize', 1e-5,             ...
    'FiniteDifferenceType',     'forward',        ...
    'Display',                  'off',            ...
    'MaxFunctionEvaluations',   2000,             ...
    'MaxIterations',            150,              ...
    'OptimalityTolerance',      1e-5,             ...
    'StepTolerance',            1e-7);

fprintf('[run_calibration] Stage 2: %d diverse restarts (no regularisation)...\n', n_restart);
restart_J   = nan(n_restart, 1);
restart_flag = nan(n_restart, 1);

for r = 1:n_restart
    fprintf('[run_calibration]   Restart %d/%d ...', r, n_restart);
    try
        [x_r, ~, flag_r] = fmincon( ...
            obj2, starts(:,r), [], [], [], [], calib.lb, calib.ub, [], opts_s2);
        f_r = obj2(x_r);   % evaluate on clean objective (no regularisation)
        restart_J(r)    = f_r;
        restart_flag(r) = flag_r;
        fprintf('  J = %.6f | exitflag = %d', f_r, flag_r);
        if f_r < fbest && flag_r >= -1
            fbest       = f_r;
            xbest       = x_r;
            best_stage  = 2;
            best_restart = r;
            fprintf('  [NEW BEST]');
        end
        fprintf('\n');
    catch ME
        fprintf('  FAILED (%s)\n', ME.message);
    end
end

fprintf('[run_calibration] Best solution: Stage %d', best_stage);
if best_stage == 2
    fprintf(', Restart %d', best_restart);
end
fprintf(' | J = %.6f\n', fbest);
fprintf('[run_calibration] Improvement vs baseline: J0=%.6f → J*=%.6f  ΔJ=%.6f\n', ...
        J0, fbest, J0 - fbest);

%% Unpack best parameters
params_best = params0;
for i = 1:numel(calib.names)
    params_best = set_param_by_name(params_best, calib.names{i}, xbest(i));
end

% Reconstruct full parameter vector for traceability (active + frozen).
xbest_all = calib.x0_all;
xbest_all(calib.mask) = xbest;

%% Collect output summary
calib_out             = struct();
calib_out.names       = calib.names;
calib_out.names_all   = calib.names_all;
calib_out.mask        = calib.mask;
calib_out.x0          = calib.x0;
calib_out.xbest       = xbest;
calib_out.x0_all      = calib.x0_all;
calib_out.xbest_all   = xbest_all;
calib_out.J0          = J0;
calib_out.fbest       = fbest;
calib_out.lb          = calib.lb;
calib_out.ub          = calib.ub;
calib_out.scenario    = scenario;
calib_out.exitflag    = exitflag_s1;   % Stage 1 exit (main convergence indicator)
calib_out.output      = output_s1;
calib_out.improvement = J0 - fbest;
calib_out.best_stage  = best_stage;
calib_out.best_restart = best_restart;
calib_out.restart_J   = restart_J;
calib_out.restart_flag = restart_flag;

end  % run_calibration

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
