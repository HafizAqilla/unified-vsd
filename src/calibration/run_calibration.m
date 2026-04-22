function [params_best, calib_out] = run_calibration(params0, clinical, scenario)
% RUN_CALIBRATION
% -----------------------------------------------------------------------
% Calibrates a scenario-specific subset of model parameters to clinical
% haemodynamic targets using fmincon (interior-point, L-BFGS Hessian).
%
% RATIONALE:
%   The baseline model (from params_from_clinical.m) already maps clinical
%   parameters to model space via Ohm's law (Phase 0, deterministic).
%   This function performs Phase 1: a bounded gradient-based polish.
%
%   fmincon (interior-point + L-BFGS) is chosen because:
%     - Hard bounds on all 11 parameters are enforced natively, avoiding
%       penalty-weight tuning required by Nelder-Mead wrappers.
%     - L-BFGS approximate Hessian is cheap (no Hessian storage) and
%       robust to noisy objectives from ODE finite-differencing.
%     - Forward-difference step 1e-5 is large enough to clear ODE noise
%       while remaining within the physiologically linear regime.
%     - Convergence at ~300–500 evals (vs 500+ for Nelder-Mead at d=11).
%
%   REQUIRES: Optimization Toolbox (fmincon).
%
% INPUTS:
%   params0   - baseline (scaled, pre-conditioned by params_from_clinical)
%   clinical  - unified clinical struct
%   scenario  - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   params_best - parameter struct with optimized parameters
%   calib_out   - struct with convergence summary, including:
%                   .exitflag   fmincon exit condition
%                   .output     fmincon output struct
%                   .mask       logical active-parameter mask
%                   .xbest_all  full parameter vector (active + frozen)
%
% REFERENCES:
%   [1] Nocedal J & Wright SJ (2006). Numerical Optimization, 2nd ed.
%       Springer. Ch. 7 (L-BFGS), Ch. 19 (interior-point).
%   [2] objective_calibration.m — weighted least-squares objective
%   [3] calibration_param_sets.m — free-parameter lists and bounds
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  3.0  (fmincon interior-point + L-BFGS)
% -----------------------------------------------------------------------

%% Retrieve scenario-specific calibration configuration
calib = calibration_param_sets(scenario, params0);
d     = numel(calib.names);

%% Objective function handle
obj = @(x) objective_calibration(x, params0, clinical, calib, scenario);

%% Run fmincon with interior-point + L-BFGS Hessian approximation
fprintf('[run_calibration] Phase 1: fmincon interior-point + L-BFGS (%d active params)...\n', d);
fprintf('[run_calibration] Starting from pre-conditioned baseline (J0 = %.6f)...\n', obj(calib.x0));

if exist('fmincon', 'file') ~= 2
    error('run_calibration:missingFmincon', ...
          ['fmincon is not available. Install Optimization Toolbox ' ...
           'to run calibration.']);
end

opts = optimoptions('fmincon', ...
    'Algorithm',                'interior-point', ...
    'HessianApproximation',     'lbfgs', ...
    'FiniteDifferenceStepSize', 1e-5,  ...
    'FiniteDifferenceType',     'forward', ...
    'Display',                  'iter-detailed', ...
    'MaxFunctionEvaluations',   8000, ...
    'MaxIterations',            600, ...
    'OptimalityTolerance',      1e-6, ...
    'StepTolerance',            1e-8);

[xbest, fbest, exitflag, output] = fmincon( ...
    obj, calib.x0, [], [], [], [], calib.lb, calib.ub, [], opts);

fprintf('[run_calibration] Phase 1 done.  J = %.6f | exitflag = %d\n', fbest, exitflag);

%% Phase 2 — refinement pass from best solution
%  Start from the Phase 1 solution with a small random perturbation on
%  parameters that may be trapped near a bound.  The perturbation is ±5%%
%  of the feasible range, clipped to [lb, ub].  This is a cheap way to
%  detect shallow local minima without the cost of MultiStart.
%
%  Only executed if Phase 1 exitflag suggests it did not fully converge
%  (exitflag == 0: iteration limit; exitflag < 0: infeasible / failure).
%  Always executed for exitflag == 1 (optimal) as a free polish step.
rng(42);   % reproducible perturbation
perturb_scale = 0.05;   % 5% of [lb, ub] range
range_x  = calib.ub - calib.lb;
perturb  = perturb_scale * range_x .* (2*rand(size(xbest)) - 1);  % ±5%
x0_phase2 = min(max(xbest + perturb, calib.lb), calib.ub);         % clip to bounds

opts2 = optimoptions(opts, ...
    'Display',                  'final', ...
    'MaxFunctionEvaluations',   3000, ...
    'MaxIterations',            300, ...
    'OptimalityTolerance',      1e-7, ...
    'StepTolerance',            1e-10);

[x2, f2, exitflag2, output2] = fmincon( ...
    obj, x0_phase2, [], [], [], [], calib.lb, calib.ub, [], opts2);

if f2 < fbest
    fprintf('[run_calibration] Phase 2 improved: J %.6f → %.6f | exitflag = %d\n', ...
            fbest, f2, exitflag2);
    xbest    = x2;
    fbest    = f2;
    exitflag = exitflag2;
    output   = output2;
else
    fprintf('[run_calibration] Phase 2 did not improve (J = %.6f). Keeping Phase 1.\n', f2);
end

fprintf('[run_calibration] Final best: J = %.6f | exitflag = %d\n', fbest, exitflag);

%% Unpack best parameters
params_best = params0;
for i = 1:numel(calib.names)
    params_best = set_param_by_name(params_best, calib.names{i}, xbest(i));
end

% Reconstruct full parameter vector for traceability (active + frozen).
% xbest_all has the same length as calib.names_all; active entries are
% replaced with the optimised values, frozen entries stay at x0_all.
xbest_all          = calib.x0_all;
xbest_all(calib.mask) = xbest;

%% Collect output summary
calib_out             = struct();
calib_out.names       = calib.names;        % active parameter names
calib_out.names_all   = calib.names_all;    % full candidate list
calib_out.mask        = calib.mask;         % active-param logical mask
calib_out.x0         = calib.x0;           % initial guess (active)
calib_out.xbest      = xbest;              % optimised values (active)
calib_out.x0_all     = calib.x0_all;       % initial guess (full)
calib_out.xbest_all  = xbest_all;          % optimised values (full)
calib_out.fbest      = fbest;              % best objective value
calib_out.lb         = calib.lb;
calib_out.ub         = calib.ub;
calib_out.exitflag   = exitflag;           % fmincon exit condition
calib_out.output     = output;             % fmincon diagnostic struct
calib_out.scenario   = scenario;
calib_out.improvement = obj(calib.x0) - fbest;

fprintf('[run_calibration] Improvement: J(x0)=%.6f → J(x*)=%.6f ΔJ=%.6f\n', ...
        obj(calib.x0), fbest, calib_out.improvement);

end  % run_calibration

% =========================================================================
%  LOCAL HELPERS
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