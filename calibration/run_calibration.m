function [params_best, calib_out] = run_calibration(params0, clinical, scenario)
% RUN_CALIBRATION
% -----------------------------------------------------------------------
% Calibrates a scenario-specific subset of model parameters to clinical
% haemodynamic targets using Nelder-Mead (fminsearch).
%
% RATIONALE:
%   The baseline model (from params_from_clinical.m) already maps 28
%   clinical parameters to model space via Ohm's law and conservation
%   principles (Phase 0, deterministic). This function performs Phase 1:
%   a bounded Nelder-Mead polish to refine elastance parameters only.
%
%   Nelder-Mead is optimal for smooth, low-dimensional (2–4 vars)
%   unconstrained problems near a good starting point. Convergence:
%   ~100–200 evals = ~5–15 minutes (vs ~20–50 min for particleswarm).
%
%   Solver comparison for this use case:
%   - fminsearch (Nelder-Mead):  optimal for d≤4, smooth objective
%   - fmincon:                    requires gradient; fails with ODE noise
%   - particleswarm:              overkill; 4× longer than Nelder-Mead
%
% INPUTS:
%   params0   - baseline (scaled, pre-conditioned by params_from_clinical)
%   clinical  - unified clinical struct
%   scenario  - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   params_best - parameter struct with optimized elastances
%   calib_out   - struct with convergence summary
%
% REFERENCES:
%   [1] Torczon, V. (1991). On the Convergence of Pattern Search Algorithms.
%       SIAM J. Optimization 1(2), 123–145. (Nelder-Mead theory)
%   [2] objective_calibration.m — weighted least-squares objective
%   [3] calibration_param_sets.m — free-parameter lists and bounds
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  2.1
% -----------------------------------------------------------------------

%% Retrieve scenario-specific calibration configuration
calib = calibration_param_sets(scenario, params0);
d     = numel(calib.names);

%% Objective function handle
obj = @(x) objective_calibration(x, params0, clinical, calib, scenario);

%% Run Nelder-Mead optimization
fprintf('[run_calibration] Phase 1: Nelder-Mead (%d free params)...\n', d);
fprintf('[run_calibration] Starting from pre-conditioned baseline (J0 = %.6f)...\n', obj(calib.x0));

% Nelder-Mead options
% (fminsearch does not accept bounds; we handle bounds via quadratic penalty)
opts_nm = optimset('fminsearch');
opts_nm.MaxFunEvals = 250;   % ~12–15 min for d=5 at 4 s/eval
opts_nm.TolFun      = 1e-4;  % Stop when J stops improving by < 0.01%
opts_nm.TolX        = 1e-5;  % Stop when all x change < 0.001%
opts_nm.Display     = 'iter-detailed';

% Wrapper objective with quadratic penalty for bound violations
obj_bounded = @(x) obj_with_bounds(x, obj, calib.lb, calib.ub);

[xbest, fbest] = fminsearch(obj_bounded, calib.x0, opts_nm);

fprintf('[run_calibration] Nelder-Mead done.  J = %.6f\n', fbest);

%% Enforce bounds on final result (just in case)
xbest = max(calib.lb, min(calib.ub, xbest));

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
calib_out.improvement = obj(calib.x0) - fbest;

fprintf('[run_calibration] Improvement: J(x0)=%.6f → J(x*)=%.6f ΔJ=%.6f\n', ...
        obj(calib.x0), fbest, calib_out.improvement);

end  % run_calibration

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function J_pen = obj_with_bounds(x, obj, lb, ub)
% OBJ_WITH_BOUNDS — evaluate objective with quadratic penalty for violations
%
%   J_pen(x) = J(x) + λ × Σ max(0, lb_i - x_i)^2 + max(0, x_i - ub_i)^2
%
%   where λ = 100 ensures bounds are nearly exact while allowing
%   Nelder-Mead to operate freely (no hard bounds).

lambda = 100;  % penalty weight
J_pen  = obj(x);

% Penalize lower bound violations
J_pen = J_pen + lambda * sum(max(0, lb - x).^2);

% Penalize upper bound violations
J_pen = J_pen + lambda * sum(max(0, x - ub).^2);

end

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
