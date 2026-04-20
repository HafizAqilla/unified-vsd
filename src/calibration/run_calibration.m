function [params_best, calib_out] = run_calibration(params0, clinical, scenario, optMask, fastMode, pce_surrogate)
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
if nargin < 5
    fastMode = false;
end

% -----------------------------------------------------------------------        Retrieve scenario-specific calibration configuration
calib = calibration_param_sets(scenario, params0, optMask);
d     = numel(calib.names);

% Capping fmincon to the PCE surrogate boundaries.
% Extract bounds from the surrogate's Internal.Input (ground truth) rather
% than cfg.lb/ub, which can be stale if a checkpoint from a prior patient
% run was loaded.  Fall back to cfg if extraction fails.
if nargin >= 6 && ~isempty(pce_surrogate)
    % Resolve the true PCE training bounds per parameter name
    surr_lb = pce_surrogate.cfg.lb(:);
    surr_ub = pce_surrogate.cfg.ub(:);
    surr_fns = fieldnames(pce_surrogate);
    for kf = 1:numel(surr_fns)
        fn = surr_fns{kf};
        if isstruct(pce_surrogate.(fn)) && isfield(pce_surrogate.(fn), 'surrogate')
            s = pce_surrogate.(fn).surrogate;
            if isfield(s, 'Internal') && isfield(s.Internal, 'Input') && ...
                    isfield(s.Internal.Input, 'Marginals') && ...
                    numel(s.Internal.Input.Marginals) == numel(surr_lb)
                margs = s.Internal.Input.Marginals;
                surr_lb = arrayfun(@(m) m.Parameters(1), margs(:));
                surr_ub = arrayfun(@(m) m.Parameters(2), margs(:));
            end
            break;
        end
    end
    for i = 1:numel(calib.names)
        idx = find(strcmp(pce_surrogate.cfg.names, calib.names{i}));
        if ~isempty(idx)
            calib.lb(i) = max(calib.lb(i), surr_lb(idx));
            calib.ub(i) = min(calib.ub(i), surr_ub(idx));
        end
    end
    % Clamp x0 to the (possibly tightened) bounds so J0 is safe to evaluate
    calib.x0 = min(max(calib.x0, calib.lb), calib.ub);
end

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
obj1 = @(x) objective_calibration(x, params0, clinical, calib1, scenario, pce_surrogate);

calib2           = calib;
calib2.regLambda = 0;       % exploration in Stage 2; unregularised evaluation
obj2 = @(x) objective_calibration(x, params0, clinical, calib2, scenario, pce_surrogate);

% Unregularised objective at baseline — used for improvement reporting.
J0 = obj2(calib.x0);
fprintf('[run_calibration] J0 (unregularised, baseline) = %.6f\n', J0);

% -----------------------------------------------------------------------
%  Stage 1 — interior-point L-BFGS from analytical baseline
% -----------------------------------------------------------------------
fprintf('[run_calibration] Stage 1: fmincon interior-point + L-BFGS (%d active params)...\n', d);
fprintf('[run_calibration]   (max %d evals / %d iters — this may take a few minutes)\n', 4000, 300);
t_s1 = tic;

opts_s1 = optimoptions('fmincon', ...
    'Algorithm',                'interior-point', ...
    'HessianApproximation',     'lbfgs',          ...
    'FiniteDifferenceStepSize', 1e-5,             ...
    'FiniteDifferenceType',     'forward',        ...
    'Display',                  'iter',           ...
    'MaxFunctionEvaluations',   6000,             ...  % raised from 4000
    'MaxIterations',            500,              ...  % raised from 300
    'OptimalityTolerance',      1e-6,             ...
    'StepTolerance',            1e-6,             ...  % loosened from 1e-8 — prevents premature exitflag=2
    'UseParallel',              false);           % serial only: avoids worker path/cache mismatch

[x_s1, ~, exitflag_s1, output_s1] = fmincon( ...
    obj1, calib.x0, [], [], [], [], calib.lb, calib.ub, [], opts_s1);

% Re-evaluate Stage 1 result on the unregularised objective for fair comparison.
f_s1 = obj2(x_s1);
fprintf('[run_calibration] Stage 1 done in %.1fs.  J (unreg) = %.6f | exitflag = %d\n', ...
        toc(t_s1), f_s1, exitflag_s1);

% Initialise best-so-far with Stage 1 result.
xbest  = x_s1;
fbest  = f_s1;
best_stage = 1;
best_restart = 0;

% -----------------------------------------------------------------------
if ~fastMode
%  Stage 2 — diverse multi-start restarts (no regularisation)
%
%  Starting points mirror the scale factors from check_initial_guess_extremes.m:
%    scales = [0.30, 0.50, 1.80, 2.50]  applied to x0 (1.00 = Stage 1 already done)
%    plus one random uniform in [lb, ub]
%
%  Why scaled-x0 instead of raw lb/ub corners:
%    lb/ub corners set ALL 15 parameters to their extremes simultaneously.
%    That is almost certainly non-physiological (e.g. R.vsd at 2% AND
%    E.LV.EA at 20% AND C.SAR at 20% all at once → ODE failure or garbage J).
%    Scaled x0 preserves the relative balance between parameters while
%    moving the search to low/high regions of the feasible space.
%
%  Each start is pre-checked: if J >= J_skip the start is skipped entirely
%  (mirrors the J >= 1e8 guard in check_initial_guess_extremes.m).
% -----------------------------------------------------------------------
rng(2026, 'combRecursive');   % fixed seed for reproducibility

% Scale factors from check_initial_guess_extremes.m (1.00 skipped = Stage 1)
clip = @(x) min(max(x, calib.lb), calib.ub);
x_random = clip(calib.lb + rand(d,1).*(calib.ub - calib.lb));
starts = [clip(calib.x0 * 0.30), ...
          clip(calib.x0 * 0.50), ...
          clip(calib.x0 * 1.80), ...
          clip(calib.x0 * 2.50), ...
          x_random];
start_labels = {'0.30*x0', '0.50*x0', '1.80*x0', '2.50*x0', 'random'};
n_restart = size(starts, 2);

J_skip = 1e6;   % pre-check threshold — skip starts that are clearly infeasible

opts_s2 = optimoptions('fmincon', ...
    'Algorithm',                'interior-point', ...
    'HessianApproximation',     'lbfgs',          ...
    'FiniteDifferenceStepSize', 1e-5,             ...
    'FiniteDifferenceType',     'forward',        ...
    'Display',                  'off',            ...
    'MaxFunctionEvaluations',   2000,             ...
    'MaxIterations',            150,              ...
    'OptimalityTolerance',      1e-5,             ...
    'StepTolerance',            1e-7,             ...
    'UseParallel',              false);

fprintf('[run_calibration] Stage 2: %d diverse restarts — scales [0.30, 0.50, 1.80, 2.50, random]\n', n_restart);
fprintf('[run_calibration]   (max %d evals per restart, skip if J0 >= %.0e)\n', 2000, J_skip);
restart_J    = nan(n_restart, 1);
restart_flag = nan(n_restart, 1);
restart_x    = nan(d, n_restart);

for r = 1:n_restart
    fprintf('[run_calibration]   Restart %d/%d [%s] pre-check...\n', r, n_restart, start_labels{r});

    % Pre-check: skip starts that are clearly infeasible
    J0_r = obj2(starts(:,r));
    if J0_r >= J_skip
        fprintf('  SKIPPED (J0=%.2e >= %.0e)\n', J0_r, J_skip);
        restart_J(r)    = J0_r;
        restart_flag(r) = -998;
        continue;
    end

    try
        [x_r, ~, flag_r] = fmincon( ...
            obj2, starts(:,r), [], [], [], [], calib.lb, calib.ub, [], opts_s2);
        f_r = obj2(x_r);   % evaluate on clean objective (no regularisation)
        restart_J(r)    = f_r;
        restart_flag(r) = flag_r;
        restart_x(:,r)  = x_r;
        fprintf('  [Finish Restart %d] J = %.6f | exitflag = %d\n', r, f_r, flag_r);
    catch ME
        fprintf('  FAILED (%s)\n', ME.message);
    end
end

% Post-process restart results to find the global best
for r = 1:n_restart
    if restart_flag(r) >= -1 && restart_J(r) < fbest
        fbest       = restart_J(r);
        xbest       = restart_x(:,r);
        best_stage  = 2;
        best_restart = r;
    end
end

end % end if ~fastMode

fprintf('[run_calibration] Best solution: Stage %d', best_stage);
if best_stage == 2
    fprintf(', Restart %d', best_restart);
end
fprintf(' | J = %.6f\n', fbest);
fprintf('[run_calibration] Improvement vs baseline: J0=%.6f → J*=%.6f  ΔJ=%.6f\n', ...
        J0, fbest, J0 - fbest);

% -----------------------------------------------------------------------        Unpack best parameters
params_best = params0;
for i = 1:numel(calib.names)
    params_best = set_param_by_name(params_best, calib.names{i}, xbest(i));
end

% Reconstruct full parameter vector for traceability (active + frozen).
xbest_all = calib.x0_all;
xbest_all(calib.mask) = xbest;

% -----------------------------------------------------------------------        Collect output summary
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
if ~fastMode
    calib_out.restart_J   = restart_J;
    calib_out.restart_flag = restart_flag;
end
if ~fastMode
    calib_out.restart_J   = restart_J;
    calib_out.restart_flag = restart_flag;
end

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


