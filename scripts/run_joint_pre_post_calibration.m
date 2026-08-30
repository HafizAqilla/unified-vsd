function out = run_joint_pre_post_calibration(clinical, varargin)
% RUN_JOINT_PRE_POST_CALIBRATION
% -----------------------------------------------------------------------
% Driver for the joint pre/post-closure calibration (PRD
% reyna_statistical_calibration_v1 Phase 4).
%
% WHY THIS EXISTS
% ---------------
% The single-scenario fit for this patient is underdetermined: 9 governed
% observations against 12 active parameters, i.e. dof = -3. A low chi2/N is
% then guaranteed regardless of whether the model is right, so no statistical
% claim can rest on it. Fitting the pre- and post-closure states jointly from
% one shared parameter vector raises the observation count to 16 without
% adding parameters, which is the only route to positive dof available from
% this patient's record.
%
% Both states come from the SAME catheterisation session (device released
% 12.06.23, post readings 12.15-12.32, unchanged anaesthesia), so they are
% directly comparable rather than two separate studies.
%
% USAGE:
%   out = run_joint_pre_post_calibration(patient_reyna())
%   out = run_joint_pre_post_calibration(patient_reyna(), 'MaxFunEvals', 50)
%
% NAME-VALUE OPTIONS:
%   'ScalingMode'   default 'zhang'
%   'MaxFunEvals'   default 300   (set small for a smoke check)
%   'MaxIterations' default 40
%   'Verbose'       default true
%
% OUTPUTS:
%   out.x0, out.xbest          shared parameter vector, before and after   [-]
%   out.names                  parameter names for that vector             [-]
%   out.info0, out.info_best   objective breakdowns (see
%                              objective_joint_pre_post)                   [-]
%   out.dof                    N_total - p, the number this exists to make
%                              positive                              [count]
%   out.exitflag               fmincon exit flag                           [-]
%
% ASSUMPTIONS:
%   - Per-scenario target tiers are built explicitly here rather than reusing
%     the pre-surgery tiers for both. Reusing them would work by accident
%     (the governed metric names overlap) but would silently apply
%     pre-surgery governance to post-surgery rows.
%   - The shunt is the only difference between the two simulations; closure
%     is handled inside objective_joint_pre_post and is mode-aware.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md §7
%   [2] docs/reyna_statistical_calibration_results_20260829.md §6.4
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-30
% VERSION:  1.0
% -----------------------------------------------------------------------

opt = parse_options(varargin{:});

if nargin < 1 || isempty(clinical)
    clinical = patient_reyna();
end

%% ---- Baseline parameters, mirroring main_run's construction path -------
params_ref = default_parameters();
patient = build_patient_struct(clinical, opt.ScalingMode);
params_scaled = apply_scaling(params_ref, patient);

case_profile = build_case_calibration_profile(clinical, 'pre_surgery');

params_pre = params_from_clinical(params_scaled, clinical, 'pre_surgery', ...
    params_ref, case_profile);
params_post = params_from_clinical(params_scaled, clinical, 'post_surgery', ...
    params_ref, case_profile);

%% ---- Calibration configuration ----------------------------------------
registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
calib = calibration_param_sets('pre_surgery', params_pre, [], {}, ...
    case_profile, registry_context);

% Explicit per-scenario governance, taken from each scenario's own CASE
% PROFILE rather than from a bare build_target_tiers call.
%
% This distinction is not cosmetic. build_target_tiers(clinical, scenario)
% without the recipe's tier config does not honour
% recipe.primary_rmse_holdout, so it governs 10 pre-surgery rows where the
% production path governs 9 -- it silently readmits Q_shunt_Lmin, the
% algebraically-derived metric deliberately excluded from the governed RMSE.
% chi2_pre would then be computed over a different set than the reported
% governed RMSE, and the two would not be comparable. This is the code-path
% disagreement recorded in reyna_zhang_scientific_assessment_20260828.md §3;
% going through the case profile is what keeps both paths in step.
post_profile = build_case_calibration_profile(clinical, 'post_surgery');
calib.targetTiersByScenario = struct( ...
    'pre_surgery',  case_profile.targetTiers, ...
    'post_surgery', post_profile.targetTiers);

x0 = calib.x0(:);
lb = calib.lb(:);
ub = calib.ub(:);

%% ---- Baseline evaluation ----------------------------------------------
[J0, info0] = objective_joint_pre_post(x0, params_pre, params_post, ...
    clinical, calib);

p = numel(calib.names);
dof = info0.n_total - p;

if opt.Verbose
    fprintf('\n=== JOINT PRE/POST CALIBRATION ===\n');
    fprintf('  N_pre  = %d\n', info0.n_pre);
    fprintf('  N_post = %d\n', info0.n_post);
    fprintf('  N      = %d\n', info0.n_total);
    fprintf('  p      = %d\n', p);
    fprintf('  dof    = %d%s\n', dof, dof_annotation(dof));
    fprintf('  J0     = %.4f  (chi2_pre %.3f + chi2_post %.3f)\n', ...
        J0, info0.chi2_pre, info0.chi2_post);
    if ~info0.post_available
        fprintf(['  [WARNING] no post targets: this has reduced to a ', ...
            'pre-only fit and the dof benefit does NOT apply.\n']);
    end
end

%% ---- Optimise ----------------------------------------------------------
obj = @(x) objective_joint_pre_post(x, params_pre, params_post, clinical, calib);

% FiniteDifferenceStepSize and StepTolerance are set EXPLICITLY and must not
% be left at fmincon's defaults here.
%
% This objective is built on an ODE steady-state solve, so it carries the
% integrator's own noise floor. fmincon's default forward-difference step is
% about sqrt(eps) ~ 1.5e-8 relative, which is far below that floor: the
% differences then measure integration noise rather than the gradient. The
% observed symptom was an enormous reported first-order optimality (2.1e6)
% together with steps of 5e-8, and the solver halting after ONE iteration
% having moved J by 0.01% -- an optimiser returning its starting point, which
% reyna_zhang_scientific_assessment_20260828.md §2.2 rightly refuses to treat
% as a result.
%
% 1e-5 and 1e-6 match run_calibration.m:728-734, the path already proven to
% converge on this model.
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', ternary(opt.Verbose, 'iter-detailed', 'off'), ...
    'MaxFunctionEvaluations', opt.MaxFunEvals, ...
    'MaxIterations', opt.MaxIterations, ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize', opt.FiniteDifferenceStepSize, ...
    'OptimalityTolerance', 1e-5, ...
    'StepTolerance', 1e-6);

% Multi-start. The single-scenario result this is compared against came from
% 6 multi-starts through a 6-stage pipeline, so a single fmincon call from x0
% cannot distinguish "the shared-parameter assumption fails" from "this fit
% is simply under-converged". Starts are built to make that comparison fair:
% x0, an optional warm start from the pre-only calibrated solution (the
% sharpest test -- can the joint fit hold a good pre-op fit while also
% explaining post?), then bound-interior perturbations.
starts = build_starts(x0, lb, ub, opt);

Jbest = Inf; xbest = x0; exitflag = NaN;
for s = 1:size(starts, 2)
    try
        [xs, Js, ef] = fmincon(obj, starts(:, s), [], [], [], [], lb, ub, [], opts);
    catch ME
        if opt.Verbose
            fprintf('  start %d failed: %s\n', s, ME.message);
        end
        continue;
    end
    if opt.Verbose
        fprintf('  start %d/%d: J = %.4f\n', s, size(starts, 2), Js);
    end
    if isfinite(Js) && Js < Jbest
        Jbest = Js; xbest = xs(:); exitflag = ef;
    end
end

[~, info_best] = objective_joint_pre_post(xbest, params_pre, params_post, ...
    clinical, calib);

if opt.Verbose
    fprintf('\n  J: %.4f -> %.4f\n', J0, Jbest);
    fprintf('  chi2_pre : %.3f -> %.3f\n', info0.chi2_pre, info_best.chi2_pre);
    fprintf('  chi2_post: %.3f -> %.3f\n', info0.chi2_post, info_best.chi2_post);
    fprintf('  chi2/N   : %.3f -> %.3f\n', ...
        info0.chi2_total / max(info0.n_total, 1), ...
        info_best.chi2_total / max(info_best.n_total, 1));
    if dof > 0
        fprintf('  chi2/dof : %.3f  (meaningful: dof > 0)\n', ...
            info_best.chi2_total / dof);
    else
        fprintf(['  chi2/dof : n/a -- dof <= 0, so a low chi2 is ', ...
            'guaranteed and proves nothing.\n']);
    end
end

% OPTIMIZER_DID_NOT_MOVE guard.
%
% reyna_zhang_scientific_assessment_20260828.md §2.2 dissects a published
% comparison whose headline was, in substance, an optimiser that returned its
% starting point: identical baseline and calibrated RMSE to six decimals,
% reported as a physiological finding. That must never be reportable from
% this driver silently, so the condition is detected and carried on the
% result rather than left for a reader to notice.
rel_improvement = (J0 - Jbest) / max(abs(J0), eps);
rel_step = max(abs(xbest(:) - x0)) / max(max(abs(x0)), eps);
out_did_not_move = rel_improvement < 1e-3 || rel_step < 1e-6;

if out_did_not_move
    warning('run_joint_pre_post_calibration:optimizerDidNotMove', ...
        ['OPTIMIZER_DID_NOT_MOVE: J improved %.4g%% and the largest ', ...
         'relative parameter step was %.3g. This is NOT a calibration ', ...
         'result -- it is the starting point. Check the finite-difference ', ...
         'step against the ODE solver noise floor before interpreting.'], ...
        100 * rel_improvement, rel_step);
end
if opt.Verbose && out_did_not_move
    fprintf('\n  [OPTIMIZER_DID_NOT_MOVE] do not report this as a fit.\n');
end

out = struct('x0', x0, 'xbest', xbest(:), 'names', {calib.names(:)}, ...
    'J0', J0, 'Jbest', Jbest, 'info0', info0, 'info_best', info_best, ...
    'n_parameters', p, 'dof', dof, 'exitflag', exitflag, ...
    'relative_improvement', rel_improvement, ...
    'relative_step', rel_step, ...
    'optimizer_did_not_move', out_did_not_move);
end

% =======================================================================
function starts = build_starts(x0, lb, ub, opt)
% BUILD_STARTS - deterministic start set for the joint fit.
starts = x0(:);

% A warm start from the pre-only calibrated vector is the most informative
% single start available: it asks whether a parameter set that demonstrably
% explains the pre state can be held while also explaining the post state.
if ~isempty(opt.StartFrom)
    ws = opt.StartFrom(:);
    if numel(ws) == numel(x0)
        starts(:, end + 1) = min(max(ws, lb(:)), ub(:));
    end
end

n_extra = max(0, opt.NumStarts - size(starts, 2));
if n_extra > 0
    rng(20260830, 'twister');   % deterministic: reproducible start set
    span = ub(:) - lb(:);
    for k = 1:n_extra
        frac = 0.25 + 0.5 * rand(numel(x0), 1);
        starts(:, end + 1) = lb(:) + frac .* span; %#ok<AGROW>
    end
end
end

% =======================================================================
function note = dof_annotation(dof)
if dof > 2
    note = '  [positive: reduced chi2 is meaningful]';
elseif dof > 0
    note = '  [positive but small: reduced chi2 is unstable]';
else
    note = '  [NOT POSITIVE: chi2 cannot support any fit-quality claim]';
end
end

% =======================================================================
function patient = build_patient_struct(clinical, scaling_mode)
patient = struct();
patient.age_years  = clinical.common.age_years;
patient.age_days   = clinical.common.age_years * 365.25;
patient.weight_kg  = clinical.common.weight_kg;
patient.height_cm  = clinical.common.height_cm;
patient.sex        = clinical.common.sex;
patient.maturation_mode = 'normal';
patient.run_mode   = '';
patient.scaling_mode = scaling_mode;
if isfield(clinical.common, 'BSA') && isfinite(clinical.common.BSA)
    patient.BSA = clinical.common.BSA;
end
end

% =======================================================================
function opt = parse_options(varargin)
opt = struct('ScalingMode', 'zhang', 'MaxFunEvals', 300, ...
    'MaxIterations', 40, 'Verbose', true, ...
    'FiniteDifferenceStepSize', 1e-5, ...
    'StartFrom', [], ...        % optional warm start (e.g. the pre-only fit)
    'NumStarts', 1);
for i = 1:2:numel(varargin)
    name = varargin{i};
    if isfield(opt, name)
        opt.(name) = varargin{i + 1};
    else
        error('run_joint_pre_post_calibration:unknownOption', ...
            'Unknown option: %s', name);
    end
end
end

% =======================================================================
function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
