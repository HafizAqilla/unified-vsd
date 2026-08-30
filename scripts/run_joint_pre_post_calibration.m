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

opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', ternary(opt.Verbose, 'iter-detailed', 'off'), ...
    'MaxFunctionEvaluations', opt.MaxFunEvals, ...
    'MaxIterations', opt.MaxIterations, ...
    'FiniteDifferenceType', 'forward');

[xbest, Jbest, exitflag] = fmincon(obj, x0, [], [], [], [], lb, ub, [], opts);

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

out = struct('x0', x0, 'xbest', xbest(:), 'names', {calib.names(:)}, ...
    'J0', J0, 'Jbest', Jbest, 'info0', info0, 'info_best', info_best, ...
    'n_parameters', p, 'dof', dof, 'exitflag', exitflag);
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
    'MaxIterations', 40, 'Verbose', true);
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
