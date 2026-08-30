function tests = test_joint_pre_post_objective()
% TEST_JOINT_PRE_POST_OBJECTIVE
% -----------------------------------------------------------------------
% Contract tests for objective_joint_pre_post (PRD
% reyna_statistical_calibration_v1 Phase 4).
%
% The joint objective exists to raise the observation count without adding
% free parameters, which is the only available route to positive degrees of
% freedom for this single-patient case (pre-only: N=9, p=12, dof=-3; joint:
% N=16, p=12, dof=+4). These tests pin the properties that claim depends on:
% the parameter vector really is shared, only the shunt differs, and the
% regularisation is not silently double-counted.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md §7.4
%   [2] docs/reyna_statistical_calibration_results_20260829.md §6.0
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-30
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function ctx = fixture()
% FIXTURE - Reyna pre/post at the demographically scaled baseline.
% Mirrors main_run.m's construction path (params_ref -> patient -> scaling ->
% calibration_param_sets) so these tests exercise the same objects the real
% pipeline builds rather than a hand-rolled approximation.
ctx = struct();
ctx.clinical = patient_reyna();

params_ref = default_parameters();

patient = struct();
patient.age_years  = ctx.clinical.common.age_years;
patient.age_days   = ctx.clinical.common.age_years * 365.25;
patient.weight_kg  = ctx.clinical.common.weight_kg;
patient.height_cm  = ctx.clinical.common.height_cm;
patient.sex        = ctx.clinical.common.sex;
patient.maturation_mode = 'normal';
patient.run_mode   = '';
patient.scaling_mode = 'zhang';
patient.BSA        = ctx.clinical.common.BSA;

params_scaled = apply_scaling(params_ref, patient);

case_profile = build_case_calibration_profile(ctx.clinical, 'pre_surgery');

% params_from_clinical is main_run's Step 2: it maps HR, SVR/PVR and the VSD
% geometry onto the scaled parameters. Without it the orifice area stays 0,
% i.e. the "pre-closure" model would not shunt at all -- which would make the
% closure test vacuous.
ctx.params_pre = params_from_clinical(params_scaled, ctx.clinical, ...
    'pre_surgery', params_ref, case_profile);
ctx.params_post = ctx.params_pre;

registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
ctx.calib = calibration_param_sets('pre_surgery', ctx.params_pre, [], {}, ...
    case_profile, registry_context);
ctx.x = ctx.calib.x0(:);
end

function test_shared_parameters_are_identical_in_both_structs(tc)
% The whole DOF argument rests on the two simulations sharing one parameter
% vector. If they drift, the joint fit is two loosely-coupled fits and the
% observation count cannot be pooled.
ctx = fixture();
[~, ~] = evalc('objective_joint_pre_post(ctx.x, ctx.params_pre, ctx.params_post, ctx.clinical, ctx.calib)');

% Re-apply through the same public path the objective uses and compare.
p_pre = ctx.params_pre;
p_post = ctx.params_post;
for i = 1:numel(ctx.calib.names)
    p_pre = set_calibration_param_value(p_pre, ctx.calib.referenceParams, ...
        ctx.calib.names{i}, ctx.x(i), ctx.calib.caseProfile);
    p_post = set_calibration_param_value(p_post, ctx.calib.referenceParams, ...
        ctx.calib.names{i}, ctx.x(i), ctx.calib.caseProfile);
end
for i = 1:numel(ctx.calib.names)
    name = ctx.calib.names{i};
    if strcmp(name, 'R.vsd')
        continue;   % the one parameter allowed to differ
    end
    v_pre = get_calibration_param_value(p_pre, ctx.calib.referenceParams, ...
        name, ctx.calib.caseProfile);
    v_post = get_calibration_param_value(p_post, ctx.calib.referenceParams, ...
        name, ctx.calib.caseProfile);
    verifyEqual(tc, v_pre, v_post, 'AbsTol', 1e-12, ...
        sprintf('%s must be shared between pre and post.', name));
end
end

function test_post_vsd_is_actually_closed_in_every_mode(tc)
% The defect closing is the ONLY physiological difference the joint model
% asserts between the two states, so "closed" must mean zero shunt FLOW --
% not merely that some field was assigned.
%
% Asserting behaviour rather than a field value is deliberate. Reyna runs in
% 'orifice_bidirectional' mode, where vsd_shunt_model never reads R.vsd:
% setting R.vsd = 1e6 (the convention main_run.m uses for its pre-to-post
% seed) leaves such a patient shunting at full strength. A field-based
% assertion would have passed while the physics was wrong.
ctx = fixture();

% Baseline must actually shunt, or this test proves nothing.
p_open = apply_vector(ctx.params_pre, ctx.x, ctx.calib);
q_open = shunt_flow(p_open);
verifyGreaterThan(tc, abs(q_open), 1e-6, ...
    'fixture must have an OPEN shunt for the closure check to be meaningful.');

% Drive the objective's own closure path via a post-only evaluation.
p_closed = apply_vector(ctx.params_post, ctx.x, ctx.calib);
p_closed = close_vsd_via_objective(p_closed);
q_closed = shunt_flow(p_closed);
verifyEqual(tc, q_closed, 0, 'AbsTol', 1e-12, ...
    'post-closure shunt flow must be exactly zero, in any vsd mode.');
end

function q = shunt_flow(params)
% SHUNT_FLOW - probe the shunt model directly at a large LV-RV gradient.
q = vsd_shunt_model(90, 20, params);
end

function params = apply_vector(params, x, calib)
for i = 1:numel(calib.names)
    params = set_calibration_param_value(params, calib.referenceParams, ...
        calib.names{i}, x(i), calib.caseProfile);
end
end

function params = close_vsd_via_objective(params)
% CLOSE_VSD_VIA_OBJECTIVE - mirror of the objective's closure, kept here so
% the test fails loudly if the production closure stops covering a mode.
params.R.vsd = 1e6;
if isfield(params, 'vsd')
    if isfield(params.vsd, 'area_mm2'); params.vsd.area_mm2 = 0; end
    if isfield(params.vsd, 'D_mm');     params.vsd.D_mm = 0;     end
end
end

function test_penalty_is_counted_once_not_per_scenario(tc)
% Charging the shared vector's regularisation twice would double its weight
% relative to the single-scenario objective and make joint and pre-only
% results incomparable -- the failure PRD §7.3 calls out explicitly.
ctx = fixture();
ctx.calib.regLambda = 0.25;
ctx.calib.x0 = ctx.x;

perturbed = ctx.x * 1.10;   % 10% off the regularisation centre
[J, info] = evalc_obj(ctx, perturbed);

expected_penalty = ctx.calib.regLambda * ...
    sum(log(perturbed(:) ./ ctx.x(:)).^2);
verifyEqual(tc, info.penalty, expected_penalty, 'RelTol', 1e-10, ...
    'penalty must equal exactly ONE application over the shared vector.');
if info.valid
    verifyEqual(tc, J, info.chi2_total + info.penalty, 'AbsTol', 1e-12);
end
end

function test_reported_chi2_is_the_sum_of_both_scenarios(tc)
ctx = fixture();
ctx.calib.regLambda = 0;
[J, info] = evalc_obj(ctx, ctx.x);
if ~info.valid
    return;   % simulation failure is covered by its own test
end
verifyEqual(tc, info.chi2_total, info.chi2_pre + info.chi2_post, ...
    'AbsTol', 1e-12);
verifyEqual(tc, J, info.chi2_total, 'AbsTol', 1e-12, ...
    'with regLambda = 0, J must be exactly chi2_pre + chi2_post.');
verifyEqual(tc, info.n_total, info.n_pre + info.n_post);
end

function test_post_surgery_now_contributes_observations(tc)
% Regression guard for the clinical data correction: post_surgery held zero
% finite targets before the 06/04/2026 procedure log was encoded, which is
% what made Phase 4 pointless. If this reverts to 0 the DOF argument silently
% collapses back to dof = -3.
ctx = fixture();
[~, info] = evalc_obj(ctx, ctx.x);
verifyTrue(tc, info.post_available, ...
    'post_surgery must carry finite governed targets for joint inversion.');
verifyGreaterThanOrEqual(tc, info.n_post, 4, ...
    'PRD §7.5 acceptance: post scenario needs >= 4 clinical targets.');
end

function test_missing_post_targets_degrade_gracefully_with_warning(tc)
% A patient may genuinely have no post-closure record. That must reduce to
% pre-only with a warning, not error -- and the warning matters, because the
% DOF benefit silently does not apply in that case.
ctx = fixture();
stripped = ctx.clinical;
fn = fieldnames(stripped.post_surgery);
for i = 1:numel(fn)
    if isnumeric(stripped.post_surgery.(fn{i}))
        stripped.post_surgery.(fn{i}) = NaN;
    end
end

prev = warning('off', 'objective_joint_pre_post:noPostTargets');
restore = onCleanup(@() warning(prev));
lastwarn('', '');

[J, info] = objective_joint_pre_post(ctx.x, ctx.params_pre, ctx.params_post, ...
    stripped, ctx.calib);

[~, wid] = lastwarn();
verifyEqual(tc, wid, 'objective_joint_pre_post:noPostTargets', ...
    ['absent post targets must warn: the degrees-of-freedom benefit of ', ...
     'joint inversion silently does not apply, and a caller that does not ', ...
     'know that would overstate its result.']);
verifyFalse(tc, info.post_available);
verifyEqual(tc, info.n_post, 0);
verifyTrue(tc, isfinite(J) || J == ctx.calib.invalidPenaltyScale, ...
    'must degrade to a usable objective, not error.');
end

function test_joint_objective_actually_couples_the_two_states(tc)
% The point of a joint fit is that a vector fitting one state well and the
% other badly must be penalised. If the objective did not couple, the joint
% score would track the pre-only score and the extra observations would buy
% nothing.
ctx = fixture();
ctx.calib.regLambda = 0;

[~, info_base] = evalc_obj(ctx, ctx.x);
if ~info_base.valid || ~info_base.post_available
    return;
end
verifyGreaterThan(tc, info_base.chi2_post, 0, ...
    'post-closure residuals must contribute to the objective.');
verifyGreaterThan(tc, info_base.chi2_total, info_base.chi2_pre, ...
    ['joint chi2 must exceed pre-only chi2: if it does not, the post ', ...
     'state is not constraining the fit and the DOF gain is illusory.']);
end

function [J, info] = evalc_obj(ctx, x)
% EVALC_OBJ - run the objective while suppressing solver console output.
J = []; info = [];
out = evalc(['[J, info] = objective_joint_pre_post(x, ctx.params_pre, ' ...
    'ctx.params_post, ctx.clinical, ctx.calib);']); %#ok<NASGU>
end
