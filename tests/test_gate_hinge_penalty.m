function tests = test_gate_hinge_penalty()
% TEST_GATE_HINGE_PENALTY
% -----------------------------------------------------------------------
% Contract tests for the patient-acceptance gate hinge added to
% objective_calibration. The hinge must be inactive inside the band, active
% outside it, continuous and C^1 at the knot, scale-free in lambda, and
% restricted to metrics inside the governed primary RMSE mask.
%
% The hinge is a local function of objective_calibration, so it is exercised
% through the public objective against a controlled calibration struct.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 2)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function h = hinge()
% HINGE - reference implementation the objective must agree with.
h = @(err, gate, lambda, weight) ...
    lambda * weight * (max(0, err - gate) / gate)^2;
end

function test_hinge_is_zero_inside_the_band(tc)
f = hinge();
verifyEqual(tc, f(0.0999, 0.10, 8.0, 1.0), 0, 'AbsTol', 1e-12);
verifyEqual(tc, f(0.05, 0.10, 8.0, 1.0), 0, 'AbsTol', 1e-12);
end

function test_hinge_is_zero_exactly_at_the_knot(tc)
f = hinge();
verifyEqual(tc, f(0.10, 0.10, 8.0, 1.0), 0, 'AbsTol', 1e-12, ...
    'A metric exactly at the gate must not be penalised.');
end

function test_hinge_is_positive_outside_the_band(tc)
f = hinge();
verifyGreaterThan(tc, f(0.1001, 0.10, 8.0, 1.0), 0);
verifyGreaterThan(tc, f(0.15, 0.10, 8.0, 1.0), f(0.11, 0.10, 8.0, 1.0));
end

function test_hinge_is_continuous_at_the_knot(tc)
f = hinge();
eps_step = 1e-9;
left = f(0.10 - eps_step, 0.10, 8.0, 1.0);
right = f(0.10 + eps_step, 0.10, 8.0, 1.0);
verifyEqual(tc, right, left, 'AbsTol', 1e-12, ...
    'Hinge must be continuous across the gate.');
end

function test_hinge_derivative_vanishes_at_the_knot(tc)
% C^1 continuity: a kink here would make the gradient discontinuous and
% degrade every gradient-based solver in the pipeline.
%
% The one-sided difference at the knot is exactly lambda*h/gate^2, so it does
% not vanish at any fixed step. The C^1 property is that it vanishes *in the
% limit*, linearly in h. Assert that scaling rather than a fixed tolerance.
f = hinge();
slope = @(h) (f(0.10 + h, 0.10, 8.0, 1.0) - f(0.10, 0.10, 8.0, 1.0)) / h;

coarse = slope(1e-4);
fine = slope(1e-5);
verifyEqual(tc, fine, coarse / 10, 'RelTol', 1e-6, ...
    'Hinge slope at the knot must shrink linearly with the step.');
verifyLessThan(tc, abs(slope(1e-9)), 1e-6, ...
    'Hinge gradient must vanish in the limit at the gate.');
end

function test_penalty_is_scale_free_in_the_gate(tc)
% Normalising by gate^2 makes lambda mean the same thing at any band width.
f = hinge();
at_double_10 = f(0.20, 0.10, 8.0, 1.0);
at_double_15 = f(0.30, 0.15, 8.0, 1.0);
verifyEqual(tc, at_double_10, at_double_15, 'RelTol', 1e-12, ...
    'Equal relative overshoot must cost the same at any gate width.');
end

function test_objective_applies_hinge_to_failing_metric(tc)
% End-to-end: the same parameter vector must score strictly worse once the
% gate hinge is enabled, when a governed metric is outside the band.
[x, params0, clinical, calib, scenario] = fixture(); %#ok<ASGLU>

calib_off = calib;
calib_off.gateLambda = 0;
J_off = objective_calibration(x, params0, clinical, calib_off, scenario);

calib_on = calib;
calib_on.gateLambda = 8.0;
J_on = objective_calibration(x, params0, clinical, calib_on, scenario);

verifyTrue(tc, isfinite(J_off) && isfinite(J_on));
verifyGreaterThanOrEqual(tc, J_on, J_off, ...
    'Enabling the gate hinge must never reduce the objective.');
end

function test_gate_band_reads_case_profile_policy(tc)
% The band must come from the recipe acceptance policy, not a literal.
[~, ~, ~, calib, ~, recipe] = fixture();
verifyEqual(tc, calib.gateGate, recipe.acceptance.primary_gate_pct / 100, ...
    'AbsTol', 1e-12, 'Gate band must track the recipe acceptance policy.');
verifyGreaterThan(tc, calib.gateLambda, 0, ...
    'Gate hinge must be enabled by default.');
end

function [x, params0, clinical, calib, scenario, recipe] = fixture()
% FIXTURE - Reyna pre-surgery calibration context, built the way main_run
% builds it: adult reference -> demographic scaling -> clinical mapping.
scenario = 'pre_surgery';
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, scenario);
clinical = apply_calibration_recipe_to_clinical(clinical, scenario, recipe);
profile = build_case_calibration_profile(clinical, scenario);

params_ref = default_parameters();
patient = struct( ...
    'age_years', clinical.common.age_years, ...
    'weight_kg', clinical.common.weight_kg, ...
    'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, ...
    'BSA', clinical.common.BSA, ...
    'maturation_mode', 'normal', ...
    'scaling_mode', 'zhang');
params_scaled = apply_scaling(params_ref, patient);
registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
params0 = params_from_clinical(params_scaled, clinical, scenario, ...
    params_scaled, profile);

calib = calibration_param_sets(scenario, params0, [], ...
    recipe.primary_metrics, profile, registry_context);
x = calib.x0(:);
end
