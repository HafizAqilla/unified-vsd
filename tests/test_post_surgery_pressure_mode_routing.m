function tests = test_post_surgery_pressure_mode_routing()
% TEST_POST_SURGERY_PRESSURE_MODE_ROUTING
% -----------------------------------------------------------------------
% Regression test for defect D2 (publication-readiness cleanup, Phase 2):
% build_case_calibration_profile.m's has_pressures check read
% SAP_mean_mmHg unconditionally, but the post_surgery clinical struct
% stores mean systemic pressure under MAP_mmHg instead (SAP_mean_mmHg is
% always NaN there). This made has_pressures always false for
% post_surgery cases with no calibration recipe, silently routing any
% such case to 'adaptive_patient' instead of 'sparse_cath' regardless of
% what pressure data was actually available.
%
% REFERENCES:
%   [1] src/calibration/build_case_calibration_profile.m
%   [2] src/calibration/build_case_calibration_profile.m:766-774
%       (has_systemic_bundle), the existing scenario-aware pattern this
%       fix now follows.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-03
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function clinical = sparse_post_surgery_fixture()
% A post-surgery case with pulmonary/systemic pressures and QpQs, but no
% cardiac output or chamber volumes, and no matching calibration recipe
% (so recipe_found is false and the has_pressures branch actually runs).
clinical = patient_template();
clinical.common.patient_name = 'd2_regression_case';
clinical.common.patient_id = '';

clinical.post_surgery.MAP_mmHg = 75;
clinical.post_surgery.PAP_mean_mmHg = 15;
clinical.post_surgery.QpQs = 1.02;

clinical.post_surgery.CO_Lmin = NaN;
clinical.post_surgery.LVEDV_mL = NaN;
clinical.post_surgery.LVESV_mL = NaN;
clinical.post_surgery.RVEDV_mL = NaN;
clinical.post_surgery.RVESV_mL = NaN;
end

function test_post_surgery_MAP_is_recognized_as_a_pressure(tc)
clinical = sparse_post_surgery_fixture();
[~, recipe_found] = load_calibration_recipe(clinical, 'post_surgery');
verifyFalse(tc, recipe_found, ...
    'Fixture must have no matching recipe for this check to be meaningful.');

profile = build_case_calibration_profile(clinical, 'post_surgery');

verifyEqual(tc, profile.mode, 'sparse_cath', ...
    ['A post-surgery case with PAP_mean_mmHg + MAP_mmHg + QpQs but no ', ...
     'CO/volumes must be routed to sparse_cath, not silently fall ', ...
     'through to adaptive_patient because MAP_mmHg was not recognized ', ...
     'as a systemic pressure.']);
end

function test_pre_surgery_still_uses_SAP_mean_mmHg(tc)
% Guard the guard: the fix must not have broken the pre_surgery path,
% which genuinely uses SAP_mean_mmHg.
clinical = sparse_post_surgery_fixture();
clinical.pre_surgery.PAP_mean_mmHg = 15;
clinical.pre_surgery.SAP_mean_mmHg = 75;
clinical.pre_surgery.QpQs = 1.02;
clinical.pre_surgery.CO_Lmin = NaN;
clinical.pre_surgery.LVEDV_mL = NaN;
clinical.pre_surgery.LVESV_mL = NaN;
clinical.pre_surgery.RVEDV_mL = NaN;
clinical.pre_surgery.RVESV_mL = NaN;

[~, recipe_found] = load_calibration_recipe(clinical, 'pre_surgery');
verifyFalse(tc, recipe_found, ...
    'Fixture must have no matching recipe for this check to be meaningful.');

profile = build_case_calibration_profile(clinical, 'pre_surgery');
verifyEqual(tc, profile.mode, 'sparse_cath');
end
