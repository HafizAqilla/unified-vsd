%% test_vascular_rc_coupling.m
% =========================================================================
% Unit tests for vascular resistance-compliance coupling.
%
% PURPOSE:
%   Verifies that vascular compliance is derived from resistance using the
%   Kung/Pennati-Fumero resting-condition relation rather than optimized
%   as an independent free parameter.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'src', 'calibration'), '-begin');
addpath(fullfile(project_root, 'config'), '-begin');

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Vascular R-C Coupling Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

params_ref = default_parameters();
profile = build_case_calibration_profile(patient_reyna(), 'pre_surgery');

%% Test 1: Systemic venous compliance follows resistance scaling.
fprintf('--- Test 1: Systemic venous R-C coupling ---\n');
params_test = params_ref;
params_test.R.SVEN = 2.0 * params_ref.R.SVEN;
params_test = enforce_vascular_rc_coupling(params_test, params_ref, profile);
expected_C = params_ref.C.SVEN * (2.0)^(-4/3);
if abs(params_test.C.SVEN - expected_C) < 1e-10
    fprintf('  [PASS] C.SVEN follows Kung/Pennati-Fumero resistance-compliance coupling.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] C.SVEN did not follow expected R-C coupling.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Grouped pulmonary resistance scale updates derived venous compliance.
fprintf('--- Test 2: Grouped pulmonary resistance scale updates derived compliance ---\n');
params_group = params_ref;
params_group = set_calibration_param_value(params_group, params_ref, 'group.R_pul_scale', 1.5, profile);
expected_C_pven = params_ref.C.PVEN * (1.5)^(-4/3);
if abs(params_group.C.PVEN - expected_C_pven) < 1e-10
    fprintf('  [PASS] group.R_pul_scale propagates to derived pulmonary venous compliance.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] group.R_pul_scale did not propagate to derived pulmonary venous compliance.\n');
    n_fail = n_fail + 1;
end

%% Test 3: Systemic arterial compliance is not overwritten by grouped resistance scaling.
fprintf('--- Test 3: Arterial compliance remains separately anchorable ---\n');
params_art = params_ref;
params_art = set_calibration_param_value(params_art, params_ref, 'group.R_sys_scale', 1.5, profile);
if abs(params_art.C.SAR - params_ref.C.SAR) < 1e-10
    fprintf('  [PASS] C.SAR is preserved for direct pulse-pressure anchoring.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] C.SAR was unexpectedly overwritten by grouped resistance scaling.\n');
    n_fail = n_fail + 1;
end

fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL VASCULAR R-C COUPLING TESTS PASSED\n');
else
    fprintf('  VASCULAR R-C COUPLING TESTS FAILED\n');
end
fprintf('==========================================\n');

assert(n_fail == 0, 'test_vascular_rc_coupling:failed', ...
    '%d vascular R-C coupling test(s) failed.', n_fail);
