%% test_scaling_modes.m
% =========================================================================
% Unit tests for selectable pediatric scaling modes.
%
% PURPOSE:
%   Verifies that Zhang weight allometry and Lundquist BSA allometry are
%   selectable, traceable, and produce deterministic scaled baselines.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'src', 'utils'), '-begin');
addpath(fullfile(project_root, 'config'), '-begin');

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Scaling Mode Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

params_ref = default_parameters();
clinical = patient_reyna();
patient = struct( ...
    'age_years', clinical.common.age_years, ...
    'age_days', clinical.common.age_years * 365.25, ...
    'weight_kg', clinical.common.weight_kg, ...
    'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, ...
    'BSA', clinical.common.BSA);

%% Test 1: Explicit modes are selectable.
fprintf('--- Test 1: Explicit scaling modes ---\n');
params_zhang = apply_physiological_scaling(params_ref, patient, 'zhang');
params_lundquist = apply_physiological_scaling(params_ref, patient, 'lundquist_bsa');
if strcmp(params_zhang.scaling.mode, 'zhang') && strcmp(params_lundquist.scaling.mode, 'lundquist_bsa')
    fprintf('  [PASS] Zhang and Lundquist BSA modes are selectable.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Scaling mode metadata mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Lundquist BSA exponents are applied as declared.
fprintf('--- Test 2: Lundquist BSA exponent checks ---\n');
s = patient.BSA / 1.73;
expected_R_SAR = params_ref.R.SAR * s^-1.00;
expected_C_SAR = params_ref.C.SAR * s^+1.00;
expected_E_LV_EA = params_ref.E.LV.EA * s^-1.00;
tolerance = 1e-10;
if abs(params_lundquist.R.SAR - expected_R_SAR) < tolerance && ...
        abs(params_lundquist.C.SAR - expected_C_SAR) < tolerance && ...
        abs(params_lundquist.E.LV.EA - expected_E_LV_EA) < tolerance
    fprintf('  [PASS] Lundquist BSA resistance, compliance, and elastance scaling match declared exponents.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Lundquist BSA scaling does not match declared exponents.\n');
    n_fail = n_fail + 1;
end

%% Test 3: Wrapper accepts environment override.
fprintf('--- Test 3: apply_scaling environment override ---\n');
old_mode = getenv('UNIFIED_VSD_SCALING_MODE');
setenv('UNIFIED_VSD_SCALING_MODE', 'lundquist_bsa');
params_wrapped = apply_scaling(params_ref, patient);
setenv('UNIFIED_VSD_SCALING_MODE', old_mode);
if strcmp(params_wrapped.scaling.requested_mode, 'lundquist_bsa')
    fprintf('  [PASS] apply_scaling honors UNIFIED_VSD_SCALING_MODE.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] apply_scaling ignored UNIFIED_VSD_SCALING_MODE.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Wrapper default now prefers Lundquist BSA.
fprintf('--- Test 4: apply_scaling default mode ---\n');
setenv('UNIFIED_VSD_SCALING_MODE', '');
params_default = apply_scaling(params_ref, patient);
setenv('UNIFIED_VSD_SCALING_MODE', old_mode);
if strcmp(params_default.scaling.requested_mode, 'lundquist_bsa')
    fprintf('  [PASS] apply_scaling defaults to Lundquist BSA.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] apply_scaling default mode is not Lundquist BSA.\n');
    n_fail = n_fail + 1;
end

fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL SCALING MODE TESTS PASSED\n');
else
    fprintf('  SCALING MODE TESTS FAILED\n');
end
fprintf('==========================================\n');

assert(n_fail == 0, 'test_scaling_modes:failed', ...
    '%d scaling mode test(s) failed.', n_fail);
