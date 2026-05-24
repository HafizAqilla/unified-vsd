%% test_scaling_modes.m
% =========================================================================
% Unit tests for the Zhang 2019 weight-based allometric scaling.
%
% PURPOSE:
%   Verifies that the Zhang weight allometry is applied correctly and that
%   apply_scaling produces a deterministic, physiologically traceable baseline.
%   NOTE: apply_scaling now uses Zhang-only monolithic scaling (v3.1).
%   The Lundquist BSA path remains available via apply_physiological_scaling
%   directly (used by compare_scaling_methods.m) but is no longer invoked
%   from apply_scaling.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-22
% VERSION:  2.0  (updated for Zhang-only apply_scaling)
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

%% Test 1: Explicit modes are selectable via apply_physiological_scaling.
fprintf('--- Test 1: Explicit scaling modes (apply_physiological_scaling) ---\n');
params_zhang = apply_physiological_scaling(params_ref, patient, 'zhang');
params_lundquist = apply_physiological_scaling(params_ref, patient, 'lundquist_bsa');
if strcmp(params_zhang.scaling.mode, 'zhang') && strcmp(params_lundquist.scaling.mode, 'lundquist_bsa')
    fprintf('  [PASS] Zhang and Lundquist BSA modes are selectable via apply_physiological_scaling.\n');
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

%% Test 3: apply_scaling always uses Zhang mode (monolithic — mode not configurable).
fprintf('--- Test 3: apply_scaling uses Zhang 2019 monolithic scaling ---\n');
old_mode = getenv('UNIFIED_VSD_SCALING_MODE');
setenv('UNIFIED_VSD_SCALING_MODE', 'lundquist_bsa');   % environment override ignored
params_wrapped = apply_scaling(params_ref, patient);
setenv('UNIFIED_VSD_SCALING_MODE', old_mode);
if strcmp(params_wrapped.scaling.mode, 'zhang') && strcmp(params_wrapped.scaling.requested_mode, 'zhang')
    fprintf('  [PASS] apply_scaling uses Zhang 2019 monolithic scaling (environment override ignored).\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] apply_scaling scaling mode metadata is incorrect.\n');
    n_fail = n_fail + 1;
end

%% Test 4: apply_scaling default always produces Zhang output.
fprintf('--- Test 4: apply_scaling default mode is Zhang ---\n');
setenv('UNIFIED_VSD_SCALING_MODE', '');
params_default = apply_scaling(params_ref, patient);
setenv('UNIFIED_VSD_SCALING_MODE', old_mode);
if strcmp(params_default.scaling.mode, 'zhang')
    fprintf('  [PASS] apply_scaling default mode is Zhang 2019.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] apply_scaling default mode is not Zhang 2019.\n');
    n_fail = n_fail + 1;
end

%% Test 5: Zhang exponents in apply_scaling match declared values.
fprintf('--- Test 5: Zhang exponent verification in apply_scaling output ---\n');
W_ref = 70;
w = patient.weight_kg / W_ref;
expected_HR    = params_ref.HR    * w^(-0.300);
expected_R_SAR_z = params_ref.R.SAR * w^(-0.475);
expected_C_SAR_z = params_ref.C.SAR * w^(+1.000);
expected_E_LV_EA_z = params_ref.E.LV.EA * w^(-0.500);
expected_Rvalve_open = params_ref.Rvalve.open * w^(-0.500);
tol = 1e-10;
if abs(params_default.HR - expected_HR) < tol && ...
        abs(params_default.R.SAR - expected_R_SAR_z) < tol && ...
        abs(params_default.C.SAR - expected_C_SAR_z) < tol && ...
        abs(params_default.E.LV.EA - expected_E_LV_EA_z) < tol && ...
        abs(params_default.Rvalve.open - expected_Rvalve_open) < tol
    fprintf('  [PASS] Zhang exponents (HR, R.SAR, C.SAR, E.LV.EA, Rvalve.open) match.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Zhang exponent mismatch in apply_scaling output.\n');
    fprintf('         HR:       got %.6f, expected %.6f\n', params_default.HR, expected_HR);
    fprintf('         R.SAR:    got %.6f, expected %.6f\n', params_default.R.SAR, expected_R_SAR_z);
    fprintf('         C.SAR:    got %.6f, expected %.6f\n', params_default.C.SAR, expected_C_SAR_z);
    fprintf('         E.LV.EA:  got %.6f, expected %.6f\n', params_default.E.LV.EA, expected_E_LV_EA_z);
    fprintf('         Rvalve:   got %.6f, expected %.6f\n', params_default.Rvalve.open, expected_Rvalve_open);
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
