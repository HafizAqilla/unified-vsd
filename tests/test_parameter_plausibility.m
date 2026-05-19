%% test_parameter_plausibility.m
% =========================================================================
% Unit tests for evaluate_parameter_plausibility.m
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
addpath(fullfile(project_root, 'src', 'utils'), '-begin');
addpath(fullfile(project_root, 'config'), '-begin');

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Parameter Plausibility Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% Build a real active registry from Reyna
params_ref = default_parameters();
clinical = patient_reyna();
patient = struct( ...
    'age_years', clinical.common.age_years, ...
    'age_days', clinical.common.age_years * 365.25, ...
    'weight_kg', clinical.common.weight_kg, ...
    'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, ...
    'BSA', clinical.common.BSA, ...
    'maturation_mode', 'normal');
params_scaled = apply_scaling(params_ref, patient);
params_seeded = params_from_clinical(params_scaled, clinical, 'pre_surgery', params_scaled, struct());
profile = build_case_calibration_profile(clinical, 'pre_surgery');
registry_context = struct('params_adult', params_ref, 'params_scaled', params_scaled);
calib = calibration_param_sets('pre_surgery', params_seeded, [], {}, profile, registry_context);

%% Test 1: seeded active vector is plausible by construction
fprintf('--- Test 1: Seeded vector plausibility ---\n');
res_1 = evaluate_parameter_plausibility(calib.x0, calib.parameterRegistryActive);
if res_1.n_fail == 0 && height(res_1.table) == numel(calib.x0)
    fprintf('  [PASS] Seeded active vector is inside registry bounds.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Seeded active vector unexpectedly failed plausibility.\n');
    n_fail = n_fail + 1;
end

%% Test 2: obvious out-of-bounds value is flagged FAIL
fprintf('--- Test 2: Out-of-bounds detection ---\n');
x_bad = calib.x0;
x_bad(1) = 10 * calib.parameterRegistryActive.ub(1);
res_2 = evaluate_parameter_plausibility(x_bad, calib.parameterRegistryActive);
if any(res_2.table.PlausibilityFlag == "FAIL")
    fprintf('  [PASS] Out-of-bounds fitted value is flagged FAIL.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Out-of-bounds fitted value was not flagged.\n');
    n_fail = n_fail + 1;
end

%% Test 3: near-boundary value is flagged WARNING
fprintf('--- Test 3: Near-boundary warning ---\n');
x_warn = calib.x0;
ub0 = calib.parameterRegistryActive.ub(1);
lb0 = calib.parameterRegistryActive.lb(1);
x_warn(1) = ub0 - 0.05 * (ub0 - lb0);
res_3 = evaluate_parameter_plausibility(x_warn, calib.parameterRegistryActive);
if res_3.table.PlausibilityFlag(1) == "WARNING"
    fprintf('  [PASS] Near-boundary fitted value is flagged WARNING.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Near-boundary fitted value did not trigger WARNING.\n');
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL PARAMETER PLAUSIBILITY TESTS PASSED\n');
else
    fprintf('  ONE OR MORE PARAMETER PLAUSIBILITY TESTS FAILED\n');
end
fprintf('==========================================\n');
