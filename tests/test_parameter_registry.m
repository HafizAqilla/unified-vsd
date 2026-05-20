%% test_parameter_registry.m
% =========================================================================
% Unit tests for parameter-registry-driven calibration bounds.
%
% PURPOSE:
%   Verifies that Batch 1 and Batch 2 wiring produces a centralized
%   registry, registry-driven vectors, and valid bounds for both full-data
%   and sparse-cath governance modes.
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
fprintf('  UNIFIED VSD MODEL - Parameter Registry Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% Common setup
params_ref = default_parameters();
clinical = patient_profile_A();
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
registry_context = struct('params_adult', params_ref, 'params_scaled', params_scaled);

%% Test 1: Full-data registry builds core fields
fprintf('--- Test 1: Full-data registry columns ---\n');
profile_full = build_case_calibration_profile(clinical, 'pre_surgery');
calib_full = calibration_param_sets('pre_surgery', params_seeded, [], {}, profile_full, registry_context);
required_cols = {'name','baseline_adult','baseline_scaled','seeded_value','lb','ub','bound_type'};
if all(ismember(required_cols, calib_full.parameterRegistry.Properties.VariableNames))
    fprintf('  [PASS] Registry contains the required metadata columns.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Registry is missing one or more required columns.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Build vector preserves requested order
fprintf('--- Test 2: Registry vector ordering ---\n');
if ismember('vsd.Cd', calib_full.parameterRegistry.name)
    vsd_param_name = 'vsd.Cd';
else
    vsd_param_name = 'R.vsd';
end
requested = {vsd_param_name,'group.R_sys_scale','E.LV.EA'};
[x0, lb, ub, names] = build_calibration_vector(calib_full.parameterRegistry, requested); %#ok<ASGLU>
if isequal(names, requested(:)) && numel(x0) == numel(requested)
    fprintf('  [PASS] build_calibration_vector preserves requested parameter order.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] build_calibration_vector ordering mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 3: Bounds validate for a real full-data case
fprintf('--- Test 3: Full-data bounds validation ---\n');
threw_3 = false;
try
    validate_bounds(calib_full.parameterRegistry, 'pre_surgery');
catch
    threw_3 = true;
end
if ~threw_3
    fprintf('  [PASS] Full-data registry passes bounds validation.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Full-data registry failed bounds validation.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Sparse cath registry exposes grouped scales
fprintf('--- Test 4: Sparse cath grouped parameter rows ---\n');
clinical_sparse = patient_profile_Razka();
patient_sparse = struct( ...
    'age_years', clinical_sparse.common.age_years, ...
    'age_days', clinical_sparse.common.age_years * 365.25, ...
    'weight_kg', clinical_sparse.common.weight_kg, ...
    'height_cm', clinical_sparse.common.height_cm, ...
    'sex', clinical_sparse.common.sex, ...
    'BSA', clinical_sparse.common.BSA, ...
    'maturation_mode', 'normal');
params_scaled_sparse = apply_scaling(params_ref, patient_sparse);
params_seeded_sparse = params_from_clinical(params_scaled_sparse, clinical_sparse, 'pre_surgery', params_scaled_sparse, struct());
profile_sparse = build_case_calibration_profile(clinical_sparse, 'pre_surgery');
registry_context_sparse = struct('params_adult', params_ref, 'params_scaled', params_scaled_sparse);
calib_sparse = calibration_param_sets('pre_surgery', params_seeded_sparse, [], {}, profile_sparse, registry_context_sparse);

expected_groups = {'group.R_sys_scale','group.R_pul_scale'};
if all(ismember(expected_groups, calib_sparse.parameterRegistry.name))
    fprintf('  [PASS] Sparse cath registry includes grouped calibration scales.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Sparse cath registry is missing grouped calibration scales.\n');
    n_fail = n_fail + 1;
end

%% Test 5: Full-data systemic scale can reduce load more than sparse cases
fprintf('--- Test 5: Full-data systemic load bound is relaxed ---\n');
idx_full_sys = find(strcmp(calib_full.parameterRegistry.name, 'group.R_sys_scale'), 1, 'first');
if ~isempty(idx_full_sys) && calib_full.parameterRegistry.lb(idx_full_sys) <= 0.25
    fprintf('  [PASS] Full-data registry allows lower systemic resistance scale.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Full-data registry did not relax group.R_sys_scale lower bound.\n');
    n_fail = n_fail + 1;
end

%% Test 6: Sparse cath bounds validate and remain stricter
fprintf('--- Test 6: Sparse cath bounds validation ---\n');
threw_5 = false;
try
    validate_bounds(calib_sparse.parameterRegistry, 'pre_surgery');
catch
    threw_5 = true;
end
idx_sparse_sys = find(strcmp(calib_sparse.parameterRegistry.name, 'group.R_sys_scale'), 1, 'first');
if ~threw_5 && ~isempty(idx_sparse_sys) && calib_sparse.parameterRegistry.lb(idx_sparse_sys) >= 0.50
    fprintf('  [PASS] Sparse cath registry passes bounds validation.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Sparse cath registry failed bounds validation or became too loose.\n');
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL PARAMETER REGISTRY TESTS PASSED\n');
else
    fprintf('  ONE OR MORE PARAMETER REGISTRY TESTS FAILED\n');
end
fprintf('==========================================\n');
