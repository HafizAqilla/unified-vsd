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

function value = nested_value(s, path_parts)
% NESTED_VALUE - read a validated nested scalar from a structure.
value = s;
for part_idx = 1:numel(path_parts)
    value = value.(path_parts{part_idx});
end
end
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

%% Test 7: Reyna Zhang post-surgery registry keeps C.SAR baseline feasible
fprintf('--- Test 7: Reyna Zhang post-surgery C.SAR bounds ---\n');
old_scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
setenv('UNIFIED_VSD_SCALING_MODE', 'zhang');
threw_7 = false;
try
    clinical_reyna = patient_reyna();
    clinical_reyna.common.scaling_mode = 'zhang';
    clinical_reyna.post_surgery.HR = 108;              % [bpm]
    clinical_reyna.post_surgery.QpQs = 0.9;            % [-]
    clinical_reyna.post_surgery.PAP_sys_mmHg = 17.0;   % [mmHg]
    clinical_reyna.post_surgery.PAP_dia_mmHg = 9.0;    % [mmHg]
    clinical_reyna.post_surgery.PAP_mean_mmHg = 13.0;  % [mmHg]
    clinical_reyna.post_surgery.SAP_sys_mmHg = 89.0;   % [mmHg]
    clinical_reyna.post_surgery.SAP_dia_mmHg = 68.0;   % [mmHg]
    clinical_reyna.post_surgery.MAP_mmHg = 75.0;       % [mmHg]
    clinical_reyna.post_surgery.RAP_mean_mmHg = 5.0;   % [mmHg]
    clinical_reyna.post_surgery.CO_Lmin = 3.423;       % [L/min]

    patient_reyna_zhang = struct( ...
        'age_years', clinical_reyna.common.age_years, ...
        'age_days', clinical_reyna.common.age_years * 365.25, ...
        'weight_kg', clinical_reyna.common.weight_kg, ...
        'height_cm', clinical_reyna.common.height_cm, ...
        'sex', clinical_reyna.common.sex, ...
        'BSA', clinical_reyna.common.BSA, ...
        'maturation_mode', 'normal', ...
        'scaling_mode', 'zhang');
    params_scaled_reyna = apply_scaling(params_ref, patient_reyna_zhang);
    params_seeded_reyna = params_from_clinical( ...
        params_scaled_reyna, clinical_reyna, 'post_surgery', ...
        params_scaled_reyna, struct());
    profile_reyna = build_case_calibration_profile(clinical_reyna, 'post_surgery');
    registry_context_reyna = struct( ...
        'params_adult', params_ref, ...
        'params_scaled', params_scaled_reyna);
    calib_reyna = calibration_param_sets( ...
        'post_surgery', params_seeded_reyna, [], {}, ...
        profile_reyna, registry_context_reyna);
    idx_c_sar = find(strcmp(calib_reyna.parameterRegistry.name, 'C.SAR'), 1, 'first');
    idx_c_par = find(strcmp(calib_reyna.parameterRegistry.name, 'C.PAR'), 1, 'first');
    row_c_sar = calib_reyna.parameterRegistry(idx_c_sar, :);
    row_c_par = calib_reyna.parameterRegistry(idx_c_par, :);
    c_sar_ok = ~isempty(idx_c_sar) && ...
        row_c_sar.lb <= row_c_sar.baseline_scaled && ...
        row_c_sar.baseline_scaled <= row_c_sar.ub && ...
        ~strcmp(row_c_sar.bound_type{1}, 'recipe_absolute_final_bounds');
    c_par_ok = ~isempty(idx_c_par) && ...
        abs(row_c_par.lb - 0.8) < 1e-12 && ...
        abs(row_c_par.ub - 2.6825) < 1e-12 && ...
        strcmp(row_c_par.bound_type{1}, 'recipe_absolute_final_bounds');
catch
    threw_7 = true;
    c_sar_ok = false;
    c_par_ok = false;
end
setenv('UNIFIED_VSD_SCALING_MODE', old_scaling_mode);

if ~threw_7 && c_sar_ok && c_par_ok
    fprintf('  [PASS] Zhang post-surgery C.SAR is feasible and C.PAR uses Zhang-specific bounds.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Zhang post-surgery registry bounds are not as expected.\n');
    n_fail = n_fail + 1;
end

%% Test 8: Zhang post baseline retains selected fields from scaled baseline
fprintf('--- Test 8: Zhang post baseline warm-start exclusions ---\n');
old_scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
setenv('UNIFIED_VSD_SCALING_MODE', 'zhang');
threw_8 = false;
try
    clinical_reyna = patient_reyna();
    clinical_reyna.common.scaling_mode = 'zhang';
    clinical_reyna.post_surgery.HR = 108;              % [bpm]
    clinical_reyna.post_surgery.QpQs = 0.9;            % [-]
    clinical_reyna.post_surgery.PAP_sys_mmHg = 17.0;   % [mmHg]
    clinical_reyna.post_surgery.PAP_dia_mmHg = 9.0;    % [mmHg]
    clinical_reyna.post_surgery.PAP_mean_mmHg = 13.0;  % [mmHg]
    clinical_reyna.post_surgery.SAP_sys_mmHg = 89.0;   % [mmHg]
    clinical_reyna.post_surgery.SAP_dia_mmHg = 68.0;   % [mmHg]
    clinical_reyna.post_surgery.MAP_mmHg = 75.0;       % [mmHg]
    clinical_reyna.post_surgery.RAP_mean_mmHg = 5.0;   % [mmHg]

    patient_reyna_zhang = struct( ...
        'age_years', clinical_reyna.common.age_years, ...
        'age_days', clinical_reyna.common.age_years * 365.25, ...
        'weight_kg', clinical_reyna.common.weight_kg, ...
        'height_cm', clinical_reyna.common.height_cm, ...
        'sex', clinical_reyna.common.sex, ...
        'BSA', clinical_reyna.common.BSA, ...
        'maturation_mode', 'normal', ...
        'scaling_mode', 'zhang');
    params_scaled_reyna = apply_scaling(params_ref, patient_reyna_zhang);
    params_pre_reyna = params_scaled_reyna;
    params_pre_reyna.C.SVEN = 999;       % [mL/mmHg]
    params_pre_reyna.E.LV.EB = 999;      % [mmHg/mL]
    params_pre_reyna.E.RA.EA = 999;      % [mmHg/mL]
    params_pre_reyna.E.RV.EB = 999;      % [mmHg/mL]
    params_pre_reyna.R.SAR = 999;        % [mmHg*s/mL]
    params_pre_reyna.R.SVEN = 999;       % [mmHg*s/mL]
    clinical_reyna.pre_surgery.CalibParams = params_pre_reyna;

    [params_post_reyna, warm_start_reyna] = apply_post_surgery_warm_start( ...
        params_scaled_reyna, clinical_reyna, 'post_surgery');
    profile_reyna = build_case_calibration_profile(clinical_reyna, 'post_surgery');
    params_post_reyna = params_from_clinical( ...
        params_post_reyna, clinical_reyna, 'post_surgery', ...
        params_post_reyna, profile_reyna, 'prediction_baseline', true);

    excluded_paths = {'C.SVEN','E.LV.EB','E.RA.EA','E.RV.EB','R.SAR','R.SVEN'};
    warm_start_ok = true;
    for excluded_idx = 1:numel(excluded_paths)
        path_parts = strsplit(excluded_paths{excluded_idx}, '.');
        post_value = nested_value(params_post_reyna, path_parts);
        scaled_value = nested_value(params_scaled_reyna, path_parts);
        warm_start_ok = warm_start_ok && ...
            abs(post_value - scaled_value) < 1e-12 && ...
            ~ismember(excluded_paths{excluded_idx}, warm_start_reyna.copiedFields);
    end
catch
    threw_8 = true;
    warm_start_ok = false;
end
setenv('UNIFIED_VSD_SCALING_MODE', old_scaling_mode);

if ~threw_8 && warm_start_ok
    fprintf('  [PASS] Zhang post baseline keeps selected fields at scaled baseline values.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Zhang post baseline copied one or more excluded pre-op fields.\n');
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
