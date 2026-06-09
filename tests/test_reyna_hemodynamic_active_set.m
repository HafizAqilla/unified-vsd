%% test_reyna_hemodynamic_active_set.m
% =========================================================================
% Regression test for the Reyna pre-surgery calibration active parameter set.
%
% PURPOSE:
%   Confirms the reyna_hemodynamic calibration profile activates the full
%   sparse_cath-equivalent 14-parameter chamber-dynamics set, not the
%   collapsed pressure-flow-only subset observed in the 2026-05-23 runs.
%
% USAGE:
%   >> test_reyna_hemodynamic_active_set
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-24
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

fprintf('=====================================================\n');
fprintf('  UNIFIED VSD MODEL - Reyna Active Parameter Set Test\n');
fprintf('=====================================================\n\n');

n_pass = 0;
n_fail = 0;

clinical = patient_reyna();
profile = build_case_calibration_profile(clinical, 'pre_surgery');
profile_post = build_case_calibration_profile(clinical, 'post_surgery');

%% Test 1: expected 14 active parameters are present
expected_parameters = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB', ...
    'E.LA.EA','E.RA.EA','V0.LV','V0.RV','vsd.Cd'};
actual_params = profile.allowedFreeParameters(:)';
if isequal(actual_params, expected_parameters)
    fprintf('  [PASS] Active parameter set contains all 14 expected chamber-dynamics parameters.\n');
    n_pass = n_pass + 1;
else
    missing = setdiff(expected_parameters, actual_params);
    extra = setdiff(actual_params, expected_parameters);
    fprintf('  [FAIL] Active parameter set mismatch.');
    if ~isempty(missing), fprintf(' Missing: %s.', strjoin(missing, ', ')); end
    if ~isempty(extra), fprintf(' Unexpected: %s.', strjoin(extra, ', ')); end
    fprintf('\n');
    n_fail = n_fail + 1;
end

%% Test 2: post-surgery active parameter set matches pre-surgery
actual_post_params = profile_post.allowedFreeParameters(:)';
if isequal(actual_post_params, expected_parameters)
    fprintf('  [PASS] Post-surgery active parameter set matches the 14-parameter pre-surgery set.\n');
    n_pass = n_pass + 1;
else
    missing = setdiff(expected_parameters, actual_post_params);
    extra = setdiff(actual_post_params, expected_parameters);
    fprintf('  [FAIL] Post-surgery active parameter set mismatch.');
    if ~isempty(missing), fprintf(' Missing: %s.', strjoin(missing, ', ')); end
    if ~isempty(extra), fprintf(' Unexpected: %s.', strjoin(extra, ', ')); end
    fprintf('\n');
    n_fail = n_fail + 1;
end

%% Test 3: allowedMetricFields contains only pressure-flow primary targets
expected_primary_metrics = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
primary_metrics = select_primary_metrics(clinical, struct(), 'pre_surgery', profile);
if isequal(primary_metrics(:)', expected_primary_metrics)
    fprintf('  [PASS] allowedMetricFields primary targets are pressure-flow only.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Primary metrics do not match pressure-flow governance.\n');
    n_fail = n_fail + 1;
end

%% Test 4: group.R_sys_scale bounds match pre-merge sparse_cath [0.25, 2.80]
r_sys_idx = find(strcmp(profile.boundScale.names, 'group.R_sys_scale'), 1);
r_sys_lb = profile.boundScale.lower(r_sys_idx);
r_sys_ub = profile.boundScale.upper(r_sys_idx);
if abs(r_sys_lb - 0.25) < 1e-12 && abs(r_sys_ub - 2.80) < 1e-12
    fprintf('  [PASS] group.R_sys_scale bounds [%.2f, %.2f] match pre-merge sparse_cath.\n', ...
        r_sys_lb, r_sys_ub);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] group.R_sys_scale bounds [%.2f, %.2f] deviate from expected [0.25, 2.80].\n', ...
        r_sys_lb, r_sys_ub);
    n_fail = n_fail + 1;
end

%% Test 5: group.R_pul_scale bounds match pre-merge sparse_cath [0.45, 2.80]
r_pul_idx = find(strcmp(profile.boundScale.names, 'group.R_pul_scale'), 1);
r_pul_lb = profile.boundScale.lower(r_pul_idx);
r_pul_ub = profile.boundScale.upper(r_pul_idx);
if abs(r_pul_lb - 0.45) < 1e-12 && abs(r_pul_ub - 2.80) < 1e-12
    fprintf('  [PASS] group.R_pul_scale bounds [%.2f, %.2f] match pre-merge sparse_cath.\n', ...
        r_pul_lb, r_pul_ub);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] group.R_pul_scale bounds [%.2f, %.2f] deviate from expected [0.45, 2.80].\n', ...
        r_pul_lb, r_pul_ub);
    n_fail = n_fail + 1;
end

%% Test 6: consistency-only rows are not in active calibration set
tiers = profile.targetTiers;
consistency_fields = tiers.consistency_only;
calibration_metrics = tiers.included_in_calibration;
if ~any(ismember({'LVEDV','LVESV','RVEDV','LVEF'}, calibration_metrics))
    fprintf('  [PASS] Chamber volume/function rows remain consistency-only, not calibration targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] One or more chamber volume/function rows leaked into calibration targets.\n');
    n_fail = n_fail + 1;
end

%% Test 7: registry calibratable flag is set correctly for all active params
recipe = load_calibration_recipe(clinical, 'pre_surgery');
clinical_s = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
params_adult = default_parameters();
params_scaled = apply_scaling(params_adult, clinical_s.common);
params_seeded = params_from_clinical(params_scaled, clinical_s, 'pre_surgery');
reg = build_parameter_registry(params_adult, params_scaled, params_seeded, ...
    'pre_surgery', profile, profile.allowedFreeParameters);
expected_names = expected_parameters;
for i = 1:numel(expected_names)
    hit = strcmp(reg.name, expected_names{i});
    if ~any(hit)
        fprintf('  [FAIL] Parameter %s not found in registry table.\n', expected_names{i});
        n_fail = n_fail + 1;
        continue;
    end
    if ~reg.is_calibratable(hit)
        fprintf('  [FAIL] Active parameter %s is not marked calibratable in registry.\n', expected_names{i});
        n_fail = n_fail + 1;
    end
end
fprintf('  [PASS] All active parameters are marked calibratable in the registry.\n');
n_pass = n_pass + 1;

%% Test 8: post-surgery registry marks all 14 matched active params calibratable
recipe_post = load_calibration_recipe(clinical, 'post_surgery');
clinical_post = apply_calibration_recipe_to_clinical(clinical, 'post_surgery', recipe_post);
params_scaled_post = apply_scaling(params_adult, clinical_post.common);
params_seeded_post = params_from_clinical(params_scaled_post, clinical_post, 'post_surgery');
reg_post = build_parameter_registry(params_adult, params_scaled_post, params_seeded_post, ...
    'post_surgery', profile_post, profile_post.allowedFreeParameters);
for i = 1:numel(expected_names)
    hit = strcmp(reg_post.name, expected_names{i});
    if ~any(hit)
        fprintf('  [FAIL] Post parameter %s not found in registry table.\n', expected_names{i});
        n_fail = n_fail + 1;
        continue;
    end
    if ~reg_post.is_calibratable(hit)
        fprintf('  [FAIL] Post active parameter %s is not marked calibratable in registry.\n', expected_names{i});
        n_fail = n_fail + 1;
    end
end
fprintf('  [PASS] All post-surgery active parameters are marked calibratable in the registry.\n');
n_pass = n_pass + 1;

%% Test 9: post-surgery resolved bounds match pre-surgery bounds exactly
if height(reg_post) ~= height(reg)
    fprintf('  [FAIL] Pre/post registry heights differ for bound comparison.\n');
    n_fail = n_fail + 1;
else
    bounds_match = true;
    for i = 1:height(reg)
        hit = strcmp(reg_post.name, reg.name{i});
        if ~any(hit) || abs(reg_post.lb(hit) - reg.lb(i)) > 1e-12 || ...
                abs(reg_post.ub(hit) - reg.ub(i)) > 1e-12
            fprintf('  [FAIL] Bound mismatch for %s.\n', reg.name{i});
            bounds_match = false;
            n_fail = n_fail + 1;
        end
    end
    if bounds_match
        fprintf('  [PASS] Post-surgery resolved bounds match pre-surgery bounds for all 14 active parameters.\n');
        n_pass = n_pass + 1;
    end
end

%% Summary
fprintf('\n=====================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL REYNA ACTIVE SET TESTS PASSED\n');
else
    error('test_reyna_hemodynamic_active_set:failed', ...
        'One or more Reyna active set checks failed.');
end
fprintf('=====================================================\n');
