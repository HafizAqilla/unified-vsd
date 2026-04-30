%% test_scaling_ic_separation.m
% Verifies that physiological scaling and initial-condition construction
% are separated: IC generation may change params.ic.V but not params.V0.*.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Scaling/IC Separation Test\n');
fprintf('==========================================\n\n');

params_ref = default_parameters();
patient = struct('age_years', 0.5, 'age_days', 183, 'weight_kg', 7.0, ...
    'height_cm', 65, 'sex', 'M', 'maturation_mode', 'normal');

params_phys = apply_physiological_scaling(params_ref, patient, 'zhang');
V0_before = params_phys.V0;
ic = build_initial_conditions(params_phys, patient);

assert(isequaln(V0_before, params_phys.V0), 'V0 changed during IC construction.');
assert(numel(ic) == 14, 'IC vector must have 14 states.');

fprintf('[PASS] V0 is unchanged by build_initial_conditions().\n');
fprintf('[PASS] IC vector constructed with 14 states.\n');
