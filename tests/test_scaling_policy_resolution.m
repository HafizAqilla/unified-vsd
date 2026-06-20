%% test_scaling_policy_resolution.m
% Unit tests for publication/exploratory scaling role policy.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'config'), '-begin');
addpath(fullfile(project_root, 'src', 'utils'), '-begin');

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Scaling Policy Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% Test 1: Publication mode defaults to Zhang primary prior.
policy = resolve_scaling_policy('', 'publication');
if strcmp(policy.ScalingMode, 'zhang') && strcmp(policy.ScalingRole, 'primary_prior')
    fprintf('  [PASS] Publication mode defaults to Zhang primary prior.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Publication default policy mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Exploratory default preserves legacy Lundquist path.
policy = resolve_scaling_policy('', 'exploratory');
if strcmp(policy.ScalingMode, 'lundquist_bsa') && strcmp(policy.ScalingRole, 'exploratory')
    fprintf('  [PASS] Exploratory mode defaults to Lundquist exploratory.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Exploratory default policy mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 3: Lundquist publication primary requires override rationale.
did_error = false;
try
    resolve_scaling_policy('lundquist_bsa', 'publication', ...
        'ScalingRole', 'primary_prior');
catch ME
    did_error = strcmp(ME.identifier, ...
        'resolve_scaling_policy:missingOverrideRationale');
end
if did_error
    fprintf('  [PASS] Lundquist primary publication run requires rationale.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Lundquist primary publication run did not require rationale.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Explicit Lundquist comparator is allowed in publication mode.
policy = resolve_scaling_policy('lundquist', 'publication');
if strcmp(policy.ScalingMode, 'lundquist_bsa') && strcmp(policy.ScalingRole, 'comparator')
    fprintf('  [PASS] Lundquist publication comparator is allowed.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Lundquist comparator policy mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 5: apply_scaling honors publication default through env.
old_run_mode = getenv('UNIFIED_VSD_RUN_MODE');
old_scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
setenv('UNIFIED_VSD_RUN_MODE', 'publication');
setenv('UNIFIED_VSD_SCALING_MODE', '');
params_ref = default_parameters();
clinical = patient_reyna();
patient = clinical.common;
params_pub = apply_scaling(params_ref, patient);
setenv('UNIFIED_VSD_RUN_MODE', old_run_mode);
setenv('UNIFIED_VSD_SCALING_MODE', old_scaling_mode);

if strcmp(params_pub.scaling.mode, 'zhang') && ...
        strcmp(params_pub.scaling.role, 'primary_prior')
    fprintf('  [PASS] apply_scaling uses Zhang primary in publication mode.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] apply_scaling publication policy mismatch.\n');
    n_fail = n_fail + 1;
end

fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
fprintf('==========================================\n');

assert(n_fail == 0, 'test_scaling_policy_resolution:failed', ...
    '%d scaling policy test(s) failed.', n_fail);

