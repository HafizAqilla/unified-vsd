%% test_optimization_mask.m
% =========================================================================
% Unit tests for create_optimization_mask.m and mask integration contract.
%
% PURPOSE:
%   Verifies that Sobol-based screening produces a stable logical mask for
%   optimisation and protects against an empty active-parameter set.
%
% PASS CRITERIA:
%   1. Vector ST: thresholding behaves as expected.
%   2. Matrix ST: per-parameter max aggregation across metrics is used.
%   3. All-below-threshold input still activates one parameter (safety).
%   4. Invalid threshold raises an error.
%   5. Mask-size mismatch in calibration_param_sets raises an error.
%
% USAGE:
%   >> test_optimization_mask
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-16
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, '..')));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Optimization Mask Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% Test 1: Vector ST thresholding
fprintf('--- Test 1: Vector ST thresholding ---\n');
ST_vec = [0.02; 0.11; 0.09; 0.20];
mask_1 = create_optimization_mask(ST_vec, 0.10);
exp_1  = logical([0; 1; 0; 1]);
if isequal(mask_1, exp_1)
    fprintf('  [PASS] Vector thresholding matches expected mask.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Vector thresholding mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Matrix ST aggregation by max across metrics
fprintf('--- Test 2: Matrix ST aggregation ---\n');
ST_mat = [0.05 0.08; 0.02 0.20; 0.09 0.03];
mask_2 = create_optimization_mask(ST_mat, 0.10);
exp_2  = logical([0; 1; 0]);
if isequal(mask_2, exp_2)
    fprintf('  [PASS] Matrix aggregation uses max(ST_i,j) as intended.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Matrix aggregation mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 3: Empty active-set safety guard
fprintf('--- Test 3: Empty active-set safety guard ---\n');
ST_low = [0.01; 0.02; 0.03];
mask_3 = create_optimization_mask(ST_low, 0.10);
if sum(mask_3) == 1 && mask_3(3)
    fprintf('  [PASS] Safety guard activates max-ST parameter when all are below threshold.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Safety guard did not activate expected parameter.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Invalid threshold should throw
fprintf('--- Test 4: Invalid threshold input ---\n');
threw_4 = false;
try
    create_optimization_mask(ST_vec, 1.5);
catch
    threw_4 = true;
end
if threw_4
    fprintf('  [PASS] Invalid threshold correctly throws an error.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Invalid threshold did not throw an error.\n');
    n_fail = n_fail + 1;
end

%% Test 5: calibration_param_sets mask-size mismatch should throw
fprintf('--- Test 5: calibration_param_sets mask-size mismatch ---\n');
threw_5 = false;
try
    params0 = default_parameters();
    badMask = true(3, 1);
    calibration_param_sets('pre_surgery', params0, badMask);
catch
    threw_5 = true;
end
if threw_5
    fprintf('  [PASS] Mask-size mismatch correctly throws an error.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Mask-size mismatch did not throw an error.\n');
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL OPTIMIZATION MASK TESTS PASSED\n');
else
    fprintf('  ONE OR MORE OPTIMIZATION MASK TESTS FAILED\n');
end
fprintf('==========================================\n');
