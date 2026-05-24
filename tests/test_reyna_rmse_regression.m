%% test_reyna_rmse_regression.m
% =========================================================================
% End-to-end RMSE regression test for Reyna pre-surgery calibration.
%
% PURPOSE:
%   Runs a full main_run with no-GSA fast mode and asserts single-digit
%   primary governed RMSE for the Reyna pre-surgery pipeline.
%
%   SKIP BEHAVIOR: the test skips when UNIFIED_VSD_RUN_HEAVY_TESTS is not
%   set to '1', keeping the default test suite cheap and fast.
%
% USAGE:
%   >> test_reyna_rmse_regression           % skips without env var
%   >> setenv('UNIFIED_VSD_RUN_HEAVY_TESTS','1');
%   >> test_reyna_rmse_regression           % runs full pipeline
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
fprintf('  UNIFIED VSD MODEL - Reyna RMSE Regression Test\n');
fprintf('=====================================================\n\n');

run_heavy = strcmp(getenv('UNIFIED_VSD_RUN_HEAVY_TESTS'), '1');
if ~run_heavy
    fprintf('  [SKIP] Set UNIFIED_VSD_RUN_HEAVY_TESTS=1 to run this end-to-end RMSE regression.\n');
    fprintf('=====================================================\n');
    return;
end

n_pass = 0;
n_fail = 0;

%% Run Reyna pre-surgery with Lundquist scaling (no GSA, fast mode)
old_gsa = getenv('UNIFIED_VSD_DO_GSA');
old_scaling = getenv('UNIFIED_VSD_SCALING_MODE');
old_fast = getenv('UNIFIED_VSD_FAST_CALIBRATION');

setenv('UNIFIED_VSD_DO_GSA', '0');
setenv('UNIFIED_VSD_FAST_CALIBRATION', '1');
setenv('UNIFIED_VSD_SCALING_MODE', 'lundquist_bsa');

scan_before = dir(fullfile(project_root, 'results', 'runs', '*_reyna_pre_surgery'));
main_run('pre_surgery', patient_reyna());
scan_after = dir(fullfile(project_root, 'results', 'runs', '*_reyna_pre_surgery'));

new_dirs = setdiff({scan_after.name}, {scan_before.name}, 'stable');
if isempty(new_dirs)
    fprintf('  [FAIL] No new run folder found after main_run.\n');
    n_fail = n_fail + 1;
else
    run_dir = fullfile(project_root, 'results', 'runs', new_dirs{end});
    rmse_file = fullfile(run_dir, 'tables', 'validation_rmse_summary_pre_surgery.csv');
    if exist(rmse_file, 'file')
        opts = detectImportOptions(rmse_file);
        opts.VariableNamesLine = 1;
        rmse_table = readtable(rmse_file, opts);
        row = strcmp(rmse_table.RMSE_Type, 'primary_governed');
        if any(row)
            pg_value = rmse_table.Calibrated(row);
            fprintf('  Primary governed RMSE (calibrated): %.6f\n', pg_value);
            if pg_value < 0.10
                fprintf('  [PASS] Primary governed RMSE %.6f < 0.10.\n', pg_value);
                n_pass = n_pass + 1;
            else
                fprintf('  [FAIL] Primary governed RMSE %.6f >= 0.10.\n', pg_value);
                n_fail = n_fail + 1;
            end
        else
            fprintf('  [FAIL] primary_governed row not found in RMSE summary.\n');
            n_fail = n_fail + 1;
        end
    else
        fprintf('  [FAIL] RMSE summary CSV not found at %s.\n', rmse_file);
        n_fail = n_fail + 1;
    end
end

setenv('UNIFIED_VSD_DO_GSA', old_gsa);
setenv('UNIFIED_VSD_FAST_CALIBRATION', old_fast);
setenv('UNIFIED_VSD_SCALING_MODE', old_scaling);

fprintf('\n=====================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL REYNA RMSE REGRESSION TESTS PASSED\n');
else
    fprintf('  REYNA RMSE REGRESSION TESTS FAILED\n');
end
fprintf('=====================================================\n');
