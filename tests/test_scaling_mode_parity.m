%% test_scaling_mode_parity.m
% =========================================================================
% Scaling mode parity test for Reyna pre-surgery calibration.
%
% PURPOSE:
%   Runs main_run under both Lundquist and Zhang scaling and asserts that
%   both produce single-digit primary governed RMSE (< 0.10) with the
%   absolute RMSE difference < 0.05.
%
%   SKIP BEHAVIOR: the test skips when UNIFIED_VSD_RUN_HEAVY_TESTS is not
%   set to '1', keeping the default test suite cheap and fast.
%
% USAGE:
%   >> test_scaling_mode_parity              % skips without env var
%   >> setenv('UNIFIED_VSD_RUN_HEAVY_TESTS','1');
%   >> test_scaling_mode_parity              % runs both scaling modes
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
fprintf('  UNIFIED VSD MODEL - Scaling Mode Parity Test\n');
fprintf('=====================================================\n\n');

run_heavy = strcmp(getenv('UNIFIED_VSD_RUN_HEAVY_TESTS'), '1');
if ~run_heavy
    fprintf('  [SKIP] Set UNIFIED_VSD_RUN_HEAVY_TESTS=1 to run this scaling-mode parity check.\n');
    fprintf('=====================================================\n');
    return;
end

n_pass = 0;
n_fail = 0;

old_gsa = getenv('UNIFIED_VSD_DO_GSA');
old_scaling = getenv('UNIFIED_VSD_SCALING_MODE');
old_fast = getenv('UNIFIED_VSD_FAST_CALIBRATION');

setenv('UNIFIED_VSD_DO_GSA', '0');
setenv('UNIFIED_VSD_FAST_CALIBRATION', '1');

modes = {'lundquist_bsa', 'zhang'};
pg_values = zeros(1, numel(modes));

for m = 1:numel(modes)
    mode = modes{m};
    setenv('UNIFIED_VSD_SCALING_MODE', mode);
    
    scan_before = dir(fullfile(project_root, 'results', 'runs', '*_reyna_pre_surgery'));
    main_run('pre_surgery', patient_reyna());
    scan_after = dir(fullfile(project_root, 'results', 'runs', '*_reyna_pre_surgery'));
    
    new_dirs = setdiff({scan_after.name}, {scan_before.name}, 'stable');
    if isempty(new_dirs)
        fprintf('  [FAIL] %s: no new run folder found.\n', mode);
        n_fail = n_fail + 1;
        continue;
    end
    
    run_dir = fullfile(project_root, 'results', 'runs', new_dirs{end});
    rmse_file = fullfile(run_dir, 'tables', 'validation_rmse_summary_pre_surgery.csv');
    if ~exist(rmse_file, 'file')
        fprintf('  [FAIL] %s: RMSE summary CSV not found.\n', mode);
        n_fail = n_fail + 1;
        continue;
    end
    
    opts = detectImportOptions(rmse_file);
    opts.VariableNamesLine = 1;
    rmse_table = readtable(rmse_file, opts);
    row = strcmp(rmse_table.RMSE_Type, 'primary_governed');
    
    if ~any(row)
        fprintf('  [FAIL] %s: primary_governed row not found.\n', mode);
        n_fail = n_fail + 1;
        continue;
    end
    
    pg = rmse_table.Calibrated(row);
    pg_values(m) = pg;
    fprintf('  %s primary governed RMSE (calibrated): %.6f\n', mode, pg);
    
    if pg < 0.10
        fprintf('  [PASS] %s primary governed RMSE %.6f < 0.10.\n', mode, pg);
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] %s primary governed RMSE %.6f >= 0.10.\n', mode, pg);
        n_fail = n_fail + 1;
    end
end

%% Parity check: RMSE difference between modes must be < 0.05
if all(pg_values > 0)
    rmse_diff = abs(pg_values(1) - pg_values(2));
    fprintf('\n  Scaling mode RMSE difference: %.6f\n', rmse_diff);
    if rmse_diff < 0.05
        fprintf('  [PASS] Lundquist/Zhang RMSE parity within tolerance (%.6f < 0.05).\n', rmse_diff);
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] Lundquist/Zhang RMSE parity outside tolerance.\n');
        n_fail = n_fail + 1;
    end
end

setenv('UNIFIED_VSD_DO_GSA', old_gsa);
setenv('UNIFIED_VSD_FAST_CALIBRATION', old_fast);
setenv('UNIFIED_VSD_SCALING_MODE', old_scaling);

fprintf('\n=====================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL SCALING MODE PARITY TESTS PASSED\n');
else
    fprintf('  SCALING MODE PARITY TESTS FAILED\n');
end
fprintf('=====================================================\n');
