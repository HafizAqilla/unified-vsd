%% test_scaling_mode_parity.m
% =========================================================================
% Scaling mode comparison test for Reyna pre-surgery calibration.
%
% PURPOSE:
%   Runs main_run under both Lundquist and Zhang scaling through the same
%   Reyna recipe/governance path. Lundquist is the accepted production path
%   and must keep single-digit primary governed RMSE (< 0.10). Zhang is an
%   independent comparator from its own scaled baseline; this test reports
%   its RMSE and the Lundquist/Zhang difference without requiring Zhang to
%   borrow the Lundquist disease seed or match the Lundquist acceptance band.
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
fprintf('  UNIFIED VSD MODEL - Scaling Mode Comparison Test\n');
fprintf('=====================================================\n\n');

run_heavy = strcmp(getenv('UNIFIED_VSD_RUN_HEAVY_TESTS'), '1');
if ~run_heavy
    fprintf('  [SKIP] Set UNIFIED_VSD_RUN_HEAVY_TESTS=1 to run this scaling-mode comparison.\n');
    fprintf('=====================================================\n');
    return;
end

n_pass = 0;
n_fail = 0;

old_gsa = getenv('UNIFIED_VSD_DO_GSA');
old_scaling = getenv('UNIFIED_VSD_SCALING_MODE');
old_fast = getenv('UNIFIED_VSD_FAST_CALIBRATION');
cleanup_env = onCleanup(@() restore_env(old_gsa, old_scaling, old_fast));

setenv('UNIFIED_VSD_DO_GSA', '0');
setenv('UNIFIED_VSD_FAST_CALIBRATION', '1');

modes = {'lundquist_bsa', 'zhang'};
pg_values = zeros(1, numel(modes));
recipe_ids = cell(1, numel(modes));
calibration_status = cell(1, numel(modes));

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
    manifest_file = fullfile(run_dir, 'run_manifest.txt');
    rmse_file = fullfile(run_dir, 'tables', 'validation_rmse_summary_pre_surgery.csv');
    if exist(manifest_file, 'file')
        manifest = fileread(manifest_file);
        recipe_ids{m} = lookup_manifest_value(manifest, 'CalibrationRecipeId');
        calibration_status{m} = lookup_manifest_value(manifest, 'CalibrationStatus');
    end
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

    if strcmp(mode, 'lundquist_bsa')
        if pg < 0.10
            fprintf('  [PASS] Lundquist primary governed RMSE %.6f < 0.10.\n', pg);
            n_pass = n_pass + 1;
        else
            fprintf('  [FAIL] Lundquist primary governed RMSE %.6f >= 0.10.\n', pg);
            n_fail = n_fail + 1;
        end
    elseif isfinite(pg)
        fprintf('  [INFO] Zhang independent-comparator RMSE %.6f; status=%s.\n', ...
            pg, string_or_unknown(calibration_status{m}));
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] Zhang primary governed RMSE is not finite.\n');
        n_fail = n_fail + 1;
    end
end

%% Comparison check: both modes must use the same explicit Reyna recipe.
if all(~cellfun(@isempty, recipe_ids)) && ...
        all(strcmp(recipe_ids, 'reyna_pre_surgery'))
    fprintf('  [PASS] Both scaling modes used the Reyna recipe governance.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Scaling modes did not report the same Reyna recipe governance.\n');
    n_fail = n_fail + 1;
end

%% Diagnostic only: RMSE difference between independent scaling modes.
if all(pg_values > 0)
    rmse_diff = abs(pg_values(1) - pg_values(2));
    fprintf('\n  Scaling mode RMSE difference: %.6f\n', rmse_diff);
    fprintf('  [INFO] Difference is diagnostic only; Zhang no longer reuses Lundquist seed.\n');
end

fprintf('\n=====================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL SCALING MODE COMPARISON CHECKS PASSED\n');
else
    fprintf('  SCALING MODE COMPARISON CHECKS FAILED\n');
end
fprintf('=====================================================\n');

function restore_env(old_gsa, old_scaling, old_fast)
setenv('UNIFIED_VSD_DO_GSA', old_gsa);
setenv('UNIFIED_VSD_FAST_CALIBRATION', old_fast);
setenv('UNIFIED_VSD_SCALING_MODE', old_scaling);
end

function value = lookup_manifest_value(manifest, key)
value = '';
pattern = sprintf('%s:\\s*([^\\r\\n]+)', key);
tokens = regexp(manifest, pattern, 'tokens', 'once');
if ~isempty(tokens)
    value = strtrim(tokens{1});
end
end

function text = string_or_unknown(value)
if isempty(value)
    text = 'unknown';
else
    text = value;
end
end
