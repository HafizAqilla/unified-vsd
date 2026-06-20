%% test_baseline_plausibility_report.m
% Unit tests for baseline plausibility gate/reporting.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'config'), '-begin');
addpath(fullfile(project_root, 'src', 'utils'), '-begin');

fprintf('====================================================\n');
fprintf('  UNIFIED VSD MODEL - Baseline Plausibility Test\n');
fprintf('====================================================\n\n');

n_pass = 0;
n_fail = 0;
clinical = patient_reyna();
case_profile = build_case_calibration_profile(clinical, 'pre_surgery');

%% Test 1: A broad plausible baseline passes.
metrics = struct( ...
    'HR', 95, ...
    'RAP_mean', 5, ...
    'PAP_mean', 18, ...
    'SAP_mean', 70, ...
    'CO_Lmin', 3.2, ...
    'QpQs', 1.4, ...
    'LVEF', 0.62, ...
    'RVEF', 0.55);
policy = resolve_scaling_policy('', 'publication');
report = baseline_plausibility_report(metrics, clinical, 'pre_surgery', case_profile, policy);
if strcmp(report.status, 'pass') && report.allowedDownstream
    fprintf('  [PASS] Plausible publication baseline passes.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Plausible baseline did not pass.\n');
    n_fail = n_fail + 1;
end

%% Test 2: A hard-failing publication primary baseline blocks downstream.
metrics_bad = metrics;
metrics_bad.CO_Lmin = 0.05;
report_bad = baseline_plausibility_report(metrics_bad, clinical, 'pre_surgery', case_profile, policy);
if strcmp(report_bad.status, 'hard_fail') && ~report_bad.allowedDownstream && ...
        any(strcmp(report_bad.failedMetrics, 'CO_Lmin'))
    fprintf('  [PASS] Hard-failing publication baseline blocks downstream.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Hard-failing baseline gate mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 3: A hard-failing comparator is labeled comparator-only.
policy_cmp = resolve_scaling_policy('lundquist_bsa', 'publication');
report_cmp = baseline_plausibility_report(metrics_bad, clinical, 'pre_surgery', case_profile, policy_cmp);
if strcmp(report_cmp.status, 'comparator_only') && ~report_cmp.allowedDownstream
    fprintf('  [PASS] Hard-failing comparator is labeled comparator-only.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Comparator-only gate mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Report exports CSV when ResultsDir is supplied.
out_dir = fullfile(tempdir(), 'unified_vsd_baseline_plausibility_test');
if exist(out_dir, 'dir')
    rmdir(out_dir, 's');
end
mkdir(out_dir);
report_file = baseline_plausibility_report(metrics, clinical, 'pre_surgery', ...
    case_profile, policy, 'ResultsDir', out_dir);
if isfield(report_file, 'tableFile') && isfile(report_file.tableFile)
    fprintf('  [PASS] Baseline plausibility report writes CSV.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Baseline plausibility CSV was not written.\n');
    n_fail = n_fail + 1;
end

fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
fprintf('==========================================\n');

assert(n_fail == 0, 'test_baseline_plausibility_report:failed', ...
    '%d baseline plausibility test(s) failed.', n_fail);

