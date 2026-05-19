%% test_calibration_classification.m
% =========================================================================
% Unit tests for classify_calibration_run.m
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

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Calibration Classification Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% Shared mock report
report = struct();
report.primary_metrics = {'QpQs','PAP_mean'};
report.primary_gate = table( ...
    {'QpQs';'PAP_mean'}, [1.1; 15], [1.1; 15], [2; 4], [true; true], ...
    'VariableNames', {'Metric','Clinical','Model','AbsError_pct','Pass_5pct'});
report.table_cal = table( ...
    {'QpQs';'PAP_mean';'SAP_max'}, [2; 4; 8], ...
    'VariableNames', {'Metric','Error_pct'});
report.rmse_baseline = 0.20;
report.rmse_cal = 0.12;

%% Test 1: ACCEPT
fprintf('--- Test 1: ACCEPT classification ---\n');
plaus_ok = struct('n_warning', 0, 'n_fail', 0, 'warning_fraction', 0);
status_1 = classify_calibration_run(report, plaus_ok);
if strcmp(status_1.label, 'ACCEPT')
    fprintf('  [PASS] Good fit + good plausibility -> ACCEPT.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Expected ACCEPT, got %s.\n', status_1.label);
    n_fail = n_fail + 1;
end

%% Test 2: OUTPUT_FIT_ONLY
fprintf('--- Test 2: OUTPUT_FIT_ONLY classification ---\n');
plaus_warn = struct('n_warning', 2, 'n_fail', 1, 'warning_fraction', 0.25);
status_2 = classify_calibration_run(report, plaus_warn);
if strcmp(status_2.label, 'OUTPUT_FIT_ONLY')
    fprintf('  [PASS] Good fit + poor plausibility -> OUTPUT_FIT_ONLY.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Expected OUTPUT_FIT_ONLY, got %s.\n', status_2.label);
    n_fail = n_fail + 1;
end

%% Test 3: PROMISING_NEAR_MISS
fprintf('--- Test 3: PROMISING_NEAR_MISS classification ---\n');
report_poor = report;
report_poor.primary_gate.Pass_5pct(2) = false;
report_poor.primary_gate.AbsError_pct(2) = 7;
plaus_near = struct('n_warning', 4, 'n_fail', 0, 'warning_fraction', 0.35);
status_3 = classify_calibration_run(report_poor, plaus_near);
if strcmp(status_3.label, 'PROMISING_NEAR_MISS')
    fprintf('  [PASS] Near-miss fit + strong improvement -> PROMISING_NEAR_MISS.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Expected PROMISING_NEAR_MISS, got %s.\n', status_3.label);
    n_fail = n_fail + 1;
end

%% Test 4: PHYSIOLOGICAL_BUT_POOR_FIT
fprintf('--- Test 4: PHYSIOLOGICAL_BUT_POOR_FIT classification ---\n');
report_structural = report_poor;
report_structural.rmse_cal = 0.19;
report_structural.table_cal.Error_pct(3) = 28;
status_4 = classify_calibration_run(report_structural, plaus_ok);
if strcmp(status_4.label, 'PHYSIOLOGICAL_BUT_POOR_FIT')
    fprintf('  [PASS] Poor fit + plausible but weak improvement -> PHYSIOLOGICAL_BUT_POOR_FIT.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Expected PHYSIOLOGICAL_BUT_POOR_FIT, got %s.\n', status_4.label);
    n_fail = n_fail + 1;
end

%% Test 5: REJECT
fprintf('--- Test 5: REJECT classification ---\n');
status_5 = classify_calibration_run(report_structural, plaus_warn);
if strcmp(status_5.label, 'REJECT')
    fprintf('  [PASS] Poor fit + poor plausibility -> REJECT.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Expected REJECT, got %s.\n', status_5.label);
    n_fail = n_fail + 1;
end

%% Test 6: Large miss remains PHYSIOLOGICAL_BUT_POOR_FIT, not near miss
fprintf('--- Test 6: Large miss is not near miss ---\n');
report_large_miss = report;
report_large_miss.primary_gate.Pass_5pct = [false; true];
report_large_miss.primary_gate.AbsError_pct = [22; 2];
report_large_miss.rmse_cal = 0.13;
status_6 = classify_calibration_run(report_large_miss, plaus_ok);
if strcmp(status_6.label, 'PHYSIOLOGICAL_BUT_POOR_FIT')
    fprintf('  [PASS] Large primary miss stays out of PROMISING_NEAR_MISS.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Expected PHYSIOLOGICAL_BUT_POOR_FIT, got %s.\n', status_6.label);
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL CALIBRATION CLASSIFICATION TESTS PASSED\n');
else
    fprintf('  ONE OR MORE CALIBRATION CLASSIFICATION TESTS FAILED\n');
end
fprintf('==========================================\n');
