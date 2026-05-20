%% test_baseline_reference_metrics.m
% =========================================================================
% Regression test for the Keisya/Hafiz adult and pediatric baseline snapshot.
%
% PURPOSE:
%   Confirms that the downloaded default-parameter baseline and the
%   screenshot/SHEET.xlsx metric snapshot are available to the codebase.
%
% USAGE:
%   >> test_baseline_reference_metrics
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-20
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

fprintf('=====================================================\n');
fprintf('  UNIFIED VSD MODEL - Baseline Reference Test\n');
fprintf('=====================================================\n\n');

n_pass = 0;
n_fail = 0;

params = default_parameters();

%% Test 1: default_parameters carries the Keisya baseline values
value_checks = {
    'V0.RA',       params.V0.RA,       3.5385
    'V0.RV',       params.V0.RV,       8.4067
    'V0.LA',       params.V0.LA,       2.3085
    'V0.LV',       params.V0.LV,       3.5385
    'E.LV.EA',     params.E.LV.EA,     3.5
    'E.LV.EB',     params.E.LV.EB,     0.08
    'E.RV.EA',     params.E.RV.EA,     0.5
    'E.RV.EB',     params.E.RV.EB,     0.042
    'E.LA.EA',     params.E.LA.EA,     0.35
    'E.RA.EA',     params.E.RA.EA,     0.30
    'Tc_LV_frac',  params.Tc_LV_frac,  0.265
    'Tr_LV_frac',  params.Tr_LV_frac,  0.40
    'Tr_RV_frac',  params.Tr_RV_frac,  0.40
    'C.SAR',       params.C.SAR,       1.33
    'R.SC',        params.R.SC,        0.80
    'ic.V_RV',     params.ic.V(params.idx.V_RV), 120.0
    };

tolerance = 1e-12;
all_values_match = true;
for idx = 1:size(value_checks, 1)
    if abs(value_checks{idx, 2} - value_checks{idx, 3}) > tolerance
        all_values_match = false;
        fprintf('  [FAIL] %s = %.12g, expected %.12g\n', ...
            value_checks{idx, 1}, value_checks{idx, 2}, value_checks{idx, 3});
    end
end
if all_values_match
    fprintf('  [PASS] default_parameters contains the Keisya baseline parameter values.\n');
    n_pass = n_pass + 1;
else
    n_fail = n_fail + 1;
end

%% Test 2: screenshot/SHEET.xlsx baseline snapshot is available
reference = baseline_reference_metrics();
expected_columns = {'Metric','Unit','Norm_Adult','Adult_ref', ...
    'Norm_Pediatric','Zhang_2019','Lundquist_2025'};
if height(reference) == 26 && all(ismember(expected_columns, reference.Properties.VariableNames))
    fprintf('  [PASS] baseline_reference_metrics exposes the expected snapshot table.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] baseline_reference_metrics table shape changed unexpectedly.\n');
    n_fail = n_fail + 1;
end

%% Test 3: key screenshot values are preserved exactly enough for audit
metric_checks = {
    'Cardiac Output', 'Adult_ref',      5.38797099834
    'Cardiac Output', 'Zhang_2019',     2.15426137243
    'Cardiac Output', 'Lundquist_2025', 1.29231633427
    'RVEF',           'Adult_ref',      60.3230976656
    'RVEDV',          'Zhang_2019',     31.885535161
    'MAP',            'Lundquist_2025', 58.2849727943
    'Qp_Qs',          'Adult_ref',      1.00052468422
    };

all_metrics_match = true;
for idx = 1:size(metric_checks, 1)
    row = strcmp(reference.Metric, metric_checks{idx, 1});
    value = reference.(metric_checks{idx, 2})(row);
    if isempty(value) || abs(value - metric_checks{idx, 3}) > 1e-10
        all_metrics_match = false;
        fprintf('  [FAIL] %s.%s = %.12g, expected %.12g\n', ...
            metric_checks{idx, 1}, metric_checks{idx, 2}, ...
            value, metric_checks{idx, 3});
    end
end
if all_metrics_match
    fprintf('  [PASS] key screenshot baseline values are preserved.\n');
    n_pass = n_pass + 1;
else
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n=====================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL BASELINE REFERENCE TESTS PASSED\n');
else
    error('test_baseline_reference_metrics:failed', ...
        'One or more baseline reference checks failed.');
end
fprintf('=====================================================\n');
