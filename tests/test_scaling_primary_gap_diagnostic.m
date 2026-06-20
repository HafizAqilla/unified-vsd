%% test_scaling_primary_gap_diagnostic.m
% Smoke test for Zhang-vs-Lundquist primary gap diagnostic dry run.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

result = compare_scaling_primary_gap('PatientID', 'P1', 'Mode', 'dry_run');

assert(strcmp(result.PatientID, 'P1'), 'PatientID mismatch.');
assert(strcmp(result.Mode, 'dry_run'), 'Mode mismatch.');
assert(strcmp(result.Classification, 'mixed_or_inconclusive'), ...
    'Dry-run classification should be inconclusive.');

fprintf('  [PASS] scaling primary gap diagnostic dry run is safe.\n');

