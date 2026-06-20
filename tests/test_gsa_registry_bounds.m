%% test_gsa_registry_bounds.m
% Unit tests for registry-backed GSA bounds.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

fprintf('===============================================\n');
fprintf('  UNIFIED VSD MODEL - GSA Registry Bounds Test\n');
fprintf('===============================================\n\n');

params_ref = default_parameters();
clinical = patient_reyna();
patient = clinical.common;
patient.scaling_mode = 'zhang';
params0 = apply_scaling(params_ref, patient);
params0 = params_from_clinical(params0, clinical, 'pre_surgery', params0);

names = {'R.SAR','R.SC','C.SAR','E.LV.EA','V0.LV','R.vsd'};
bounds = build_gsa_registry_bounds(params0, 'pre_surgery', names);

assert(numel(bounds.names) == numel(names), 'Bounds name count mismatch.');
assert(all(isfinite(bounds.x0)), 'All x0 values should be finite.');
assert(all(isfinite(bounds.lb)), 'All lower bounds should be finite.');
assert(all(isfinite(bounds.ub)), 'All upper bounds should be finite.');
assert(all(bounds.lb < bounds.ub), 'Each lower bound must be less than upper bound.');
assert(all(ismember({'name','x0','lb','ub','source','note'}, ...
    bounds.table.Properties.VariableNames)), 'Bounds table schema mismatch.');
assert(any(bounds.table.source == "registry"), ...
    'Expected at least one registry-backed GSA bound.');

fprintf('  [PASS] registry-backed GSA bounds are finite and traceable.\n');

