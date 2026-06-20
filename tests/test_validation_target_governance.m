%% test_validation_target_governance.m
% Ensures derived/holdout targets are not silently promoted to primary.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

fprintf('================================================\n');
fprintf('  UNIFIED VSD MODEL - Validation Governance Test\n');
fprintf('================================================\n\n');

clinical = patient_reyna();
scenario = 'pre_surgery';
profile = build_case_calibration_profile(clinical, scenario);

metrics = struct();
metrics.RAP_mean = 5;
metrics.PAP_mean = 15;
metrics.SAP_mean = 71;
metrics.QpQs = 1.2;
metrics.CO_Lmin = 3.2;
metrics.PVR = 6;
metrics.SVR = 20;

report = validation_report(clinical, metrics, metrics, scenario, ...
    'PrimaryMetrics', {'QpQs','SAP_mean','PVR'}, ...
    'CaseProfile', profile, ...
    'TargetTiers', profile_target_tiers_for_test(profile), ...
    'ClinicalConsistencyAudit', profile_audit_for_test(profile));

assert(any(strcmp(report.primary_metrics, 'QpQs')), 'QpQs should remain primary.');
assert(any(strcmp(report.primary_metrics, 'SAP_mean')), 'SAP_mean should remain primary.');
assert(~any(strcmp(report.primary_metrics, 'PVR')), ...
    'PVR should not be silently promoted to primary.');
assert(any(strcmp(report.primary_metric_governance.blockedMetrics, 'PVR')), ...
    'PVR should be listed as blocked.');

fprintf('  [PASS] validation_report blocks derived/holdout primary promotion.\n');

function target_tiers = profile_target_tiers_for_test(profile)
if isfield(profile, 'targetTiers')
    target_tiers = profile.targetTiers;
else
    target_tiers = [];
end
end

function audit = profile_audit_for_test(profile)
if isfield(profile, 'clinicalConsistencyAudit')
    audit = profile.clinicalConsistencyAudit;
else
    audit = [];
end
end

