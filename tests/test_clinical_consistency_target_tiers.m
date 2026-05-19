%% test_clinical_consistency_target_tiers.m
% =========================================================================
% Regression checks for clinical consistency audit and target-tier RMSE.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-14
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Target Governance Test\n');
fprintf('==========================================\n\n');

clinical = patient_reyna();
scenario = 'pre_surgery';
audit = audit_clinical_consistency(clinical, scenario);
tiers = build_target_tiers(clinical, scenario, audit);

assert(abs(audit.SV_Qs - clinical.pre_surgery.CO_Lmin * 1000 / clinical.common.HR) < 1e-10, ...
    'SV_Qs calculation mismatch.');
assert(abs(audit.SV_Qp - clinical.pre_surgery.CO_Lmin * clinical.pre_surgery.QpQs * 1000 / clinical.common.HR) < 1e-10, ...
    'SV_Qp calculation mismatch.');
assert(strcmp(audit.severity, 'critical'), ...
    'Reyna stroke-volume inconsistency should be critical.');
assert(ismember('RVEDV', tiers.consistency_only), ...
    'RVEDV must be downgraded to consistency-check-only for Reyna.');
assert(~ismember('RVEDV', tiers.included_in_calibration), ...
    'RVEDV must not be included in calibration targets when flagged.');
assert(~ismember('RVEDV', tiers.included_in_primary_rmse), ...
    'RVEDV must not be included in governed primary RMSE when flagged.');
assert(ismember('LVEDV', tiers.consistency_only), ...
    'LVEDV must be downgraded to consistency-check-only for Reyna when LV stroke volume conflicts.');
assert(ismember('LVESV', tiers.consistency_only), ...
    'LVESV must be downgraded to consistency-check-only for Reyna when LV stroke volume conflicts.');
assert(~ismember('LVEDV', tiers.included_in_calibration), ...
    'LVEDV must not be included in calibration targets when flagged.');
assert(~ismember('LVESV', tiers.included_in_calibration), ...
    'LVESV must not be included in calibration targets when flagged.');
assert(~ismember('LVEDV', tiers.included_in_primary_rmse), ...
    'LVEDV must not be included in governed primary RMSE when flagged.');
assert(~ismember('LVESV', tiers.included_in_primary_rmse), ...
    'LVESV must not be included in governed primary RMSE when flagged.');
assert(ismember('LVEF', tiers.consistency_only), ...
    'EF/LVEF must be consistency-check-only when LVEDV and LVESV are already available.');
assert(~ismember('LVEF', tiers.included_in_calibration), ...
    'Derived EF must not be independently fitted with its source volumes.');

targets = get_calibration_targets(scenario, clinical);
metrics = struct();
for idx = 1:numel(targets)
    if isfinite(targets(idx).ClinicalValue)
        metrics.(targets(idx).Metric) = targets(idx).ClinicalValue;
    end
end
metrics.RVEDV = clinical.pre_surgery.RVEDV_mL * 2;

report = validation_report(clinical, metrics, metrics, scenario, ...
    'PrimaryMetrics', {'CO_Lmin','QpQs','PAP_mean','SAP_mean','RAP_mean'}, ...
    'TargetTiers', tiers, ...
    'ClinicalConsistencyAudit', audit);

assert(report.rmse_full_cal > 0, ...
    'Full RMSE must include the deliberately perturbed RVEDV.');
assert(report.rmse_primary_cal < 1e-12, ...
    'Primary RMSE must exclude deliberately perturbed consistency-only targets.');
assert(isfield(report, 'clinical_data_rank') && ~isempty(report.clinical_data_rank), ...
    'Validation report must expose clinical data availability ranking.');
assert(any(strcmp(report.clinical_data_rank.Field, 'EF') & ...
    strcmp(report.clinical_data_rank.EvidenceClass, 'directly_derived')), ...
    'EF should be ranked as directly derived clinical evidence.');

fprintf('\n[PASS] Clinical audit, target tiers, and RMSE masks behave as expected.\n');
