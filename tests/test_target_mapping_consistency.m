%% test_target_mapping_consistency.m
% =========================================================================
% Unit tests for shared calibration/validation target mapping.
%
% PURPOSE:
%   Confirms that calibration targets are defined in one shared registry and
%   that validation consumes the same clinical field definitions.
%
% PASS CRITERIA:
%   1. Post-surgery SAP_min maps to SAP_dia_mmHg.
%   2. Post-surgery SAP_max maps to SAP_sys_mmHg.
%   3. Validation clinical values match get_calibration_targets().
%   4. GSA-guided primary selection returns up to five available metrics.
%
% USAGE:
%   >> test_target_mapping_consistency
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Target Mapping Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

clinical = patient_template();
clinical.post_surgery.SAP_dia_mmHg = 55;    % [mmHg]
clinical.post_surgery.SAP_sys_mmHg = 92;    % [mmHg]
clinical.post_surgery.MAP_mmHg     = 67;    % [mmHg]
clinical.post_surgery.QpQs         = 1.0;   % [-]
clinical.post_surgery.PVR_WU       = 2.1;   % [WU]
clinical.post_surgery.SVR_WU       = 12.0;  % [WU]
clinical.post_surgery.EF           = 0.62;  % [-]
clinical.post_surgery.RVEF         = 0.58;  % [-]
clinical.post_surgery.LVEDV_mL     = 45;    % [mL]
clinical.post_surgery.RVEDV_mL     = 50;    % [mL]
clinical.post_surgery.RAP_mean_mmHg = 4;    % [mmHg]
clinical.post_surgery.PAP_mean_mmHg = 16;   % [mmHg]
clinical.post_surgery.CO_Lmin      = 3.5;   % [L/min]

targets = get_calibration_targets('post_surgery', clinical);
target_names = {targets.Metric};

%% Test 1: Post-op SAP_min maps to diastolic pressure
idx_sap_min = find(strcmp(target_names, 'SAP_min'), 1, 'first');
if strcmp(targets(idx_sap_min).ClinicalField, 'SAP_dia_mmHg') && targets(idx_sap_min).ClinicalValue == 55
    fprintf('  [PASS] SAP_min maps to post.SAP_dia_mmHg.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] SAP_min mapping is incorrect.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Post-op SAP_max maps to systolic pressure
idx_sap_max = find(strcmp(target_names, 'SAP_max'), 1, 'first');
if strcmp(targets(idx_sap_max).ClinicalField, 'SAP_sys_mmHg') && targets(idx_sap_max).ClinicalValue == 92
    fprintf('  [PASS] SAP_max maps to post.SAP_sys_mmHg.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] SAP_max mapping is incorrect.\n');
    n_fail = n_fail + 1;
end

%% Test 3: validation_report consumes the shared target values
metrics = struct();
for k = 1:numel(targets)
    metrics.(targets(k).Metric) = targets(k).ClinicalValue;
end
report = validation_report(clinical, metrics, [], 'post_surgery', ...
    'PrimaryMetrics', {'QpQs', 'SAP_mean', 'PVR'});

row_min = find(strcmp(report.table_baseline.Metric, 'SAP_min'), 1, 'first');
row_max = find(strcmp(report.table_baseline.Metric, 'SAP_max'), 1, 'first');
if report.table_baseline.Clinical(row_min) == 55 && report.table_baseline.Clinical(row_max) == 92
    fprintf('  [PASS] validation_report uses shared SAP clinical values.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] validation_report clinical values do not match target registry.\n');
    n_fail = n_fail + 1;
end

%% Test 4: PVR/SVR are derived validation-only targets
idx_pvr = find(strcmp(target_names, 'PVR'), 1, 'first');
idx_svr = find(strcmp(target_names, 'SVR'), 1, 'first');
tiers = build_target_tiers(clinical, 'post_surgery');
tier_tbl = tiers.table;
pvr_tier_idx = find(strcmp(tier_tbl.Metric, 'PVR'), 1, 'first');
svr_tier_idx = find(strcmp(tier_tbl.Metric, 'SVR'), 1, 'first');
if ~targets(idx_pvr).UseForCalibration && ~targets(idx_svr).UseForCalibration && ...
        strcmp(targets(idx_pvr).Reliability, 'Derived') && ...
        strcmp(targets(idx_svr).Reliability, 'Derived') && ...
        strcmp(tier_tbl.Tier{pvr_tier_idx}, 'derived_validation') && ...
        strcmp(tier_tbl.Tier{svr_tier_idx}, 'derived_validation') && ...
        ~tier_tbl.IncludedInCalibration(pvr_tier_idx) && ...
        ~tier_tbl.IncludedInCalibration(svr_tier_idx) && ...
        ~tier_tbl.IncludedInPrimaryRMSE(pvr_tier_idx) && ...
        ~tier_tbl.IncludedInPrimaryRMSE(svr_tier_idx)
    fprintf('  [PASS] PVR/SVR are derived validation-only targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] PVR/SVR target governance is not derived validation-only.\n');
    n_fail = n_fail + 1;
end

%% Test 5: GSA-guided primary selector is deterministic and bounded
[primary_metrics, selection_table] = select_primary_metrics(clinical, [], 'post_surgery');
if numel(primary_metrics) <= 5 && all(ismember(primary_metrics, selection_table.Metric)) && ...
        ~any(ismember(primary_metrics, {'PVR','SVR'}))
    fprintf('  [PASS] select_primary_metrics returns a bounded reproducible set.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] select_primary_metrics returned an invalid set.\n');
    n_fail = n_fail + 1;
end

%% Test 6: Post-op BP echo derivatives do not double-count source evidence
clinical_echo = patient_template();
clinical_echo.common.HR = 108;                  % [bpm]
clinical_echo.post_surgery.LVESV_mL = 13.8;     % [mL] direct BP echo ESV
clinical_echo.post_surgery.EF = 0.618;          % [-] direct BP echo EF
clinical_echo.post_surgery.LVEDV_mL = ...
    clinical_echo.post_surgery.LVESV_mL / ...
    (1 - clinical_echo.post_surgery.EF);         % [mL] derived from ESV/EF
clinical_echo.post_surgery.CO_Lmin = ...
    (clinical_echo.post_surgery.LVEDV_mL - ...
    clinical_echo.post_surgery.LVESV_mL) * ...
    clinical_echo.common.HR / 1000;              % [L/min] echo-derived LVCO

tiers_echo = build_target_tiers(clinical_echo, 'post_surgery');
tier_tbl_echo = tiers_echo.table;
lvesv_tier_idx = find(strcmp(tier_tbl_echo.Metric, 'LVESV'), 1, 'first');
lvef_tier_idx = find(strcmp(tier_tbl_echo.Metric, 'LVEF'), 1, 'first');
lvedv_tier_idx = find(strcmp(tier_tbl_echo.Metric, 'LVEDV'), 1, 'first');
co_tier_idx = find(strcmp(tier_tbl_echo.Metric, 'CO_Lmin'), 1, 'first');

if strcmp(tier_tbl_echo.Tier{lvesv_tier_idx}, 'hard') && ...
        strcmp(tier_tbl_echo.Tier{lvef_tier_idx}, 'hard') && ...
        strcmp(tier_tbl_echo.Tier{lvedv_tier_idx}, 'derived_validation') && ...
        strcmp(tier_tbl_echo.Tier{co_tier_idx}, 'derived_validation') && ...
        tier_tbl_echo.IncludedInCalibration(lvesv_tier_idx) && ...
        tier_tbl_echo.IncludedInCalibration(lvef_tier_idx) && ...
        ~tier_tbl_echo.IncludedInCalibration(lvedv_tier_idx) && ...
        ~tier_tbl_echo.IncludedInCalibration(co_tier_idx)
    fprintf('  [PASS] Post-op BP echo source targets and derived checks are tiered correctly.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Post-op BP echo derived target governance is incorrect.\n');
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL TARGET MAPPING TESTS PASSED\n');
else
    error('test_target_mapping_consistency:failed', ...
          'One or more target mapping checks failed.');
end
fprintf('==========================================\n');
