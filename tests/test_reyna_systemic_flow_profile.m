%% test_reyna_systemic_flow_profile.m
% =========================================================================
% Regression test for Reyna raw-protocol systemic-flow calibration profile.
%
% PURPOSE:
%   Confirms that the active Reyna case uses protocol anthropometry,
%   direct PAP waveform evidence, PAP_mean-centered fitting, and the systemic preload-flow-volume
%   polish profile.
%
% USAGE:
%   >> test_reyna_systemic_flow_profile
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-11
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

fprintf('=====================================================\n');
fprintf('  UNIFIED VSD MODEL - Reyna Systemic Profile Test\n');
fprintf('=====================================================\n\n');

n_pass = 0;
n_fail = 0;

clinical = patient_reyna();
profile = build_case_calibration_profile(clinical, 'pre_surgery');
[primary_metrics, ~] = select_primary_metrics(clinical, struct(), 'pre_surgery', profile);

%% Test 1: raw protocol anthropometry remains active
if clinical.common.weight_kg == 13.4 && clinical.common.height_cm == 95.0 && ...
        abs(clinical.common.BSA - 0.588) < 1e-12
    fprintf('  [PASS] Reyna active anthropometry uses raw protocol values.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Reyna active anthropometry changed unexpectedly.\n');
    n_fail = n_fail + 1;
end

%% Test 2: pulmonary calibration keeps waveform evidence but fits PAP_mean
tiers = profile.targetTiers;
pap_min_idx = find(strcmp(tiers.table.Metric, 'PAP_min'), 1);
pap_max_idx = find(strcmp(tiers.table.Metric, 'PAP_max'), 1);
if clinical.pre_surgery.PAP_sys_mmHg == 20 && ...
        clinical.pre_surgery.PAP_dia_mmHg == 10 && ...
        clinical.pre_surgery.PAP_mean_mmHg == 15 && ...
        strcmp(tiers.table.Tier{pap_min_idx}, 'validation_only') && ...
        strcmp(tiers.table.Tier{pap_max_idx}, 'validation_only')
    fprintf('  [PASS] PAP_sys/PAP_dia are retained as validation evidence while PAP_mean anchors fitting.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] PAP waveform evidence/tiering is not configured as expected.\n');
    n_fail = n_fail + 1;
end

%% Test 3: primary metrics target systemic flow consistency
expected_primary = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
if isequal(primary_metrics(:)', expected_primary)
    fprintf('  [PASS] Primary metrics include CO_Lmin with systemic pressure targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Primary metrics do not match the Reyna systemic profile.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Reyna volume-flow guard respects target-tier governance
guard_ok = isfield(profile, 'volumeFlowGuard') && profile.volumeFlowGuard.enabled && ...
    isinf(profile.volumeFlowGuard.RVEDV_upper_ratio) && ...
    profile.volumeFlowGuard.RVEDV_upper_scale == 0 && ...
    profile.volumeFlowGuard.LVEF_upper_ratio == 1.10 && ...
    profile.volumeFlowGuard.CO_low_rel_thresh == 0.12;
if guard_ok
    fprintf('  [PASS] Reyna guard disables consistency-only RVEDV penalty while retaining EF/CO checks.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Reyna volume-flow guard governance is missing or changed.\n');
    n_fail = n_fail + 1;
end

%% Summary
fprintf('\n=====================================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL REYNA SYSTEMIC PROFILE TESTS PASSED\n');
else
    error('test_reyna_systemic_flow_profile:failed', ...
          'One or more Reyna systemic profile checks failed.');
end
fprintf('=====================================================\n');
