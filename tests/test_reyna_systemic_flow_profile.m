%% test_reyna_systemic_flow_profile.m
% =========================================================================
% Regression test for Reyna hemodynamic-only pre-surgery calibration profile.
%
% PURPOSE:
%   Confirms that the active Reyna case uses protocol anthropometry,
%   direct PAP waveform evidence, PAP_mean-centered fitting, and excludes
%   H+1 post-operative volume/function targets from pre-surgery fitting.
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

%% Test 3: H+1 post-operative echo volumes are excluded from pre-surgery
volume_values = [clinical.pre_surgery.LVEDV_mL, clinical.pre_surgery.LVESV_mL, ...
    clinical.pre_surgery.RVEDV_mL, clinical.pre_surgery.RVESV_mL, clinical.pre_surgery.EF];
if all(isnan(volume_values)) && isequal(clinical.pre_surgery.override_IC, false)
    fprintf('  [PASS] Pre-surgery Reyna excludes H+1 post-operative volume/function targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Pre-surgery Reyna still contains volume/function targets.\n');
    n_fail = n_fail + 1;
end

%% Test 4: primary metrics target hemodynamics only
expected_primary = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
if isequal(primary_metrics(:)', expected_primary)
    fprintf('  [PASS] Primary metrics include only direct hemodynamic targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Primary metrics do not match the Reyna hemodynamic-only profile.\n');
    n_fail = n_fail + 1;
end

%% Test 5: sparse cath profile keeps CO in the fit surface
if strcmp(profile.mode, 'sparse_cath') && ...
        ismember('CO_Lmin', profile.allowedMetricFields) && ...
        ismember('Q_shunt_Lmin', profile.allowedMetricFields) && ...
        ~ismember('LVEDV', profile.allowedMetricFields) && ...
        ~ismember('LVEF', profile.allowedMetricFields)
    fprintf('  [PASS] Sparse cath profile fits CO and derived shunt flow without volume/function targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Sparse cath profile target surface is not hemodynamic-only.\n');
    n_fail = n_fail + 1;
end

%% Test 6: sparse cath profile excludes unsupported chamber knobs
chamber_parameters = {'E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB','V0.LV','V0.RV'};
if ~any(ismember(chamber_parameters, profile.allowedFreeParameters)) && ...
        all(ismember({'group.R_sys_scale','R.SVEN','C.SAR', ...
        'group.R_pul_scale','C.PAR','vsd.Cd'}, profile.allowedFreeParameters))
    fprintf('  [PASS] Sparse cath profile does not tune chamber mechanics without pre-op chamber evidence.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Sparse cath profile allows unsupported chamber parameters.\n');
    n_fail = n_fail + 1;
end

%% Test 7: derived shunt flow is visible but not selected over direct anchors
shunt_idx = find(strcmp(tiers.table.Metric, 'Q_shunt_Lmin'), 1);
if isfinite(clinical.pre_surgery.Q_shunt_Lmin) && ...
        ~isempty(shunt_idx) && ...
        strcmp(tiers.table.Tier{shunt_idx}, 'soft') && ...
        ~ismember('Q_shunt_Lmin', primary_metrics)
    fprintf('  [PASS] Derived Q_shunt target is a soft guard behind direct hemodynamic anchors.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Derived Q_shunt target governance is not configured as expected.\n');
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
