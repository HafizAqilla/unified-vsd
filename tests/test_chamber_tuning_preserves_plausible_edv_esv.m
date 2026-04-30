%% test_chamber_tuning_preserves_plausible_edv_esv.m
% Ensures chamber tuning preserves plausible EDV/ESV with echo-informed overrides.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

clinical = patient_template();
clinical.common.age_years = 0.5;
clinical.common.weight_kg = 7.0;
clinical.common.height_cm = 65.0;
clinical.common.sex = 'M';
clinical.common.HR = 140;

pre = clinical.pre_surgery;
pre.LVEDV_mL = 40;
pre.LVESV_mL = 15;
pre.RVEDV_mL = 45;
pre.RVESV_mL = 18;
pre.LAP_mean_mmHg = 6;
pre.RAP_mean_mmHg = 4;
pre.SAP_sys_mmHg = 90;
pre.SAP_dia_mmHg = 55;
pre.PAP_sys_mmHg = 30;
pre.PAP_dia_mmHg = 12;
pre.override_IC = true;
clinical.pre_surgery = pre;

params_ref = default_parameters();
patient = struct('age_years', clinical.common.age_years, ...
    'age_days', clinical.common.age_years * 365.25, ...
    'weight_kg', clinical.common.weight_kg, ...
    'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, ...
    'maturation_mode', 'normal');
params0 = apply_scaling(params_ref, patient);
params1 = params_from_clinical(params0, clinical, 'pre_surgery');

assert(params1.V0.LV > 0 && params1.V0.LV < pre.LVESV_mL, 'LV V0 outside plausible range.');
assert(params1.V0.RV > 0 && params1.V0.RV < pre.RVESV_mL, 'RV V0 outside plausible range.');
assert(params1.E.LV.EB > 0 && params1.E.RV.EB > 0, 'Passive elastance must be positive.');

P_lv_ed = pre.LAP_mean_mmHg;
P_lv_es = pre.SAP_sys_mmHg;
P_rv_ed = pre.RAP_mean_mmHg;
P_rv_es = pre.PAP_sys_mmHg;

Emax_LV = P_lv_es / max(pre.LVESV_mL - params1.V0.LV, 0.1);
Emax_RV = P_rv_es / max(pre.RVESV_mL - params1.V0.RV, 0.1);

ratio_LV = params1.E.LV.EB / max(Emax_LV, 1e-6);
ratio_RV = params1.E.RV.EB / max(Emax_RV, 1e-6);
assert(ratio_LV >= 0.01 && ratio_LV <= 0.08, 'LV EB/Emax ratio out of expected range.');
assert(ratio_RV >= 0.01 && ratio_RV <= 0.10, 'RV EB/Emax ratio out of expected range.');

LVEDV_pred = params1.V0.LV + P_lv_ed / params1.E.LV.EB;
RVEDV_pred = params1.V0.RV + P_rv_ed / params1.E.RV.EB;
assert(abs(LVEDV_pred - pre.LVEDV_mL) <= 0.10 * pre.LVEDV_mL, 'LVEDV prediction drift too large.');
assert(abs(RVEDV_pred - pre.RVEDV_mL) <= 0.10 * pre.RVEDV_mL, 'RVEDV prediction drift too large.');

fprintf('[PASS] Chamber tuning preserves plausible EDV/ESV relationships.\n');
