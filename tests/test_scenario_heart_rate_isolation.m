%% test_scenario_heart_rate_isolation.m
% Pre- and post-closure states must each run at the heart rate recorded on
% their own protocol form, and the mean arterial pressure targets must follow
% the protocol's MAP definition applied to the catheter readings.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Scenario HR Isolation Test\n');
fprintf('==========================================\n\n');

clinical = patient_reyna();

% ---- 1. Reyna carries no shared HR; each scenario has its own --------------
assert(~isfield(clinical.common, 'HR'), ...
    'A shared common.HR would let one state inherit the other''s rate.');
assert(clinical_heart_rate(clinical, 'pre_surgery') == 119, ...
    'Pre-closure HR must be 119 bpm (pre protocol row 6).');
assert(clinical_heart_rate(clinical, 'post_surgery') == 108, ...
    'Post-closure HR must be 108 bpm (post protocol row 6).');

% ---- 2. Patients with a single recorded HR still resolve through common ----
single = struct('common', struct('HR', 100), ...
    'pre_surgery', struct(), 'post_surgery', struct());
assert(clinical_heart_rate(single, 'pre_surgery') == 100);
assert(clinical_heart_rate(single, 'post_surgery') == 100);
single.post_surgery.HR = 90;
assert(clinical_heart_rate(single, 'post_surgery') == 90, ...
    'A scenario HR must take precedence over common.HR.');
assert(clinical_heart_rate(single, 'pre_surgery') == 100, ...
    'Setting the post HR must not change the pre HR.');

% ---- 3. params_from_clinical applies the scenario HR and cycle timing -------
patient = struct('age_years', clinical.common.age_years, ...
    'age_days', clinical.common.age_years * 365.25, ...
    'weight_kg', clinical.common.weight_kg, 'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, 'BSA', clinical.common.BSA, ...
    'maturation_mode', 'normal', 'scaling_mode', 'zhang');
params0 = apply_scaling(default_parameters(), patient);
for sc = {'pre_surgery', 'post_surgery'}
    p = params_from_clinical(params0, clinical, sc{1});
    HR_expected = clinical_heart_rate(clinical, sc{1});
    assert(p.HR == HR_expected, sprintf('%s: params.HR not set from scenario.', sc{1}));
    T_HB = 60 / HR_expected;
    assert(abs(p.Tc_LV - p.Tc_LV_frac * T_HB) < 1e-12, ...
        sprintf('%s: ventricular timing not rebuilt for the scenario HR.', sc{1}));
    assert(abs(p.t_ac_RA - p.t_ac_RA_frac * T_HB) < 1e-12, ...
        sprintf('%s: atrial timing not rebuilt for the scenario HR.', sc{1}));
end

% ---- 4. MAP targets follow protocol MAP = (sys + 2*dia)/3 on catheter RFA --
pre = clinical.pre_surgery;
post = clinical.post_surgery;
assert(abs(pre.SAP_mean_mmHg - (pre.SAP_sys_mmHg + 2*pre.SAP_dia_mmHg)/3) < 1e-12, ...
    'Pre-closure MAP must be the protocol formula on RFA 100/57.');
assert(abs(pre.SAP_mean_mmHg - 71.3333333333) < 1e-6);
assert(abs(post.MAP_mmHg - (post.SAP_sys_mmHg + 2*post.SAP_dia_mmHg)/3) < 1e-12, ...
    'Post-closure MAP must be the protocol formula on RFA 89/68.');
assert(post.MAP_mmHg == 75);

fprintf('  PASS: pre 119 bpm / post 108 bpm isolated; MAP 71.33 / 75 from protocol formula.\n');
