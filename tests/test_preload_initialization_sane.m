%% test_preload_initialization_sane.m
% Checks that preload initialization is pressure-aware and blood-volume consistent.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

params_ref = default_parameters();
patient = struct('age_years', 0.5, 'age_days', 183, 'weight_kg', 7.0, ...
    'height_cm', 65, 'sex', 'M', 'maturation_mode', 'normal');

params = apply_scaling(params_ref, patient);
ic = params.ic.V;
idx = params.idx;

P_SAR = (ic(idx.V_SAR) - params.V0.SAR) / params.C.SAR;
P_SC = (ic(idx.V_SC) - params.V0.SC) / params.C.SC;
P_SVEN = (ic(idx.V_SVEN) - params.V0.SVEN) / params.C.SVEN;
P_PAR = (ic(idx.V_PAR) - params.V0.PAR) / params.C.PAR;
P_PVEN = (ic(idx.V_PVEN) - params.V0.PVEN) / params.C.PVEN;
P_RA = params.E.RA.EB * (ic(idx.V_RA) - params.V0.RA);
P_LA = params.E.LA.EB * (ic(idx.V_LA) - params.V0.LA);
P_LV_ED = params.E.LV.EB * (ic(idx.V_LV) - params.V0.LV);
P_RV_ED = params.E.RV.EB * (ic(idx.V_RV) - params.V0.RV);

BV_target = patient.weight_kg * blood_volume_per_kg(patient.age_years);
BV_total = sum(ic([idx.V_RA idx.V_RV idx.V_LA idx.V_LV idx.V_SAR idx.V_SC idx.V_SVEN idx.V_PAR idx.V_PVEN]));

assert(P_SVEN >= 0.5 && P_SVEN <= 8.0, 'Systemic venous pressure out of range.');
assert(P_PVEN >= 2.0 && P_PVEN <= 12.0, 'Pulmonary venous pressure out of range.');
assert(P_SAR >= 40.0 && P_SAR <= 90.0, 'Systemic arterial pressure out of range.');
assert(P_PAR >= 8.0 && P_PAR <= 25.0, 'Pulmonary arterial pressure out of range.');
assert(P_RA >= 0.5 && P_RA <= 8.0, 'Right atrial preload out of range.');
assert(P_LA >= 1.0 && P_LA <= 12.0, 'Left atrial preload out of range.');
assert(P_LV_ED >= 3.0 && P_LV_ED <= 12.0, 'LV end-diastolic pressure out of range.');
assert(P_RV_ED >= 1.0 && P_RV_ED <= 8.0, 'RV end-diastolic pressure out of range.');
assert(abs(BV_total - BV_target) <= 0.02 * BV_target, 'Total blood volume mismatch > 2%%.');

fprintf('[PASS] Preload initialization is pressure-aware and volume-consistent.\n');

function BV_per_kg = blood_volume_per_kg(age_years)
if age_years < 1
    BV_per_kg = 82;
else
    BV_per_kg = 70;
end
end
