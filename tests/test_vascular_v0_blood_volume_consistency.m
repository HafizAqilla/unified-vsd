%% test_vascular_v0_blood_volume_consistency.m
% Ensures vascular V0 reconciliation aligns blood volume with nominal pressures.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

params_ref = default_parameters();
patient = struct('age_years', 0.5, 'age_days', 183, 'weight_kg', 7.0, ...
    'height_cm', 65, 'sex', 'M', 'maturation_mode', 'normal');

params_scaled = apply_scaling(params_ref, patient);

BV_target = patient.weight_kg * blood_volume_per_kg(patient.age_years);
P_nom_SAR = 65;
P_nom_SC = max(0.5 * P_nom_SAR, 15);
P_nom_SVEN = 2;
P_nom_PAR = 15;
P_nom_PVEN = 6;
P_nom_LV_ED = 8;
P_nom_RV_ED = 4;
P_nom_LA_ED = 6;
P_nom_RA_ED = 4;

V_RA = params_scaled.V0.RA + P_nom_RA_ED / max(params_scaled.E.RA.EB, 1e-6);
V_RV = params_scaled.V0.RV + P_nom_RV_ED / max(params_scaled.E.RV.EB, 1e-6);
V_LA = params_scaled.V0.LA + P_nom_LA_ED / max(params_scaled.E.LA.EB, 1e-6);
V_LV = params_scaled.V0.LV + P_nom_LV_ED / max(params_scaled.E.LV.EB, 1e-6);
V_SAR = params_scaled.V0.SAR + P_nom_SAR * params_scaled.C.SAR;
V_SC = params_scaled.V0.SC + P_nom_SC * params_scaled.C.SC;
V_PAR = params_scaled.V0.PAR + P_nom_PAR * params_scaled.C.PAR;
V_PVEN = params_scaled.V0.PVEN + P_nom_PVEN * params_scaled.C.PVEN;

V_other = V_RA + V_RV + V_LA + V_LV + V_SAR + V_SC + V_PAR + V_PVEN;
V0_required = BV_target - V_other - P_nom_SVEN * params_scaled.C.SVEN;

assert(isfield(params_scaled.scaling, 'V0_SVEN_reconciled') && params_scaled.scaling.V0_SVEN_reconciled, ...
    'V0.SVEN reconciliation flag missing.');
assert(abs(params_scaled.V0.SVEN - V0_required) <= 0.02 * BV_target, ...
    'V0.SVEN does not reconcile blood volume within tolerance.');

fprintf('[PASS] Vascular V0 reconciliation matches blood volume target.\n');

function BV_per_kg = blood_volume_per_kg(age_years)
if age_years < 1
    BV_per_kg = 82;
else
    BV_per_kg = 70;
end
end
