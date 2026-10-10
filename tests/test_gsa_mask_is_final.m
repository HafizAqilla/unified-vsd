%% test_gsa_mask_is_final.m
% Parameters screened out by the global sensitivity analysis must stay fixed:
% no case-profile list (e.g. the Reyna systemic-polish knobs) may re-activate
% them. Previously the polish list was OR-ed into the mask, so a 7-parameter
% GSA selection silently became a 12-parameter calibration.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - GSA Mask Is Final Test\n');
fprintf('==========================================\n\n');

scenario = 'pre_surgery';
clinical = patient_reyna();
profile = build_case_calibration_profile(clinical, scenario);
assert(isfield(profile, 'systemicPolishNames') && ...
    any(ismember({'E.LV.EA','V0.LV'}, profile.systemicPolishNames)), ...
    'Precondition: the Reyna polish list includes chamber knobs the GSA screens out.');

patient = struct('age_years', clinical.common.age_years, ...
    'age_days', clinical.common.age_years * 365.25, ...
    'weight_kg', clinical.common.weight_kg, 'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, 'BSA', clinical.common.BSA, ...
    'maturation_mode', 'normal', 'scaling_mode', 'zhang');
params_ref = default_parameters();
params_scaled = apply_scaling(params_ref, patient);
params0 = params_from_clinical(params_scaled, clinical, scenario, params_scaled, profile);
registry = struct('params_adult', params_ref, 'params_scaled', params_scaled);

% Resolve the candidate order, then mask exactly the 7 GSA-selected parameters.
calib_all = calibration_param_sets(scenario, params0, [], {}, profile, registry);
gsa_selected = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.RV.EB','vsd.Cd'};
assert(all(ismember(gsa_selected, calib_all.names_all)), ...
    'Precondition: every GSA-selected name must be a calibration candidate.');
mask = ismember(calib_all.names_all(:), gsa_selected);

calib = calibration_param_sets(scenario, params0, mask, {}, profile, registry);

assert(isequal(sort(calib.names(:)'), sort(gsa_selected)), sprintf( ...
    'Active set must equal the GSA selection; got: %s', strjoin(calib.names, ', ')));
screened_out = setdiff(calib_all.names_all(:)', gsa_selected);
assert(~any(ismember(screened_out, calib.names)), ...
    'A GSA-screened parameter was re-activated.');
assert(isequal(calib.mask(:), mask(:)), 'calib.mask must equal the GSA mask.');

fprintf('  PASS: %d GSA-selected parameters stay the active set; %d screened-out stay fixed.\n', ...
    numel(gsa_selected), numel(screened_out));
