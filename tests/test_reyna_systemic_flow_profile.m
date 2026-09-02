%% test_reyna_systemic_flow_profile.m
% =========================================================================
% Regression test for Reyna hemodynamic-only pre-surgery calibration profile.
%
% PURPOSE:
%   Confirms that the active Reyna case resolves through the explicit
%   pre-surgery recipe instead of generic sparse-cath inference.
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

%% Test 1: explicit recipe is selected
if strcmp(profile.mode, 'reyna_recipe') && ...
        isfield(profile, 'recipe_id') && strcmp(profile.recipe_id, 'reyna_pre_surgery')
    fprintf('  [PASS] Reyna pre-surgery resolves to the explicit recipe.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Reyna pre-surgery did not resolve to the explicit recipe.\n');
    n_fail = n_fail + 1;
end

%% Test 2: recipe demographics are applied before scaling
[recipe, recipe_found] = load_calibration_recipe(clinical, 'pre_surgery');
clinical_recipe = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
% Measurement-day anthropometry from the study "reyna" procedure log
% (06/04/2026 07.53; full provenance in
% config/private/patient_provenance.local.m) — the session that produced
% the pressures this recipe fits. Superseded a 2026-05-11 revision
% (14.0 kg / 98.0 cm /
% Mosteller BSA 0.6173) recorded five weeks later, which scaled the model to
% a larger child than the one who was measured. This assertion exists because
% recipe.demographics is merged OVER clinical.common, so a mismatch between
% the recipe and patient_reyna silently wins here rather than erroring.
if recipe_found && clinical_recipe.common.weight_kg == 13.4 && ...
        clinical_recipe.common.height_cm == 95.0 && ...
        abs(clinical_recipe.common.BSA - 0.588) < 1e-12
    fprintf('  [PASS] Recipe demographics are deterministic.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Recipe demographics were not applied as expected.\n');
    n_fail = n_fail + 1;
end

%% Test 3: pulmonary calibration keeps waveform evidence but fits PAP_mean
tiers = profile.targetTiers;
pap_min_idx = find(strcmp(tiers.table.Metric, 'PAP_min'), 1);
pap_max_idx = find(strcmp(tiers.table.Metric, 'PAP_max'), 1);
if clinical_recipe.pre_surgery.PAP_sys_mmHg == 20 && ...
        clinical_recipe.pre_surgery.PAP_dia_mmHg == 10 && ...
        clinical_recipe.pre_surgery.PAP_mean_mmHg == 15 && ...
        strcmp(tiers.table.Tier{pap_min_idx}, 'validation_only') && ...
        strcmp(tiers.table.Tier{pap_max_idx}, 'validation_only')
    fprintf('  [PASS] PAP_sys/PAP_dia are retained as validation evidence while PAP_mean anchors fitting.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] PAP waveform evidence/tiering is not configured as expected.\n');
    n_fail = n_fail + 1;
end

%% Test 4: H+1 echo rows are retained only as consistency checks
volume_values = [clinical_recipe.pre_surgery.LVEDV_mL, ...
    clinical_recipe.pre_surgery.LVESV_mL, ...
    clinical_recipe.pre_surgery.RVEDV_mL, ...
    clinical_recipe.pre_surgery.RVESV_mL, clinical_recipe.pre_surgery.EF];
lvedv_idx = find(strcmp(tiers.table.Metric, 'LVEDV'), 1);
lvesv_idx = find(strcmp(tiers.table.Metric, 'LVESV'), 1);
rvedv_idx = find(strcmp(tiers.table.Metric, 'RVEDV'), 1);
lvef_idx = find(strcmp(tiers.table.Metric, 'LVEF'), 1);
if all(isfinite(volume_values)) && clinical_recipe.pre_surgery.override_IC && ...
        strcmp(tiers.table.Tier{lvedv_idx}, 'consistency_check_only') && ...
        strcmp(tiers.table.Tier{lvesv_idx}, 'consistency_check_only') && ...
        strcmp(tiers.table.Tier{rvedv_idx}, 'consistency_check_only') && ...
        strcmp(tiers.table.Tier{lvef_idx}, 'consistency_check_only') && ...
        ~tiers.table.IncludedInPrimaryRMSE(lvedv_idx)
    fprintf('  [PASS] H+1 volume/function rows are consistency-only evidence.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] H+1 volume/function governance is not recipe-stable.\n');
    n_fail = n_fail + 1;
end

%% Test 5: primary metrics target hemodynamics only
expected_primary = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
if isequal(primary_metrics(:)', expected_primary)
    fprintf('  [PASS] Primary metrics include only direct hemodynamic targets.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Primary metrics do not match the Reyna hemodynamic-only profile.\n');
    n_fail = n_fail + 1;
end

%% Test 6: active parameter list is recipe-pinned
expected_parameters = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB', ...
    'E.LA.EA','E.RA.EA','V0.LV','V0.RV','vsd.Cd'};
if isequal(profile.allowedFreeParameters(:)', expected_parameters)
    fprintf('  [PASS] Active parameter list is pinned by the recipe.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Active parameter list drifted from the recipe.\n');
    n_fail = n_fail + 1;
end

%% Test 7: derived validation and manifest recipe fields are protected
if isfield(profile, 'initialParameterValues') && ...
        isequal(profile.initialParameterValues.names(:)', expected_parameters) && ...
        abs(profile.initialParameterValues.values(1) - 0.457982333986) < 1e-12 && ...
        isfield(profile.initialParameterValues, 'scaling_modes') && ...
        isequal(profile.initialParameterValues.scaling_modes(:)', {'lundquist_bsa'}) && ...
        isfield(profile, 'acceptInitialSeedIfPass') && profile.acceptInitialSeedIfPass && ...
        all(ismember({'SVR','PVR'}, tiers.derived_validation)) && ...
        contains(fileread(fullfile(project_root, 'main_run.m')), 'CalibrationRecipeId')
    fprintf('  [PASS] Lundquist recipe seed, derived validation, and manifest reporting are protected.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Recipe seed, derived validation, or manifest reporting drifted.\n');
    n_fail = n_fail + 1;
end

%% Test 8: derived shunt flow is visible but not selected over direct anchors
shunt_idx = find(strcmp(tiers.table.Metric, 'Q_shunt_Lmin'), 1);
if isfinite(clinical_recipe.pre_surgery.Q_shunt_Lmin) && ...
        ~isempty(shunt_idx) && ...
        strcmp(tiers.table.Tier{shunt_idx}, 'soft') && ...
        ~tiers.table.IncludedInPrimaryRMSE(shunt_idx) && ...
        ~ismember('Q_shunt_Lmin', primary_metrics)
    fprintf('  [PASS] Derived Q_shunt target is a soft guard outside primary RMSE.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Derived Q_shunt target governance is not configured as expected.\n');
    n_fail = n_fail + 1;
end

%% Test 9: accepted Lundquist disease seed is not universal
seed_modes_ok = isfield(recipe.initial_parameter_values, 'scaling_modes') && ...
    isequal(recipe.initial_parameter_values.scaling_modes(:)', {'lundquist_bsa'}) && ...
    isfield(recipe.fixed_parameter_values, 'scaling_modes') && ...
    isequal(recipe.fixed_parameter_values.scaling_modes(:)', {'lundquist_bsa'}) && ...
    isfield(recipe.initial_conditions, 'scaling_modes') && ...
    isequal(recipe.initial_conditions.scaling_modes(:)', {'lundquist_bsa'}) && ...
    isfield(profile, 'acceptInitialSeedScalingModes') && ...
    isequal(profile.acceptInitialSeedScalingModes(:)', {'lundquist_bsa'});
if seed_modes_ok
    fprintf('  [PASS] Accepted disease seed is explicitly Lundquist-only.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Accepted disease seed is not protected from Zhang scaling.\n');
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
