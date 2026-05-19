function scan_tbl = scan_reyna_targeted_sensitivity(start_package_path)
% SCAN_REYNA_TARGETED_SENSITIVITY
% -----------------------------------------------------------------------
% Performs a targeted one-at-a-time sensitivity scan around a saved Reyna
% best candidate for residuals that dominate the 10 percent gate.
%
% INPUTS:
%   start_package_path - optional params_best_candidate MAT file          [-]
%
% OUTPUTS:
%   scan_tbl           - local sensitivity table                          [-]
%
% ASSUMPTIONS:
%   - Perturbations are diagnostic only; they do not define accepted bounds.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-15
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

if nargin < 1 || isempty(start_package_path)
    start_package_path = fullfile(root, 'results', 'runs', ...
        '20260514_155253_reyna_pre_surgery', 'mat', ...
        'params_best_candidate_pre_surgery.mat');
end

scenario = 'pre_surgery';
clinical = patient_reyna();
case_profile = build_case_calibration_profile(clinical, scenario);
[params_scaled, params_seeded, registry_context] = ...
    build_reyna_reference_context(clinical, scenario, case_profile);
[primary_metrics, ~] = select_primary_metrics(clinical, [], scenario, case_profile);
calib_seed = calibration_param_sets( ...
    scenario, params_seeded, [], primary_metrics, case_profile, registry_context);

loaded = load(start_package_path, 'best_candidate');
params_start = loaded.best_candidate.params;
base_metrics = compute_clinical_indices(integrate_system(params_start), params_start);
targets = get_calibration_targets(scenario, clinical);

names = {'E.LA.EA','E.RA.EA','E.LV.EA','E.LV.EB','V0.LV', ...
    'group.R_sys_scale','R.SVEN','C.SAR','group.R_pul_scale','vsd.Cd'};
names = names(ismember(names, calib_seed.names_all));
metric_names = {'LAP_mean','LVEF','LVESV','LVEDV','CO_Lmin', ...
    'SAP_mean','SAP_max','QpQs','PAP_mean','RAP_mean','SVR'};
factors = [0.70, 0.85, 1.15, 1.30];

rows = {};
for name_idx = 1:numel(names)
    name = names{name_idx};
    registry_idx = find(strcmp(calib_seed.parameterRegistry.name, name), 1, 'first');
    base_value = get_calibration_param_value(params_start, params_scaled, name, case_profile);
    lb = calib_seed.parameterRegistry.lb(registry_idx);
    ub = calib_seed.parameterRegistry.ub(registry_idx);
    for factor_idx = 1:numel(factors)
        factor = factors(factor_idx);
        trial_value = min(max(base_value * factor, lb), ub);
        if abs(trial_value - base_value) < 1e-12
            continue;
        end
        params_trial = set_calibration_param_value( ...
            params_start, params_scaled, name, trial_value, case_profile);
        try
            trial_metrics = compute_clinical_indices(integrate_system(params_trial), params_trial);
            row = {name, factor, base_value, trial_value};
            for metric_idx = 1:numel(metric_names)
                metric_name = metric_names{metric_idx};
                row{end + 1} = metric_error_delta_pct( ...
                    base_metrics, trial_metrics, targets, metric_name); %#ok<AGROW>
            end
            row{end + 1} = primary_rmse(trial_metrics, targets, case_profile); %#ok<AGROW>
            rows(end + 1, :) = row; %#ok<AGROW>
        catch ME
            fprintf('[scan_reyna_targeted_sensitivity] %s factor %.2f failed: %s\n', ...
                name, factor, ME.message);
        end
    end
end

var_names = [{'Parameter','Factor','BaseValue','TrialValue'}, ...
    strcat(metric_names, '_DeltaErrorPct'), {'PrimaryRMSE'}];
scan_tbl = cell2table(rows, 'VariableNames', var_names);

out_dir = fullfile(root, 'results', 'tables');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
out_csv = fullfile(out_dir, sprintf('reyna_targeted_sensitivity_%s.csv', ...
    datestr(now, 'yyyymmdd_HHMMSS')));
writetable(scan_tbl, out_csv);
fprintf('[scan_reyna_targeted_sensitivity] Exported: %s\n', out_csv);
disp(scan_tbl);
end

function delta_error_pct = metric_error_delta_pct(base_metrics, trial_metrics, targets, metric_name)
idx = find(strcmp({targets.Metric}, metric_name), 1, 'first');
if isempty(idx) || ~isfield(base_metrics, metric_name) || ~isfield(trial_metrics, metric_name)
    delta_error_pct = NaN;
    return;
end
clinical_value = targets(idx).ClinicalValue;
if ~isfinite(clinical_value)
    delta_error_pct = NaN;
    return;
end
base_error = 100 * (base_metrics.(metric_name) - clinical_value) / ...
    max(abs(clinical_value), 1e-9);
trial_error = 100 * (trial_metrics.(metric_name) - clinical_value) / ...
    max(abs(clinical_value), 1e-9);
delta_error_pct = trial_error - base_error;
end

function value = primary_rmse(metrics, targets, case_profile)
tier_tbl = case_profile.targetTiers.table;
metric_names = {targets.Metric};
clinical_values = [targets.ClinicalValue];
err = nan(numel(targets), 1);
included = false(numel(targets), 1);
for idx = 1:numel(targets)
    metric_name = metric_names{idx};
    hit = find(strcmp(tier_tbl.Metric, metric_name), 1, 'first');
    if ~isempty(hit)
        included(idx) = logical(tier_tbl.IncludedInPrimaryRMSE(hit));
    end
    if included(idx) && isfield(metrics, metric_name) && ...
            isfinite(metrics.(metric_name)) && isfinite(clinical_values(idx))
        err(idx) = (metrics.(metric_name) - clinical_values(idx)) / ...
            max(abs(clinical_values(idx)), 1e-9);
    end
end
mask = included & isfinite(err);
if any(mask)
    value = sqrt(mean(err(mask).^2));
else
    value = NaN;
end
end

function [params_scaled, params_seeded, registry_context] = ...
    build_reyna_reference_context(clinical, scenario, case_profile)
params_ref = default_parameters();
patient = struct();
patient.age_years = clinical.common.age_years;
patient.age_days = clinical.common.age_years * 365.25;
patient.weight_kg = clinical.common.weight_kg;
patient.height_cm = clinical.common.height_cm;
patient.sex = clinical.common.sex;
patient.maturation_mode = 'normal';
patient.scaling_mode = case_profile.preferredScalingMode;
if isfield(clinical.common, 'BSA') && isfinite(clinical.common.BSA)
    patient.BSA = clinical.common.BSA;
end
params_scaled = apply_scaling(params_ref, patient);
params_seeded = params_from_clinical(params_scaled, clinical, scenario, ...
    params_scaled, case_profile);
registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
end

function project_path = build_clean_project_path(root)
project_paths = strsplit(genpath(root), pathsep);
project_paths = project_paths(~cellfun('isempty', project_paths));
is_shadow = contains(project_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
    contains(project_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
    contains(project_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, project_paths);
project_path = strjoin(project_paths(~is_shadow & is_existing), pathsep);
end
