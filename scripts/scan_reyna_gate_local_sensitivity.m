function scan_tbl = scan_reyna_gate_local_sensitivity(start_package_path)
% SCAN_REYNA_GATE_LOCAL_SENSITIVITY
% -----------------------------------------------------------------------
% One-at-a-time local sensitivity scan around a saved Reyna best candidate.
% The scan is intentionally small and direct-ODE based; it is used to
% identify whether the remaining 15% gate errors have a local feasible
% direction before running another constrained polish.
%
% INPUTS:
%   start_package_path - optional params_best_candidate MAT file          [-]
%
% OUTPUTS:
%   scan_tbl           - local perturbation table                         [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-12
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

if nargin < 1 || isempty(start_package_path)
    start_package_path = find_latest_reyna_best_candidate(root);
end

clinical = patient_reyna();
scenario = 'pre_surgery';
case_profile = build_case_calibration_profile(clinical, scenario);
[params_scaled, params_seeded, registry_context] = build_reyna_reference_context( ...
    clinical, scenario, case_profile);
[primary_metrics, ~] = select_primary_metrics(clinical, [], scenario, case_profile);
calib_seed = calibration_param_sets( ...
    scenario, params_seeded, [], primary_metrics, case_profile, registry_context);

loaded = load(start_package_path, 'best_candidate');
params_start = loaded.best_candidate.params;
base_metrics = compute_clinical_indices(integrate_system(params_start), params_start);
targets = get_calibration_targets(scenario, clinical);

names = {'group.R_sys_scale','R.SVEN','C.SAR', ...
    'E.LV.EA','E.LV.EB','V0.LV', ...
    'E.RV.EA','E.RV.EB','V0.RV','vsd.Cd'};
names = names(ismember(names, calib_seed.names_all));
metric_names = {'CO_Lmin','SAP_max','RVEDV','LVEF','LVESV','SAP_mean','QpQs','PAP_mean'};

rows = {};
for name_idx = 1:numel(names)
    name = names{name_idx};
    registry_idx = find(strcmp(calib_seed.parameterRegistry.name, name), 1, 'first');
    base_value = get_calibration_param_value(params_start, params_scaled, name, case_profile);
    lb = calib_seed.parameterRegistry.lb(registry_idx);
    ub = calib_seed.parameterRegistry.ub(registry_idx);
    step = max(0.05 * abs(base_value), 0.02 * (ub - lb));

    for sign_value = [-1, 1]
        trial_value = min(max(base_value + sign_value * step, lb), ub);
        if abs(trial_value - base_value) < 1e-12
            continue;
        end
        params_trial = set_calibration_param_value( ...
            params_start, params_scaled, name, trial_value, case_profile);
        try
            trial_metrics = compute_clinical_indices(integrate_system(params_trial), params_trial);
            row = {name, sign_value, base_value, trial_value};
            for metric_idx = 1:numel(metric_names)
                metric_name = metric_names{metric_idx};
                row{end + 1} = metric_error_delta_pct( ...
                    base_metrics, trial_metrics, targets, metric_name); %#ok<AGROW>
            end
            rows(end + 1, :) = row; %#ok<AGROW>
        catch ME
            fprintf('[scan_reyna_gate_local_sensitivity] %s sign %+d failed: %s\n', ...
                name, sign_value, ME.message);
        end
    end
end

var_names = [{'Parameter','Direction','BaseValue','TrialValue'}, ...
    strcat(metric_names, '_DeltaErrorPct')];
scan_tbl = cell2table(rows, 'VariableNames', var_names);

out_dir = fullfile(root, 'results', 'tables');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
out_csv = fullfile(out_dir, sprintf('reyna_gate_local_sensitivity_%s.csv', ...
    datestr(now, 'yyyymmdd_HHMMSS')));
writetable(scan_tbl, out_csv);
fprintf('[scan_reyna_gate_local_sensitivity] Exported: %s\n', out_csv);
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
params_seeded = params_from_clinical(params_scaled, clinical, scenario, params_scaled, case_profile);
registry_context = struct('params_adult', params_ref, 'params_scaled', params_scaled);
end

function best_candidate_path = find_latest_reyna_best_candidate(root)
listing = dir(fullfile(root, 'results', 'runs', '*reyna_pre_surgery', ...
    'mat', 'params_best_candidate_pre_surgery.mat'));
if isempty(listing)
    error('scan_reyna_gate_local_sensitivity:noBestCandidate', ...
        'No Reyna best-candidate package found under results/runs.');
end
[~, order] = sort([listing.datenum], 'descend');
best_candidate_path = fullfile(listing(order(1)).folder, listing(order(1)).name);
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
