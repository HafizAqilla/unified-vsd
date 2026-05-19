function run_reyna_uncertainty_ensemble(n_samples, rng_seed)
% RUN_REYNA_UNCERTAINTY_ENSEMBLE
% -----------------------------------------------------------------------
% Propagates baseline-physiology uncertainty around the Reyna pre-surgery
% case without invoking calibration. The goal is to quantify how much the
% predicted metrics move when blood volume and vascular maturation priors
% are perturbed inside defensible ranges.
%
% INPUTS:
%   n_samples   - number of ensemble members                           [-]
%   rng_seed    - random seed for reproducibility                      [-]
%
% OUTPUTS:
%   results/tables/reyna_uncertainty_ensemble_*.csv
%   results/tables/reyna_uncertainty_summary_*.txt
%
% UNCERTAINTY KNOBS:
%   - BV_total_mL multiplier:   0.85 to 1.15
%   - systemic vascular scale: 0.90 to 1.10
%   - pulmonary vascular scale:0.85 to 1.15
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-08
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1 || isempty(n_samples)
    n_samples = 24;
end
if nargin < 2 || isempty(rng_seed)
    rng_seed = 20260508;
end

rng(rng_seed, 'twister');

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));
results_dir = fullfile(root, 'results', 'tables');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

clinical = patient_reyna();
scenario = 'pre_surgery';
case_profile = build_case_calibration_profile(clinical, scenario);

patient.age_years = clinical.common.age_years;
patient.age_days = clinical.common.age_years * 365.25;
patient.weight_kg = clinical.common.weight_kg;
patient.height_cm = clinical.common.height_cm;
patient.sex = clinical.common.sex;
patient.BSA = clinical.common.BSA;
patient.scaling_mode = 'lundquist_bsa';
patient.maturation_mode = 'normal';

params_ref = default_parameters();
params_scaled_nominal = apply_scaling(params_ref, patient);
params_seed_nominal = params_from_clinical(params_scaled_nominal, clinical, scenario, params_scaled_nominal, case_profile);
sim_nominal = integrate_system(params_seed_nominal);
metrics_nominal = compute_clinical_indices(sim_nominal, params_seed_nominal);
rmse_nominal = compute_case_rmse(metrics_nominal, clinical, scenario);

sample_ids = (1:n_samples)';
bv_scale = nan(n_samples, 1);
svr_scale = nan(n_samples, 1);
pvr_scale = nan(n_samples, 1);
rmse_col = nan(n_samples, 1);
qpqs_col = nan(n_samples, 1);
pap_min_col = nan(n_samples, 1);
pap_mean_col = nan(n_samples, 1);
svr_col = nan(n_samples, 1);
sap_mean_col = nan(n_samples, 1);
co_col = nan(n_samples, 1);
lvedv_col = nan(n_samples, 1);
lvef_col = nan(n_samples, 1);

for idx = 1:n_samples
    bv_scale(idx) = 0.85 + 0.30 * rand();
    svr_scale(idx) = 0.90 + 0.20 * rand();
    pvr_scale(idx) = 0.85 + 0.30 * rand();

    clinical_i = clinical;
    base_bv_mL = resolve_nominal_blood_volume_mL(clinical.common.age_years, clinical.common.weight_kg);
    clinical_i.pre_surgery.BV_total_mL = base_bv_mL * bv_scale(idx);

    params_scaled = apply_scaling(params_ref, patient);
    params_seeded = params_from_clinical(params_scaled, clinical_i, scenario, params_scaled, case_profile);
    params_perturbed = apply_uncertainty_scales(params_seeded, svr_scale(idx), pvr_scale(idx));

    sim_i = integrate_system(params_perturbed);
    metrics_i = compute_clinical_indices(sim_i, params_perturbed);

    rmse_col(idx) = compute_case_rmse(metrics_i, clinical_i, scenario);
    qpqs_col(idx) = metrics_i.QpQs;
    pap_min_col(idx) = metrics_i.PAP_min;
    pap_mean_col(idx) = metrics_i.PAP_mean;
    svr_col(idx) = metrics_i.SVR;
    sap_mean_col(idx) = metrics_i.SAP_mean;
    co_col(idx) = metrics_i.CO_Lmin;
    lvedv_col(idx) = metrics_i.LVEDV;
    lvef_col(idx) = metrics_i.LVEF;
end

ensemble_tbl = table(sample_ids, bv_scale, svr_scale, pvr_scale, rmse_col, ...
    qpqs_col, pap_min_col, pap_mean_col, svr_col, sap_mean_col, ...
    co_col, lvedv_col, lvef_col, ...
    'VariableNames', {'SampleId','BVScale','SVRScale','PVRScale','RMSE', ...
    'QpQs','PAP_min','PAP_mean','SVR','SAP_mean','CO_Lmin','LVEDV','LVEF'});

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_path = fullfile(results_dir, sprintf('reyna_uncertainty_ensemble_%s.csv', timestamp));
mat_path = fullfile(results_dir, sprintf('reyna_uncertainty_ensemble_%s.mat', timestamp));
txt_path = fullfile(results_dir, sprintf('reyna_uncertainty_summary_%s.txt', timestamp));
writetable(ensemble_tbl, csv_path);

summary_stats = struct();
summary_stats.rmse_nominal = rmse_nominal;
summary_stats.rmse_mean = mean(rmse_col, 'omitnan');
summary_stats.rmse_min = min(rmse_col);
summary_stats.rmse_max = max(rmse_col);
summary_stats.qpqs_range = [min(qpqs_col) max(qpqs_col)];
summary_stats.pap_min_range = [min(pap_min_col) max(pap_min_col)];
summary_stats.pap_mean_range = [min(pap_mean_col) max(pap_mean_col)];
summary_stats.svr_range = [min(svr_col) max(svr_col)];
summary_stats.sap_mean_range = [min(sap_mean_col) max(sap_mean_col)];
summary_stats.co_range = [min(co_col) max(co_col)];
summary_stats.lvedv_range = [min(lvedv_col) max(lvedv_col)];
summary_stats.lvef_range = [min(lvef_col) max(lvef_col)];
summary_stats.age_validity = build_age_validity_annotation(clinical.common.age_years, patient.scaling_mode, patient.maturation_mode);
save(mat_path, 'ensemble_tbl', 'summary_stats', 'clinical');

fid = fopen(txt_path, 'w');
if fid < 0
    error('run_reyna_uncertainty_ensemble:openFailed', ...
        'Unable to write uncertainty summary.');
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Reyna Uncertainty Ensemble Summary\n');
fprintf(fid, '=================================\n');
fprintf(fid, 'Samples: %d\n', n_samples);
fprintf(fid, 'Seed: %d\n', rng_seed);
fprintf(fid, 'Nominal RMSE: %.4f\n', rmse_nominal);
fprintf(fid, 'Ensemble RMSE mean/min/max: %.4f / %.4f / %.4f\n', ...
    summary_stats.rmse_mean, summary_stats.rmse_min, summary_stats.rmse_max);
fprintf(fid, 'QpQs range: %.4f to %.4f\n', summary_stats.qpqs_range(1), summary_stats.qpqs_range(2));
fprintf(fid, 'PAP_min range: %.4f to %.4f mmHg\n', ...
    summary_stats.pap_min_range(1), summary_stats.pap_min_range(2));
fprintf(fid, 'PAP_mean range: %.4f to %.4f mmHg\n', ...
    summary_stats.pap_mean_range(1), summary_stats.pap_mean_range(2));
fprintf(fid, 'SVR range: %.4f to %.4f WU\n', summary_stats.svr_range(1), summary_stats.svr_range(2));
fprintf(fid, 'SAP_mean range: %.4f to %.4f mmHg\n', ...
    summary_stats.sap_mean_range(1), summary_stats.sap_mean_range(2));
fprintf(fid, 'CO range: %.4f to %.4f L/min\n', summary_stats.co_range(1), summary_stats.co_range(2));
fprintf(fid, 'LVEDV range: %.4f to %.4f mL\n', summary_stats.lvedv_range(1), summary_stats.lvedv_range(2));
fprintf(fid, 'LVEF range: %.4f to %.4f\n', summary_stats.lvef_range(1), summary_stats.lvef_range(2));
fprintf(fid, 'Age validity regime: %s\n', summary_stats.age_validity.regime);
fprintf(fid, 'Age validity summary: %s\n', summary_stats.age_validity.summary);

fprintf('[run_reyna_uncertainty_ensemble] Ensemble exported:\n  %s\n  %s\n  %s\n', ...
    csv_path, mat_path, txt_path);
disp(ensemble_tbl(1:min(10, height(ensemble_tbl)), :));
end

function params = apply_uncertainty_scales(params, svr_scale, pvr_scale)
% APPLY_UNCERTAINTY_SCALES - perturb vascular priors around the seeded case.

params.R.SAR = params.R.SAR * svr_scale;
params.R.SC = params.R.SC * svr_scale;
params.R.SVEN = params.R.SVEN * svr_scale;
params.C.SAR = params.C.SAR / max(svr_scale, 1e-6);
params.C.SC = params.C.SC / max(svr_scale, 1e-6);
params.C.SVEN = params.C.SVEN / max(sqrt(svr_scale), 1e-6);

params.R.PAR = params.R.PAR * pvr_scale;
params.R.PCOX = params.R.PCOX * pvr_scale;
params.R.PCNO = params.R.PCNO * pvr_scale;
params.R.PVEN = params.R.PVEN * pvr_scale;
params.C.PAR = params.C.PAR / max(sqrt(pvr_scale), 1e-6);
params.C.PCOX = params.C.PCOX / max(sqrt(pvr_scale), 1e-6);
params.C.PCNO = params.C.PCNO / max(sqrt(pvr_scale), 1e-6);
params.C.PVEN = params.C.PVEN / max(sqrt(pvr_scale), 1e-6);
end

function rmse = compute_case_rmse(metrics, clinical, scenario)
targets = get_calibration_targets(scenario, clinical);
clinical_values = nan(numel(targets), 1);
model_values = nan(numel(targets), 1);
for idx = 1:numel(targets)
    clinical_values(idx) = targets(idx).ClinicalValue;
    if isfield(metrics, targets(idx).Metric)
        model_values(idx) = metrics.(targets(idx).Metric);
    end
end
valid = isfinite(clinical_values) & isfinite(model_values);
if ~any(valid)
    rmse = NaN;
    return;
end
rel_err = (model_values(valid) - clinical_values(valid)) ./ max(abs(clinical_values(valid)), 1e-9);
rmse = sqrt(mean(rel_err .^ 2));
end

function bv_mL = resolve_nominal_blood_volume_mL(age_years, weight_kg)
if age_years < 1
    bv_per_kg = 82;
else
    bv_per_kg = 70;
end
bv_mL = bv_per_kg * weight_kg;
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
