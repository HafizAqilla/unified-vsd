function result = run_patient_cohort_fast_pipeline(options)
% RUN_PATIENT_COHORT_FAST_PIPELINE
% -----------------------------------------------------------------------
% Runs the available VSD cohort through scaling, clinical consistency
% audit, baseline simulation, fast staged calibration, and validation
% statistics.
%
% INPUTS:
%   options - optional struct controlling calibration and output          [-]
%
% OUTPUTS:
%   result  - struct containing cohort tables and output paths            [-]
%
% ASSUMPTIONS:
%   - Fast calibration is used for cohort triage. Final publication runs
%     should rerun selected cases with the full Reyna-grade pipeline.
%   - Cases without anthropometry are audited but not simulated.
%
% REFERENCES:
%   [1] config/patient_cohort_cases.m
%   [2] src/validation/audit_clinical_consistency.m
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-15
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1 || isempty(options)
    options = struct();
end
options = default_options(options);

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));
setenv('UNIFIED_VSD_FMINCON_PARALLEL', '0');

cases = patient_cohort_cases();
cases = filter_cases(cases, options.case_labels);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
run_dir = fullfile(root, 'results', 'runs', ...
    sprintf('%s_patient_cohort_fast_pipeline', timestamp));
tables_dir = fullfile(run_dir, 'tables');
mat_dir = fullfile(run_dir, 'mat');
if ~exist(tables_dir, 'dir'), mkdir(tables_dir); end
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end

summary_rows = {};
error_tables = {};
consistency_rows = {};
case_results = cell(numel(cases), 1);

for case_idx = 1:numel(cases)
    case_item = cases(case_idx);
    clinical = case_item.clinical;
    scenario = case_item.scenario;
    label = case_item.label;
    fprintf('\n[cohort] Case %d/%d: %s\n', case_idx, numel(cases), label);

    audit = audit_clinical_consistency(clinical, scenario);
    consistency_rows(end + 1, :) = consistency_row(label, audit, clinical, scenario); %#ok<AGROW>

    case_result = struct('label', label, 'scenario', scenario, ...
        'clinical', clinical, 'audit', audit, 'status', 'not_started');

    if ~case_item.can_simulate
        case_result.status = 'skipped';
        case_result.skip_reason = case_item.skip_reason;
        summary_rows(end + 1, :) = skipped_summary_row(label, scenario, case_item.skip_reason, audit); %#ok<AGROW>
        case_results{case_idx} = case_result;
        fprintf('[cohort] Skipped %s: %s\n', label, case_item.skip_reason);
        continue;
    end

    try
        case_profile = build_case_calibration_profile(clinical, scenario);
        [params_scaled, params0, registry_context] = ...
            build_case_reference_context(clinical, scenario, case_profile);
        sim_base = integrate_system(params0);
        metrics_base = compute_clinical_indices(sim_base, params0);
        validity_base = evaluate_simulation_validity(sim_base, params0, ...
            metrics_base, scenario, clinical);

        [primary_metrics, primary_selection_table] = ...
            select_primary_metrics(clinical, [], scenario, case_profile); %#ok<ASGLU>
        [params_cal, calib_out] = run_calibration(params0, clinical, scenario, ...
            [], options.fast_calibration, [], primary_metrics, ...
            case_profile, registry_context);
        sim_cal = integrate_system(params_cal);
        metrics_cal = compute_clinical_indices(sim_cal, params_cal);
        validity_cal = evaluate_simulation_validity(sim_cal, params_cal, ...
            metrics_cal, scenario, clinical);

        targets = build_target_table(clinical, scenario, case_profile);
        error_table = build_error_table(label, targets, metrics_base, metrics_cal);
        stats = summarize_errors(error_table);
        summary_rows(end + 1, :) = summary_row(label, scenario, 'calibrated', ...
            '', case_profile, audit, validity_base, validity_cal, stats); %#ok<AGROW>
        error_tables{end + 1} = error_table; %#ok<AGROW>

        case_result.status = 'calibrated';
        case_result.case_profile = case_profile;
        case_result.params_scaled = params_scaled;
        case_result.params0 = params0;
        case_result.params_cal = params_cal;
        case_result.metrics_base = metrics_base;
        case_result.metrics_cal = metrics_cal;
        case_result.validity_base = validity_base;
        case_result.validity_cal = validity_cal;
        case_result.calib_out = calib_out;
        case_result.error_table = error_table;
        case_result.stats = stats;
        fprintf('[cohort] %s primary RMSE %.4f, calibration pass10 %d/%d\n', ...
            label, stats.primary_rmse_cal, stats.calibration_pass10, ...
            stats.calibration_n);
    catch ME
        case_result.status = 'failed';
        case_result.error_message = ME.message;
        summary_rows(end + 1, :) = failed_summary_row(label, scenario, ME.message, audit); %#ok<AGROW>
        fprintf('[cohort] FAILED %s: %s\n', label, ME.message);
    end
    case_results{case_idx} = case_result;
end

summary_table = cell2table(summary_rows, 'VariableNames', summary_variable_names());
if isempty(error_tables)
    error_table_all = table();
else
    error_table_all = vertcat(error_tables{:});
end
consistency_table = cell2table(consistency_rows, 'VariableNames', consistency_variable_names());

summary_csv = fullfile(tables_dir, 'patient_cohort_summary.csv');
errors_csv = fullfile(tables_dir, 'patient_cohort_errors.csv');
consistency_csv = fullfile(tables_dir, 'patient_cohort_consistency.csv');
markdown_path = fullfile(run_dir, 'patient_cohort_report.md');
mat_path = fullfile(mat_dir, 'patient_cohort_fast_pipeline_result.mat');

writetable(summary_table, summary_csv);
writetable(error_table_all, errors_csv);
writetable(consistency_table, consistency_csv);
result = struct('run_dir', run_dir, 'summary_table', summary_table, ...
    'error_table', error_table_all, 'consistency_table', consistency_table, ...
    'case_results', case_results, 'options', options, ...
    'paths', struct('summary_csv', summary_csv, 'errors_csv', errors_csv, ...
    'consistency_csv', consistency_csv, 'markdown', markdown_path, 'mat', mat_path));
write_markdown_report(markdown_path, result);
save(mat_path, 'result');

fprintf('\n[cohort] Exported:\n  %s\n  %s\n  %s\n', summary_csv, errors_csv, markdown_path);
end

function options = default_options(options)
defaults = struct('fast_calibration', true);
defaults.case_labels = {};
fields = fieldnames(defaults);
for idx = 1:numel(fields)
    field_name = fields{idx};
    if ~isfield(options, field_name) || isempty(options.(field_name))
        options.(field_name) = defaults.(field_name);
    end
end
end

function cases = filter_cases(cases, labels)
if isempty(labels)
    return;
end
if ischar(labels) || isstring(labels)
    labels = cellstr(labels);
end
keep = ismember({cases.label}, labels);
cases = cases(keep);
end

function [params_scaled, params0, registry_context] = ...
    build_case_reference_context(clinical, scenario, case_profile)
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
params0 = params_from_clinical(params_scaled, clinical, scenario, ...
    params_scaled, case_profile);
registry_context = struct('params_adult', params_ref, ...
    'params_scaled', params_scaled);
end

function targets = build_target_table(clinical, scenario, case_profile)
target_struct = get_calibration_targets(scenario, clinical);
tier_tbl = case_profile.targetTiers.table;
n_targets = numel(target_struct);
metric_col = {target_struct.Metric}';
clinical_col = [target_struct.ClinicalValue]';
unit_col = {target_struct.Unit}';
tier_col = repmat({'unavailable'}, n_targets, 1);
included_cal = false(n_targets, 1);
included_primary = false(n_targets, 1);
for idx = 1:n_targets
    hit = find(strcmp(tier_tbl.Metric, metric_col{idx}), 1, 'first');
    if ~isempty(hit)
        tier_col{idx} = tier_tbl.Tier{hit};
        included_cal(idx) = logical(tier_tbl.IncludedInCalibration(hit));
        included_primary(idx) = logical(tier_tbl.IncludedInPrimaryRMSE(hit));
    end
end
targets = table(metric_col, clinical_col, unit_col, tier_col, ...
    included_cal, included_primary, ...
    'VariableNames', {'Metric','Clinical','Unit','Tier', ...
    'IncludedInCalibration','IncludedInPrimaryRMSE'});
end

function tbl = build_error_table(label, targets, metrics_base, metrics_cal)
n_targets = height(targets);
base_col = nan(n_targets, 1);
cal_col = nan(n_targets, 1);
for idx = 1:n_targets
    metric_name = targets.Metric{idx};
    if isfield(metrics_base, metric_name)
        base_col(idx) = metrics_base.(metric_name);
    end
    if isfield(metrics_cal, metric_name)
        cal_col(idx) = metrics_cal.(metric_name);
    end
end
base_err = 100 * (base_col - targets.Clinical) ./ max(abs(targets.Clinical), 1e-9);
cal_err = 100 * (cal_col - targets.Clinical) ./ max(abs(targets.Clinical), 1e-9);
base_err(~isfinite(targets.Clinical)) = NaN;
cal_err(~isfinite(targets.Clinical)) = NaN;
case_col = repmat({label}, n_targets, 1);
tbl = table(case_col, targets.Metric, targets.Unit, targets.Tier, ...
    targets.IncludedInCalibration, targets.IncludedInPrimaryRMSE, ...
    targets.Clinical, base_col, base_err, cal_col, cal_err, ...
    abs(cal_err) <= 5, abs(cal_err) <= 10, ...
    'VariableNames', {'Case','Metric','Unit','Tier', ...
    'IncludedInCalibration','IncludedInPrimaryRMSE','Clinical', ...
    'BaselineModel','BaselineError_pct','CalibratedModel', ...
    'CalibratedError_pct','Pass_5pct','Pass_10pct'});
end

function stats = summarize_errors(tbl)
stats = struct();
all_mask = isfinite(tbl.Clinical) & isfinite(tbl.CalibratedModel);
cal_mask = all_mask & tbl.IncludedInCalibration;
primary_mask = all_mask & tbl.IncludedInPrimaryRMSE;
pressure_flow_names = {'RAP_mean','PAP_mean','SAP_min','SAP_max', ...
    'SAP_mean','QpQs','PVR','SVR','CO_Lmin'};
pf_mask = all_mask & ismember(tbl.Metric, pressure_flow_names);
stats.all_n = nnz(all_mask);
stats.all_pass5 = nnz(abs(tbl.CalibratedError_pct(all_mask)) <= 5);
stats.all_pass10 = nnz(abs(tbl.CalibratedError_pct(all_mask)) <= 10);
stats.calibration_n = nnz(cal_mask);
stats.calibration_pass5 = nnz(abs(tbl.CalibratedError_pct(cal_mask)) <= 5);
stats.calibration_pass10 = nnz(abs(tbl.CalibratedError_pct(cal_mask)) <= 10);
stats.primary_n = nnz(primary_mask);
stats.primary_pass5 = nnz(abs(tbl.CalibratedError_pct(primary_mask)) <= 5);
stats.primary_pass10 = nnz(abs(tbl.CalibratedError_pct(primary_mask)) <= 10);
stats.pressure_flow_n = nnz(pf_mask);
stats.pressure_flow_pass5 = nnz(abs(tbl.CalibratedError_pct(pf_mask)) <= 5);
stats.pressure_flow_pass10 = nnz(abs(tbl.CalibratedError_pct(pf_mask)) <= 10);
stats.primary_rmse_cal = rmse_pct(tbl.CalibratedError_pct(primary_mask));
stats.full_rmse_cal = rmse_pct(tbl.CalibratedError_pct(all_mask));
stats.calibration_mae_pct = mean(abs(tbl.CalibratedError_pct(cal_mask)), 'omitnan');
stats.primary_mae_pct = mean(abs(tbl.CalibratedError_pct(primary_mask)), 'omitnan');
stats.max_abs_error_pct = max(abs(tbl.CalibratedError_pct(all_mask)), [], 'omitnan');
end

function value = rmse_pct(error_pct)
mask = isfinite(error_pct);
if any(mask)
    value = sqrt(mean((error_pct(mask) / 100).^2));
else
    value = NaN;
end
end

function row = summary_row(label, scenario, status, message, profile, audit, validity_base, validity_cal, stats)
row = {label, scenario, status, message, profile.mode, audit.severity, ...
    audit.max_relative_difference, validity_base.is_valid, validity_cal.is_valid, ...
    stats.primary_rmse_cal, stats.full_rmse_cal, stats.primary_mae_pct, ...
    stats.calibration_mae_pct, stats.max_abs_error_pct, stats.all_n, ...
    stats.all_pass5, stats.all_pass10, stats.calibration_n, ...
    stats.calibration_pass5, stats.calibration_pass10, stats.primary_n, ...
    stats.primary_pass5, stats.primary_pass10, stats.pressure_flow_n, ...
    stats.pressure_flow_pass5, stats.pressure_flow_pass10};
end

function row = skipped_summary_row(label, scenario, reason, audit)
row = {label, scenario, 'skipped', reason, 'not_runnable', audit.severity, ...
    audit.max_relative_difference, false, false, NaN, NaN, NaN, NaN, ...
    NaN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
end

function row = failed_summary_row(label, scenario, message, audit)
row = {label, scenario, 'failed', message, 'failed', audit.severity, ...
    audit.max_relative_difference, false, false, NaN, NaN, NaN, NaN, ...
    NaN, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
end

function names = summary_variable_names()
names = {'Case','Scenario','Status','Message','CaseMode','ConsistencySeverity', ...
    'MaxStrokeVolumeRelDiff','BaselinePhysiologyValid','CalibratedPhysiologyValid', ...
    'PrimaryRMSE','FullRMSE','PrimaryMAE_pct','CalibrationMAE_pct', ...
    'MaxAbsError_pct','AllN','AllPass5','AllPass10', ...
    'CalibrationN','CalibrationPass5','CalibrationPass10', ...
    'PrimaryN','PrimaryPass5','PrimaryPass10', ...
    'PressureFlowN','PressureFlowPass5','PressureFlowPass10'};
end

function row = consistency_row(label, audit, clinical, scenario)
src = clinical.(scenario);
row = {label, audit.severity, audit.max_relative_difference, ...
    field_or_nan(clinical.common, 'HR'), field_or_nan(src, 'CO_Lmin'), ...
    field_or_nan(src, 'QpQs'), audit.SV_Qs, audit.SV_Qp, ...
    audit.SV_LV, audit.SV_RV, strjoin(audit.flags, '; ')};
end

function names = consistency_variable_names()
names = {'Case','Severity','MaxStrokeVolumeRelDiff','HR_bpm', ...
    'CO_Lmin','QpQs','SV_Qs_mL','SV_Qp_mL','SV_LV_mL', ...
    'SV_RV_mL','Flags'};
end

function value = field_or_nan(s, field_name)
if isstruct(s) && isfield(s, field_name) && isnumeric(s.(field_name)) && ...
        isscalar(s.(field_name)) && isfinite(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

function write_markdown_report(path, result)
fid = fopen(path, 'w');
if fid < 0
    warning('run_patient_cohort_fast_pipeline:reportOpenFailed', ...
        'Unable to write markdown report.');
    return;
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# Patient Cohort Fast Pipeline Report\n\n');
fprintf(fid, 'Run folder: `%s`\n\n', result.run_dir);
fprintf(fid, '## Summary\n\n');
fprintf(fid, '| Case | Status | Mode | Severity | Primary RMSE | Calib <10%% | Pressure-flow <10%% |\n');
fprintf(fid, '|---|---|---|---|---:|---:|---:|\n');
tbl = result.summary_table;
for idx = 1:height(tbl)
    fprintf(fid, '| %s | %s | %s | %s | %.4f | %d/%d | %d/%d |\n', ...
        tbl.Case{idx}, tbl.Status{idx}, tbl.CaseMode{idx}, ...
        tbl.ConsistencySeverity{idx}, tbl.PrimaryRMSE(idx), ...
        tbl.CalibrationPass10(idx), tbl.CalibrationN(idx), ...
        tbl.PressureFlowPass10(idx), tbl.PressureFlowN(idx));
end
fprintf(fid, '\n## Interpretation Notes\n\n');
fprintf(fid, '- Fast calibration is cohort triage, not the final publication-grade run.\n');
fprintf(fid, '- Cases with `QpQs < 1` are retained as raw catheter evidence and should be reviewed.\n');
fprintf(fid, '- Missing anthropometry blocks scaling and simulation.\n');
fprintf(fid, '- Use the error table to identify target classes that repeatedly fail.\n');
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
