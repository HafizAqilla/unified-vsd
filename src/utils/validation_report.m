function report = validation_report(clinical, metrics_baseline, metrics_cal, scenario, varargin)
% VALIDATION_REPORT
% -----------------------------------------------------------------------
% Produces a scenario-aware comparison table, RMSE summary, primary target
% gate checks, sorted error ranking, and baseline-to-calibrated delta table.
% Combines Hafiz-style formatted output with Keisya-style RMSE computation.
%
% INPUTS:
%   clinical         - unified clinical struct (from patient_template.m)
%   metrics_baseline - struct from compute_clinical_indices (pre-calibration)
%   metrics_cal      - struct from compute_clinical_indices (post-calibration)
%                      pass [] to skip calibrated column
%   scenario         - 'pre_surgery' | 'post_surgery'
%   varargin         - optional name/value inputs:
%       'ResultsDir'   output directory for exported tables [char/string]
%       'GsaInitOut'   initial Sobol output struct from gsa_run_sobol
%       'GsaFinalOut'  final Sobol output struct from gsa_run_sobol
%
% OUTPUTS:
%   report           - struct with:
%       .table_baseline   MATLAB table: metric, clinical, baseline, error%
%       .table_cal        MATLAB table: metric, clinical, calibrated, error%
%       .rmse_baseline    scalar overall RMSE (dimensionless, normalised)
%       .rmse_cal         scalar overall RMSE after calibration
%       .primary_gate     table for selected primary metrics at 5% error
%       .table_delta      per-metric delta: BaseErr_pct, CalErr_pct, Delta_pct
%       .sorted_errors    calibrated errors sorted by |Error_pct| descending
%
% METRIC ROWS:
%   pre_surgery:  RAP/LAP min/mean/max, PAP_min/max/mean, QpQs, PVR, SVR,
%                 CO_Lmin, VSD_frac_pct, LVEDV, LVESV, RVEDV, RVESV, LVEF
%   post_surgery: SAP_min/max/mean, MAP (=SAP_mean), SVR, PVR,
%                 LVEF, RVEF, QpQs, CO_Lmin, VSD_frac_pct, LVEDV, RVEDV
%   Rows with NaN clinical value are shown but excluded from RMSE.
%
% REFERENCES:
%   [1] Hafiz validation report style (VSD Model V-series).
%   [2] Keisya compare_metrics_table / compute_overall_rmse approach.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-08
% VERSION:  3.0
% -----------------------------------------------------------------------

opts = parse_options(varargin{:});

targets = get_calibration_targets(scenario, clinical);
clinical_audit = opts.ClinicalConsistencyAudit;
if isempty(clinical_audit)
    clinical_audit = audit_clinical_consistency(clinical, scenario);
end
target_tiers = opts.TargetTiers;
if isempty(target_tiers)
    target_tiers = build_target_tiers(clinical, scenario, clinical_audit);
end
n_rows = numel(targets);

%% Build output arrays
metric_names = cell(n_rows, 1);
units_col    = cell(n_rows, 1);
descr_col    = cell(n_rows, 1);
comparator_col = cell(n_rows, 1);
sigma_rel_col = nan(n_rows, 1);
sigma_abs_col = nan(n_rows, 1);
clin_col     = nan(n_rows, 1);
base_col     = nan(n_rows, 1);
cal_col      = nan(n_rows, 1);
tier_col     = cell(n_rows, 1);
included_cal_col = false(n_rows, 1);
included_primary_rmse_col = false(n_rows, 1);
flag_col     = cell(n_rows, 1);

for i = 1:n_rows
    model_field  = targets(i).Metric;   % field in metrics struct
    metric_names{i} = model_field;
    units_col{i}    = targets(i).Unit;
    descr_col{i}    = targets(i).Description;
    comparator_col{i} = targets(i).Comparator;
    sigma_rel_col(i) = targets(i).UncertaintyFraction;
    sigma_abs_col(i) = targets(i).UncertaintyAbs;
    clin_col(i)     = targets(i).ClinicalValue;

    if isfield(metrics_baseline, model_field)
        base_col(i) = metrics_baseline.(model_field);
    end
    if ~isempty(metrics_cal) && isfield(metrics_cal, model_field)
        cal_col(i) = metrics_cal.(model_field);
    end

    [tier_col{i}, included_cal_col(i), included_primary_rmse_col(i), flag_col{i}] = ...
        target_tier_metadata(target_tiers, model_field);
end

%% Percentage errors
base_err = 100 * (base_col - clin_col) ./ max(abs(clin_col), 1e-9);
cal_err  = 100 * (cal_col  - clin_col) ./ max(abs(clin_col), 1e-9);

base_err(isnan(clin_col)) = NaN;
cal_err(isnan(clin_col))  = NaN;

%% Assemble comparison tables
report.table_baseline = table(metric_names, units_col, descr_col, ...
    comparator_col, sigma_rel_col, sigma_abs_col, ...
    clin_col, base_col, base_err, tier_col, included_cal_col, ...
    included_primary_rmse_col, flag_col, ...
    'VariableNames', {'Metric','Unit','Description','Comparator', ...
    'UncertaintyFraction','UncertaintyAbs','Clinical','Baseline','Error_pct', ...
    'Tier','IncludedInCalibration','IncludedInPrimaryRMSE','Flag'});

if ~isempty(metrics_cal)
    report.table_cal = table(metric_names, units_col, descr_col, ...
        comparator_col, sigma_rel_col, sigma_abs_col, ...
        clin_col, cal_col, cal_err, tier_col, included_cal_col, ...
        included_primary_rmse_col, flag_col, ...
        'VariableNames', {'Metric','Unit','Description','Comparator', ...
        'UncertaintyFraction','UncertaintyAbs','Clinical','Calibrated','Error_pct', ...
        'Tier','IncludedInCalibration','IncludedInPrimaryRMSE','Flag'});
else
    report.table_cal = [];
end

%% RMSE  (normalised, with explicit target-governance masks)
mask_full = true(n_rows, 1);
mask_primary = included_primary_rmse_col;
mask_hard = strcmp(tier_col, 'hard');
mask_soft = strcmp(tier_col, 'soft');

report.rmse_full_baseline = compute_rmse(clin_col, base_col, mask_full);
report.rmse_full_cal      = compute_rmse(clin_col, cal_col, mask_full);
report.rmse_primary_baseline = compute_rmse(clin_col, base_col, mask_primary);
report.rmse_primary_cal      = compute_rmse(clin_col, cal_col, mask_primary);
report.rmse_hard_baseline = compute_rmse(clin_col, base_col, mask_hard);
report.rmse_hard_cal      = compute_rmse(clin_col, cal_col, mask_hard);
report.rmse_soft_baseline = compute_rmse(clin_col, base_col, mask_soft);
report.rmse_soft_cal      = compute_rmse(clin_col, cal_col, mask_soft);

% Backward-compatible fields now point to the governed primary RMSE.
% Full RMSE remains available in rmse_full_* for transparent reporting.
report.rmse_baseline = report.rmse_primary_baseline;
report.rmse_cal      = report.rmse_primary_cal;
report.rmse_summary = build_rmse_summary(report);
report.target_tiers = target_tiers;
report.target_tier_table = target_tiers.table;
report.clinical_consistency_audit = clinical_audit;
report.clinical_data_rank = audit_clinical_data_availability(clinical, scenario);

%% Per-metric delta table (baseline → calibrated)
% Delta_pct < 0 means the error shrank (improvement).
% Rows with NaN clinical or both model values NaN are included but delta = NaN.
delta_err = cal_err - base_err;   % [%] negative = improvement
report.table_delta = table(metric_names, units_col, clin_col, base_err, cal_err, delta_err, ...
    tier_col, included_cal_col, included_primary_rmse_col, flag_col, ...
    'VariableNames', {'Metric','Unit','Clinical_val', ...
                      'BaseErr_pct','CalErr_pct','Delta_pct', ...
                      'Tier','IncludedInCalibration','IncludedInPrimaryRMSE','Flag'});

%% Sorted top absolute errors (calibrated, or baseline if no calibration)
err_for_sort = cal_err;
if all(isnan(err_for_sort))
    err_for_sort = base_err;
end
valid_rows   = find(~isnan(err_for_sort));
[~, sort_ix] = sort(abs(err_for_sort(valid_rows)), 'descend');
sorted_idx   = valid_rows(sort_ix);
report.sorted_errors = table( ...
    metric_names(sorted_idx), units_col(sorted_idx), ...
    clin_col(sorted_idx), err_for_sort(sorted_idx), ...
    abs(err_for_sort(sorted_idx)), tier_col(sorted_idx), flag_col(sorted_idx), ...
    'VariableNames', {'Metric','Unit','Clinical','Error_pct','AbsError_pct', ...
    'Tier','Flag'});

%% Primary metric gate check (Batch 5: strict 5% target)
primary_metrics = opts.PrimaryMetrics;
if isempty(primary_metrics)
    [primary_metrics, ~] = select_primary_metrics(clinical, opts.GsaInitOut, scenario);
end
report.primary_metrics = primary_metrics(:)';
report.validation_holdout_metrics = opts.ValidationHoldoutMetrics(:)';
report.primary_gate = primary_metric_gate(report.table_cal, report.table_baseline, report.primary_metrics);

%% Print to console
fprintf('\n==========================================================\n');
fprintf('  VALIDATION REPORT — %s\n', upper(strrep(scenario,'_',' ')));
fprintf('==========================================================\n');
print_clinical_consistency_audit(clinical_audit);
fprintf('\n--- TARGET TIER GOVERNANCE ---\n');
disp(target_tiers.table);
fprintf('\n--- CLINICAL DATA AVAILABILITY RANK ---\n');
print_clinical_data_rank(report.clinical_data_rank);
disp(report.table_baseline);
fprintf('RMSE Baseline (primary governed): %.4f\n', report.rmse_primary_baseline);
fprintf('RMSE Baseline (full transparent): %.4f\n', report.rmse_full_baseline);
if ~isempty(metrics_cal)
    fprintf('\n--- After calibration ---\n');
    disp(report.table_cal);
    fprintf('RMSE Calibrated (primary governed): %.4f   (improvement: %.1f%%)\n', ...
        report.rmse_primary_cal, ...
        100*(report.rmse_primary_baseline - report.rmse_primary_cal)/max(report.rmse_primary_baseline, 1e-9));
    fprintf('RMSE Calibrated (full transparent): %.4f\n', report.rmse_full_cal);
    fprintf('RMSE Hard-only: %.4f | Soft-only: %.4f\n', ...
        report.rmse_hard_cal, report.rmse_soft_cal);
end

print_primary_gate(report.primary_gate);

%% Sorted top absolute errors
print_sorted_errors(report.sorted_errors);

%% Baseline → calibrated delta table
if ~isempty(metrics_cal)
    print_delta_table(report.table_delta);
end

%% Export overlay table for initial vs final GSA (Batch 5)
if ~isempty(opts.ResultsDir) && ~isempty(opts.GsaInitOut) && ~isempty(opts.GsaFinalOut)
    export_gsa_overlay_tables(opts.ResultsDir, scenario, opts.GsaInitOut, opts.GsaFinalOut);
end

end  % validation_report

% =========================================================================
function opts = parse_options(varargin)
% PARSE_OPTIONS — parse optional name/value arguments.
parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'ResultsDir', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'GsaInitOut', [], @(x) isempty(x) || isstruct(x));
addParameter(parser, 'GsaFinalOut', [], @(x) isempty(x) || isstruct(x));
addParameter(parser, 'PrimaryMetrics', {}, @(x) iscell(x) || isstring(x));
addParameter(parser, 'ValidationHoldoutMetrics', {}, @(x) iscell(x) || isstring(x));
addParameter(parser, 'TargetTiers', [], @(x) isempty(x) || isstruct(x));
addParameter(parser, 'ClinicalConsistencyAudit', [], @(x) isempty(x) || isstruct(x));
parse(parser, varargin{:});
opts = parser.Results;
if ischar(opts.PrimaryMetrics)
    opts.PrimaryMetrics = {opts.PrimaryMetrics};
elseif isstring(opts.PrimaryMetrics)
    opts.PrimaryMetrics = cellstr(opts.PrimaryMetrics);
end
if ischar(opts.ValidationHoldoutMetrics)
    opts.ValidationHoldoutMetrics = {opts.ValidationHoldoutMetrics};
elseif isstring(opts.ValidationHoldoutMetrics)
    opts.ValidationHoldoutMetrics = cellstr(opts.ValidationHoldoutMetrics);
end
end

function [tier, included_cal, included_primary_rmse, flag] = target_tier_metadata(target_tiers, metric_name)
% TARGET_TIER_METADATA - lookup target governance for one metric.
tier = 'validation_only';
included_cal = false;
included_primary_rmse = true;
flag = 'none';
if ~isstruct(target_tiers) || ~isfield(target_tiers, 'table') || isempty(target_tiers.table)
    return;
end
tier_tbl = target_tiers.table;
idx = find(strcmp(tier_tbl.Metric, metric_name), 1, 'first');
if isempty(idx)
    return;
end
tier = tier_tbl.Tier{idx};
included_cal = tier_tbl.IncludedInCalibration(idx);
included_primary_rmse = tier_tbl.IncludedInPrimaryRMSE(idx);
if ismember('Flag', tier_tbl.Properties.VariableNames)
    flag = tier_tbl.Flag{idx};
end
end

function rmse_summary = build_rmse_summary(report)
% BUILD_RMSE_SUMMARY - compact RMSE table with explicit inclusion semantics.
rmse_summary = table( ...
    {'primary_governed'; 'full_transparent'; 'hard_only'; 'soft_only'}, ...
    [report.rmse_primary_baseline; report.rmse_full_baseline; ...
     report.rmse_hard_baseline; report.rmse_soft_baseline], ...
    [report.rmse_primary_cal; report.rmse_full_cal; ...
     report.rmse_hard_cal; report.rmse_soft_cal], ...
    {'Excludes consistency-check-only targets'; ...
     'Includes all available clinical validation targets'; ...
     'Includes hard target tier only'; ...
     'Includes soft target tier only'}, ...
    'VariableNames', {'RMSE_Type','Baseline','Calibrated','Definition'});
end

function print_clinical_consistency_audit(audit)
% PRINT_CLINICAL_CONSISTENCY_AUDIT - concise audit display.
if ~isstruct(audit) || ~isfield(audit, 'summary')
    return;
end
fprintf('\n--- CLINICAL CONSISTENCY AUDIT ---\n');
fprintf('%s\n', audit.summary);
if isfield(audit, 'stroke_volume_table') && ~isempty(audit.stroke_volume_table)
    disp(audit.stroke_volume_table);
end
if isfield(audit, 'pairwise_table') && ~isempty(audit.pairwise_table)
    disp(audit.pairwise_table);
end
if isfield(audit, 'notes') && ~isempty(audit.notes)
    for idx = 1:numel(audit.notes)
        fprintf('[audit] %s\n', audit.notes{idx});
    end
end
end

function print_clinical_data_rank(evidence_table)
% PRINT_CLINICAL_DATA_RANK - compact evidence availability summary.
if isempty(evidence_table)
    return;
end
classes = unique(evidence_table.EvidenceClass, 'stable');
count_col = zeros(numel(classes), 1);
for idx = 1:numel(classes)
    count_col(idx) = nnz(strcmp(evidence_table.EvidenceClass, classes{idx}));
end
summary_tbl = table(classes, count_col, ...
    'VariableNames', {'EvidenceClass', 'Count'});
disp(summary_tbl);

visible = evidence_table(:, {'AvailabilityRank','Section','Field','Value', ...
    'EvidenceClass','RecommendedUse'});
disp(visible);
end

function gate_tbl = primary_metric_gate(table_cal, table_base, primary_names)
% PRIMARY_METRIC_GATE — evaluate strict 5% absolute error gate.
metric_col = primary_names(:);
clinical_col = nan(numel(primary_names), 1);
model_col = nan(numel(primary_names), 1);
abs_err_pct_col = nan(numel(primary_names), 1);
pass_col = false(numel(primary_names), 1);

src_tbl = table_base;
if ~isempty(table_cal)
    src_tbl = table_cal;
end

for i = 1:numel(primary_names)
    idx = find(strcmp(src_tbl.Metric, primary_names{i}), 1, 'first');
    if isempty(idx)
        continue;
    end

    clinical_col(i) = src_tbl.Clinical(idx);
    if ismember('Calibrated', src_tbl.Properties.VariableNames)
        model_col(i) = src_tbl.Calibrated(idx);
    else
        model_col(i) = src_tbl.Baseline(idx);
    end

    abs_err_pct_col(i) = abs(src_tbl.Error_pct(idx));
    pass_col(i) = ~isnan(abs_err_pct_col(i)) && abs_err_pct_col(i) <= 5.0;
end

gate_tbl = table(metric_col, clinical_col, model_col, abs_err_pct_col, pass_col, ...
    'VariableNames', {'Metric', 'Clinical', 'Model', 'AbsError_pct', 'Pass_5pct'});
end

function print_primary_gate(gate_tbl)
% PRINT_PRIMARY_GATE — print explicit pass/fail warning lines for 5% gate.
fprintf('\n--- PRIMARY 5%% TARGET GATE (GSA-guided primary metrics) ---\n');
disp(gate_tbl);

idx_fail = find(~gate_tbl.Pass_5pct | isnan(gate_tbl.Pass_5pct));
if isempty(idx_fail)
    fprintf('[PASS] All primary metrics are within 5%% absolute error.\n');
    return;
end

for i = 1:numel(idx_fail)
    k = idx_fail(i);
    fprintf(2, '[WARNING] %s exceeds 5%% absolute error (|error| = %.2f%%).\n', ...
        gate_tbl.Metric{k}, gate_tbl.AbsError_pct(k));
end
end

function export_gsa_overlay_tables(results_dir, scenario, gsa_init_out, gsa_final_out)
% EXPORT_GSA_OVERLAY_TABLES — export initial vs final Sobol ST overlay table.
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

overlay_tbl = build_overlay_table(gsa_init_out, gsa_final_out);
if isempty(overlay_tbl)
    warning('validation_report:noGsaOverlayRows', ...
        'No GSA overlay rows were generated; skipping table export.');
    return;
end

csv_path = fullfile(results_dir, sprintf('gsa_overlay_%s.csv', scenario));
mat_path = fullfile(results_dir, sprintf('gsa_overlay_%s.mat', scenario));

writetable(overlay_tbl, csv_path);
save(mat_path, 'overlay_tbl');

fprintf('[validation_report] GSA overlay table saved:\n  %s\n  %s\n', ...
    csv_path, mat_path);
end

function overlay_tbl = build_overlay_table(gsa_init_out, gsa_final_out)
% BUILD_OVERLAY_TABLE — long table with initial/final ST and delta per metric.
metrics_init = fieldnames(gsa_init_out);
metrics_final = fieldnames(gsa_final_out);
ignore_fields = {'scenario', 'cfg'};
metrics_init = metrics_init(~ismember(metrics_init, ignore_fields));
metrics_final = metrics_final(~ismember(metrics_final, ignore_fields));
metrics = intersect(metrics_init, metrics_final, 'stable');

rows = {};
for m = 1:numel(metrics)
    mf = metrics{m};
    if ~isfield(gsa_init_out.(mf), 'table') || ~isfield(gsa_final_out.(mf), 'table')
        continue;
    end

    t_init = gsa_init_out.(mf).table;
    t_final = gsa_final_out.(mf).table;
    if isempty(t_init) || isempty(t_final)
        continue;
    end

    [is_common, idx_final] = ismember(t_init.Parameter, t_final.Parameter);
    idx_init = find(is_common);
    idx_final = idx_final(is_common);
    if isempty(idx_init)
        continue;
    end

    for r = 1:numel(idx_init)
        p_name = t_init.Parameter{idx_init(r)};
        st_init = t_init.Sobol_ST(idx_init(r));
        st_final = t_final.Sobol_ST(idx_final(r));
        rows(end+1, :) = {mf, p_name, st_init, st_final, st_final - st_init}; %#ok<AGROW>
    end
end

if isempty(rows)
    overlay_tbl = table();
    return;
end

overlay_tbl = cell2table(rows, ...
    'VariableNames', {'Metric', 'Parameter', 'ST_Initial', 'ST_Final', 'Delta_ST'});
overlay_tbl = sortrows(overlay_tbl, {'Metric', 'ST_Final'}, {'ascend', 'descend'});
end

% =========================================================================
function print_sorted_errors(sorted_tbl)
% PRINT_SORTED_ERRORS — print metrics ranked by absolute error (worst first).
TOP_N = 5;
n = min(TOP_N, height(sorted_tbl));
fprintf('\n--- TOP %d ABSOLUTE ERRORS ---\n', n);
fprintf('  %-10s  %6s  %10s  %10s  %10s\n', ...
        'Metric', 'Unit', 'Clinical', 'Error(%)', '|Error|(%)');
fprintf('  %s\n', repmat('-', 1, 54));
for i = 1:n
    fprintf('  %-10s  %6s  %10.4g  %+10.2f  %10.2f\n', ...
        sorted_tbl.Metric{i}, sorted_tbl.Unit{i}, ...
        sorted_tbl.Clinical(i), sorted_tbl.Error_pct(i), ...
        sorted_tbl.AbsError_pct(i));
end
end

% =========================================================================
function print_delta_table(delta_tbl)
% PRINT_DELTA_TABLE — per-metric baseline→calibrated error delta table.
% Rows sorted by |delta| descending so largest changes appear first.
% Delta_pct < 0 is an improvement (error shrank).
valid = ~isnan(delta_tbl.Delta_pct);
if ~any(valid)
    return;
end
[~, sort_ix] = sort(abs(delta_tbl.Delta_pct(valid)), 'descend');
rows = find(valid);
rows = rows(sort_ix);

fprintf('\n--- BASELINE → CALIBRATED DELTA (sorted by |change|) ---\n');
fprintf('  %-10s  %8s  %8s  %8s  %s\n', ...
        'Metric', 'Base(%)', 'Cal(%)', 'Delta(%)', 'Dir');
fprintf('  %s\n', repmat('-', 1, 52));
for i = 1:numel(rows)
    r   = rows(i);
    dir = '';
    if ~isnan(delta_tbl.Delta_pct(r))
        if delta_tbl.Delta_pct(r) < -0.1
            dir = '[+]';   % improved
        elseif delta_tbl.Delta_pct(r) > 0.1
            dir = '[-]';   % worsened
        else
            dir = '[=]';
        end
    end
    fprintf('  %-10s  %+8.2f  %+8.2f  %+8.2f  %s\n', ...
        delta_tbl.Metric{r}, delta_tbl.BaseErr_pct(r), ...
        delta_tbl.CalErr_pct(r), delta_tbl.Delta_pct(r), dir);
end
end

% =========================================================================
function rmse = compute_rmse(y_clin, y_model, include_mask)
% COMPUTE_RMSE — normalised RMSE across valid (non-NaN) rows
if nargin < 3 || isempty(include_mask)
    include_mask = true(size(y_clin));
end
valid = include_mask(:) & ~isnan(y_clin) & ~isnan(y_model);
if ~any(valid), rmse = NaN; return; end
pct_err = (y_model(valid) - y_clin(valid)) ./ max(abs(y_clin(valid)), 1e-9);
rmse = sqrt(mean(pct_err.^2));
end
