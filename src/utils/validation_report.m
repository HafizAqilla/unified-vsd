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
n_rows = numel(targets);

%% Build output arrays
metric_names = cell(n_rows, 1);
units_col    = cell(n_rows, 1);
descr_col    = cell(n_rows, 1);
clin_col     = nan(n_rows, 1);
base_col     = nan(n_rows, 1);
cal_col      = nan(n_rows, 1);

for i = 1:n_rows
    model_field  = targets(i).Metric;   % field in metrics struct
    metric_names{i} = model_field;
    units_col{i}    = targets(i).Unit;
    descr_col{i}    = targets(i).Description;
    clin_col(i)     = targets(i).ClinicalValue;

    if isfield(metrics_baseline, model_field)
        base_col(i) = metrics_baseline.(model_field);
    end
    if ~isempty(metrics_cal) && isfield(metrics_cal, model_field)
        cal_col(i) = metrics_cal.(model_field);
    end
end

%% Percentage errors
base_err = 100 * (base_col - clin_col) ./ max(abs(clin_col), 1e-9);
cal_err  = 100 * (cal_col  - clin_col) ./ max(abs(clin_col), 1e-9);

base_err(isnan(clin_col)) = NaN;
cal_err(isnan(clin_col))  = NaN;

%% Assemble comparison tables
report.table_baseline = table(metric_names, units_col, descr_col, ...
    clin_col, base_col, base_err, ...
    'VariableNames', {'Metric','Unit','Description','Clinical','Baseline','Error_pct'});

if ~isempty(metrics_cal)
    report.table_cal = table(metric_names, units_col, descr_col, ...
        clin_col, cal_col, cal_err, ...
        'VariableNames', {'Metric','Unit','Description','Clinical','Calibrated','Error_pct'});
else
    report.table_cal = [];
end

%% RMSE  (normalised, excluding NaN rows)
report.rmse_baseline = compute_rmse(clin_col, base_col);
report.rmse_cal      = compute_rmse(clin_col, cal_col);

%% Per-metric delta table (baseline → calibrated)
% Delta_pct < 0 means the error shrank (improvement).
% Rows with NaN clinical or both model values NaN are included but delta = NaN.
delta_err = cal_err - base_err;   % [%] negative = improvement
report.table_delta = table(metric_names, units_col, clin_col, base_err, cal_err, delta_err, ...
    'VariableNames', {'Metric','Unit','Clinical_val', ...
                      'BaseErr_pct','CalErr_pct','Delta_pct'});

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
    abs(err_for_sort(sorted_idx)), ...
    'VariableNames', {'Metric','Unit','Clinical','Error_pct','AbsError_pct'});

%% Primary metric gate check (Batch 5: strict 5% target)
primary_metrics = opts.PrimaryMetrics;
if isempty(primary_metrics)
    [primary_metrics, ~] = select_primary_metrics(clinical, opts.GsaInitOut, scenario);
end
report.primary_metrics = primary_metrics(:)';
report.primary_gate = primary_metric_gate(report.table_cal, report.table_baseline, report.primary_metrics);

%% Print to console
fprintf('\n==========================================================\n');
fprintf('  VALIDATION REPORT — %s\n', upper(strrep(scenario,'_',' ')));
fprintf('==========================================================\n');
disp(report.table_baseline);
fprintf('RMSE Baseline:    %.4f\n', report.rmse_baseline);
if ~isempty(metrics_cal)
    fprintf('\n--- After calibration ---\n');
    disp(report.table_cal);
    fprintf('RMSE Calibrated:  %.4f   (improvement: %.1f%%)\n', ...
        report.rmse_cal, ...
        100*(report.rmse_baseline - report.rmse_cal)/max(report.rmse_baseline, 1e-9));
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
parse(parser, varargin{:});
opts = parser.Results;
if ischar(opts.PrimaryMetrics)
    opts.PrimaryMetrics = {opts.PrimaryMetrics};
elseif isstring(opts.PrimaryMetrics)
    opts.PrimaryMetrics = cellstr(opts.PrimaryMetrics);
end
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
function rmse = compute_rmse(y_clin, y_model)
% COMPUTE_RMSE — normalised RMSE across valid (non-NaN) rows
valid = ~isnan(y_clin) & ~isnan(y_model);
if ~any(valid), rmse = NaN; return; end
pct_err = (y_model(valid) - y_clin(valid)) ./ max(abs(y_clin(valid)), 1e-9);
rmse = sqrt(mean(pct_err.^2));
end
