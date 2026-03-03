function report = validation_report(clinical, metrics_baseline, metrics_cal, scenario)
% VALIDATION_REPORT
% -----------------------------------------------------------------------
% Produces a scenario-aware comparison table and RMSE summary.
% Combines Hafiz-style formatted output with Keisya-style RMSE computation.
%
% INPUTS:
%   clinical         - unified clinical struct (from patient_template.m)
%   metrics_baseline - struct from compute_clinical_indices (pre-calibration)
%   metrics_cal      - struct from compute_clinical_indices (post-calibration)
%                      pass [] to skip calibrated column
%   scenario         - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   report           - struct with:
%       .table_baseline   MATLAB table: metric, clinical, baseline, error%
%       .table_cal        MATLAB table: metric, clinical, calibrated, error%
%       .rmse_baseline    scalar overall RMSE (dimensionless, normalised)
%       .rmse_cal         scalar overall RMSE after calibration
%
% METRIC ROWS:
%   pre_surgery:  RAP_mean, PAP_min/max/mean, QpQs, PVR, SVR,
%                 LVEDV, LVESV, RVEDV, RVESV, LVEF
%   post_surgery: SAP_min/max/mean, MAP (=SAP_mean), SVR, PVR,
%                 LVEF, RVEF, QpQs, LVEDV, RVEDV
%   Rows with NaN clinical value are shown but excluded from RMSE.
%
% REFERENCES:
%   [1] Hafiz validation report style (VSD Model V-series).
%   [2] Keisya compare_metrics_table / compute_overall_rmse approach.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

%% Select scenario-specific metric list
switch scenario
    case 'pre_surgery'
        src = clinical.pre_surgery;
        metric_defs = {
            'RAP_mean',  'RAP_mean_mmHg',  'mmHg',   'Right atrial mean pressure'
            'PAP_min',   'PAP_sys_mmHg',   'mmHg',   'PA systolic pressure'
            'PAP_max',   'PAP_sys_mmHg',   'mmHg',   'PA systolic pressure (max)'
            'PAP_mean',  'PAP_mean_mmHg',  'mmHg',   'PA mean pressure'
            'SAP_mean',  'SAP_mean_mmHg',  'mmHg',   'Mean arterial pressure (MAP)'
            'QpQs',      'QpQs',           '—',      'Pulmonary/Systemic flow ratio'
            'PVR',       'PVR_WU',         'WU',     'Pulmonary vascular resistance'
            'SVR',       'SVR_WU',         'WU',     'Systemic vascular resistance'
            'LVEDV',     'LVEDV_mL',       'mL',     'LV end-diastolic volume'
            'LVESV',     'LVESV_mL',       'mL',     'LV end-systolic volume'
            'RVEDV',     'RVEDV_mL',       'mL',     'RV end-diastolic volume'
            'RVESV',     'RVESV_mL',       'mL',     'RV end-systolic volume'
            'LVEF',      'LVEF',           '—',      'LV ejection fraction'
            };

    case 'post_surgery'
        src = clinical.post_surgery;
        metric_defs = {
            'SAP_min',   'SAP_sys_mmHg',   'mmHg',   'Systemic arterial systolic (min)'
            'SAP_max',   'SAP_sys_mmHg',   'mmHg',   'Systemic arterial systolic (max)'
            'SAP_mean',  'MAP_mmHg',       'mmHg',   'Mean arterial pressure'
            'SVR',       'SVR_WU',         'WU',     'Systemic vascular resistance'
            'PVR',       'PVR_WU',         'WU',     'Pulmonary vascular resistance'
            'LVEF',      'LVEF',           '—',      'LV ejection fraction'
            'RVEF',      'RVEF',           '—',      'RV ejection fraction'
            'QpQs',      'QpQs',           '—',      'Qp/Qs ratio (should be ~1.0)'
            'LVEDV',     'LVEDV_mL',       'mL',     'LV end-diastolic volume'
            'RVEDV',     'RVEDV_mL',       'mL',     'RV end-diastolic volume'
            };
    otherwise
        error('validation_report:unknownScenario', ...
              'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

n_rows = size(metric_defs, 1);

%% Build output arrays
metric_names = cell(n_rows, 1);
units_col    = cell(n_rows, 1);
descr_col    = cell(n_rows, 1);
clin_col     = nan(n_rows, 1);
base_col     = nan(n_rows, 1);
cal_col      = nan(n_rows, 1);

for i = 1:n_rows
    model_field  = metric_defs{i, 1};   % field in metrics struct
    clin_field   = metric_defs{i, 2};   % field in clinical.pre/post_surgery
    metric_names{i} = model_field;
    units_col{i}    = metric_defs{i, 3};
    descr_col{i}    = metric_defs{i, 4};

    if isfield(src, clin_field) && ~isnan(src.(clin_field))
        clin_col(i) = src.(clin_field);
    end
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

end  % validation_report

% =========================================================================
function rmse = compute_rmse(y_clin, y_model)
% COMPUTE_RMSE — normalised RMSE across valid (non-NaN) rows
valid = ~isnan(y_clin) & ~isnan(y_model);
if ~any(valid), rmse = NaN; return; end
pct_err = (y_model(valid) - y_clin(valid)) ./ max(abs(y_clin(valid)), 1e-9);
rmse = sqrt(mean(pct_err.^2));
end
