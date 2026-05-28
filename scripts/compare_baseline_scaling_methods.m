%% compare_baseline_scaling_methods.m
% -------------------------------------------------------------------------
% Compare adult baseline simulation against pediatric Zhang and Lundquist
% scaling baselines without modifying the model pipeline.
%
% OUTPUTS:
%   results/baseline_scaling_comparison/tables/baseline_scaling_metrics_*.xlsx
%   results/baseline_scaling_comparison/figures/baseline_scaling_*.pdf
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-26
% VERSION:  1.0
% -------------------------------------------------------------------------

clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

output_root = fullfile(project_root, 'results', 'baseline_scaling_comparison');
figure_dir = fullfile(output_root, 'figures');
table_dir = fullfile(output_root, 'tables');
if ~exist(figure_dir, 'dir'), mkdir(figure_dir); end
if ~exist(table_dir, 'dir'), mkdir(table_dir); end

fprintf('============================================================\n');
fprintf('  BASELINE SCALING METHOD COMPARISON\n');
fprintf('============================================================\n');
fprintf('Adult reference: default_parameters(), no patient scaling.\n');
fprintf('Pediatric scaling: Zhang weight allometry vs Lundquist BSA allometry.\n\n');

%% 1. Adult reference baseline
params_adult = default_parameters();
params_adult = ensure_absolute_timing(params_adult);

fprintf('[baseline] Integrating adult reference baseline...\n');
sim_adult = integrate_system(params_adult);
metrics_adult = compute_clinical_indices(sim_adult, params_adult);
fprintf('[baseline] Adult reference complete. Steady state: %d\n\n', sim_adult.ss_reached);

%% 2. Pediatric demographics for scaled baselines
patient_baseline.age_years = 3 + 2 / 12;  % [years]
patient_baseline.weight_kg = 14;          % [kg]
patient_baseline.height_cm = 98;          % [cm]
patient_baseline.sex = 'M';               % [-]

fprintf('[baseline] Pediatric baseline subject: %.2f years | %.1f kg | %.0f cm | %s\n\n', ...
    patient_baseline.age_years, patient_baseline.weight_kg, ...
    patient_baseline.height_cm, patient_baseline.sex);

%% 3. Zhang 2019 weight-based baseline
params_ref_zhang = default_parameters();
patient_zhang = patient_baseline;
patient_zhang.scaling_mode = 'zhang';

fprintf('[baseline] Integrating Zhang-scaled pediatric baseline...\n');
params_zhang = apply_scaling(params_ref_zhang, patient_zhang);
sim_zhang = integrate_system(params_zhang);
metrics_zhang = compute_clinical_indices(sim_zhang, params_zhang);
fprintf('[baseline] Zhang baseline complete. Steady state: %d\n\n', sim_zhang.ss_reached);

%% 4. Lundquist BSA-based baseline
params_ref_lundquist = default_parameters();
patient_lundquist = patient_baseline;
patient_lundquist.scaling_mode = 'lundquist_bsa';

fprintf('[baseline] Integrating Lundquist-scaled pediatric baseline...\n');
params_lundquist = apply_scaling(params_ref_lundquist, patient_lundquist);
sim_lundquist = integrate_system(params_lundquist);
metrics_lundquist = compute_clinical_indices(sim_lundquist, params_lundquist);
fprintf('[baseline] Lundquist baseline complete. Steady state: %d\n\n', sim_lundquist.ss_reached);

%% 5. Metrics table
metrics_table = build_metrics_table(metrics_adult, params_adult, ...
    metrics_zhang, params_zhang, metrics_lundquist, params_lundquist);

disp(metrics_table);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
metrics_file = fullfile(table_dir, sprintf('baseline_scaling_metrics_%s.xls', timestamp));
writetable(metrics_table, metrics_file, 'Sheet', 'Metrics');
fprintf('[baseline] Metrics table saved: %s\n', metrics_file);

scaling_table = compare_scaling_methods(default_parameters(), patient_baseline);
scaling_file = fullfile(table_dir, sprintf('baseline_scaling_factors_%s.xlsx', timestamp));
writetable(scaling_table, scaling_file, 'Sheet', 'ScalingFactors');
fprintf('[baseline] Scaling-factor table saved: %s\n', scaling_file);

%% 6. Figures
plot_baseline_pressure_overlay(sim_adult, params_adult, ...
    sim_zhang, params_zhang, sim_lundquist, params_lundquist, figure_dir, timestamp);
plot_baseline_pv_overlay(sim_adult, params_adult, ...
    sim_zhang, params_zhang, sim_lundquist, params_lundquist, figure_dir, timestamp);

fprintf('\n============================================================\n');
fprintf('  Done. Outputs written under:\n');
fprintf('  %s\n', output_root);
fprintf('============================================================\n');

function params = ensure_absolute_timing(params)
% ENSURE_ABSOLUTE_TIMING - convert fractional chamber timing to seconds.
T_HB = 60 / params.HR;             % [s]
params.Tc_LV   = params.Tc_LV_frac   * T_HB;  % [s]
params.Tr_LV   = params.Tr_LV_frac   * T_HB;  % [s]
params.Tc_RV   = params.Tc_RV_frac   * T_HB;  % [s]
params.Tr_RV   = params.Tr_RV_frac   * T_HB;  % [s]
params.t_ac_LA = params.t_ac_LA_frac * T_HB;  % [s]
params.Tc_LA   = params.Tc_LA_frac   * T_HB;  % [s]
params.t_ar_LA = params.t_ac_LA + params.Tc_LA; % [s]
params.Tr_LA   = params.Tr_LA_frac   * T_HB;  % [s]
params.t_ac_RA = params.t_ac_RA_frac * T_HB;  % [s]
params.Tc_RA   = params.Tc_RA_frac   * T_HB;  % [s]
params.t_ar_RA = params.t_ac_RA + params.Tc_RA; % [s]
params.Tr_RA   = params.Tr_RA_frac   * T_HB;  % [s]
end

function metrics_table = build_metrics_table(metrics_adult, params_adult, ...
    metrics_zhang, params_zhang, metrics_lundquist, params_lundquist)
% BUILD_METRICS_TABLE - assemble baseline clinical metrics with units.
metric = { ...
    'Heart Rate'
    'Cardiac Output'
    'LV Stroke Volume'
    'RV Stroke Volume'
    'LVEF'
    'RVEF'
    'LVEDV'
    'LVESV'
    'RVEDV'
    'RVESV'
    'SVR'
    'PVR'
    'SBP'
    'DBP'
    'MAP'
    'PAP_sys'
    'PAP_dia'
    'PAP_mean'
    'LVESP'
    'LVEDP'
    'RVESP'
    'RVEDP'
    'LAP_mean'
    'RAP_mean'
    'PWP_mean'
    'Qp_Qs'};

unit = { ...
    'bpm'
    'L/min'
    'mL'
    'mL'
    '%'
    '%'
    'mL'
    'mL'
    'mL'
    'mL'
    'WU'
    'WU'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    'mmHg'
    '-'};

adult_ref = metric_values(metrics_adult, params_adult);
zhang_2019 = metric_values(metrics_zhang, params_zhang);
lundquist_2025 = metric_values(metrics_lundquist, params_lundquist);

[norm_adult, norm_pediatric, adult_bounds, pediatric_bounds] = clinical_reference_columns();
adult_ref_status = range_status(adult_ref, adult_bounds);
zhang_status = range_status(zhang_2019, pediatric_bounds);
lundquist_status = range_status(lundquist_2025, pediatric_bounds);

metrics_table = table(metric, unit, norm_adult, norm_pediatric, ...
    adult_ref, adult_ref_status, zhang_2019, zhang_status, ...
    lundquist_2025, lundquist_status, ...
    'VariableNames', {'Metric','Unit','Norm_Adult','Norm_Pediatric', ...
    'Adult_ref','Adult_ref_status','Zhang_2019','Zhang_status', ...
    'Lundquist_2025','Lundquist_status'});
end

function [norm_adult, norm_pediatric, adult_bounds, pediatric_bounds] = clinical_reference_columns()
% CLINICAL_REFERENCE_COLUMNS - clinical ranges from user-provided sheet.
norm_adult = { ...
    '61-87'
    '4-8'
    '60-100'
    '--'
    '52-72'
    '50-66'
    '87-137'
    '31-51'
    '64.6-160.5'
    '18.5-81.2'
    '9.1-31.5'
    '0.3-2'
    '128-156'
    '80-96'
    '70-100'
    '15-30'
    '4-12'
    '8-20'
    '90-140'
    '3-12'
    '20-30'
    '<8'
    '2-12'
    '2-6'
    '--'
    '1.0 (normal)'};

norm_pediatric = { ...
    '73-142'
    '2.2-4.85'
    '15-25'
    '--'
    '55-73'
    '45-78'
    '25-41'
    '10-16'
    '30.73-65.68'
    '7-29'
    '13.8-33'
    '<6'
    '90-103'
    '47-59'
    '61-74'
    '<32.9'
    '<14.95'
    '<21'
    '--'
    '--'
    '<35'
    '--'
    '2-10'
    '3-6'
    '<12'
    '1.0 (normal)'};

adult_bounds = [ ...
    61 87
    4 8
    60 100
    NaN NaN
    52 72
    50 66
    87 137
    31 51
    64.6 160.5
    18.5 81.2
    9.1 31.5
    0.3 2
    128 156
    80 96
    70 100
    15 30
    4 12
    8 20
    90 140
    3 12
    20 30
    -Inf 8
    2 12
    2 6
    NaN NaN
    0.95 1.05];

% Qp/Qs tolerance uses +/-5% around the nominal normal value of 1.0.
pediatric_bounds = [ ...
    73 142
    2.2 4.85
    15 25
    NaN NaN
    55 73
    45 78
    25 41
    10 16
    30.73 65.68
    7 29
    13.8 33
    -Inf 6
    90 103
    47 59
    61 74
    -Inf 32.9
    -Inf 14.95
    -Inf 21
    NaN NaN
    NaN NaN
    -Inf 35
    NaN NaN
    2 10
    3 6
    -Inf 12
    0.95 1.05];
end

function status = range_status(values, bounds)
% RANGE_STATUS - classify metric values against low/high clinical bounds.
status = cell(size(values));
for i_value = 1:numel(values)
    lower_bound = bounds(i_value, 1);  % [metric unit]
    upper_bound = bounds(i_value, 2);  % [metric unit]
    value = values(i_value);           % [metric unit]

    if isnan(lower_bound) && isnan(upper_bound)
        status{i_value} = 'Not assessed';
    elseif value >= lower_bound && value <= upper_bound
        status{i_value} = 'In range';
    elseif value < lower_bound
        status{i_value} = 'Below range';
    else
        status{i_value} = 'Above range';
    end
end
end

function values = metric_values(metrics, params)
% METRIC_VALUES - extract ordered clinical metric values for one simulation.
values = [ ...
    params.HR
    metrics.CO_Lmin
    metrics.LVSV
    metrics.RVSV
    metrics.LVEF * 100
    metrics.RVEF * 100
    metrics.LVEDV
    metrics.LVESV
    metrics.RVEDV
    metrics.RVESV
    metrics.SVR
    metrics.PVR
    metrics.SAP_max
    metrics.SAP_min
    metrics.SAP_mean
    metrics.PAP_max
    metrics.PAP_min
    metrics.PAP_mean
    metrics.LVP_max
    metrics.LVEDP
    metrics.RVP_max
    metrics.RVEDP
    metrics.LAP_mean
    metrics.RAP_mean
    metrics.PWP_mean
    metrics.QpQs];
end

function plot_baseline_pressure_overlay(sim_adult, params_adult, ...
    sim_zhang, params_zhang, sim_lundquist, params_lundquist, figure_dir, timestamp)
% PLOT_BASELINE_PRESSURE_OVERLAY - plot systemic and chamber pressures.
[t_adult, P_adult] = last_cycle_pressures(sim_adult, params_adult);
[t_zhang, P_zhang] = last_cycle_pressures(sim_zhang, params_zhang);
[t_lundquist, P_lundquist] = last_cycle_pressures(sim_lundquist, params_lundquist);

figure_handle = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 14]);
set(figure_handle, 'Toolbar', 'none');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_pressure_panel(t_adult, P_adult.LV, t_zhang, P_zhang.LV, ...
    t_lundquist, P_lundquist.LV, 'Left Ventricle', 'P_{LV} [mmHg]');
plot_pressure_panel(t_adult, P_adult.RV, t_zhang, P_zhang.RV, ...
    t_lundquist, P_lundquist.RV, 'Right Ventricle', 'P_{RV} [mmHg]');
plot_pressure_panel(t_adult, P_adult.SAR, t_zhang, P_zhang.SAR, ...
    t_lundquist, P_lundquist.SAR, 'Systemic Artery', 'P_{SAR} [mmHg]');
plot_pressure_panel(t_adult, P_adult.PAR, t_zhang, P_zhang.PAR, ...
    t_lundquist, P_lundquist.PAR, 'Pulmonary Artery', 'P_{PAR} [mmHg]');

sgtitle('Baseline Scaling Method Pressure Overlay', 'FontName', 'Arial', 'FontSize', 12);
output_file = fullfile(figure_dir, sprintf('baseline_scaling_pressure_overlay_%s.pdf', timestamp));
exportgraphics(figure_handle, output_file, 'ContentType', 'vector', 'Resolution', 300);
fprintf('[baseline] Pressure overlay saved: %s\n', output_file);
end

function plot_pressure_panel(t_adult, P_adult, t_zhang, P_zhang, ...
    t_lundquist, P_lundquist, panel_title, y_label)
% PLOT_PRESSURE_PANEL - draw one pressure panel with three scaling cases.
nexttile;
plot(t_adult, P_adult, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Adult ref'); hold on;
plot(t_zhang, P_zhang, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Zhang 2019');
plot(t_lundquist, P_lundquist, 'g:', 'LineWidth', 1.8, 'DisplayName', 'Lundquist BSA');
xlabel('Time [s]', 'FontName', 'Arial');
ylabel(y_label, 'FontName', 'Arial');
title(panel_title, 'FontName', 'Arial');
legend('Location', 'best');
grid on;
set(gca, 'FontName', 'Arial', 'Box', 'on');
end

function plot_baseline_pv_overlay(sim_adult, params_adult, ...
    sim_zhang, params_zhang, sim_lundquist, params_lundquist, figure_dir, timestamp)
% PLOT_BASELINE_PV_OVERLAY - plot LV and RV pressure-volume loops.
[~, P_adult, V_adult] = last_cycle_pressures(sim_adult, params_adult);
[~, P_zhang, V_zhang] = last_cycle_pressures(sim_zhang, params_zhang);
[~, P_lundquist, V_lundquist] = last_cycle_pressures(sim_lundquist, params_lundquist);

figure_handle = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 8]);
set(figure_handle, 'Toolbar', 'none');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

plot_pv_panel(V_adult.LV, P_adult.LV, V_zhang.LV, P_zhang.LV, ...
    V_lundquist.LV, P_lundquist.LV, 'LV Pressure-Volume Loop');
plot_pv_panel(V_adult.RV, P_adult.RV, V_zhang.RV, P_zhang.RV, ...
    V_lundquist.RV, P_lundquist.RV, 'RV Pressure-Volume Loop');

sgtitle('Baseline Scaling Method PV Loop Overlay', 'FontName', 'Arial', 'FontSize', 12);
output_file = fullfile(figure_dir, sprintf('baseline_scaling_pv_overlay_%s.pdf', timestamp));
exportgraphics(figure_handle, output_file, 'ContentType', 'vector', 'Resolution', 300);
fprintf('[baseline] PV overlay saved: %s\n', output_file);
end

function plot_pv_panel(V_adult, P_adult, V_zhang, P_zhang, ...
    V_lundquist, P_lundquist, panel_title)
% PLOT_PV_PANEL - draw one ventricular pressure-volume loop panel.
nexttile;
plot(V_adult, P_adult, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Adult ref'); hold on;
plot(V_zhang, P_zhang, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Zhang 2019');
plot(V_lundquist, P_lundquist, 'g:', 'LineWidth', 1.8, 'DisplayName', 'Lundquist BSA');
xlabel('Volume [mL]', 'FontName', 'Arial');
ylabel('Pressure [mmHg]', 'FontName', 'Arial');
title(panel_title, 'FontName', 'Arial');
legend('Location', 'best');
grid on;
set(gca, 'FontName', 'Arial', 'Box', 'on');
end

function [tc, P, V] = last_cycle_pressures(sim, params)
% LAST_CYCLE_PRESSURES - reconstruct pressure traces over the final cycle.
t = sim.t(:);                  % [s]
X = sim.V;                     % [state units]
idx = params.idx;
T_HB = 60 / params.HR;         % [s]
time_mask = t >= (t(end) - T_HB);
tc = t(time_mask);             % [s]
Xc = X(time_mask, :);          % [state units]

[P_all, ~] = reconstruct_hemodynamic_signals(tc, Xc, params);
P = P_all;
V.LV = Xc(:, idx.V_LV);        % [mL]
V.RV = Xc(:, idx.V_RV);        % [mL]
end
