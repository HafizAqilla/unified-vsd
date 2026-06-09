%% compare_baseline_scaling_methods_demo.m
% -------------------------------------------------------------------------
% Compare adult baseline simulation against pediatric Zhang and Lundquist
% scaling baselines using the default parameters from:
% Baseline analisis/default_parameters_demo_metrics_20260514_215213/
% default_parameters_demo_metrics_20260514_215213.txt
%
% OUTPUTS:
%   results/baseline_scaling_comparison_demo/tables/baseline_scaling_metrics_*.xlsx
%   results/baseline_scaling_comparison_demo/figures/baseline_scaling_*.pdf
%
% AUTHOR:   Unified VSD Model (Demo Run)
% DATE:     2026-06-08
% VERSION:  1.0
% -------------------------------------------------------------------------

clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

output_root = fullfile(project_root, 'results', 'baseline_scaling_comparison_demo');
figure_dir = fullfile(output_root, 'figures');
table_dir = fullfile(output_root, 'tables');
if ~exist(figure_dir, 'dir'), mkdir(figure_dir); end
if ~exist(table_dir, 'dir'), mkdir(table_dir); end

fprintf('============================================================\n');
fprintf('  BASELINE SCALING METHOD COMPARISON (DEMO PARAMETERS)\n');
fprintf('============================================================\n');
fprintf('Adult reference: default_parameters_demo(), no patient scaling.\n');
fprintf('Pediatric scaling: Zhang weight allometry vs Lundquist BSA allometry.\n\n');

%% 1. Adult reference baseline
params_adult = default_parameters_demo();
params_adult = ensure_absolute_timing(params_adult);

fprintf('[baseline-demo] Integrating adult reference baseline...\n');
sim_adult = integrate_system(params_adult);
metrics_adult = compute_clinical_indices(sim_adult, params_adult);
fprintf('[baseline-demo] Adult reference complete. Steady state: %d\n\n', sim_adult.ss_reached);

%% 2. Pediatric demographics for scaled baselines
patient_baseline.age_years = 3 + 2 / 12;  % [years]
patient_baseline.weight_kg = 14;          % [kg]
patient_baseline.height_cm = 98;          % [cm]
patient_baseline.sex = 'M';               % [-]

fprintf('[baseline-demo] Pediatric baseline subject: %.2f years | %.1f kg | %.0f cm | %s\n\n', ...
    patient_baseline.age_years, patient_baseline.weight_kg, ...
    patient_baseline.height_cm, patient_baseline.sex);

%% 3. Zhang 2019 weight-based baseline
params_ref_zhang = default_parameters_demo();
patient_zhang = patient_baseline;
patient_zhang.scaling_mode = 'zhang';

fprintf('[baseline-demo] Integrating Zhang-scaled pediatric baseline...\n');
params_zhang = apply_scaling(params_ref_zhang, patient_zhang);
sim_zhang = integrate_system(params_zhang);
metrics_zhang = compute_clinical_indices(sim_zhang, params_zhang);
fprintf('[baseline-demo] Zhang baseline complete. Steady state: %d\n\n', sim_zhang.ss_reached);

%% 4. Lundquist BSA-based baseline
params_ref_lundquist = default_parameters_demo();
patient_lundquist = patient_baseline;
patient_lundquist.scaling_mode = 'lundquist_bsa';

fprintf('[baseline-demo] Integrating Lundquist-scaled pediatric baseline...\n');
params_lundquist = apply_scaling(params_ref_lundquist, patient_lundquist);
sim_lundquist = integrate_system(params_lundquist);
metrics_lundquist = compute_clinical_indices(sim_lundquist, params_lundquist);
fprintf('[baseline-demo] Lundquist baseline complete. Steady state: %d\n\n', sim_lundquist.ss_reached);

%% 5. Metrics table
metrics_table = build_metrics_table(metrics_adult, params_adult, ...
    metrics_zhang, params_zhang, metrics_lundquist, params_lundquist);

disp(metrics_table);

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
metrics_file = fullfile(table_dir, sprintf('baseline_scaling_metrics_%s.xlsx', timestamp));
writetable(metrics_table, metrics_file, 'Sheet', 'Metrics');
fprintf('[baseline-demo] Metrics table saved: %s\n', metrics_file);

scaling_table = compare_scaling_methods(default_parameters_demo(), patient_baseline);
scaling_file = fullfile(table_dir, sprintf('baseline_scaling_factors_%s.xlsx', timestamp));
writetable(scaling_table, scaling_file, 'Sheet', 'ScalingFactors');
fprintf('[baseline-demo] Scaling-factor table saved: %s\n', scaling_file);

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
% Calculate LVEDP and RVEDP using the passive filling pressure formula
% to avoid the isovolumic contraction max-volume indexing artifact.
metrics.LVEDP = params.E.LV.EB * (metrics.LVEDV - params.V0.LV);
metrics.RVEDP = params.E.RV.EB * (metrics.RVEDV - params.V0.RV);

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

sgtitle('Baseline Scaling Method Pressure Overlay (Demo)', 'FontName', 'Arial', 'FontSize', 12);
output_file = fullfile(figure_dir, sprintf('baseline_scaling_pressure_overlay_%s.pdf', timestamp));
exportgraphics(figure_handle, output_file, 'ContentType', 'vector', 'Resolution', 300);
fprintf('[baseline-demo] Pressure overlay saved: %s\n', output_file);
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

sgtitle('Baseline Scaling Method PV Loop Overlay (Demo)', 'FontName', 'Arial', 'FontSize', 12);
output_file = fullfile(figure_dir, sprintf('baseline_scaling_pv_overlay_%s.pdf', timestamp));
exportgraphics(figure_handle, output_file, 'ContentType', 'vector', 'Resolution', 300);
fprintf('[baseline-demo] PV overlay saved: %s\n', output_file);
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
Xc = sim.V;                    % [state units]
idx = params.idx;
T_HB = 60 / params.HR;         % [s]
time_mask = t >= (t(end) - T_HB);
tc = t(time_mask);             % [s]
Xc = Xc(time_mask, :);         % [state units]

[P_all, ~] = reconstruct_hemodynamic_signals(tc, Xc, params);
P = P_all;
V.LV = Xc(:, idx.V_LV);        % [mL]
V.RV = Xc(:, idx.V_RV);        % [mL]
end

function params = default_parameters_demo()
% DEFAULT_PARAMETERS_DEMO
% Loaded from: Baseline analisis/default_parameters_demo_metrics_20260518_144652/default_parameters_demo_metrics_20260518_144652.txt
params = struct();

idx.V_RA   = 1;    % Right atrial volume              [mL]
idx.V_RV   = 2;    % Right ventricular volume         [mL]
idx.V_LA   = 3;    % Left atrial volume               [mL]
idx.V_LV   = 4;    % Left ventricular volume          [mL]
idx.V_SAR  = 5;    % Systemic arterial volume         [mL]
idx.Q_SAR  = 6;    % Systemic arterial flow           [mL/s]
idx.V_SC   = 7;    % Systemic capillary volume        [mL]
idx.V_SVEN = 8;    % Systemic venous volume           [mL]
idx.Q_SVEN = 9;    % Systemic venous flow             [mL/s]
idx.V_PAR  = 10;   % Pulmonary arterial volume        [mL]
idx.Q_PAR  = 11;   % Pulmonary arterial flow          [mL/s]
idx.P_PC   = 12;   % Pulmonary capillary pressure     [mmHg]
idx.V_PVEN = 13;   % Pulmonary venous volume          [mL]
idx.Q_PVEN = 14;   % Pulmonary venous flow            [mL/s]
params.idx = idx;

params.sim.nCyclesSteady = 80;
params.sim.nCyclesKeep   = 2;
params.sim.batch_size    = 5;
params.sim.rtol          = 1e-5;
params.sim.atol          = 1e-6;
params.sim.ss_tol_P      = 0.1;
params.sim.ss_tol_V      = 0.1;
params.sim.ss_rtol       = 0.005;

params.HR = 75;

params.conv.WU_to_R    = 60/1000;
params.conv.R_to_WU    = 1000/60;
params.conv.mLs_to_Lmin = 60/1000;
params.conv.Lmin_to_mLs = 1000/60;
params.conv.mmHg_to_Pa = 133.322;
params.conv.mm_to_m    = 1e-3;
params.conv.m3_to_mL   = 1e6;

params.epsilon_valve = 0.5;
params.epsilon_vsd = 0.1;

params.vsd.mode = 'linear_bidirectional';        % default calibration mode (R.vsd active)
params.vsd.area_mm2 = 0.0;                       % [mm^2]
params.vsd.diameter_mm = 0.0;                   % [mm]
params.vsd.Cd = 0.7;                            % [-] discharge coefficient
params.vsd.rho_blood = 1060;                    % [kg/m^3]
params.vsd.reference_gradient_mmHg = 20;        % [mmHg]
params.vsd.reverse_leak_fraction = 0.02;        % [-] used only in legacy mode
params.maturation.mode = 'none';
params.calibration.primary_target_pct = 5;
params.calibration.secondary_target_pct = 10;
params.calibration.secondary_lambda = 1.0;
params.calibration.invalid_penalty_scale = 1e3;

params.V0.RA = 3.5385;
params.V0.RV = 8.4067;
params.V0.LA = 2.3085;
params.V0.LV = 3.5385;

params.C.RA = 5.0;        % [mL/mmHg]  atrial placeholder (not used in elastance state eq.)
params.C.LA = 5.0;        % [mL/mmHg]  atrial placeholder (not used in elastance state eq.)

params.E.LV.EA = 3.5;
params.E.LV.EB = 0.08;
params.E.RV.EA = 0.5;
params.E.RV.EB = 0.042;
params.E.LA.EA = 0.35;
params.E.LA.EB = 0.20;
params.E.RA.EA = 0.30;
params.E.RA.EB = 0.10;

params.Tc_LV_frac  = 0.265;
params.Tr_LV_frac  = 0.40;
params.Tc_RV_frac  = 0.30;
params.Tr_RV_frac  = 0.40;

params.t_ac_LA_frac = 0.75;
params.Tc_LA_frac   = 0.10;
params.Tr_LA_frac   = 0.80;
params.t_ac_RA_frac = 0.80;
params.Tc_RA_frac   = 0.10;
params.Tr_RA_frac   = 0.70;

params.R.SAR  = 0.050;
params.C.SAR  = 1.33;
params.L.SAR  = 1.0e-5;

params.R.SC   = 0.80;
params.C.SC   = 1.000;

params.R.SVEN = 0.050;
params.C.SVEN = 30.0;
params.L.SVEN = 1.0e-5;

params.V0.SAR  = 90.9;
params.V0.SC   = 49.0;
params.V0.SVEN = 450.0;

params.R.PAR  = 0.010;
params.C.PAR  = 5.000;
params.L.PAR  = 1.0e-5;

params.R.PCOX = 0.0630;
params.C.PCOX = 0.1983;

params.R.PCNO = 1.2640;
params.C.PCNO = 0.0017;

params.R.PVEN = 0.020;
params.C.PVEN = 13.750;
params.L.PVEN = 1.0e-5;

params.V0.PAR  = 196.0;
params.V0.PCOX = 0.0;
params.V0.PCNO = 0.0;
params.V0.PVEN = 0.0;

params.Rvalve.open   = 4.0e-3;
params.Rvalve.closed = 9.4168e+4;

params.R.vsd = 1e6;

ic = zeros(14, 1);
ic(idx.V_RA)   = 25.0;
ic(idx.V_RV)   = 120.0;
ic(idx.V_LA)   = 25.0;
ic(idx.V_LV)   = 120.0;
ic(idx.V_SAR)  = 202.5;
ic(idx.Q_SAR)  = 100.0;
ic(idx.V_SC)   = 100.0;
ic(idx.V_SVEN) = 750.0;
ic(idx.Q_SVEN) = 100.0;
ic(idx.V_PAR)  = 261.0;
ic(idx.Q_PAR)  = 100.0;
ic(idx.P_PC)   = 10.0;
ic(idx.V_PVEN) = 55.0;
ic(idx.Q_PVEN) = 100.0;
params.ic.V = ic';

end
