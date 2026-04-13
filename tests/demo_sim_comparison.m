%% demo_sim_comparison.m
% =========================================================================
% Demonstration: apply_scaling → integrate_system → plot_simulation_comparison
%
% PURPOSE:
%   Runs the full ODE pipeline for (1) an adult reference and (2) a
%   healthy paediatric patient, then renders the 4-figure comparison:
%
%     Figure 1 — LV and RV PV loops
%     Figure 2 — Pressure waveforms (LV, RV, SAR, PAR, LA, RA)
%     Figure 3 — Flow waveforms (AV, MV, TV, PVv, SVEN, PVEN)
%     Figure 4 — Clinical metrics panel (EF, SV, SVR/PVR, summary)
%
%   The patient has NO structural defect.  Only allometric scaling
%   (Zhang et al. 2019 via apply_scaling) is applied — no VSD shunt,
%   no clinical overrides, no params_from_clinical.
%
% USAGE:
%   >> cd <project_root>          % e.g. C:\Users\asus\Documents\MATLAB\VSD
%   >> demo_sim_comparison
%
% To change the patient demographics, edit Section 2 below.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-10
% VERSION:  2.0  (healthy paediatric demo)
% =========================================================================

clear; clc; close all;

%% ---- 0. Path setup -----------------------------------------------------
root = fileparts(mfilename('fullpath'));   % /tests
root = fullfile(root, '..');              % project root
addpath(genpath(root));
fprintf('[demo] Project path ready.\n\n');

%% =========================================================================
%  1. ADULT REFERENCE  (Valenti Table 3.3 — raw reference, no scaling)
%
%  Uses default_parameters() directly — no allometric scaling, no clinical
%  overrides.  This guarantees the adult simulation reflects Valenti's
%  reference model exactly (all parameters straight from Table 3.3).
%
%  The only computation needed is converting fractional timing fields
%  (Tc_LV_frac, Tr_LV_frac, …) into absolute seconds at the reference
%  HR = 75 bpm.  apply_scaling Section D normally does this; we replicate
%  it here without calling apply_scaling.
% =========================================================================
params_adult = default_parameters();

% Compute absolute timing fields from fractional values  (Valenti Table 2.1)
T_HB_adult           = 60 / params_adult.HR;
params_adult.Tc_LV   = params_adult.Tc_LV_frac   * T_HB_adult;
params_adult.Tr_LV   = params_adult.Tr_LV_frac   * T_HB_adult;
params_adult.Tc_RV   = params_adult.Tc_RV_frac   * T_HB_adult;
params_adult.Tr_RV   = params_adult.Tr_RV_frac   * T_HB_adult;
params_adult.t_ac_LA = params_adult.t_ac_LA_frac  * T_HB_adult;
params_adult.Tc_LA   = params_adult.Tc_LA_frac    * T_HB_adult;
params_adult.t_ar_LA = params_adult.t_ac_LA + params_adult.Tc_LA;
params_adult.Tr_LA   = params_adult.Tr_LA_frac    * T_HB_adult;
params_adult.t_ac_RA = params_adult.t_ac_RA_frac  * T_HB_adult;
params_adult.Tc_RA   = params_adult.Tc_RA_frac    * T_HB_adult;
params_adult.t_ar_RA = params_adult.t_ac_RA + params_adult.Tc_RA;
params_adult.Tr_RA   = params_adult.Tr_RA_frac    * T_HB_adult;

fprintf('[demo] Integrating adult ODE (raw Valenti Table 3.3, no scaling)...\n');
sim_adult = integrate_system(params_adult);
fprintf('[demo] Adult simulation complete.\n\n');

%% =========================================================================
%  2. HEALTHY PAEDIATRIC PATIENT  ← edit demographics here
%
%  No VSD or other structural defect.  Only Zhang 2019 allometric
%  scaling is applied via apply_scaling().
%  params_from_clinical is NOT called — no pathology to configure.
% =========================================================================

% --- Demographics (typical healthy 6-month-old male)
patient_pt.age_years = 6 / 12;    % [yr]   6 months
patient_pt.weight_kg = 7.5;        % [kg]   50th percentile WHO growth chart
patient_pt.height_cm = 67;         % [cm]   50th percentile WHO growth chart
patient_pt.sex       = 'M';

%% =========================================================================
%  3. SCALE + SIMULATE HEALTHY PAEDIATRIC PATIENT
%
%  apply_scaling applies Zhang 2019 weight-based allometric scaling:
%    HR    ~ w^(-0.30)    V0   ~ w^(+0.80)
%    E_LV  ~ w^(-0.50)    R    ~ w^(-0.475) systemic
%    E_RV  ~ w^(-0.75)    R    ~ w^(-0.70)  pulmonary
%    C     ~ w^(+1.00)
%  params_from_clinical is intentionally skipped (no pathology).
% =========================================================================
fprintf('[demo] Scaling healthy paediatric parameters (Zhang 2019)...\n');
params_ref = default_parameters();           % adult reference baseline
params_pat = apply_scaling(params_ref, patient_pt);   % allometric scaling only

fprintf('[demo] Integrating patient ODE...\n');
sim_pat = integrate_system(params_pat);
fprintf('[demo] Patient simulation complete.\n\n');

%% =========================================================================
%  4. BUILD COMPARISON LABEL
% =========================================================================
label = sprintf('Healthy paediatric (%.0f mo, %.1f kg, M)', ...
    patient_pt.age_years * 12, patient_pt.weight_kg);

%% =========================================================================
%  5. PLOT
% =========================================================================
fprintf('[demo] Rendering 4-figure comparison...\n');
plot_simulation_comparison(sim_adult, params_adult, ...
                            sim_pat,   params_pat, ...
                            label);
fprintf('[demo] Done. 4 figures are now open.\n');

%% =========================================================================
%  6. CLINICAL METRICS — Stroke Volume and Haemodynamic Summary
%
%  compute_clinical_indices derives LVSV, RVSV, EF, CO, SVR, PVR, etc.
%  from the last complete cardiac cycle of each simulation.
% =========================================================================
fprintf('[demo] Computing clinical indices...\n');
metrics_adult = compute_clinical_indices(sim_adult, params_adult);
metrics_pat   = compute_clinical_indices(sim_pat,   params_pat);

% Cardiac output: SV [mL] × HR [bpm] → CO [L/min]
CO_adult = metrics_adult.LVSV * params_adult.HR / 1000;   % [L/min]
CO_pat   = metrics_pat.LVSV   * params_pat.HR   / 1000;   % [L/min]

fprintf('\n');
fprintf('============================================================\n');
fprintf('  CLINICAL METRICS SUMMARY\n');
fprintf('============================================================\n');
fprintf('  Metric                Adult (ref)     Patient (%s)\n', label);
fprintf('  ----------------------------------------------------------\n');
fprintf('  LV Stroke Volume      %6.1f mL       %6.1f mL\n', ...
    metrics_adult.LVSV, metrics_pat.LVSV);
fprintf('  RV Stroke Volume      %6.1f mL       %6.1f mL\n', ...
    metrics_adult.RVSV, metrics_pat.RVSV);
fprintf('  LV Ejection Fraction  %6.1f %%        %6.1f %%\n', ...
    metrics_adult.LVEF * 100, metrics_pat.LVEF * 100);
fprintf('  RV Ejection Fraction  %6.1f %%        %6.1f %%\n', ...
    metrics_adult.RVEF * 100, metrics_pat.RVEF * 100);
fprintf('  Heart Rate            %6.0f bpm      %6.0f bpm\n', ...
    params_adult.HR, params_pat.HR);
fprintf('  Cardiac Output        %6.2f L/min    %6.2f L/min\n', ...
    CO_adult, CO_pat);
fprintf('  SVR                   %6.2f WU       %6.2f WU\n', ...
    metrics_adult.SVR, metrics_pat.SVR);
fprintf('  PVR                   %6.2f WU       %6.2f WU\n', ...
    metrics_adult.PVR, metrics_pat.PVR);
fprintf('  PA sys / dia / mean   %4.0f / %4.0f / %4.0f mmHg   %4.0f / %4.0f / %4.0f mmHg\n', ...
    metrics_adult.PAP_max, metrics_adult.PAP_min, metrics_adult.PAP_mean, ...
    metrics_pat.PAP_max,   metrics_pat.PAP_min,   metrics_pat.PAP_mean);
fprintf('============================================================\n\n');
