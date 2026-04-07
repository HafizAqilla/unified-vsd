%% demo_sim_comparison.m
% =========================================================================
% Demonstration: apply_scaling → integrate_system → plot_simulation_comparison
%
% PURPOSE:
%   Runs the full ODE pipeline for (1) an adult reference and (2) a
%   paediatric VSD patient, then renders the 4-figure comparison:
%
%     Figure 1 — LV and RV PV loops
%     Figure 2 — Pressure waveforms (LV, RV, SAR, PAR, LA, RA)
%     Figure 3 — Flow waveforms (AV, MV, TV, PVv, SVEN, PVEN, VSD shunt)
%     Figure 4 — Clinical metrics panel (Qp/Qs, EF, SV, SVR/PVR, summary)
%
% USAGE:
%   >> cd <project_root>          % e.g. C:\Users\asus\Documents\MATLAB\VSD
%   >> demo_sim_comparison
%
% To change the patient, edit Section 2 below.
% To change the scenario ('pre_surgery' | 'post_surgery'), edit Section 3.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-31
% VERSION:  1.0
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
%  2. PAEDIATRIC VSD PATIENT  ← edit this block for your patient
% =========================================================================

% --- Demographics
patient_pt.age_years = 1.6 / 12;   % [yr]   1.6 months
patient_pt.weight_kg = 3.7;         % [kg]
patient_pt.height_cm = 54;          % [cm]
patient_pt.sex       = 'M';

% --- Clinical measurements (fill what you have; leave NaN if unavailable)
clinical = patient_template();
clinical.common.age_years = patient_pt.age_years;
clinical.common.weight_kg = patient_pt.weight_kg;
clinical.common.height_cm = patient_pt.height_cm;
clinical.common.sex       = patient_pt.sex;
% clinical.common.HR      = NaN;    % [bpm] — leave NaN → scaled by BSA^-0.33

% VSD shunt (pre-surgery)
clinical.pre_surgery.VSD_diameter_mm   = 6.0;   % [mm]     defect diameter
clinical.pre_surgery.VSD_gradient_mmHg = NaN;   % [mmHg]   Doppler gradient

% Pulmonary haemodynamics
clinical.pre_surgery.PAP_mean_mmHg = NaN;   % [mmHg]
clinical.pre_surgery.PVR_WU        = NaN;   % [WU]

% Systemic haemodynamics
clinical.pre_surgery.SVR_WU        = NaN;   % [WU]

%% =========================================================================
%  3. SCENARIO: 'pre_surgery' | 'post_surgery'
% =========================================================================
scenario = 'pre_surgery';

%% =========================================================================
%  4. SCALE + SIMULATE PATIENT
% =========================================================================
fprintf('[demo] Scaling patient parameters...\n');
params_ref = default_parameters();          % reference for patient allometric scaling
params_pat = apply_scaling(params_ref, patient_pt);
params_pat = params_from_clinical(params_pat, clinical, scenario);

fprintf('[demo] Integrating patient ODE...\n');
sim_pat = integrate_system(params_pat);
fprintf('[demo] Patient simulation complete.\n\n');

%% =========================================================================
%  5. BUILD COMPARISON LABEL
% =========================================================================
label = sprintf('Patient (%.0f mo, %.1f kg) – %s', ...
    patient_pt.age_years * 12, patient_pt.weight_kg, scenario);

%% =========================================================================
%  6. PLOT
% =========================================================================
fprintf('[demo] Rendering 4-figure comparison...\n');
plot_simulation_comparison(sim_adult, params_adult, ...
                            sim_pat,   params_pat, ...
                            label);
fprintf('[demo] Done. 4 figures are now open.\n');
