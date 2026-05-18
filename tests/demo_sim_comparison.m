%% 
%% demo_sim_comparison.m
% =========================================================================
% Demonstration: Comparing Two Allometric Scaling Methods
%   (1) apply_scaling   — Zhang 2019, weight-based power law
%   (2) apply_scaling_L — Lundquist 2025, BSA-based power law
%
% PURPOSE:
%   Runs the full ODE pipeline for three configurations:
%     [A] Adult reference      — Valenti Table 3.3, no scaling
%     [Z] Healthy paediatric   — Zhang 2019  (apply_scaling,   weight ratio)
%     [L] Healthy paediatric   — Lundquist25 (apply_scaling_L, BSA ratio)
%
%   Then renders:
%     Figures 1-4   — Single 4-figure set with all 3 traces overlaid
%                     (Adult=blue solid | Zhang=red dashed | Lundquist=green dotted)
%     Console       — 6-column clinical metrics table:
%                     Norm:Adult | Adult(ref) | Norm:Pediatric | Zhang | Lundquist
%
% USAGE:
%   >> cd <project_root>          % e.g. C:\Users\asus\Documents\MATLAB\VSD
%   >> demo_sim_comparison
%
% PATIENT DATA SOURCES:
%   Adult (reference)       — Valenti Table 3.3 (raw, no scaling)
%   Healthy paediatric [Z]  — Zhang 2019 allometric scaling (apply_scaling)
%   Healthy paediatric [L]  — Lundquist 2025 BSA scaling   (apply_scaling_L)
%
% NOTE ON main_run.m:
%   This demo script does NOT alter main_run.m or any established pipeline
%   functions.  It calls apply_scaling and apply_scaling_L independently.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  7.0  (6-column layout: Norm:Adult | Adult | Norm:Ped | Zhang | Lundquist)
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
%  reference model exactly.
%
%  Fractional timing fields (Tc_LV_frac, Tr_LV_frac, ...) are converted to
%  absolute seconds at HR = 75 bpm; normally done by apply_scaling Section D.
% =========================================================================
params_adult = default_parameters();

% Absolute timing from fractional values  (Valenti Table 2.1)
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
%  2. PAEDIATRIC PATIENT DEMOGRAPHICS
%
%  One patient, two scaling methods applied independently.
%  Edit demographics here to change the comparison subject.
% =========================================================================
patient_pt.age_years = 3 + (2 / 12);    % [yr]   3 years 2 months
patient_pt.weight_kg = 14;              % [kg]
patient_pt.height_cm = 98;              % [cm]
patient_pt.sex       = 'M';

fprintf('[demo] Pediatric: %.1f yr | %.1f kg | %.0f cm | %s\n\n', ...
    patient_pt.age_years, patient_pt.weight_kg, ...
    patient_pt.height_cm, patient_pt.sex);

%% =========================================================================
%  3A. ZHANG 2019 SCALING  (apply_scaling — weight-based power law)
%
%  Scale factor: w = W_patient / W_ref  (W_ref = 70 kg)
%  Exponents from Zhang et al. (2019):
%    HR    ~ w^(-0.30)    V0   ~ w^(+0.80)
%    E_LV  ~ w^(-0.50)    R    ~ w^(-0.475) systemic
%    E_RV  ~ w^(-0.75)    R    ~ w^(-0.70)  pulmonary
%    C     ~ w^(+1.00)    R_valve ~ w^(-0.50)
% =========================================================================
fprintf('[demo] ---- Scaling Method [Z]: Zhang 2019 (weight-based) ----\n');
params_ref_Z  = default_parameters();
params_zhang  = apply_scaling(params_ref_Z, patient_pt);

fprintf('[demo] Integrating Zhang-scaled ODE...\n');
sim_zhang = integrate_system(params_zhang);
fprintf('[demo] Zhang-scaled simulation complete.\n\n');

%% =========================================================================
%  3B. LUNDQUIST 2025 SCALING  (apply_scaling_L — BSA-based power law)
%
%  Scale factor: s = BSA_patient / BSA_ref  (BSA_ref = 1.73 m², Valenti ref)
%  BSA computed by Mosteller formula unless patient.BSA is supplied.
%  Exponents from Lundquist et al. (2025), Tables 2 & 3:
%    HR    ~ s^(-0.33)    V0   re-scaled by BV ratio (not s)
%    E_LV  ~ s^(-1.0)     R    ~ s^(-1.0)  (all vascular)
%    E_RV  ~ s^(-1.5)     C    ~ s^(+1.0)
%    L     ~ s^(-1.0)     R_valve_open ~ s^(-1.0)
% =========================================================================
fprintf('[demo] ---- Scaling Method [L]: Lundquist 2025 (BSA-based) ----\n');
params_ref_L    = default_parameters();
params_lundquist = apply_scaling_L(params_ref_L, patient_pt);

fprintf('[demo] Integrating Lundquist-scaled ODE...\n');
sim_lundquist = integrate_system(params_lundquist);
fprintf('[demo] Lundquist-scaled simulation complete.\n\n');

%% =========================================================================
%  4. BUILD COMPARISON LABELS
% =========================================================================
label_zhang = sprintf('Zhang-2019 (%.0f mo, %.1f kg)', ...
    patient_pt.age_years * 12, patient_pt.weight_kg);
label_lundquist = sprintf('Lundquist-2025 (%.0f mo, %.1f kg)', ...
    patient_pt.age_years * 12, patient_pt.weight_kg);

%% =========================================================================
%  5. 3-WAY OVERLAY — All scaling methods in the same 4 figures
%
%  plot_comparison_3way overlays all three traces:
%    Adult ref   — blue  solid
%    Zhang 2019  — red   dashed
%    Lundquist   — green dotted
% =========================================================================
fprintf('[demo] Rendering 3-way comparison figures (4 figures, 3 traces each)...\n');
plot_comparison_3way(sim_adult, params_adult, ...
                     sim_zhang, params_zhang, ...
                     sim_lundquist, params_lundquist, ...
                     label_zhang, label_lundquist);
fprintf('[demo] 3-way figures done.\n\n');

%% =========================================================================
%  7. CLINICAL METRICS — console table
%
%  Column order:
%    Norm: Adult  |  Adult (ref)  |  Norm: Pediatric  |  Zhang 2019  |  Lundquist 2025
%
%  The healthy reference for each simulation sits immediately to its left,
%  making it easy to judge at a glance whether each result is in range.
% =========================================================================
fprintf('[demo] Computing clinical indices for all three cases...\n');
metrics_adult = compute_clinical_indices(sim_adult,     params_adult);
metrics_zhang = compute_clinical_indices(sim_zhang,     params_zhang);
metrics_lund  = compute_clinical_indices(sim_lundquist, params_lundquist);

% Cardiac output: SV [mL] x HR [bpm] / 1000 = CO [L/min]
CO_adult = metrics_adult.LVSV * params_adult.HR      / 1000;   % [L/min]
CO_zhang = metrics_zhang.LVSV * params_zhang.HR      / 1000;   % [L/min]
CO_lund  = metrics_lund.LVSV  * params_lundquist.HR  / 1000;   % [L/min]

% ---- Load clinical reference ranges (both populations) -----------------
ref_ad  = clinical_reference_ranges('adult');       % healthy adult norms
ref_ped = clinical_reference_ranges('pediatric');   % healthy paed  norms

% Helper: format [lo hi] range as 'lo-hi' or '<hi' or '--'
fmt_r = @(r) demo_fmt_ref(r);

% ---- Print header -------------------------------------------------------
% Column order: Metric | Norm:Adult | Adult(ref) | Norm:Pediatric | Zhang | Lundquist
SEP = repmat('=', 1, 125);
sep = repmat('-', 1, 125);
fprintf('\n%s\n', SEP);
fprintf('  CLINICAL METRICS - 3-WAY SCALING COMPARISON\n');
fprintf('  Pediatric patient: %.1f yr | %.1f kg | %.0f cm | %s\n', ...
    patient_pt.age_years, patient_pt.weight_kg, patient_pt.height_cm, patient_pt.sex);
fprintf('  Column order: Norm:Adult | Adult(ref) | Norm:Pediatric | Zhang 2019 | Lundquist 2025\n');
fprintf('%s\n', SEP);
fprintf('  %-24s  %-16s  %-13s  %-16s  %-14s  %s\n', ...
    'Metric', 'Norm: Adult', 'Adult (ref)', 'Norm: Pediatric', 'Zhang 2019', 'Lundquist 2025');
fprintf('  %-24s  %-16s  %-13s  %-16s  %-14s  %s\n', ...
    '', '(healthy adult)', '(no scaling)', '(healthy paed)', '(w-based)', '(BSA-based)');
fprintf('%s\n', sep);

% ---- Heart rate ---------------------------------------------------------
fprintf('  %-24s  %-16s  %6.1f bpm     %-16s  %6.1f bpm     %6.1f bpm\n', ...
    'Heart Rate [bpm]', ...
    fmt_r(ref_ad.HR),  params_adult.HR, ...
    fmt_r(ref_ped.HR), params_zhang.HR, params_lundquist.HR);

% ---- Cardiac output -----------------------------------------------------
fprintf('  %-24s  %-16s  %6.2f L/min   %-16s  %6.2f L/min   %6.2f L/min\n', ...
    'Cardiac Output [L/min]', ...
    fmt_r(ref_ad.CO),  CO_adult, ...
    fmt_r(ref_ped.CO), CO_zhang, CO_lund);

% ---- Stroke volumes -----------------------------------------------------
fprintf('  %-24s  %-16s  %6.1f mL      %-16s  %6.1f mL      %6.1f mL\n', ...
    'LV Stroke Vol [mL]', ...
    fmt_r(ref_ad.SV),  metrics_adult.LVSV, ...
    fmt_r(ref_ped.SV), metrics_zhang.LVSV, metrics_lund.LVSV);
fprintf('  %-24s  %-16s  %6.1f mL      %-16s  %6.1f mL      %6.1f mL\n', ...
    'RV Stroke Vol [mL]', ...
    '--', metrics_adult.RVSV, '--', metrics_zhang.RVSV, metrics_lund.RVSV);

% ---- Ejection fractions -------------------------------------------------
fprintf('  %-24s  %-16s  %6.1f %%       %-16s  %6.1f %%       %6.1f %%\n', ...
    'LV Ejection Frac [%]', ...
    fmt_r(ref_ad.LVEF),  metrics_adult.LVEF*100, ...
    fmt_r(ref_ped.LVEF), metrics_zhang.LVEF*100, metrics_lund.LVEF*100);
fprintf('  %-24s  %-16s  %6.1f %%       %-16s  %6.1f %%       %6.1f %%\n', ...
    'RV Ejection Frac [%]', ...
    fmt_r(ref_ad.RVEF),  metrics_adult.RVEF*100, ...
    fmt_r(ref_ped.RVEF), metrics_zhang.RVEF*100, metrics_lund.RVEF*100);

% ---- Volumes -------------------------------------------------------------
fprintf('  %-24s  %-16s  %6.1f mL      %-16s  %6.1f mL      %6.1f mL\n', ...
    'LVEDV [mL]', ...
    fmt_r(ref_ad.LVEDV),  metrics_adult.LVEDV, ...
    fmt_r(ref_ped.LVEDV), metrics_zhang.LVEDV, metrics_lund.LVEDV);
fprintf('  %-24s  %-16s  %6.1f mL      %-16s  %6.1f mL      %6.1f mL\n', ...
    'LVESV [mL]', ...
    fmt_r(ref_ad.LVESV),  metrics_adult.LVESV, ...
    fmt_r(ref_ped.LVESV), metrics_zhang.LVESV, metrics_lund.LVESV);
fprintf('  %-24s  %-16s  %6.1f mL      %-16s  %6.1f mL      %6.1f mL\n', ...
    'RVEDV [mL]', ...
    fmt_r(ref_ad.RVEDV),  metrics_adult.RVEDV, ...
    fmt_r(ref_ped.RVEDV), metrics_zhang.RVEDV, metrics_lund.RVEDV);
fprintf('  %-24s  %-16s  %6.1f mL      %-16s  %6.1f mL      %6.1f mL\n', ...
    'RVESV [mL]', ...
    fmt_r(ref_ad.RVESV),  metrics_adult.RVESV, ...
    fmt_r(ref_ped.RVESV), metrics_zhang.RVESV, metrics_lund.RVESV);

% ---- Resistances (Wood units) -------------------------------------------
fprintf('  %-24s  %-16s  %6.2f WU      %-16s  %6.2f WU      %6.2f WU\n', ...
    'SVR [WU]', ...
    fmt_r(ref_ad.SVR),  metrics_adult.SVR, ...
    fmt_r(ref_ped.SVR), metrics_zhang.SVR, metrics_lund.SVR);
fprintf('  %-24s  %-16s  %6.2f WU      %-16s  %6.2f WU      %6.2f WU\n', ...
    'PVR [WU]', ...
    fmt_r(ref_ad.PVR),  metrics_adult.PVR, ...
    fmt_r(ref_ped.PVR), metrics_zhang.PVR, metrics_lund.PVR);

% ---- Systemic pressures -------------------------------------------------
fprintf('  %-24s  %-16s  %3.0f/%3.0f/%3.0f mmHg  %-16s  %3.0f/%3.0f/%3.0f mmHg  %3.0f/%3.0f/%3.0f mmHg\n', ...
    'SA SBP/DBP/MAP [mmHg]', ...
    [fmt_r(ref_ad.SBP) '/' fmt_r(ref_ad.DBP) '/' fmt_r(ref_ad.MAP)], ...
    metrics_adult.SAP_max, metrics_adult.SAP_min, metrics_adult.SAP_mean, ...
    [fmt_r(ref_ped.SBP) '/' fmt_r(ref_ped.DBP) '/' fmt_r(ref_ped.MAP)], ...
    metrics_zhang.SAP_max, metrics_zhang.SAP_min, metrics_zhang.SAP_mean, ...
    metrics_lund.SAP_max,  metrics_lund.SAP_min,  metrics_lund.SAP_mean);

% ---- Pulmonary pressures ------------------------------------------------
fprintf('  %-24s  %-16s  %3.0f/%3.0f/%3.0f mmHg  %-16s  %3.0f/%3.0f/%3.0f mmHg  %3.0f/%3.0f/%3.0f mmHg\n', ...
    'PA sys/dia/mean [mmHg]', ...
    ['s' fmt_r(ref_ad.PAP_max) ' d' fmt_r(ref_ad.PAP_min) ' m' fmt_r(ref_ad.PAP_mean)], ...
    metrics_adult.PAP_max, metrics_adult.PAP_min, metrics_adult.PAP_mean, ...
    ['s' fmt_r(ref_ped.PAP_max) ' d' fmt_r(ref_ped.PAP_min) ' m' fmt_r(ref_ped.PAP_mean)], ...
    metrics_zhang.PAP_max, metrics_zhang.PAP_min, metrics_zhang.PAP_mean, ...
    metrics_lund.PAP_max,  metrics_lund.PAP_min,  metrics_lund.PAP_mean);

% ---- LV pressures -------------------------------------------------------
fprintf('  %-24s  %-16s  %6.0f mmHg    %-16s  %6.0f mmHg    %6.0f mmHg\n', ...
    'LVESP [mmHg]', ...
    fmt_r(ref_ad.LVESP),  metrics_adult.LVP_max, ...
    fmt_r(ref_ped.LVESP), metrics_zhang.LVP_max, metrics_lund.LVP_max);
fprintf('  %-24s  %-16s  %6.0f mmHg    %-16s  %6.0f mmHg    %6.0f mmHg\n', ...
    'LVEDP [mmHg]', ...
    fmt_r(ref_ad.LVEDP),  metrics_adult.LVEDP, ...
    fmt_r(ref_ped.LVEDP), metrics_zhang.LVEDP, metrics_lund.LVEDP);

% ---- RV pressures -------------------------------------------------------
fprintf('  %-24s  %-16s  %6.0f mmHg    %-16s  %6.0f mmHg    %6.0f mmHg\n', ...
    'RVESP [mmHg]', ...
    fmt_r(ref_ad.RVESP),  metrics_adult.RVP_max, ...
    fmt_r(ref_ped.RVESP), metrics_zhang.RVP_max, metrics_lund.RVP_max);
fprintf('  %-24s  %-16s  %6.0f mmHg    %-16s  %6.0f mmHg    %6.0f mmHg\n', ...
    'RVEDP [mmHg]', ...
    fmt_r(ref_ad.RVEDP),  metrics_adult.RVEDP, ...
    fmt_r(ref_ped.RVEDP), metrics_zhang.RVEDP, metrics_lund.RVEDP);

% ---- Mean atrial pressures ----------------------------------------------
fprintf('  %-24s  %-16s  %6.1f mmHg    %-16s  %6.1f mmHg    %6.1f mmHg\n', ...
    'LAP mean [mmHg]', ...
    fmt_r(ref_ad.LAP_mean),  metrics_adult.LAP_mean, ...
    fmt_r(ref_ped.LAP_mean), metrics_zhang.LAP_mean, metrics_lund.LAP_mean);
fprintf('  %-24s  %-16s  %6.1f mmHg    %-16s  %6.1f mmHg    %6.1f mmHg\n', ...
    'RAP mean [mmHg]', ...
    fmt_r(ref_ad.RAP_mean),  metrics_adult.RAP_mean, ...
    fmt_r(ref_ped.RAP_mean), metrics_zhang.RAP_mean, metrics_lund.RAP_mean);

% ---- PWP mean -----------------------------------------------------------
fprintf('  %-24s  %-16s  %6.1f mmHg    %-16s  %6.1f mmHg    %6.1f mmHg\n', ...
    'PWP mean [mmHg]', ...
    '--', metrics_adult.PWP_mean, '--', metrics_zhang.PWP_mean, metrics_lund.PWP_mean);

% ---- Qp/Qs --------------------------------------------------------------
fprintf('  %-24s  %-16s  %6.2f          %-16s  %6.2f          %6.2f\n', ...
    'Qp/Qs [-]', ...
    '1.0 (normal)', metrics_adult.QpQs, '1.0 (normal)', metrics_zhang.QpQs, metrics_lund.QpQs);

% ---- Scaling parameters (informational) ---------------------------------
fprintf('%s\n', sep);
fprintf('  %-24s\n', 'Scaling parameters used:');
fprintf('  %-24s  %6s           %6.4f w        %6.4f s\n', ...
    '  Scale factor', '-', params_zhang.scaling.w, params_lundquist.scaling.s);
fprintf('  %-24s  %6s           %6.1f kg       %6.3f m2\n', ...
    '  Patient value', '-', params_zhang.scaling.W_patient, params_lundquist.scaling.BSA_patient);
fprintf('  %-24s  %6s           %6.1f kg       %6.3f m2\n', ...
    '  Reference value', '-', params_zhang.scaling.W_ref, params_lundquist.scaling.BSA_ref);
fprintf('  %-24s  %6.0f mL      %6.0f mL      %6.0f mL\n', ...
    '  BV target', ...
    params_adult.ic.V(params_adult.idx.V_SVEN) + ...     % rough adult BV estimate
    params_adult.ic.V(params_adult.idx.V_SAR)  + ...
    params_adult.ic.V(params_adult.idx.V_SC)   + ...
    params_adult.ic.V(params_adult.idx.V_PAR)  + ...
    params_adult.ic.V(params_adult.idx.V_PVEN) + ...
    params_adult.ic.V(params_adult.idx.V_LV)   + ...
    params_adult.ic.V(params_adult.idx.V_RV)   + ...
    params_adult.ic.V(params_adult.idx.V_LA)   + ...
    params_adult.ic.V(params_adult.idx.V_RA),  ...
    params_zhang.scaling.BV_patient, params_lundquist.scaling.BV_patient);
fprintf('%s\n\n', SEP);
fprintf('[demo] Done. 4 figures open with all 3 scaling methods overlaid.\n');


%% =========================================================================
%  8. CSV EXPORT — write clinical metrics table to results/tables/
%
%  Output columns (one row per metric):
%    Metric, Unit, Norm_Adult_lo, Norm_Adult_hi,
%    Adult_ref, Zhang_2019, Lundquist_2025,
%    Norm_Ped_lo, Norm_Ped_hi
% =========================================================================

% ---- Resolve output directory -------------------------------------------
out_dir = fullfile(root, 'results', 'tables');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% ---- Assemble data rows -------------------------------------------------
% Each row: {metric_label, unit, norm_adult, adult_val, norm_ped, zhang_val, lund_val}
rows = { ...
    'Heart Rate',         'bpm',   fmt_r(ref_ad.HR),       params_adult.HR,          fmt_r(ref_ped.HR),       params_zhang.HR,          params_lundquist.HR; ...
    'Cardiac Output',     'L/min', fmt_r(ref_ad.CO),       CO_adult,                 fmt_r(ref_ped.CO),       CO_zhang,                 CO_lund; ...
    'LV Stroke Volume',   'mL',    fmt_r(ref_ad.SV),       metrics_adult.LVSV,       fmt_r(ref_ped.SV),       metrics_zhang.LVSV,       metrics_lund.LVSV; ...
    'RV Stroke Volume',   'mL',    '--',                   metrics_adult.RVSV,       '--',                    metrics_zhang.RVSV,       metrics_lund.RVSV; ...
    'LVEF',               '%',     fmt_r(ref_ad.LVEF),     metrics_adult.LVEF*100,   fmt_r(ref_ped.LVEF),     metrics_zhang.LVEF*100,   metrics_lund.LVEF*100; ...
    'RVEF',               '%',     fmt_r(ref_ad.RVEF),     metrics_adult.RVEF*100,   fmt_r(ref_ped.RVEF),     metrics_zhang.RVEF*100,   metrics_lund.RVEF*100; ...
    'LVEDV',              'mL',    fmt_r(ref_ad.LVEDV),    metrics_adult.LVEDV,      fmt_r(ref_ped.LVEDV),    metrics_zhang.LVEDV,      metrics_lund.LVEDV; ...
    'LVESV',              'mL',    fmt_r(ref_ad.LVESV),    metrics_adult.LVESV,      fmt_r(ref_ped.LVESV),    metrics_zhang.LVESV,      metrics_lund.LVESV; ...
    'RVEDV',              'mL',    fmt_r(ref_ad.RVEDV),    metrics_adult.RVEDV,      fmt_r(ref_ped.RVEDV),    metrics_zhang.RVEDV,      metrics_lund.RVEDV; ...
    'RVESV',              'mL',    fmt_r(ref_ad.RVESV),    metrics_adult.RVESV,      fmt_r(ref_ped.RVESV),    metrics_zhang.RVESV,      metrics_lund.RVESV; ...
    'SVR',                'WU',    fmt_r(ref_ad.SVR),      metrics_adult.SVR,        fmt_r(ref_ped.SVR),      metrics_zhang.SVR,        metrics_lund.SVR; ...
    'PVR',                'WU',    fmt_r(ref_ad.PVR),      metrics_adult.PVR,        fmt_r(ref_ped.PVR),      metrics_zhang.PVR,        metrics_lund.PVR; ...
    'SBP',                'mmHg',  fmt_r(ref_ad.SBP),      metrics_adult.SAP_max,    fmt_r(ref_ped.SBP),      metrics_zhang.SAP_max,    metrics_lund.SAP_max; ...
    'DBP',                'mmHg',  fmt_r(ref_ad.DBP),      metrics_adult.SAP_min,    fmt_r(ref_ped.DBP),      metrics_zhang.SAP_min,    metrics_lund.SAP_min; ...
    'MAP',                'mmHg',  fmt_r(ref_ad.MAP),      metrics_adult.SAP_mean,   fmt_r(ref_ped.MAP),      metrics_zhang.SAP_mean,   metrics_lund.SAP_mean; ...
    'PAP_sys',            'mmHg',  fmt_r(ref_ad.PAP_max),  metrics_adult.PAP_max,    fmt_r(ref_ped.PAP_max),  metrics_zhang.PAP_max,    metrics_lund.PAP_max; ...
    'PAP_dia',            'mmHg',  fmt_r(ref_ad.PAP_min),  metrics_adult.PAP_min,    fmt_r(ref_ped.PAP_min),  metrics_zhang.PAP_min,    metrics_lund.PAP_min; ...
    'PAP_mean',           'mmHg',  fmt_r(ref_ad.PAP_mean), metrics_adult.PAP_mean,   fmt_r(ref_ped.PAP_mean), metrics_zhang.PAP_mean,   metrics_lund.PAP_mean; ...
    'LVESP',              'mmHg',  fmt_r(ref_ad.LVESP),    metrics_adult.LVP_max,    fmt_r(ref_ped.LVESP),    metrics_zhang.LVP_max,    metrics_lund.LVP_max; ...
    'LVEDP',              'mmHg',  fmt_r(ref_ad.LVEDP),    metrics_adult.LVEDP,    fmt_r(ref_ped.LVEDP),    metrics_zhang.LVEDP,    metrics_lund.LVEDP; ...
    'RVESP',              'mmHg',  fmt_r(ref_ad.RVESP),    metrics_adult.RVP_max,    fmt_r(ref_ped.RVESP),    metrics_zhang.RVP_max,    metrics_lund.RVP_max; ...
    'RVEDP',              'mmHg',  fmt_r(ref_ad.RVEDP),    metrics_adult.RVEDP,    fmt_r(ref_ped.RVEDP),    metrics_zhang.RVEDP,    metrics_lund.RVEDP; ...
    'LAP_mean',           'mmHg',  fmt_r(ref_ad.LAP_mean), metrics_adult.LAP_mean,   fmt_r(ref_ped.LAP_mean), metrics_zhang.LAP_mean,   metrics_lund.LAP_mean; ...
    'RAP_mean',           'mmHg',  fmt_r(ref_ad.RAP_mean), metrics_adult.RAP_mean,   fmt_r(ref_ped.RAP_mean), metrics_zhang.RAP_mean,   metrics_lund.RAP_mean; ...
    'PWP_mean',           'mmHg',  '--',                   metrics_adult.PWP_mean,   '--',                    metrics_zhang.PWP_mean,   metrics_lund.PWP_mean; ...
    'Qp_Qs',              '-',     '1.0 (normal)',         metrics_adult.QpQs,       '1.0 (normal)',          metrics_zhang.QpQs,       metrics_lund.QpQs };

% ---- Write XLSX ---------------------------------------------------------
xlsx_fname = fullfile(out_dir, ...
    sprintf('demo_metrics_%s.xlsx', datestr(now, 'yyyymmdd_HHMMSS')));

% Convert the cell array of rows into a MATLAB table
colNames = {'Metric', 'Unit', 'Norm_Adult', 'Adult_ref', 'Norm_Pediatric', 'Zhang_2019', 'Lundquist_2025'};
T = cell2table(rows, 'VariableNames', colNames);

try
    writetable(T, xlsx_fname, 'Sheet', 'Metrics');
    fprintf('[demo] XLSX saved: %s\n', xlsx_fname);
catch ME
    warning('[demo] Could not export XLSX: %s', ME.message);
end


%  LOCAL HELPER — reference range formatter
% =========================================================================
function s = demo_fmt_ref(r)
% DEMO_FMT_REF — format a [lo hi] reference pair as a compact string.
%   [NaN NaN] -> '--'
%   [0   hi]  -> '<hi'
%   [lo  hi]  -> 'lo-hi'
if all(isnan(r))
    s = '--';
elseif r(1) == 0
    s = sprintf('<%g', r(2));
else
    s = sprintf('%g-%g', r(1), r(2));
end
end
