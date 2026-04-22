%% demo_sim_comparison.m
% =========================================================================
% Demonstration: apply_scaling → integrate_system → plot_simulation_comparison
%
% PURPOSE:
%   Runs the full ODE pipeline for (1) an adult reference, (2) a healthy
%   paediatric patient, and (3) a paediatric VSD patient (pre-surgery),
%   then renders the 4-figure comparison:
%
%     Figure 1 — LV and RV PV loops
%     Figure 2 — Pressure waveforms (LV, RV, SAR, PAR, LA, RA)
%     Figure 3 — Flow waveforms (AV, MV, TV, PVv, SVEN, PVEN, VSD)
%     Figure 4 — Clinical metrics panel (EF, SV, SVR/PVR, summary)
%
%   The base comparison (adult vs healthy paediatric) uses only allometric
%   scaling.  The VSD simulation additionally calls params_from_clinical
%   to configure the shunt resistance and clinical haemodynamic targets
%   (HR, SVR, PVR, echo volumes, and IC override).
%
% USAGE:
%   >> cd <project_root>          % e.g. C:\Users\asus\Documents\MATLAB\VSD
%   >> demo_sim_comparison
%
% PATIENT DATA SOURCES:
%   Adult (reference)       — Valenti Table 3.3 (raw, no scaling)
%   Healthy paediatric      — Zhang 2019 allometric scaling only
%   Paediatric VSD patient  — run_patient_case.m clinical data
%     Age  : 3 yr 2 mo   Weight : 13.4 kg   Height : 95 cm   Sex: F
%     VSD  : 5.0 mm diameter   Gradient : 94 mmHg (peak echo Doppler)
%     PAP  : 20/10/15 mmHg     SAP : 119/83/95 mmHg
%     QpQs : 1.194   LVEDV : 32 mL   LVESV : 23.6 mL
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-20
% VERSION:  3.0  (3-way comparison: adult / healthy-paediatric / VSD-paediatric)
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

% --- Demographics (3 years 2 months)
patient_pt.age_years = 3 + (2 / 12);    % [yr]   3 years 2 months
patient_pt.weight_kg = 13.4;       % [kg]
patient_pt.height_cm = 95;         % [cm]
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
%  5. PLOT (adult vs healthy paediatric — base 2-way comparison)
% =========================================================================
fprintf('[demo] Rendering 4-figure comparison (adult vs healthy paediatric)...\n');
plot_simulation_comparison(sim_adult, params_adult, ...
                            sim_pat,   params_pat, ...
                            label);
fprintf('[demo] Base 2-way comparison done.\n\n');

%% =========================================================================
%  6. PAEDIATRIC VSD PATIENT  (pre-surgery)
%
%  Clinical data from run_patient_case.m — 3 yr 2 mo, 13.4 kg, Female.
%  Pipeline:
%    apply_scaling        — allometric scaling (Zhang 2019)
%    params_from_clinical — HR, SVR, PVR, R_VSD from gradient + diameter,
%                           echo volume IC override (override_IC = true)
%
%  VSD geometry (Gorlin fallback — only diameter provided without
%  simultaneous Q_shunt flow; peak Doppler gradient 94 mmHg used to
%  compute Option A mean gradient directly):
%    D_VSD = 5.0 mm   VSD_gradient = 94 mmHg (peak echo Doppler)
%    Q_shunt_Lmin  = (QpQs - 1) × Qs
%      where Qs ≈ (SAP_mean - RAP) / SVR  and SVR from clinical data
%
%  Source: run_patient_case.m (this repository) & clinical catheter report.
% =========================================================================
fprintf('[demo] Configuring paediatric VSD patient (pre-surgery)...\n');

% ---- Demographics (same age/weight/height as healthy-paediatric; sex = F)
patient_vsd.age_years = 3 + (2 / 12);    % [yr]   3 years 2 months
patient_vsd.weight_kg = 13.4;            % [kg]
patient_vsd.height_cm = 95;              % [cm]
patient_vsd.sex       = 'F';            % Female — per clinical record

% ---- Build a minimal clinical struct using patient_template scaffold
clinical_vsd = patient_template();

clinical_vsd.common.age_years = patient_vsd.age_years;   % [yr]
clinical_vsd.common.weight_kg = patient_vsd.weight_kg;   % [kg]
clinical_vsd.common.height_cm = patient_vsd.height_cm;   % [cm]
clinical_vsd.common.sex       = patient_vsd.sex;
clinical_vsd.common.HR        = 110;    % [bpm]  slightly elevated for 3-yr VSD child
%   Estimated: healthy toddler HR 100–110 bpm (Fleming S et al. 2011 PLoS Med;
%   VSD children show mild sympathetic activation even with small-to-moderate shunt).

% ---- Pre-surgery clinical measurements (catheter + echo) ----------------
% VSD shunt geometry
clinical_vsd.pre_surgery.VSD_diameter_mm   = 5.0;    % [mm]   clinical echo
clinical_vsd.pre_surgery.VSD_gradient_mmHg = 94;     % [mmHg] peak systolic VSD gradient (echo Doppler)

% Haemodynamic targets
clinical_vsd.pre_surgery.PAP_sys_mmHg  = 20;         % [mmHg]  PAP systolic
clinical_vsd.pre_surgery.PAP_dia_mmHg  = 10;         % [mmHg]  PAP diastolic
clinical_vsd.pre_surgery.PAP_mean_mmHg = 15;         % [mmHg]  mean PAP

clinical_vsd.pre_surgery.SAP_sys_mmHg  = 119;        % [mmHg]  systolic SAP
clinical_vsd.pre_surgery.SAP_dia_mmHg  = 83;         % [mmHg]  diastolic SAP
clinical_vsd.pre_surgery.SAP_mean_mmHg = 95;         % [mmHg]  MAP

clinical_vsd.pre_surgery.RAP_mean_mmHg = 5;          % [mmHg]  right atrial pressure
clinical_vsd.pre_surgery.LAP_mean_mmHg = 9;          % [mmHg]  estimated LAP
%   LAP ≈ LVEDP + mitral inflow gradient (<2 mmHg for moderate shunt).
%   LVEDP assumed ~7 mmHg for a 3-yr child with QpQs 1.194 (mild overload).

clinical_vsd.pre_surgery.QpQs          = 1.194;      % [-]  Qp/Qs ratio from cath

% Echo volumes — activate IC override for correct chamber sizing
clinical_vsd.pre_surgery.LVEDV_mL      = 32;         % [mL]  LV end-diastolic volume
clinical_vsd.pre_surgery.LVESV_mL      = 23.6;       % [mL]  LV end-systolic volume
clinical_vsd.pre_surgery.RVEDV_mL      = 30.5;       % [mL]  RV end-diastolic volume
clinical_vsd.pre_surgery.RVESV_mL      = 12;         % [mL]  RV end-systolic volume
clinical_vsd.pre_surgery.LVEF          = (32 - 23.6) / 32;  % [-]  LVEF fraction

clinical_vsd.pre_surgery.LVEDP_mmHg    = 7;          % [mmHg] assumed clinical LVEDP
%   Q_shunt = (QpQs - 1) * Qs  assumes  Qs = (SAP_mean - RAP)/SVR
%   SVR not directly available; SVR_WU omitted — calibration will refine.

% Activate IC and elastance override (clinical echo volumes replace allometric ICs)
clinical_vsd.pre_surgery.override_IC   = true;       % [-]   use echo volumes for IC
clinical_vsd.pre_surgery.BV_total_mL   = 13.4 * 75; % [mL]  75 mL/kg (3-yr child;
%   Source: Lundquist 2025, Table 1 — BV 70-80 mL/kg at age 2-4 yr)

% Q_shunt for Option A R_VSD calculation (Ohm's law instead of Gorlin alone)
%   Qs  ≈ LVSV × HR × 1e-3 = (32-23.6) × 110 / 1000 = 0.924 L/min
%   Qp  = QpQs × Qs         = 1.194 × 0.924         = 1.103 L/min
%   Q_shunt = Qp - Qs        = 1.103 - 0.924         = 0.179 L/min
LVSV_vsd   = clinical_vsd.pre_surgery.LVEDV_mL - clinical_vsd.pre_surgery.LVESV_mL; % [mL]
Qs_vsd     = LVSV_vsd * clinical_vsd.common.HR / 1000;                              % [L/min]
Qp_vsd     = clinical_vsd.pre_surgery.QpQs * Qs_vsd;                               % [L/min]
clinical_vsd.pre_surgery.Q_shunt_Lmin = Qp_vsd - Qs_vsd;   % [L/min]  L→R shunt
fprintf('[demo]   Estimated Q_shunt = %.3f L/min (Qp=%.3f, Qs=%.3f)\n', ...
        clinical_vsd.pre_surgery.Q_shunt_Lmin, Qp_vsd, Qs_vsd);

% ---- Scale → clinical override → simulate
params_ref_vsd = default_parameters();
params_vsd     = apply_scaling(params_ref_vsd, patient_vsd);
params_vsd     = params_from_clinical(params_vsd, clinical_vsd, 'pre_surgery');

fprintf('[demo] Integrating paediatric VSD ODE...\n');
sim_vsd = integrate_system(params_vsd);
fprintf('[demo] Paediatric VSD simulation complete.\n\n');

label_vsd = sprintf('Paediatric VSD pre-surg (%.0f mo, %.1f kg, F, VSD 5 mm)', ...
    patient_vsd.age_years * 12, patient_vsd.weight_kg);

%% =========================================================================
%  7. OVERLAY VSD PATIENT TRACES ON THE 4 OPEN FIGURES
%
%  The existing 4 figures (from plot_simulation_comparison) show adult (blue)
%  and healthy paediatric (red).  We now add the VSD patient as a third
%  trace in green (dashed).  Axes are addressed by figure handle; no new
%  figures are opened.
% =========================================================================
fprintf('[demo] Overlaying VSD patient traces on 4 figures...\n');

% Shared style token for VSD trace
clr_vsd  = [0.13 0.60 0.30];   % forest green
lw_vsd   = 1.8;
ls_vsd   = '-.';                % dash-dot to distinguish from patient (dashed)

% ---- Reconstruct P/Q signals for VSD patient ----------------------------
%   This replicates the logic of helper_reconstruct in plot_simulation_comparison.
%   It rebuilds pressures/flows from the ODE state vector and extracts the last
%   cardiac cycle for consistent comparison.
[P_vsd, Q_vsd, tc_vsd, Pc_vsd, Qc_vsd] = local_reconstruct(sim_vsd, params_vsd);
m_vsd = compute_clinical_indices(sim_vsd, params_vsd);

T_vsd  = 60 / params_vsd.HR;          % [s]
tn_vsd = tc_vsd / T_vsd;              % normalised cycle fraction

sidx_vsd = params_vsd.idx;

% =========================================================================
%  FIGURE 1 — PV LOOPS: add VSD patient to LV and RV subplots
% =========================================================================
fh1 = findobj('Type','figure','Name','PV Loops');
if ~isempty(fh1)
    figure(fh1(1));
    allAxes = findobj(fh1(1), 'Type', 'axes');
    % subplot(1,2,1) = LV;  subplot(1,2,2) = RV
    % MATLAB stacks subplots in reverse order in findobj results
    ax_lv = allAxes(end);   % subplot 1 — LV (first created = last in list)
    ax_rv = allAxes(end-1); % subplot 2 — RV

    mask_vsd = local_last_mask(sim_vsd, params_vsd);
    V_LV_vsd = sim_vsd.V(mask_vsd, sidx_vsd.V_LV);   % [mL]
    V_RV_vsd = sim_vsd.V(mask_vsd, sidx_vsd.V_RV);   % [mL]

    hold(ax_lv, 'on');
    plot(ax_lv, V_LV_vsd, Pc_vsd.LV, 'Color', clr_vsd, ...
         'LineWidth', lw_vsd, 'LineStyle', ls_vsd);
    text(ax_lv, 0.98, 0.83, sprintf('VSD  LVEF = %.0f%%', m_vsd.LVEF*100), ...
        'Units','normalized','HorizontalAlignment','right', ...
        'FontSize',7,'Color',clr_vsd);
    % Update legend
    ch = get(ax_lv,'Children');
    lg = findobj(ax_lv, 'Type','Legend');
    if ~isempty(lg)
        lg.String{end+1} = label_vsd;
    else
        legend(ax_lv, [label_vsd], 'Location','northwest', 'FontSize',8);
    end
    hold(ax_lv, 'off');

    hold(ax_rv, 'on');
    plot(ax_rv, V_RV_vsd, Pc_vsd.RV, 'Color', clr_vsd, ...
         'LineWidth', lw_vsd, 'LineStyle', ls_vsd);
    text(ax_rv, 0.98, 0.83, sprintf('VSD  RVEF = %.0f%%', m_vsd.RVEF*100), ...
        'Units','normalized','HorizontalAlignment','right', ...
        'FontSize',7,'Color',clr_vsd);
    lg2 = findobj(ax_rv, 'Type','Legend');
    if ~isempty(lg2)
        lg2.String{end+1} = label_vsd;
    else
        legend(ax_rv, [label_vsd], 'Location','northwest', 'FontSize',8);
    end
    hold(ax_rv, 'off');

    % Update sgtitle
    sgtitle(fh1(1), sprintf('PV Loops — Adult / Healthy Paediatric / %s', label_vsd), ...
        'FontSize', 11, 'FontWeight', 'bold');
end

% =========================================================================
%  FIGURE 2 — PRESSURE WAVEFORMS: overlay VSD patient on all 6 subplots
% =========================================================================
fh2 = findobj('Type','figure','Name','Pressure Waveforms');
if ~isempty(fh2)
    figure(fh2(1));
    allAx2 = findobj(fh2(1), 'Type','axes');
    allAx2 = flipud(allAx2);   % restore creation order (subplot 1 … 6)

    pressure_fields = {'LV','RV','SAR','PAR','LA','RA'};
    for k = 1:numel(pressure_fields)
        fld = pressure_fields{k};
        hold(allAx2(k), 'on');
        plot(allAx2(k), tn_vsd, Pc_vsd.(fld), ...
             'Color', clr_vsd, 'LineWidth', lw_vsd, 'LineStyle', ls_vsd);
        hold(allAx2(k), 'off');
    end
    % Update legend on first subplot only
    lg3 = findobj(allAx2(1), 'Type','Legend');
    if ~isempty(lg3)
        lg3.String{end+1} = label_vsd;
    end
    sgtitle(fh2(1), sprintf('Pressure Waveforms — Adult / Healthy Paediatric / %s', label_vsd), ...
        'FontSize', 11, 'FontWeight', 'bold');
end

% =========================================================================
%  FIGURE 3 — FLOW WAVEFORMS: overlay VSD patient on all subplots
% =========================================================================
fh3 = findobj('Type','figure','Name','Flow Waveforms');
if ~isempty(fh3)
    figure(fh3(1));
    allAx3 = findobj(fh3(1), 'Type','axes');
    allAx3 = flipud(allAx3);   % restore creation order

    flow_fields = {'AV','PVv','MV','TV','SVEN','PVEN','VSD'};
    % subplots 1–6 = flow panels; subplot 7 = VSD; subplot 8 = text box
    for k = 1:numel(flow_fields)
        if k > numel(allAx3), break; end
        fld = flow_fields{k};
        hold(allAx3(k), 'on');
        plot(allAx3(k), tn_vsd, Qc_vsd.(fld), ...
             'Color', clr_vsd, 'LineWidth', lw_vsd, 'LineStyle', ls_vsd);
        hold(allAx3(k), 'off');
    end
    % Re-annotate the flow summary text box (subplot 8) with 3-column data
    if numel(allAx3) >= 8
        ax_txt = allAx3(8);
        delete(findobj(ax_txt, 'Type', 'Text'));
        m_ad_f = compute_clinical_indices(sim_adult, params_adult);
        m_pt_f = compute_clinical_indices(sim_pat,   params_pat);
        txt = { ...
            sprintf('  %-14s  %6s %6s %6s', 'Metric', 'Adult', 'Hlthy', 'VSD'), ...
            sprintf('  %-14s  %6.2f %6.2f %6.2f', 'Qp/Qs', ...
                m_ad_f.QpQs, m_pt_f.QpQs, m_vsd.QpQs), ...
            sprintf('  %-14s  %6.1f %6.1f %6.1f', 'Q_VSD [mL/s]', ...
                m_ad_f.Q_shunt_mean_mLs, m_pt_f.Q_shunt_mean_mLs, m_vsd.Q_shunt_mean_mLs) };
        text(ax_txt, 0.05, 0.5, strjoin(txt, '\n'), ...
            'FontSize', 8.5, 'FontName', 'Courier New', ...
            'VerticalAlignment', 'middle', 'Units', 'normalized');
    end
    % Refresh legend on first flow subplot
    lg4 = findobj(allAx3(1), 'Type','Legend');
    if ~isempty(lg4)
        lg4.String{end+1} = label_vsd;
    end
    sgtitle(fh3(1), sprintf('Flow Waveforms — Adult / Healthy Paediatric / %s', label_vsd), ...
        'FontSize', 11, 'FontWeight', 'bold');
end

% =========================================================================
%  FIGURE 4 — CLINICAL METRICS: overlay VSD patient bar and update table
% =========================================================================
fh4 = findobj('Type','figure','Name','Clinical Metrics');
if ~isempty(fh4)
    figure(fh4(1));
    allAx4 = findobj(fh4(1), 'Type','axes');
    allAx4 = flipud(allAx4);   % subplots 1–6

    m_ad3 = compute_clinical_indices(sim_adult, params_adult);
    m_pt3 = compute_clinical_indices(sim_pat,   params_pat);

    clr_adult = [0.18 0.44 0.70];
    clr_pat   = [0.84 0.22 0.22];

    % ---- 4a: Qp/Qs — add 3rd bar
    if numel(allAx4) >= 1
        ax4a = allAx4(1);
        bh_a = bar(ax4a, [m_ad3.QpQs, m_pt3.QpQs, m_vsd.QpQs], 0.55, ...
                   'FaceColor', 'flat', ...
                   'CData', [clr_adult; clr_pat; clr_vsd], 'EdgeColor', 'none');
        yline(ax4a, 1.0, 'k--', 'LineWidth', 1.2, 'Label', 'Balanced (1.0)');
        xticks(ax4a, [1 2 3]);
        xticklabels(ax4a, {'Adult ref', 'Healthy paed.', 'VSD pre-surg'});
        xtickangle(ax4a, 20);
        ylabel(ax4a, 'Qp/Qs  [-]', 'FontSize', 9);
        title(ax4a, 'Qp/Qs Ratio', 'FontSize', 10, 'FontWeight', 'bold');
        text(ax4a, 3, m_vsd.QpQs*1.04, sprintf('%.2f', m_vsd.QpQs), ...
            'HorizontalAlignment','center','FontSize',9,'Color',clr_vsd);
        grid(ax4a,'on'); box(ax4a,'off');
    end

    % ---- 4b: Ejection fractions — add VSD column
    if numel(allAx4) >= 2
        ax4b = allAx4(2);
        cla(ax4b);
        EF_data = [m_ad3.LVEF, m_pt3.LVEF, m_vsd.LVEF; ...
                   m_ad3.RVEF, m_pt3.RVEF, m_vsd.RVEF] * 100;
        bh2 = bar(ax4b, EF_data, 0.65, 'grouped');
        bh2(1).FaceColor = clr_adult;  bh2(1).EdgeColor = 'none';
        bh2(2).FaceColor = clr_pat;    bh2(2).EdgeColor = 'none';
        bh2(3).FaceColor = clr_vsd;    bh2(3).EdgeColor = 'none';
        yline(ax4b, 55, 'k--', 'LineWidth', 0.8, 'Label', 'Normal LV EF (55%)');
        xticks(ax4b,[1 2]); xticklabels(ax4b,{'LV','RV'});
        ylabel(ax4b,'EF  [%]','FontSize',9);
        title(ax4b,'Ejection Fractions','FontSize',10,'FontWeight','bold');
        legend(ax4b,{'Adult ref', 'Healthy paed.', 'VSD pre-surg'}, ...
               'Location','southeast','FontSize',7);
        ylim(ax4b,[0 100]); grid(ax4b,'on'); box(ax4b,'off');
    end

    % ---- 4c: Stroke volumes — add VSD column
    if numel(allAx4) >= 3
        ax4c = allAx4(3);
        cla(ax4c);
        SV_data = [m_ad3.LVSV, m_pt3.LVSV, m_vsd.LVSV; ...
                   m_ad3.RVSV, m_pt3.RVSV, m_vsd.RVSV];
        bh3 = bar(ax4c, SV_data, 0.65, 'grouped');
        bh3(1).FaceColor = clr_adult;  bh3(1).EdgeColor = 'none';
        bh3(2).FaceColor = clr_pat;    bh3(2).EdgeColor = 'none';
        bh3(3).FaceColor = clr_vsd;    bh3(3).EdgeColor = 'none';
        xticks(ax4c,[1 2]); xticklabels(ax4c,{'LV','RV'});
        ylabel(ax4c,'SV  [mL]','FontSize',9);
        title(ax4c,'Stroke Volumes','FontSize',10,'FontWeight','bold');
        legend(ax4c,{'Adult ref', 'Healthy paed.', 'VSD pre-surg'}, ...
               'Location','northeast','FontSize',7);
        grid(ax4c,'on'); box(ax4c,'off');
    end

    % ---- 4d: Mean pressures — add VSD bar
    if numel(allAx4) >= 4
        ax4d = allAx4(4);
        cla(ax4d);
        P_fields = {'LVP_mean','RVP_mean','SAP_mean','PAP_mean','LAP_mean','RAP_mean'};
        P_labels  = {'LV','RV','SAR','PAR','LA','RA'};
        P_ad_vals = cellfun(@(f) m_ad3.(f), P_fields);
        P_pt_vals = cellfun(@(f) m_pt3.(f), P_fields);
        P_vd_vals = cellfun(@(f) m_vsd.(f), P_fields);
        bh4 = bar(ax4d, [P_ad_vals; P_pt_vals; P_vd_vals]', 0.65, 'grouped');
        bh4(1).FaceColor = clr_adult;  bh4(1).EdgeColor = 'none';
        bh4(2).FaceColor = clr_pat;    bh4(2).EdgeColor = 'none';
        bh4(3).FaceColor = clr_vsd;    bh4(3).EdgeColor = 'none';
        xticks(ax4d, 1:numel(P_labels)); xticklabels(ax4d, P_labels);
        ylabel(ax4d,'Mean Pressure  [mmHg]','FontSize',9);
        title(ax4d,'Mean Chamber/Vessel Pressures','FontSize',10,'FontWeight','bold');
        legend(ax4d,{'Adult ref', 'Healthy paed.', 'VSD pre-surg'}, ...
               'Location','northeast','FontSize',7);
        grid(ax4d,'on'); box(ax4d,'off');
    end

    % ---- 4e: SVR / PVR — add VSD bar
    if numel(allAx4) >= 5
        ax4e = allAx4(5);
        cla(ax4e);
        R_data = [m_ad3.SVR, m_pt3.SVR, m_vsd.SVR; ...
                  m_ad3.PVR, m_pt3.PVR, m_vsd.PVR];
        bh5 = bar(ax4e, R_data, 0.65, 'grouped');
        bh5(1).FaceColor = clr_adult;  bh5(1).EdgeColor = 'none';
        bh5(2).FaceColor = clr_pat;    bh5(2).EdgeColor = 'none';
        bh5(3).FaceColor = clr_vsd;    bh5(3).EdgeColor = 'none';
        xticks(ax4e,[1 2]); xticklabels(ax4e,{'SVR','PVR'});
        ylabel(ax4e,'Resistance  [Wood Units]','FontSize',9);
        title(ax4e,'SVR and PVR','FontSize',10,'FontWeight','bold');
        legend(ax4e,{'Adult ref', 'Healthy paed.', 'VSD pre-surg'}, ...
               'Location','northeast','FontSize',7);
        grid(ax4e,'on'); box(ax4e,'off');
    end

    % ---- 4f: Summary table — extend to 3 columns
    if numel(allAx4) >= 6
        ax4f = allAx4(6);
        cla(ax4f); axis(ax4f,'off');
        CO_ad3 = m_ad3.LVSV * params_adult.HR / 1000;  % [L/min]
        CO_pt3 = m_pt3.LVSV * params_pat.HR   / 1000;  % [L/min]
        CO_vsd = m_vsd.LVSV  * params_vsd.HR  / 1000;  % [L/min]
        tbl = { ...
            sprintf('%-17s %6s %6s %6s', 'Metric', 'Adult', 'Hlthy', 'VSD'); ...
            repmat('-',1,44); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'HR [bpm]', ...
                params_adult.HR, params_pat.HR, params_vsd.HR); ...
            sprintf('%-17s %6.2f %6.2f %6.2f', 'CO [L/min]', CO_ad3, CO_pt3, CO_vsd); ...
            sprintf('%-17s %6.2f %6.2f %6.2f', 'Qp/Qs', ...
                m_ad3.QpQs, m_pt3.QpQs, m_vsd.QpQs); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'LVEDV [mL]', ...
                m_ad3.LVEDV, m_pt3.LVEDV, m_vsd.LVEDV); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'LVESV [mL]', ...
                m_ad3.LVESV, m_pt3.LVESV, m_vsd.LVESV); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'LVEF [%%]', ...
                m_ad3.LVEF*100, m_pt3.LVEF*100, m_vsd.LVEF*100); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'RVEF [%%]', ...
                m_ad3.RVEF*100, m_pt3.RVEF*100, m_vsd.RVEF*100); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'PAP_mean [mmHg]', ...
                m_ad3.PAP_mean, m_pt3.PAP_mean, m_vsd.PAP_mean); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'SVR [WU]', ...
                m_ad3.SVR, m_pt3.SVR, m_vsd.SVR); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'PVR [WU]', ...
                m_ad3.PVR, m_pt3.PVR, m_vsd.PVR); ...
            sprintf('%-17s %6.1f %6.1f %6.1f', 'Q_VSD [mL/s]', ...
                m_ad3.Q_shunt_mean_mLs, m_pt3.Q_shunt_mean_mLs, m_vsd.Q_shunt_mean_mLs) };
        text(ax4f, 0.02, 0.97, strjoin(tbl, '\n'), ...
            'FontSize', 7, 'FontName', 'Courier New', ...
            'VerticalAlignment', 'top', 'Units', 'normalized');
        title(ax4f, 'Summary Table', 'FontSize', 10, 'FontWeight', 'bold');
    end

    sgtitle(fh4(1), ...
        sprintf('Clinical Metrics — Adult / Healthy Paediatric / %s', label_vsd), ...
        'FontSize', 11, 'FontWeight', 'bold');
end

fprintf('[demo] VSD patient overlay complete. All 4 figures updated.\n\n');

%% =========================================================================
%  8. CLINICAL METRICS — 3-column console summary
%
%  compute_clinical_indices derives LVSV, RVSV, EF, CO, SVR, PVR, etc.
%  from the last complete cardiac cycle of each simulation.
% =========================================================================
fprintf('[demo] Computing clinical indices...\n');
metrics_adult = compute_clinical_indices(sim_adult, params_adult);
metrics_pat   = compute_clinical_indices(sim_pat,   params_pat);
metrics_vsd   = compute_clinical_indices(sim_vsd,   params_vsd);

% Cardiac output: SV [mL] × HR [bpm] → CO [L/min]
CO_adult = metrics_adult.LVSV * params_adult.HR / 1000;   % [L/min]
CO_pat   = metrics_pat.LVSV   * params_pat.HR   / 1000;   % [L/min]
CO_vsd_c = metrics_vsd.LVSV   * params_vsd.HR   / 1000;   % [L/min]

fprintf('\n');
fprintf('====================================================================\n');
fprintf('  CLINICAL METRICS SUMMARY — 3-WAY COMPARISON\n');
fprintf('====================================================================\n');
fprintf('  %-22s  %-12s  %-12s  %s\n', 'Metric', 'Adult (ref)', ...
        'Healthy paed.', label_vsd);
fprintf('  ------------------------------------------------------------------\n');
fprintf('  %-22s  %6.1f mL       %6.1f mL       %6.1f mL\n', ...
    'LV Stroke Volume', metrics_adult.LVSV, metrics_pat.LVSV, metrics_vsd.LVSV);
fprintf('  %-22s  %6.1f mL       %6.1f mL       %6.1f mL\n', ...
    'RV Stroke Volume', metrics_adult.RVSV, metrics_pat.RVSV, metrics_vsd.RVSV);
fprintf('  %-22s  %6.1f %%        %6.1f %%        %6.1f %%\n', ...
    'LV Ejection Fraction', ...
    metrics_adult.LVEF * 100, metrics_pat.LVEF * 100, metrics_vsd.LVEF * 100);
fprintf('  %-22s  %6.1f %%        %6.1f %%        %6.1f %%\n', ...
    'RV Ejection Fraction', ...
    metrics_adult.RVEF * 100, metrics_pat.RVEF * 100, metrics_vsd.RVEF * 100);
fprintf('  %-22s  %6.0f bpm      %6.0f bpm      %6.0f bpm\n', ...
    'Heart Rate', params_adult.HR, params_pat.HR, params_vsd.HR);
fprintf('  %-22s  %6.2f L/min    %6.2f L/min    %6.2f L/min\n', ...
    'Cardiac Output', CO_adult, CO_pat, CO_vsd_c);
fprintf('  %-22s  %6.2f WU       %6.2f WU       %6.2f WU\n', ...
    'SVR', metrics_adult.SVR, metrics_pat.SVR, metrics_vsd.SVR);
fprintf('  %-22s  %6.2f WU       %6.2f WU       %6.2f WU\n', ...
    'PVR', metrics_adult.PVR, metrics_pat.PVR, metrics_vsd.PVR);
fprintf('  %-22s  %4.0f/%4.0f/%4.0f     %4.0f/%4.0f/%4.0f     %4.0f/%4.0f/%4.0f  mmHg\n', ...
    'PA sys/dia/mean', ...
    metrics_adult.PAP_max, metrics_adult.PAP_min, metrics_adult.PAP_mean, ...
    metrics_pat.PAP_max,   metrics_pat.PAP_min,   metrics_pat.PAP_mean, ...
    metrics_vsd.PAP_max,   metrics_vsd.PAP_min,   metrics_vsd.PAP_mean);
fprintf('  %-22s  %6.2f           %6.2f           %6.2f\n', ...
    'Qp/Qs', metrics_adult.QpQs, metrics_pat.QpQs, metrics_vsd.QpQs);
fprintf('  %-22s  %6.1f mL/s      %6.1f mL/s      %6.1f mL/s\n', ...
    'Mean VSD shunt flow', ...
    metrics_adult.Q_shunt_mean_mLs, metrics_pat.Q_shunt_mean_mLs, metrics_vsd.Q_shunt_mean_mLs);
fprintf('====================================================================\n\n');
fprintf('[demo] Done. 4 figures open with 3-way comparison.\n');

% =========================================================================
%  LOCAL HELPERS  (replicate plot_simulation_comparison logic without
%  depending on its private helper functions)
% =========================================================================

function [P, Q, tc, Pc, Qc] = local_reconstruct(sim, params)
% LOCAL_RECONSTRUCT — Rebuild P/Q signals from ODE solution; extract last cycle.
%
% Replicates helper_reconstruct from plot_simulation_comparison.m but exposed
% here as a script-level function so demo_sim_comparison can call it for the
% third (VSD) simulation without modifying the shared utility.
%
% INPUTS:
%   sim    - struct from integrate_system  (fields: .t, .V)
%   params - parameter struct with .idx, .HR, .V0, .C, .E, .R.vsd
%
% OUTPUTS:
%   P      - struct of full pressure vectors        [mmHg]
%   Q      - struct of full flow vectors             [mL/s]
%   tc     - time vector for last cycle (re-zeroed)  [s]
%   Pc     - struct of last-cycle pressure vectors   [mmHg]
%   Qc     - struct of last-cycle flow vectors       [mL/s]

t  = sim.t(:);
XV = sim.V;
n  = numel(t);
sidx = params.idx;

P.RA   = zeros(n,1); P.RV   = zeros(n,1);
P.LA   = zeros(n,1); P.LV   = zeros(n,1);
P.SAR  = zeros(n,1); P.SVEN = zeros(n,1);
P.PAR  = zeros(n,1); P.PVEN = zeros(n,1);

Q.TV   = zeros(n,1); Q.PVv  = zeros(n,1);
Q.MV   = zeros(n,1); Q.AV   = zeros(n,1);
Q.SVEN = zeros(n,1); Q.PVEN = zeros(n,1);
Q.VSD  = zeros(n,1);

for i = 1:n
    xi     = XV(i,:)';
    V_RA   = xi(sidx.V_RA);
    V_RV   = xi(sidx.V_RV);
    V_LA   = xi(sidx.V_LA);
    V_LV   = xi(sidx.V_LV);
    V_SAR  = xi(sidx.V_SAR);
    V_SVEN = xi(sidx.V_SVEN);
    V_PAR  = xi(sidx.V_PAR);
    V_PVEN = xi(sidx.V_PVEN);

    [E_LV_i, E_RV_i, E_LA_i, E_RA_i] = elastance_model(t(i), params);

    P_RA_i   = max(E_RA_i * (V_RA   - params.V0.RA),   -5);  % [mmHg]
    P_RV_i   = E_RV_i     * (V_RV   - params.V0.RV);         % [mmHg]
    P_LA_i   = max(E_LA_i * (V_LA   - params.V0.LA),   -5);  % [mmHg]
    P_LV_i   = E_LV_i     * (V_LV   - params.V0.LV);         % [mmHg]
    P_SAR_i  = (V_SAR  - params.V0.SAR)  / params.C.SAR;     % [mmHg]
    P_SVEN_i = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);  % [mmHg]
    P_PAR_i  = (V_PAR  - params.V0.PAR)  / params.C.PAR;     % [mmHg]
    P_PVEN_i = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);  % [mmHg]

    P.RA(i)   = P_RA_i;   P.RV(i)   = P_RV_i;
    P.LA(i)   = P_LA_i;   P.LV(i)   = P_LV_i;
    P.SAR(i)  = P_SAR_i;  P.SVEN(i) = P_SVEN_i;
    P.PAR(i)  = P_PAR_i;  P.PVEN(i) = P_PVEN_i;

    Q.TV(i)   = valve_model(P_RA_i,  P_RV_i,  params);   % [mL/s]
    Q.PVv(i)  = valve_model(P_RV_i,  P_PAR_i, params);   % [mL/s]
    Q.MV(i)   = valve_model(P_LA_i,  P_LV_i,  params);   % [mL/s]
    Q.AV(i)   = valve_model(P_LV_i,  P_SAR_i, params);   % [mL/s]
    Q.SVEN(i) = xi(sidx.Q_SVEN);                         % [mL/s]
    Q.PVEN(i) = xi(sidx.Q_PVEN);                         % [mL/s]
    Q.VSD(i)  = (P_LV_i - P_RV_i) / params.R.vsd;       % [mL/s]  L→R positive
end

% Extract last complete cardiac cycle
mask = local_last_mask(sim, params);
tc   = t(mask) - t(find(mask, 1, 'first'));   % re-zero to [0, T_HB]

flds_P = fieldnames(P);
flds_Q = fieldnames(Q);
Pc = struct(); Qc = struct();
for k = 1:numel(flds_P), Pc.(flds_P{k}) = P.(flds_P{k})(mask); end
for k = 1:numel(flds_Q), Qc.(flds_Q{k}) = Q.(flds_Q{k})(mask); end
end  % local_reconstruct


function mask = local_last_mask(sim, params)
% LOCAL_LAST_MASK — logical mask for the last cardiac cycle in sim.t
t    = sim.t(:);
T_HB = 60 / params.HR;   % [s]
T1   = t(end);
T0   = T1 - T_HB;
mask = (t >= T0) & (t <= T1);
end  % local_last_mask
