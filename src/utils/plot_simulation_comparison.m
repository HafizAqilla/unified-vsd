function plot_simulation_comparison(sim_adult, params_adult, sim_pat, params_pat, label_pat)
% PLOT_SIMULATION_COMPARISON
% -----------------------------------------------------------------------
% Compares a patient simulation (post apply_scaling) against the adult
% reference baseline across four figure windows:
%
%   Figure 1 — PV Loops  (LV and RV)
%   Figure 2 — Pressure Waveforms  (LV, RV, SAR, PAR, LA, RA)
%   Figure 3 — Flow Waveforms  (AV, PVv, MV, TV, Q_SAR, Q_PAR, VSD)
%   Figure 4 — Clinical Metrics Panel  (Qp/Qs bar, EF, volumes, SVR/PVR,
%               CO summary, mean pressures)
%
% INPUTS:
%   sim_adult   - struct from integrate_system(params_adult)
%   params_adult- adult reference params (from default_parameters)
%   sim_pat     - struct from integrate_system(params_pat)
%   params_pat  - patient-scaled params (from apply_scaling + params_from_clinical)
%   label_pat   - (optional) string label for the patient, e.g. 'Patient A (2yr, 12kg)'
%
% USAGE EXAMPLE:
%   params_ref  = default_parameters();
%   params_adult = params_from_clinical(apply_scaling(params_ref, ...
%       struct('age_years',30,'weight_kg',70,'height_cm',175,'sex','M')), ...
%       patient_template(), 'pre_surgery');
%   sim_adult   = integrate_system(params_adult);
%
%   patient.age_years = 2; patient.weight_kg = 12;
%   patient.height_cm = 87; patient.sex = 'M';
%   params_pat  = apply_scaling(params_ref, patient);
%   params_pat  = params_from_clinical(params_pat, clinical, 'pre_surgery');
%   sim_pat     = integrate_system(params_pat);
%
%   plot_simulation_comparison(sim_adult, params_adult, sim_pat, params_pat, ...
%       'Patient (2 yr, 12 kg)');
%
% SIGN CONVENTIONS:
%   Q_VSD > 0  →  left-to-right shunt  (net L→R)
%   Flows in [mL/s] for waveforms; [L/min] for clinical metrics.
%
% REFERENCES:
%   [1] compute_clinical_indices.m — metric definitions.
%   [2] reconstruct_signals (local copy below) — P/Q signal rebuild.
%   [3] VIBECODING_GUARDRAILS.md §11 — Plotting Standards.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-31
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 5 || isempty(label_pat)
    if isfield(params_pat, 'scaling') && isfield(params_pat.scaling, 'patient')
        p = params_pat.scaling.patient;
        label_pat = sprintf('Patient (%.1f yr, %.0f kg)', p.age_years, p.weight_kg);
    else
        label_pat = 'Patient';
    end
end

%% ---- shared colour / style tokens ------------------------------------
clr_adult = [0.18 0.44 0.70];   % blue
clr_pat   = [0.84 0.22 0.22];   % red
lw_adult  = 1.5;
lw_pat    = 1.8;
font_ax   = 9;
font_ti   = 10;

%% ---- reconstruct P and Q signals for both simulations ----------------
fprintf('[plot_simulation_comparison] Reconstructing signals...\n');
[P_ad, Q_ad, tc_ad, Pc_ad, Qc_ad] = helper_reconstruct(sim_adult, params_adult);
[P_pt, Q_pt, tc_pt, Pc_pt, Qc_pt] = helper_reconstruct(sim_pat,   params_pat);

%% ---- clinical metrics for both --------------------------------------
m_ad = compute_clinical_indices(sim_adult, params_adult);
m_pt = compute_clinical_indices(sim_pat,   params_pat);

mLs_Lmin_ad = params_adult.conv.mLs_to_Lmin;
mLs_Lmin_pt = params_pat.conv.mLs_to_Lmin;

%% ====================================================================
%  FIGURE 1 — PV LOOPS
%% ====================================================================
fh1 = figure('Name', 'PV Loops', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [40 500 820 420]);
sgtitle(sprintf('PV Loops — Adult vs %s', label_pat), ...
    'FontSize', font_ti+1, 'FontWeight', 'bold');

sidx_ad = params_adult.idx;
sidx_pt = params_pat.idx;

% --- LV PV loop -------------------------------------------------------
ax_lv = subplot(1,2,1);
V_LV_ad = sim_adult.V(helper_last_mask(sim_adult, params_adult), sidx_ad.V_LV);
V_LV_pt = sim_pat.V(helper_last_mask(sim_pat, params_pat),       sidx_pt.V_LV);

plot(ax_lv, V_LV_ad, Pc_ad.LV, 'Color', clr_adult, 'LineWidth', lw_adult);
hold(ax_lv, 'on');
plot(ax_lv, V_LV_pt, Pc_pt.LV, 'Color', clr_pat,   'LineWidth', lw_pat);
hold(ax_lv, 'off');
xlabel(ax_lv, 'V_{LV}  [mL]',  'FontSize', font_ax);
ylabel(ax_lv, 'P_{LV}  [mmHg]','FontSize', font_ax);
title(ax_lv,  'Left Ventricle PV Loop', 'FontSize', font_ti, 'FontWeight', 'bold');
legend(ax_lv, {'Adult ref', label_pat}, 'Location', 'northwest', 'FontSize', 8);
grid(ax_lv, 'on'); box(ax_lv, 'off');

% Annotate EF on each loop
text(ax_lv, 0.98, 0.97, sprintf('Adult  LVEF = %.0f%%', m_ad.LVEF*100), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color',clr_adult);
text(ax_lv, 0.98, 0.90, sprintf('Patient LVEF = %.0f%%', m_pt.LVEF*100), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color',clr_pat);

% --- RV PV loop -------------------------------------------------------
ax_rv = subplot(1,2,2);
V_RV_ad = sim_adult.V(helper_last_mask(sim_adult, params_adult), sidx_ad.V_RV);
V_RV_pt = sim_pat.V(helper_last_mask(sim_pat, params_pat),       sidx_pt.V_RV);

plot(ax_rv, V_RV_ad, Pc_ad.RV, 'Color', clr_adult, 'LineWidth', lw_adult);
hold(ax_rv, 'on');
plot(ax_rv, V_RV_pt, Pc_pt.RV, 'Color', clr_pat,   'LineWidth', lw_pat);
hold(ax_rv, 'off');
xlabel(ax_rv, 'V_{RV}  [mL]',  'FontSize', font_ax);
ylabel(ax_rv, 'P_{RV}  [mmHg]','FontSize', font_ax);
title(ax_rv,  'Right Ventricle PV Loop', 'FontSize', font_ti, 'FontWeight', 'bold');
legend(ax_rv, {'Adult ref', label_pat}, 'Location', 'northwest', 'FontSize', 8);
grid(ax_rv, 'on'); box(ax_rv, 'off');

text(ax_rv, 0.98, 0.97, sprintf('Adult  RVEF = %.0f%%', m_ad.RVEF*100), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color',clr_adult);
text(ax_rv, 0.98, 0.90, sprintf('Patient RVEF = %.0f%%', m_pt.RVEF*100), ...
    'Units','normalized','HorizontalAlignment','right','FontSize',7,'Color',clr_pat);

%% ====================================================================
%  FIGURE 2 — PRESSURE WAVEFORMS  (last 1 cardiac cycle, normalised time)
%% ====================================================================
fh2 = figure('Name', 'Pressure Waveforms', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [60 300 1200 680]);
sgtitle(sprintf('Pressure Waveforms — Adult vs %s', label_pat), ...
    'FontSize', font_ti+1, 'FontWeight', 'bold');

T_ad = 60 / params_adult.HR;   % [s]
T_pt = 60 / params_pat.HR;     % [s]

% Normalise time axes to [0 1] (fraction of cardiac cycle)
tn_ad = tc_ad / T_ad;
tn_pt = tc_pt / T_pt;

pressure_panels = { ...
    'LV',   'P_{LV}  [mmHg]',   'Left Ventricular Pressure'; ...
    'RV',   'P_{RV}  [mmHg]',   'Right Ventricular Pressure'; ...
    'SAR',  'P_{SAR}  [mmHg]',  'Systemic Arterial Pressure'; ...
    'PAR',  'P_{PAR}  [mmHg]',  'Pulmonary Arterial Pressure'; ...
    'LA',   'P_{LA}  [mmHg]',   'Left Atrial Pressure'; ...
    'RA',   'P_{RA}  [mmHg]',   'Right Atrial Pressure' };

for k = 1:size(pressure_panels,1)
    fld  = pressure_panels{k,1};
    ylbl = pressure_panels{k,2};
    ttl  = pressure_panels{k,3};

    ax = subplot(3,2,k);
    plot(ax, tn_ad, Pc_ad.(fld), 'Color', clr_adult, 'LineWidth', lw_adult);
    hold(ax,'on');
    plot(ax, tn_pt, Pc_pt.(fld), 'Color', clr_pat,   'LineWidth', lw_pat, 'LineStyle','--');
    hold(ax,'off');
    xlabel(ax, 'Cardiac cycle fraction', 'FontSize', font_ax);
    ylabel(ax, ylbl, 'FontSize', font_ax);
    title(ax,  ttl,  'FontSize', font_ti, 'FontWeight', 'bold');
    if k == 1
        legend(ax, {'Adult ref', label_pat}, 'Location', 'northeast', 'FontSize', 8);
    end
    grid(ax,'on'); box(ax,'off');
end

%% ====================================================================
%  FIGURE 3 — FLOW WAVEFORMS
%% ====================================================================
fh3 = figure('Name', 'Flow Waveforms', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [80 100 1200 720]);
sgtitle(sprintf('Flow Waveforms — Adult vs %s', label_pat), ...
    'FontSize', font_ti+1, 'FontWeight', 'bold');

flow_panels = { ...
    'AV',   'Q_{AV}  [mL/s]',   'Aortic Valve Flow'; ...
    'PVv',  'Q_{PVv}  [mL/s]',  'Pulmonary Valve Flow'; ...
    'MV',   'Q_{MV}  [mL/s]',   'Mitral Valve Flow'; ...
    'TV',   'Q_{TV}  [mL/s]',   'Tricuspid Valve Flow'; ...
    'SVEN', 'Q_{SVEN}  [mL/s]', 'Systemic Venous Flow  (Q_{SAR} state)'; ...
    'PVEN', 'Q_{PVEN}  [mL/s]', 'Pulmonary Venous Flow  (Q_{PVEN} state)' };

for k = 1:size(flow_panels,1)
    fld  = flow_panels{k,1};
    ylbl = flow_panels{k,2};
    ttl  = flow_panels{k,3};

    ax = subplot(4,2,k);
    plot(ax, tn_ad, Qc_ad.(fld), 'Color', clr_adult, 'LineWidth', lw_adult);
    hold(ax,'on');
    plot(ax, tn_pt, Qc_pt.(fld), 'Color', clr_pat,   'LineWidth', lw_pat, 'LineStyle','--');
    yline(ax, 0, 'k:', 'LineWidth', 0.8);   % zero-flow reference
    hold(ax,'off');
    xlabel(ax, 'Cardiac cycle fraction', 'FontSize', font_ax);
    ylabel(ax, ylbl, 'FontSize', font_ax);
    title(ax,  ttl,  'FontSize', font_ti, 'FontWeight', 'bold');
    if k == 1
        legend(ax, {'Adult ref', label_pat}, 'Location', 'northeast', 'FontSize', 8);
    end
    grid(ax,'on'); box(ax,'off');
end

% --- VSD shunt (subplot 7) ----------------------------------------
ax_vsd = subplot(4,2,7);
plot(ax_vsd, tn_ad, Qc_ad.VSD, 'Color', clr_adult, 'LineWidth', lw_adult);
hold(ax_vsd,'on');
plot(ax_vsd, tn_pt, Qc_pt.VSD, 'Color', clr_pat,   'LineWidth', lw_pat, 'LineStyle','--');
yline(ax_vsd, 0, 'k:', 'LineWidth', 0.8);
hold(ax_vsd,'off');
xlabel(ax_vsd, 'Cardiac cycle fraction', 'FontSize', font_ax);
ylabel(ax_vsd, 'Q_{VSD}  [mL/s]',        'FontSize', font_ax);
title(ax_vsd,  'VSD Shunt Flow  (+ve = L→R)', 'FontSize', font_ti, 'FontWeight', 'bold');
grid(ax_vsd,'on'); box(ax_vsd,'off');

% --- Flow legend / annotation box (subplot 8) ----------------------
ax_leg = subplot(4,2,8);
axis(ax_leg,'off');
txt = { ...
    sprintf('  %-18s  Adult     Patient', 'Metric'), ...
    sprintf('  %-18s  %-8.2f  %-8.2f',  'Qp [L/min]',  ...
        mean(Qc_ad.PVEN)*mLs_Lmin_ad, mean(Qc_pt.PVEN)*mLs_Lmin_pt), ...
    sprintf('  %-18s  %-8.2f  %-8.2f',  'Qs [L/min]',  ...
        mean(Qc_ad.SVEN)*mLs_Lmin_ad, mean(Qc_pt.SVEN)*mLs_Lmin_pt), ...
    sprintf('  %-18s  %-8.2f  %-8.2f',  'Qp/Qs',       m_ad.QpQs, m_pt.QpQs), ...
    sprintf('  %-18s  %-8.1f  %-8.1f',  'Q_VSD_mean [mL/s]', m_ad.Q_shunt_mean_mLs, m_pt.Q_shunt_mean_mLs) };
text(ax_leg, 0.05, 0.5, strjoin(txt, '\n'), ...
    'FontSize', 8.5, 'FontName', 'Courier New', ...
    'VerticalAlignment', 'middle', 'Units', 'normalized');
title(ax_leg, 'Flow Summary', 'FontSize', font_ti, 'FontWeight', 'bold');

%% ====================================================================
%  FIGURE 4 — CLINICAL METRICS PANEL
%% ====================================================================
fh4 = figure('Name', 'Clinical Metrics', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [100 60 1050 680]);
sgtitle(sprintf('Clinical Metrics — Adult vs %s', label_pat), ...
    'FontSize', font_ti+1, 'FontWeight', 'bold');

% --- 4a: Qp/Qs bar ---------------------------------------------------
ax4a = subplot(2,3,1);
bar(ax4a, [m_ad.QpQs, m_pt.QpQs], 0.55, 'FaceColor', 'flat', ...
    'CData', [clr_adult; clr_pat], 'EdgeColor', 'none');
yline(ax4a, 1.0, 'k--', 'LineWidth', 1.2, 'Label', 'Balanced (1.0)');
xticks(ax4a, [1 2]);  xticklabels(ax4a, {'Adult ref', label_pat});
ylabel(ax4a, 'Qp/Qs  [-]', 'FontSize', font_ax);
title(ax4a, 'Qp/Qs Ratio', 'FontSize', font_ti, 'FontWeight', 'bold');
text(ax4a, 1, m_ad.QpQs*1.02, sprintf('%.2f', m_ad.QpQs), ...
    'HorizontalAlignment','center','FontSize',9,'Color',clr_adult);
text(ax4a, 2, m_pt.QpQs*1.02, sprintf('%.2f', m_pt.QpQs), ...
    'HorizontalAlignment','center','FontSize',9,'Color',clr_pat);
grid(ax4a,'on'); box(ax4a,'off');

% --- 4b: Ejection fractions (LV + RV) --------------------------------
ax4b = subplot(2,3,2);
EF_data = [m_ad.LVEF, m_pt.LVEF; m_ad.RVEF, m_pt.RVEF] * 100;
bh = bar(ax4b, EF_data, 0.65, 'grouped');
bh(1).FaceColor = clr_adult;  bh(1).EdgeColor = 'none';
bh(2).FaceColor = clr_pat;    bh(2).EdgeColor = 'none';
yline(ax4b, 55, 'k--', 'LineWidth', 0.8, 'Label', 'Normal LV EF (55%)');
xticks(ax4b,[1 2]); xticklabels(ax4b,{'LV','RV'});
ylabel(ax4b,'EF  [%]','FontSize',font_ax);
title(ax4b,'Ejection Fractions','FontSize',font_ti,'FontWeight','bold');
legend(ax4b,{'Adult ref', label_pat},'Location','southeast','FontSize',7);
ylim(ax4b,[0 100]); grid(ax4b,'on'); box(ax4b,'off');

% --- 4c: Stroke volumes (LV + RV) ------------------------------------
ax4c = subplot(2,3,3);
SV_data = [m_ad.LVSV, m_pt.LVSV; m_ad.RVSV, m_pt.RVSV];
bh2 = bar(ax4c, SV_data, 0.65, 'grouped');
bh2(1).FaceColor = clr_adult;  bh2(1).EdgeColor = 'none';
bh2(2).FaceColor = clr_pat;    bh2(2).EdgeColor = 'none';
xticks(ax4c,[1 2]); xticklabels(ax4c,{'LV','RV'});
ylabel(ax4c,'SV  [mL]','FontSize',font_ax);
title(ax4c,'Stroke Volumes','FontSize',font_ti,'FontWeight','bold');
legend(ax4c,{'Adult ref', label_pat},'Location','northeast','FontSize',7);
grid(ax4c,'on'); box(ax4c,'off');

% --- 4d: Mean pressures (grouped bar) --------------------------------
ax4d = subplot(2,3,4);
P_fields = {'LVP_mean','RVP_mean','SAP_mean','PAP_mean','LAP_mean','RAP_mean'};
P_labels  = {'LV','RV','SAR','PAR','LA','RA'};
P_ad_vals = cellfun(@(f) m_ad.(f), P_fields);
P_pt_vals = cellfun(@(f) m_pt.(f), P_fields);
bh3 = bar(ax4d, [P_ad_vals; P_pt_vals]', 0.65, 'grouped');
bh3(1).FaceColor = clr_adult;  bh3(1).EdgeColor = 'none';
bh3(2).FaceColor = clr_pat;    bh3(2).EdgeColor = 'none';
xticks(ax4d, 1:numel(P_labels)); xticklabels(ax4d, P_labels);
ylabel(ax4d,'Mean Pressure  [mmHg]','FontSize',font_ax);
title(ax4d,'Mean Chamber/Vessel Pressures','FontSize',font_ti,'FontWeight','bold');
legend(ax4d,{'Adult ref', label_pat},'Location','northeast','FontSize',7);
grid(ax4d,'on'); box(ax4d,'off');

% --- 4e: SVR / PVR (Wood units) --------------------------------------
ax4e = subplot(2,3,5);
R_data = [m_ad.SVR, m_pt.SVR; m_ad.PVR, m_pt.PVR];
bh4 = bar(ax4e, R_data, 0.65, 'grouped');
bh4(1).FaceColor = clr_adult;  bh4(1).EdgeColor = 'none';
bh4(2).FaceColor = clr_pat;    bh4(2).EdgeColor = 'none';
xticks(ax4e,[1 2]); xticklabels(ax4e,{'SVR','PVR'});
ylabel(ax4e,'Resistance  [Wood Units]','FontSize',font_ax);
title(ax4e,'SVR and PVR','FontSize',font_ti,'FontWeight','bold');
legend(ax4e,{'Adult ref', label_pat},'Location','northeast','FontSize',7);
grid(ax4e,'on'); box(ax4e,'off');

% --- 4f: Summary text table ------------------------------------------
ax4f = subplot(2,3,6);
axis(ax4f,'off');
CO_ad = m_ad.LVSV * params_adult.HR / 1000;    % [L/min]  SV[mL] × HR[bpm] / 1000
CO_pt = m_pt.LVSV * params_pat.HR   / 1000;    % [L/min]

tbl = { ...
    sprintf('%-20s  %6s  %6s',  'Metric',           'Adult',   'Patient'); ...
    repmat('-',1,42); ...
    sprintf('%-20s  %6.1f  %6.1f',  'HR [bpm]',   params_adult.HR,         params_pat.HR); ...
    sprintf('%-20s  %6.2f  %6.2f',  'CO [L/min]', CO_ad,                   CO_pt); ...
    sprintf('%-20s  %6.2f  %6.2f',  'Qp/Qs',      m_ad.QpQs,               m_pt.QpQs); ...
    sprintf('%-20s  %6.1f  %6.1f',  'LVEDV [mL]', m_ad.LVEDV,              m_pt.LVEDV); ...
    sprintf('%-20s  %6.1f  %6.1f',  'LVESV [mL]', m_ad.LVESV,              m_pt.LVESV); ...
    sprintf('%-20s  %6.1f  %6.1f',  'LVEF [%%]',  m_ad.LVEF*100,           m_pt.LVEF*100); ...
    sprintf('%-20s  %6.1f  %6.1f',  'RVEF [%%]',  m_ad.RVEF*100,           m_pt.RVEF*100); ...
    sprintf('%-20s  %6.1f  %6.1f',  'PAP_mean [mmHg]', m_ad.PAP_mean,      m_pt.PAP_mean); ...
    sprintf('%-20s  %6.1f  %6.1f',  'SVR [WU]',   m_ad.SVR,               m_pt.SVR); ...
    sprintf('%-20s  %6.1f  %6.1f',  'PVR [WU]',   m_ad.PVR,               m_pt.PVR); ...
    sprintf('%-20s  %6.1f  %6.1f',  'Q_VSD [mL/s]',m_ad.Q_shunt_mean_mLs, m_pt.Q_shunt_mean_mLs) };
text(ax4f, 0.02, 0.97, strjoin(tbl, '\n'), ...
    'FontSize', 7.5, 'FontName', 'Courier New', ...
    'VerticalAlignment', 'top', 'Units', 'normalized');
title(ax4f, 'Summary Table', 'FontSize', font_ti, 'FontWeight', 'bold');

fprintf('[plot_simulation_comparison] 4 figures rendered for: %s\n', label_pat);

end  % plot_simulation_comparison


% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function [P, Q, tc, Pc, Qc] = helper_reconstruct(sim, params)
% HELPER_RECONSTRUCT — rebuild P/Q from full sim, then extract last cycle.
t  = sim.t(:);
XV = sim.V;

% ---- full-signal reconstruction ----
n    = numel(t);
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
    P_RV_i   = E_RV_i     * (V_RV   - params.V0.RV);
    P_LA_i   = max(E_LA_i * (V_LA   - params.V0.LA),   -5);
    P_LV_i   = E_LV_i     * (V_LV   - params.V0.LV);
    P_SAR_i  = (V_SAR  - params.V0.SAR)  / params.C.SAR;
    P_SVEN_i = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);
    P_PAR_i  = (V_PAR  - params.V0.PAR)  / params.C.PAR;
    P_PVEN_i = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);

    P.RA(i)   = P_RA_i;   P.RV(i)   = P_RV_i;
    P.LA(i)   = P_LA_i;   P.LV(i)   = P_LV_i;
    P.SAR(i)  = P_SAR_i;  P.SVEN(i) = P_SVEN_i;
    P.PAR(i)  = P_PAR_i;  P.PVEN(i) = P_PVEN_i;

    Q.TV(i)   = valve_model(P_RA_i,  P_RV_i,  params);
    Q.PVv(i)  = valve_model(P_RV_i,  P_PAR_i, params);
    Q.MV(i)   = valve_model(P_LA_i,  P_LV_i,  params);
    Q.AV(i)   = valve_model(P_LV_i,  P_SAR_i, params);
    Q.SVEN(i) = xi(sidx.Q_SVEN);
    Q.PVEN(i) = xi(sidx.Q_PVEN);
    Q.VSD(i)  = (P_LV_i - P_RV_i) / params.R.vsd;
end

% ---- extract last complete cardiac cycle ----
T_HB = 60 / params.HR;
mask = helper_last_mask(sim, params);
tc = t(mask) - t(find(mask,1,'first'));   % re-zero to [0, T_HB]

flds_P = fieldnames(P);
flds_Q = fieldnames(Q);
Pc = struct();  Qc = struct();
for k = 1:numel(flds_P), Pc.(flds_P{k}) = P.(flds_P{k})(mask); end
for k = 1:numel(flds_Q), Qc.(flds_Q{k}) = Q.(flds_Q{k})(mask); end
end


function mask = helper_last_mask(sim, params)
% HELPER_LAST_MASK — logical mask for the last cardiac cycle in sim.t
t    = sim.t(:);
T_HB = 60 / params.HR;   % [s]
T1   = t(end);
T0   = T1 - T_HB;
mask = (t >= T0) & (t <= T1);
end
