function plotting_tools(sim, params, tag, scenario)
% PLOTTING_TOOLS
% -----------------------------------------------------------------------
% Generate publication-quality haemodynamic plots from a steady-state
% simulation.  Plots are scenario-aware (pre/post surgery annotations).
%
% INPUTS:
%   sim      - struct from integrate_system.m  (.t, .V)
%   params   - parameter struct (from apply_scaling + params_from_clinical)
%   tag      - string label appended to figure titles, e.g. 'Baseline'
%   scenario - 'pre_surgery' | 'post_surgery'   (controls annotations)
%
% PLOTS GENERATED (all exported as vector PDF, Guardrail §11.3):
%   1. Four-chamber pressure traces  (RA, RV, LA, LV)
%   2. Four-vessel pressure traces   (SAR, SVEN, PAR, PVEN)
%   3. LV and RV pressure-volume loops
%   4. VSD shunt flow  Q_VSD(t)  (pre-surgery only)
%   5. Pulmonary valve flow  Q_PVv(t)
%
% OUTPUT FILES:
%   results/figures/ChamberPressures_<tag>_<scenario>.pdf
%   results/figures/VascularPressures_<tag>_<scenario>.pdf
%   results/figures/PVLoops_<tag>_<scenario>.pdf
%   results/figures/VSDShuntFlow_<tag>_<scenario>.pdf   (pre-surgery only)
%   results/figures/PulmValveFlow_<tag>_<scenario>.pdf
%
% ASSUMPTIONS:
%   - All figures are publication-size: 16×12 cm, Arial 11 pt (Guardrail §11.1).
%   - Axes labelled with units in parentheses.
%   - Exported as vector PDF (never PNG, Guardrail §11.3).
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.1
% -----------------------------------------------------------------------

t  = sim.t(:);
XV = sim.V;

%% Locate results/figures directory (create if absent)
root_dir   = fileparts(fileparts(mfilename('fullpath')));   % project root
fig_dir    = fullfile(root_dir, 'results', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Extract last cycle
T_HB = 60 / params.HR;     % [s]  cardiac period
T1   = t(end);
T0   = T1 - T_HB;
time_mask = (t >= T0) & (t <= T1);
tc   = t(time_mask);        % [s]  trimmed time axis

%% Reconstruct signals (inline helper below)
[P, Q] = signals_for_plot(tc, XV(time_mask,:), params);

lw = 1.8;   % [pt]  line width — publication standard

%% Helper: apply publication axes style (Guardrail §11.1)
% (Called after each plot command)
function style_axes(ax)
    set(ax, 'FontSize', 10, 'FontName', 'Arial', 'Box', 'on', ...
            'TickDir', 'out', 'LineWidth', 0.8);
end

function h = std_fig(fig_name)
    h = figure('Name', fig_name, 'Color', 'w', 'NumberTitle', 'off', ...
               'Units', 'centimeters', 'Position', [2 2 16 12]);
end

function export_fig(fig_handle, fig_dir, filename)
    out_path = fullfile(fig_dir, [filename '.pdf']);
    exportgraphics(fig_handle, out_path, 'ContentType', 'vector', 'Resolution', 300);
    fprintf('[plotting_tools] Saved: %s\n', out_path);
end

%% ---- 1. Chamber pressures -------------------------------------------
fig1_name = sprintf('ChamberPressures_%s_%s', tag, scenario);
fig1 = std_fig(fig1_name);

subplot(2,2,1);
plot(tc, P.RA, 'b', 'LineWidth', lw);
xlabel('Time (s)',        'FontSize', 11, 'FontName', 'Arial');
ylabel('P_{RA} (mmHg)',  'FontSize', 11, 'FontName', 'Arial');
title('Right Atrium',     'FontSize', 12, 'FontName', 'Arial');
style_axes(gca); grid on;

subplot(2,2,2);
plot(tc, P.RV, 'r', 'LineWidth', lw);
xlabel('Time (s)',        'FontSize', 11, 'FontName', 'Arial');
ylabel('P_{RV} (mmHg)',  'FontSize', 11, 'FontName', 'Arial');
title('Right Ventricle',  'FontSize', 12, 'FontName', 'Arial');
style_axes(gca); grid on;

subplot(2,2,3);
plot(tc, P.LA, 'b--', 'LineWidth', lw);
xlabel('Time (s)',        'FontSize', 11, 'FontName', 'Arial');
ylabel('P_{LA} (mmHg)',  'FontSize', 11, 'FontName', 'Arial');
title('Left Atrium',      'FontSize', 12, 'FontName', 'Arial');
style_axes(gca); grid on;

subplot(2,2,4);
plot(tc, P.LV, 'r--', 'LineWidth', lw);
xlabel('Time (s)',        'FontSize', 11, 'FontName', 'Arial');
ylabel('P_{LV} (mmHg)',  'FontSize', 11, 'FontName', 'Arial');
title('Left Ventricle',   'FontSize', 12, 'FontName', 'Arial');
style_axes(gca); grid on;

sgtitle(sprintf('Chamber Pressures — %s | %s', tag, scenario), ...
        'FontSize', 13, 'FontName', 'Arial', 'FontWeight', 'bold');
export_fig(fig1, fig_dir, fig1_name);

%% ---- 2. Vascular pressures -----------------------------------------
fig2_name = sprintf('VascularPressures_%s_%s', tag, scenario);
fig2 = std_fig(fig2_name);

plot(tc, P.SAR,  'k',   'LineWidth', lw, 'DisplayName', 'SAR (systemic artery)');
hold on;
plot(tc, P.SVEN, 'b',   'LineWidth', lw, 'DisplayName', 'SVEN (systemic vein)');
plot(tc, P.PAR,  'r',   'LineWidth', lw, 'DisplayName', 'PAR (pulmonary artery)');
plot(tc, P.PVEN, 'm--', 'LineWidth', lw, 'DisplayName', 'PVEN (pulmonary vein)');
xlabel('Time (s)',          'FontSize', 11, 'FontName', 'Arial');
ylabel('Pressure (mmHg)',   'FontSize', 11, 'FontName', 'Arial');
title(sprintf('Vascular Pressures — %s | %s', tag, scenario), ...
      'FontSize', 12, 'FontName', 'Arial');
legend('FontSize', 10, 'FontName', 'Arial', 'Location', 'best');
style_axes(gca); grid on;
export_fig(fig2, fig_dir, fig2_name);

%% ---- 3. PV loops ---------------------------------------------------
sidx    = params.idx;
V_LV_c  = XV(time_mask, sidx.V_LV);   % [mL]  (Guardrail §7.1)
V_RV_c  = XV(time_mask, sidx.V_RV);   % [mL]

fig3_name = sprintf('PVLoops_%s_%s', tag, scenario);
fig3 = std_fig(fig3_name);

subplot(1,2,1);
plot(V_LV_c, P.LV, 'r', 'LineWidth', lw);
xlabel('V_{LV} (mL)',                'FontSize', 11, 'FontName', 'Arial');
ylabel('P_{LV} (mmHg)',              'FontSize', 11, 'FontName', 'Arial');
title('LV Pressure-Volume Loop',     'FontSize', 12, 'FontName', 'Arial');
style_axes(gca); grid on;

subplot(1,2,2);
plot(V_RV_c, P.RV, 'b', 'LineWidth', lw);
xlabel('V_{RV} (mL)',                'FontSize', 11, 'FontName', 'Arial');
ylabel('P_{RV} (mmHg)',              'FontSize', 11, 'FontName', 'Arial');
title('RV Pressure-Volume Loop',     'FontSize', 12, 'FontName', 'Arial');
style_axes(gca); grid on;

sgtitle(sprintf('PV Loops — %s | %s', tag, scenario), ...
        'FontSize', 13, 'FontName', 'Arial', 'FontWeight', 'bold');
export_fig(fig3, fig_dir, fig3_name);

%% ---- 4. VSD shunt flow (pre-surgery only) ---------------------------
if strcmp(scenario, 'pre_surgery')
    fig4_name = sprintf('VSDShuntFlow_%s_%s', tag, scenario);
    fig4 = std_fig(fig4_name);
    plot(tc, Q.VSD, 'g', 'LineWidth', lw);
    xlabel('Time (s)',      'FontSize', 11, 'FontName', 'Arial');
    ylabel('Q_{VSD} (mL/s)', 'FontSize', 11, 'FontName', 'Arial');
    title(sprintf('VSD Shunt Flow (positive = L\\rightarrowR) — %s', tag), ...
          'FontSize', 12, 'FontName', 'Arial');
    yline(0, 'k--', 'LineWidth', 1);
    style_axes(gca); grid on;
    export_fig(fig4, fig_dir, fig4_name);
end

%% ---- 5. Pulmonary valve flow ----------------------------------------
fig5_name = sprintf('PulmValveFlow_%s_%s', tag, scenario);
fig5 = std_fig(fig5_name);
plot(tc, Q.PVv, 'b', 'LineWidth', lw);
xlabel('Time (s)',        'FontSize', 11, 'FontName', 'Arial');
ylabel('Q_{PVv} (mL/s)', 'FontSize', 11, 'FontName', 'Arial');
title(sprintf('Pulmonary Valve Flow — %s | %s', tag, scenario), ...
      'FontSize', 12, 'FontName', 'Arial');
yline(0, 'k--', 'LineWidth', 1);
style_axes(gca); grid on;
export_fig(fig5, fig_dir, fig5_name);

end  % plotting_tools

% =========================================================================
%  LOCAL HELPER
% =========================================================================

function [P, Q] = signals_for_plot(tc, XV, params)
% SIGNALS_FOR_PLOT — rebuild P and Q arrays from a trimmed state window
%   Returns P and Q structs for an already-trimmed time+state window.
%   Uses params.idx for state access (Guardrail §7.1).
%   P fields [mmHg]: RA, RV, LA, LV, SAR, SVEN, PAR, PVEN
%   Q fields [mL/s]: TV, PVv, MV, AV, SVEN, PVEN, VSD

n    = numel(tc);
sidx = params.idx;   % state index struct

P.RA = zeros(n,1); P.RV = zeros(n,1); P.LA = zeros(n,1); P.LV = zeros(n,1);
P.SAR= zeros(n,1); P.SVEN=zeros(n,1); P.PAR=zeros(n,1); P.PVEN=zeros(n,1);
Q.TV = zeros(n,1); Q.PVv=zeros(n,1); Q.MV=zeros(n,1); Q.AV=zeros(n,1);
Q.VSD= zeros(n,1); Q.SVEN=zeros(n,1); Q.PVEN=zeros(n,1);

for i = 1:n
    xi       = XV(i,:)';
    V_RA     = xi(sidx.V_RA);     % [mL]
    V_RV     = xi(sidx.V_RV);     % [mL]
    V_LA     = xi(sidx.V_LA);     % [mL]
    V_LV     = xi(sidx.V_LV);     % [mL]
    V_SAR    = xi(sidx.V_SAR);    % [mL]
    V_SVEN   = xi(sidx.V_SVEN);   % [mL]
    Q_SVEN_i = xi(sidx.Q_SVEN);   % [mL/s]
    V_PAR    = xi(sidx.V_PAR);    % [mL]
    V_PVEN   = xi(sidx.V_PVEN);   % [mL]
    Q_PVEN_i = xi(sidx.Q_PVEN);   % [mL/s]

    [E_LV_i, E_RV_i, E_LA_i, E_RA_i] = elastance_model(tc(i), params);

    P_RA_i = max(E_RA_i*(V_RA - params.V0.RA), -5);   % [mmHg]
    P_RV_i = E_RV_i*(V_RV - params.V0.RV);             % [mmHg]
    P_LA_i = max(E_LA_i*(V_LA - params.V0.LA), -5);   % [mmHg]
    P_LV_i = E_LV_i*(V_LV - params.V0.LV);             % [mmHg]

    P_SAR_i  = (V_SAR  - params.V0.SAR)  / params.C.SAR;              % [mmHg]
    P_SVEN_i = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);    % [mmHg]
    P_PAR_i  = (V_PAR  - params.V0.PAR)  / params.C.PAR;              % [mmHg]
    P_PVEN_i = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);    % [mmHg]

    P.RA(i) = P_RA_i; P.RV(i) = P_RV_i;
    P.LA(i) = P_LA_i; P.LV(i) = P_LV_i;
    P.SAR(i) = P_SAR_i; P.SVEN(i) = P_SVEN_i;
    P.PAR(i) = P_PAR_i; P.PVEN(i) = P_PVEN_i;

    Q.TV(i)  = valve_model(P_RA_i, P_RV_i,  params);   % [mL/s]
    Q.PVv(i) = valve_model(P_RV_i, P_PAR_i, params);   % [mL/s]
    Q.MV(i)  = valve_model(P_LA_i, P_LV_i,  params);   % [mL/s]
    Q.AV(i)  = valve_model(P_LV_i, P_SAR_i, params);   % [mL/s]
    Q.SVEN(i)= Q_SVEN_i;                                % [mL/s]
    Q.PVEN(i)= Q_PVEN_i;                                % [mL/s]
    Q.VSD(i) = (P_LV_i - P_RV_i) / params.R.vsd;       % [mL/s] positive = L→R

    V_PAR; %#ok<VUNUS>
end
end  % signals_for_plot
