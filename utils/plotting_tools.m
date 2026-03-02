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
% PLOTS GENERATED:
%   1. Four-chamber pressure traces  (RA, RV, LA, LV)
%   2. Four-vessel pressure traces   (SAR, SVEN, PAR, PVEN)
%   3. LV and RV pressure-volume loops
%   4. VSD shunt flow  Q_VSD(t)  (pre-surgery only; omitted post-surgery)
%   5. Pulmonary valve flow  Q_PVv(t)
%
% ASSUMPTIONS:
%   - All figures use print-ready colours and 1.8 pt line width.
%   - Axes labelled in mmHg, mL, or mL/s as appropriate.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

t  = sim.t(:);
XV = sim.V;

%% Extract last cycle
T_HB = 60 / params.HR;
T1   = t(end);
T0   = T1 - T_HB;
idx  = (t >= T0) & (t <= T1);
tc   = t(idx);

%% Reconstruct signals (inline helper below)
[P, Q] = signals_for_plot(tc, XV(idx,:), params);

lw = 1.8;   % line width for publication figures

%% ---- 1. Chamber pressures -------------------------------------------
figure('Name', sprintf('Chamber Pressures — %s (%s)', tag, scenario), ...
       'Color', 'w', 'NumberTitle', 'off');

subplot(2,2,1);
plot(tc, P.RA, 'b', 'LineWidth', lw);
xlabel('Time (s)'); ylabel('P_{RA} (mmHg)');
title('Right Atrium'); grid on;

subplot(2,2,2);
plot(tc, P.RV, 'r', 'LineWidth', lw);
xlabel('Time (s)'); ylabel('P_{RV} (mmHg)');
title('Right Ventricle'); grid on;

subplot(2,2,3);
plot(tc, P.LA, 'b--', 'LineWidth', lw);
xlabel('Time (s)'); ylabel('P_{LA} (mmHg)');
title('Left Atrium'); grid on;

subplot(2,2,4);
plot(tc, P.LV, 'r--', 'LineWidth', lw);
xlabel('Time (s)'); ylabel('P_{LV} (mmHg)');
title('Left Ventricle'); grid on;

sgtitle(sprintf('Chamber Pressures — %s | %s', tag, scenario));

%% ---- 2. Vascular pressures -----------------------------------------
figure('Name', sprintf('Vascular Pressures — %s', tag), ...
       'Color', 'w', 'NumberTitle', 'off');

plot(tc, P.SAR,  'k',   'LineWidth', lw); hold on;
plot(tc, P.SVEN, 'b',   'LineWidth', lw);
plot(tc, P.PAR,  'r',   'LineWidth', lw);
plot(tc, P.PVEN, 'm--', 'LineWidth', lw);
xlabel('Time (s)'); ylabel('Pressure (mmHg)');
legend('SAR','SVEN','PAR','PVEN', 'Location','best');
title(sprintf('Vascular Pressures — %s | %s', tag, scenario));
grid on;

%% ---- 3. PV loops ---------------------------------------------------
figure('Name', sprintf('PV Loops — %s', tag), ...
       'Color', 'w', 'NumberTitle', 'off');

V_LV_c = XV(idx, 4);
V_RV_c = XV(idx, 2);

subplot(1,2,1);
plot(V_LV_c, P.LV, 'r', 'LineWidth', lw);
xlabel('V_{LV} (mL)'); ylabel('P_{LV} (mmHg)');
title('LV Pressure-Volume Loop'); grid on;

subplot(1,2,2);
plot(V_RV_c, P.RV, 'b', 'LineWidth', lw);
xlabel('V_{RV} (mL)'); ylabel('P_{RV} (mmHg)');
title('RV Pressure-Volume Loop'); grid on;

sgtitle(sprintf('PV Loops — %s | %s', tag, scenario));

%% ---- 4. VSD shunt flow (pre-surgery only) ---------------------------
if strcmp(scenario, 'pre_surgery')
    figure('Name', sprintf('VSD Shunt Flow — %s', tag), ...
           'Color', 'w', 'NumberTitle', 'off');
    plot(tc, Q.VSD, 'g', 'LineWidth', lw);
    xlabel('Time (s)'); ylabel('Q_{VSD} (mL/s)');
    title(sprintf('VSD Shunt Flow (positive = L\\rightarrowR) — %s', tag));
    yline(0, 'k--', 'LineWidth', 1);
    grid on;
end

%% ---- 5. Pulmonary valve flow ----------------------------------------
figure('Name', sprintf('Pulmonary Valve Flow — %s', tag), ...
       'Color', 'w', 'NumberTitle', 'off');
plot(tc, Q.PVv, 'b', 'LineWidth', lw);
xlabel('Time (s)'); ylabel('Q_{PVv} (mL/s)');
title(sprintf('Pulmonary Valve Flow — %s | %s', tag, scenario));
yline(0, 'k--', 'LineWidth', 1);
grid on;

end  % plotting_tools

% =========================================================================
%  LOCAL HELPER
% =========================================================================

function [P, Q] = signals_for_plot(tc, XV, params)
% SIGNALS_FOR_PLOT — thin wrapper around compute_clinical_indices internals
%   Returns P and Q structs for an already-trimmed time+state window.
%   (Duplicates reconstruct_signals from compute_clinical_indices for
%    standalone use; keep in sync if elastance / valve logic changes.)

n = numel(tc);
P.RA = zeros(n,1); P.RV = zeros(n,1); P.LA = zeros(n,1); P.LV = zeros(n,1);
P.SAR= zeros(n,1); P.SVEN=zeros(n,1); P.PAR=zeros(n,1); P.PVEN=zeros(n,1);
Q.TV = zeros(n,1); Q.PVv=zeros(n,1); Q.MV=zeros(n,1); Q.AV=zeros(n,1);
Q.VSD= zeros(n,1); Q.SVEN=zeros(n,1); Q.PVEN=zeros(n,1);

for i = 1:n
    xi     = XV(i,:)';
    V_RA   = xi(1); V_RV = xi(2); V_LA = xi(3); V_LV = xi(4);
    V_SAR  = xi(5);
    V_SVEN = xi(8);  Q_SVEN_i = xi(9);
    V_PAR  = xi(10);
    V_PVEN = xi(13); Q_PVEN_i = xi(14);

    [E_LV_i, E_RV_i, E_LA_i, E_RA_i] = elastance_model(tc(i), params);

    P_RA_i = max(E_RA_i*(V_RA - params.V0.RA), -5);
    P_RV_i = E_RV_i*(V_RV - params.V0.RV);
    P_LA_i = max(E_LA_i*(V_LA - params.V0.LA), -5);
    P_LV_i = E_LV_i*(V_LV - params.V0.LV);

    P_SAR_i  = (V_SAR  - params.V0.SAR)  / params.C.SAR;
    P_SVEN_i = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);
    P_PAR_i  = (V_PAR  - params.V0.PAR)  / params.C.PAR;
    P_PVEN_i = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);

    P.RA(i) = P_RA_i; P.RV(i) = P_RV_i;
    P.LA(i) = P_LA_i; P.LV(i) = P_LV_i;
    P.SAR(i) = P_SAR_i; P.SVEN(i) = P_SVEN_i;
    P.PAR(i) = P_PAR_i; P.PVEN(i) = P_PVEN_i;

    Q.TV(i)  = valve_model(P_RA_i, P_RV_i,  params);
    Q.PVv(i) = valve_model(P_RV_i, P_PAR_i, params);
    Q.MV(i)  = valve_model(P_LA_i, P_LV_i,  params);
    Q.AV(i)  = valve_model(P_LV_i, P_SAR_i, params);
    Q.SVEN(i)= Q_SVEN_i;
    Q.PVEN(i)= Q_PVEN_i;
    Q.VSD(i) = (P_LV_i - P_RV_i) / params.R.vsd;

    % Suppress unused variable warnings
    V_PAR; %#ok<VUNUS>
end
end
