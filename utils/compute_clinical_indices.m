function metrics = compute_clinical_indices(sim, params)
% COMPUTE_CLINICAL_INDICES
% -----------------------------------------------------------------------
% Derive haemodynamic metrics from a steady-state simulation output.
%
% INPUTS:
%   sim     - struct from integrate_system.m  (.t, .V)
%   params  - parameter struct (post-scaling/calibration)
%
% OUTPUTS:
%   metrics - struct of computed indices (all in clinical reporting units)
%
% FIELDS RETURNED:
%   RAP_mean   [mmHg]   right atrial mean pressure
%   LAP_mean   [mmHg]   left atrial mean pressure
%   PAP_min    [mmHg]   PA pressure minimum
%   PAP_max    [mmHg]   PA pressure maximum
%   PAP_mean   [mmHg]   PA time-averaged mean
%   SAP_min    [mmHg]   systemic arterial minimum
%   SAP_max    [mmHg]   systemic arterial maximum
%   SAP_mean   [mmHg]   systemic arterial time-averaged mean
%   SVR        [WU]     systemic vascular resistance (Wood units)
%   PVR        [WU]     pulmonary vascular resistance (Wood units)
%   QpQs       [-]      pulmonary-to-systemic flow ratio
%   LVEDV      [mL]     LV end-diastolic volume
%   LVESV      [mL]     LV end-systolic volume
%   RVEDV      [mL]     RV end-diastolic volume
%   RVESV      [mL]     RV end-systolic volume
%   LVEF       [-]      LV ejection fraction (fraction)
%   RVEF       [-]      RV ejection fraction (fraction)
%   LVSV       [mL]     LV stroke volume
%   RVSV       [mL]     RV stroke volume
%   Q_shunt_mean_mLs  [mL/s]   mean VSD shunt flow
%
% REFERENCES:
%   [1] system_rhs.m — state layout and ODE formulation.
%   [2] Valenti (2023). Thesis. Eq. (2.5)–(2.7).
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

t  = sim.t(:);
XV = sim.V;           % [n × 14]

%% Reconstruct pressure and flow signals at each time point
[P, Q] = reconstruct_signals(t, XV, params);

%% Use last complete cardiac cycle
T_HB = 60 / params.HR;
T1   = t(end);
T0   = T1 - T_HB;
idx  = (t >= T0) & (t <= T1);

tc = t(idx);
Pc = struct();
Qc = struct();
flds_P = fieldnames(P);
flds_Q = fieldnames(Q);
for i = 1:numel(flds_P), Pc.(flds_P{i}) = P.(flds_P{i})(idx); end
for i = 1:numel(flds_Q), Qc.(flds_Q{i}) = Q.(flds_Q{i})(idx); end

mean_t = @(y) trapz(tc, y) / max(tc(end) - tc(1), 1e-9);

metrics = struct();

%% Atrial pressures
metrics.RAP_mean = mean_t(Pc.RA);
metrics.LAP_mean = mean_t(Pc.LA);

%% Pulmonary artery  (P_PAR)
metrics.PAP_min  = min(Pc.PAR);
metrics.PAP_max  = max(Pc.PAR);
metrics.PAP_mean = mean_t(Pc.PAR);

%% Systemic artery  (P_SAR)
metrics.SAP_min  = min(Pc.SAR);
metrics.SAP_max  = max(Pc.SAR);
metrics.SAP_mean = mean_t(Pc.SAR);

%% Flows  (mL/s → L/min)
Qsys_mLs  = mean_t(Qc.SVEN);    % systemic: venous return to RA
Qpul_mLs  = mean_t(Qc.PVEN);    % pulmonary: venous return to LA
Qsys_Lmin = Qsys_mLs  * 60/1000;
Qpul_Lmin = Qpul_mLs  * 60/1000;

%% Resistances  (Wood units = mmHg / [L/min])
metrics.SVR  = (metrics.SAP_mean - metrics.RAP_mean) / max(Qsys_Lmin, 1e-6);
metrics.PVR  = (metrics.PAP_mean - metrics.LAP_mean) / max(Qpul_Lmin, 1e-6);
metrics.QpQs = Qpul_Lmin / max(Qsys_Lmin, 1e-6);

%% Ventricular volumes and ejection fractions
V_LV_c = XV(idx, 4);
V_RV_c = XV(idx, 2);

metrics.LVEDV = max(V_LV_c);
metrics.LVESV = min(V_LV_c);
metrics.RVEDV = max(V_RV_c);
metrics.RVESV = min(V_RV_c);

metrics.LVEF  = (metrics.LVEDV - metrics.LVESV) / max(metrics.LVEDV, 1e-6);
metrics.RVEF  = (metrics.RVEDV - metrics.RVESV) / max(metrics.RVEDV, 1e-6);
metrics.LVSV  = metrics.LVEDV - metrics.LVESV;
metrics.RVSV  = metrics.RVEDV - metrics.RVESV;

%% VSD shunt
metrics.Q_shunt_mean_mLs = mean_t(Qc.VSD);

end  % compute_clinical_indices

% =========================================================================
%  LOCAL HELPER — reconstruct pressure and flow signals
% =========================================================================

function [P, Q] = reconstruct_signals(t, XV, params)
% RECONSTRUCT_SIGNALS — rebuild P and Q arrays from the 14-state vector
%   P.RA, P.RV, P.LA, P.LV, P.SAR, P.SVEN, P.PAR, P.PVEN [mmHg]
%   Q.TV, Q.PVv, Q.MV, Q.AV, Q.SVEN, Q.PVEN, Q.VSD       [mL/s]

n = numel(t);

P.RA   = zeros(n,1); P.RV  = zeros(n,1);
P.LA   = zeros(n,1); P.LV  = zeros(n,1);
P.SAR  = zeros(n,1); P.SVEN= zeros(n,1);
P.PAR  = zeros(n,1); P.PVEN= zeros(n,1);

Q.TV   = zeros(n,1); Q.PVv  = zeros(n,1);
Q.MV   = zeros(n,1); Q.AV   = zeros(n,1);
Q.SVEN = zeros(n,1); Q.PVEN = zeros(n,1);
Q.VSD  = zeros(n,1);

for i = 1:n
    xi     = XV(i,:)';
    V_RA   = xi(1); V_RV  = xi(2); V_LA = xi(3); V_LV = xi(4);
    V_SAR  = xi(5);
    V_SVEN = xi(8);  Q_SVEN_i = xi(9);
    V_PAR  = xi(10);
    P_PC_i = xi(12);
    V_PVEN = xi(13); Q_PVEN_i = xi(14);

    [E_LV_i, E_RV_i, E_LA_i, E_RA_i] = elastance_model(t(i), params);

    P_RA_i = max(E_RA_i * (V_RA - params.V0.RA), -5);
    P_RV_i = E_RV_i * (V_RV - params.V0.RV);
    P_LA_i = max(E_LA_i * (V_LA - params.V0.LA), -5);
    P_LV_i = E_LV_i * (V_LV - params.V0.LV);

    P_SAR_i  = (V_SAR  - params.V0.SAR)  / params.C.SAR;
    P_SVEN_i = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);
    P_PAR_i  = (V_PAR  - params.V0.PAR)  / params.C.PAR;
    P_PVEN_i = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);

    P.RA(i)   = P_RA_i; P.RV(i) = P_RV_i;
    P.LA(i)   = P_LA_i; P.LV(i) = P_LV_i;
    P.SAR(i)  = P_SAR_i; P.SVEN(i) = P_SVEN_i;
    P.PAR(i)  = P_PAR_i; P.PVEN(i) = P_PVEN_i;

    Q.TV(i)   = valve_model(P_RA_i,  P_RV_i,  params);
    Q.PVv(i)  = valve_model(P_RV_i,  P_PAR_i, params);
    Q.MV(i)   = valve_model(P_LA_i,  P_LV_i,  params);
    Q.AV(i)   = valve_model(P_LV_i,  P_SAR_i, params);
    Q.SVEN(i) = Q_SVEN_i;
    Q.PVEN(i) = Q_PVEN_i;
    Q.VSD(i)  = (P_LV_i - P_RV_i) / params.R.vsd;

    % Suppress unused variable warnings
    P_PC_i; %#ok<VUNUS>
end
end  % reconstruct_signals
