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
%   RAP_max    [mmHg]   right atrial maximal pressure  (max P_RA over last cycle)
%   RAP_min    [mmHg]   right atrial minimal pressure  (min P_RA over last cycle)
%   RAP_mean   [mmHg]   right atrial mean pressure
%   LAP_mean   [mmHg]   left atrial mean pressure
%   PAP_min    [mmHg]   PA pressure minimum
%   PAP_max    [mmHg]   PA pressure maximum
%   PAP_mean   [mmHg]   PA time-averaged mean
%   PVP_mean   [mmHg]   pulmonary venous mean pressure
%   RVP_min    [mmHg]   RV pressure minimum
%   RVP_max    [mmHg]   RV pressure maximum
%   RVP_mean   [mmHg]   RV time-averaged mean pressure
%   LVP_min    [mmHg]   LV pressure minimum
%   LVP_max    [mmHg]   LV pressure maximum
%   LVP_mean   [mmHg]   LV time-averaged mean pressure
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
%   Q_AV_mean  [mL/s]  mean aortic valve flow (forward = LV→Ao)
%   Q_PVv_mean [mL/s]  mean pulmonary valve flow (forward = RV→PA)
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

%% Unpack state index struct (Guardrail §7.1)
sidx = params.idx;

%% Unpack unit conversion factors (Guardrail §6.3)
mLs_to_Lmin = params.conv.mLs_to_Lmin;   % [L/min per mL/s]

%% Reconstruct pressure and flow signals at each time point
[P, Q] = reconstruct_signals(t, XV, params);

%% Use last complete cardiac cycle
T_HB = 60 / params.HR;
T1   = t(end);
T0   = T1 - T_HB;
time_mask = (t >= T0) & (t <= T1);   % logical mask — renamed from 'idx' to avoid clash

tc = t(time_mask);
Pc = struct();
Qc = struct();
flds_P = fieldnames(P);
flds_Q = fieldnames(Q);
for i = 1:numel(flds_P), Pc.(flds_P{i}) = P.(flds_P{i})(time_mask); end
for i = 1:numel(flds_Q), Qc.(flds_Q{i}) = Q.(flds_Q{i})(time_mask); end

mean_t = @(y) trapz(tc, y) / max(tc(end) - tc(1), 1e-9);

metrics = struct();

%% Atrial pressures
metrics.RAP_max  = max(Pc.RA);      % [mmHg]  Right Atrial Maximal Pressure
metrics.RAP_min  = min(Pc.RA);      % [mmHg]  Right Atrial Minimal Pressure
metrics.RAP_mean = mean_t(Pc.RA);   % [mmHg]  Right Atrial Mean Pressure
metrics.LAP_mean = mean_t(Pc.LA);   % [mmHg]

%% Pulmonary artery  (P_PAR)
metrics.PAP_min  = min(Pc.PAR);     % [mmHg]
metrics.PAP_max  = max(Pc.PAR);     % [mmHg]
metrics.PAP_mean = mean_t(Pc.PAR);  % [mmHg]

%% Pulmonary venous pressure (P_PVEN)
metrics.PVP_mean = mean_t(Pc.PVEN); % [mmHg]

%% Right ventricular pressure (P_RV)
metrics.RVP_min  = min(Pc.RV);      % [mmHg]
metrics.RVP_max  = max(Pc.RV);      % [mmHg]
metrics.RVP_mean = mean_t(Pc.RV);   % [mmHg]

%% Left ventricular pressure (P_LV)
metrics.LVP_min  = min(Pc.LV);      % [mmHg]
metrics.LVP_max  = max(Pc.LV);      % [mmHg]
metrics.LVP_mean = mean_t(Pc.LV);   % [mmHg]

%% Systemic artery  (P_SAR)
metrics.SAP_min  = min(Pc.SAR);     % [mmHg]
metrics.SAP_max  = max(Pc.SAR);     % [mmHg]
metrics.SAP_mean = mean_t(Pc.SAR);  % [mmHg]

%% Flows: convert mL/s → L/min for resistance and output reporting
Qsys_mLs  = mean_t(Qc.SVEN);                   % [mL/s]  systemic venous return
Qpul_mLs  = mean_t(Qc.PVEN);                   % [mL/s]  pulmonary venous return
Qsys_Lmin = Qsys_mLs  * mLs_to_Lmin;          % [L/min]
Qpul_Lmin = Qpul_mLs  * mLs_to_Lmin;          % [L/min]

%% Resistances  (Wood units = mmHg / [L/min])
metrics.SVR  = (metrics.SAP_mean - metrics.RAP_mean) / max(Qsys_Lmin, 1e-6);   % [WU]
metrics.PVR  = (metrics.PAP_mean - metrics.LAP_mean) / max(Qpul_Lmin, 1e-6);   % [WU]
metrics.QpQs = Qpul_Lmin / max(Qsys_Lmin, 1e-6);                               % [-]

%% Ventricular volumes and ejection fractions
%  Use params.idx to access state columns — never hardcode column numbers (Guardrail §7.1)
V_LV_c = XV(time_mask, sidx.V_LV);   % [mL]
V_RV_c = XV(time_mask, sidx.V_RV);   % [mL]

metrics.LVEDV = max(V_LV_c);     % [mL]
metrics.LVESV = min(V_LV_c);     % [mL]
metrics.RVEDV = max(V_RV_c);     % [mL]
metrics.RVESV = min(V_RV_c);     % [mL]

% EF stored as fraction (0–1), NOT percentage (Guardrail §6.4)
metrics.LVEF  = (metrics.LVEDV - metrics.LVESV) / max(metrics.LVEDV, 1e-6);   % [-]
metrics.RVEF  = (metrics.RVEDV - metrics.RVESV) / max(metrics.RVEDV, 1e-6);   % [-]
metrics.LVSV  = metrics.LVEDV - metrics.LVESV;   % [mL]
metrics.RVSV  = metrics.RVEDV - metrics.RVESV;   % [mL]

%% VSD shunt
metrics.Q_shunt_mean_mLs = mean_t(Qc.VSD);   % [mL/s]  positive = L→R

%% Valve flows  (time-averaged over last cardiac cycle)
%   Positive = forward physiological direction (Guardrail §3.4 sign convention)
metrics.Q_AV_mean  = mean_t(Qc.AV);    % [mL/s]  aortic valve flow  (LV→Ao)
metrics.Q_PVv_mean = mean_t(Qc.PVv);   % [mL/s]  pulmonary valve flow (RV→PA)

end  % compute_clinical_indices

% =========================================================================
%  LOCAL HELPER — reconstruct pressure and flow signals
% =========================================================================

function [P, Q] = reconstruct_signals(t, XV, params)
% RECONSTRUCT_SIGNALS — rebuild P and Q arrays from the 14-state vector
%   P.RA, P.RV, P.LA, P.LV, P.SAR, P.SVEN, P.PAR, P.PVEN [mmHg]
%   Q.TV, Q.PVv, Q.MV, Q.AV, Q.SVEN, Q.PVEN, Q.VSD       [mL/s]

n    = numel(t);
sidx = params.idx;   % state index struct (Guardrail §7.1)

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
    V_RA   = xi(sidx.V_RA);    % [mL]
    V_RV   = xi(sidx.V_RV);    % [mL]
    V_LA   = xi(sidx.V_LA);    % [mL]
    V_LV   = xi(sidx.V_LV);    % [mL]
    V_SAR  = xi(sidx.V_SAR);   % [mL]
    V_SVEN = xi(sidx.V_SVEN);  % [mL]
    Q_SVEN_i = xi(sidx.Q_SVEN); % [mL/s]
    V_PAR  = xi(sidx.V_PAR);   % [mL]
    P_PC_i = xi(sidx.P_PC);    % [mmHg]  (unused in pressure reconstruction; suppress warning)
    V_PVEN = xi(sidx.V_PVEN);  % [mL]
    Q_PVEN_i = xi(sidx.Q_PVEN); % [mL/s]

    [E_LV_i, E_RV_i, E_LA_i, E_RA_i] = elastance_model(t(i), params);

    P_RA_i = max(E_RA_i * (V_RA - params.V0.RA), -5);    % [mmHg]
    P_RV_i = E_RV_i * (V_RV - params.V0.RV);              % [mmHg]
    P_LA_i = max(E_LA_i * (V_LA - params.V0.LA), -5);    % [mmHg]
    P_LV_i = E_LV_i * (V_LV - params.V0.LV);              % [mmHg]

    P_SAR_i  = (V_SAR  - params.V0.SAR)  / params.C.SAR;             % [mmHg]
    P_SVEN_i = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);   % [mmHg]
    P_PAR_i  = (V_PAR  - params.V0.PAR)  / params.C.PAR;             % [mmHg]
    P_PVEN_i = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);   % [mmHg]

    P.RA(i)   = P_RA_i; P.RV(i) = P_RV_i;
    P.LA(i)   = P_LA_i; P.LV(i) = P_LV_i;
    P.SAR(i)  = P_SAR_i; P.SVEN(i) = P_SVEN_i;
    P.PAR(i)  = P_PAR_i; P.PVEN(i) = P_PVEN_i;

    Q.TV(i)   = valve_model(P_RA_i,  P_RV_i,  params);   % [mL/s]
    Q.PVv(i)  = valve_model(P_RV_i,  P_PAR_i, params);   % [mL/s]
    Q.MV(i)   = valve_model(P_LA_i,  P_LV_i,  params);   % [mL/s]
    Q.AV(i)   = valve_model(P_LV_i,  P_SAR_i, params);   % [mL/s]
    Q.SVEN(i) = Q_SVEN_i;                                 % [mL/s]
    Q.PVEN(i) = Q_PVEN_i;                                 % [mL/s]
    Q.VSD(i)  = (P_LV_i - P_RV_i) / params.R.vsd;        % [mL/s] positive = L→R

    P_PC_i; %#ok<VUNUS>  % P_PC is a state but not re-derived from volumes here
end
end  % reconstruct_signals
