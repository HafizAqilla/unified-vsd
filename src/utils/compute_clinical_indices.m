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
%   RAP_min    [mmHg]   right atrial diastolic-like minimum
%   RAP_max    [mmHg]   right atrial systolic-like maximum
%   LAP_mean   [mmHg]   left atrial mean pressure
%   LAP_min    [mmHg]   left atrial diastolic-like minimum
%   LAP_max    [mmHg]   left atrial systolic-like maximum
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
%   CO_Lmin    [L/min]  cardiac output from LV stroke volume and HR
%   Qp_mean_mLs [mL/s]  mean pulmonary flow proxy (Q_PVEN)
%   Qs_mean_mLs [mL/s]  mean systemic flow proxy (Q_SVEN)
%   Q_shunt_mean_mLs  [mL/s]   mean VSD shunt flow
%   VSD_frac_pct [%]   shunt fraction relative to pulmonary flow
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
metrics.RAP_mean = mean_t(Pc.RA);   % [mmHg]
metrics.RAP_min  = min(Pc.RA);      % [mmHg]
metrics.RAP_max  = max(Pc.RA);      % [mmHg]
metrics.LAP_mean = mean_t(Pc.LA);   % [mmHg]
metrics.LAP_min  = min(Pc.LA);      % [mmHg]
metrics.LAP_max  = max(Pc.LA);      % [mmHg]

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
metrics.Qs_mean_mLs = Qsys_mLs;               % [mL/s]
metrics.Qp_mean_mLs = Qpul_mLs;               % [mL/s]

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
metrics.VSD_frac_pct = 100 * metrics.Q_shunt_mean_mLs / max(metrics.Qp_mean_mLs, 1e-6); % [%]

%% Cardiac output (LV-derived)
metrics.CO_Lmin = metrics.LVSV * params.HR / 1000;   % [L/min]

end  % compute_clinical_indices

% =========================================================================
%  LOCAL HELPER — reconstruct pressure and flow signals
% =========================================================================

function [P, Q] = reconstruct_signals(t, XV, params)
% RECONSTRUCT_SIGNALS — rebuild P and Q arrays from the 14-state vector
%   P.RA, P.RV, P.LA, P.LV, P.SAR, P.SVEN, P.PAR, P.PVEN [mmHg]
%   Q.TV, Q.PVv, Q.MV, Q.AV, Q.SVEN, Q.PVEN, Q.VSD       [mL/s]

sidx = params.idx;   % state index struct (Guardrail §7.1)

V_RA   = XV(:, sidx.V_RA);      % [mL]
V_RV   = XV(:, sidx.V_RV);      % [mL]
V_LA   = XV(:, sidx.V_LA);      % [mL]
V_LV   = XV(:, sidx.V_LV);      % [mL]
V_SAR  = XV(:, sidx.V_SAR);     % [mL]
V_SVEN = XV(:, sidx.V_SVEN);    % [mL]
Q_SVEN = XV(:, sidx.Q_SVEN);    % [mL/s]
V_PAR  = XV(:, sidx.V_PAR);     % [mL]
V_PVEN = XV(:, sidx.V_PVEN);    % [mL]
Q_PVEN = XV(:, sidx.Q_PVEN);    % [mL/s]

[E_LV, E_RV, E_LA, E_RA] = elastance_model(t, params);

P_RA   = max(E_RA .* (V_RA - params.V0.RA), -5);            % [mmHg]
P_RV   = E_RV .* (V_RV - params.V0.RV);                     % [mmHg]
P_LA   = max(E_LA .* (V_LA - params.V0.LA), -5);            % [mmHg]
P_LV   = E_LV .* (V_LV - params.V0.LV);                     % [mmHg]

P_SAR  = (V_SAR  - params.V0.SAR)  ./ params.C.SAR;         % [mmHg]
P_SVEN = max((V_SVEN - params.V0.SVEN) ./ params.C.SVEN, -5); % [mmHg]
P_PAR  = (V_PAR  - params.V0.PAR)  ./ params.C.PAR;         % [mmHg]
P_PVEN = max((V_PVEN - params.V0.PVEN) ./ params.C.PVEN, -5); % [mmHg]

P.RA   = P_RA;   P.RV = P_RV;
P.LA   = P_LA;   P.LV = P_LV;
P.SAR  = P_SAR;  P.SVEN = P_SVEN;
P.PAR  = P_PAR;  P.PVEN = P_PVEN;

Q.TV   = valve_model(P_RA,  P_RV,  params);                 % [mL/s]
Q.PVv  = valve_model(P_RV,  P_PAR, params);                 % [mL/s]
Q.MV   = valve_model(P_LA,  P_LV,  params);                 % [mL/s]
Q.AV   = valve_model(P_LV,  P_SAR, params);                 % [mL/s]
Q.SVEN = Q_SVEN;                                             % [mL/s]
Q.PVEN = Q_PVEN;                                             % [mL/s]
Q.VSD  = vsd_shunt_model(P_LV, P_RV, params);               % [mL/s] positive = L→R
end  % reconstruct_signals
