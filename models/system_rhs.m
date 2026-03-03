function dXdt = system_rhs(t, X, params)
% SYSTEM_RHS
% -----------------------------------------------------------------------
% Right-hand side of the 14-state Valenti RLC cardiovascular ODE.
%
% Implements the full lumped-parameter model for a 4-chamber heart with
% biventricular function, RLC systemic and pulmonary circuits, and an
% optional ventricular septal defect (VSD) shunt.
%
% INPUTS:
%   t       - current time                                  [s]
%   X       - state vector, 14×1                            (see below)
%   params  - parameter struct (from default_parameters.m + apply_scaling.m)
%
% OUTPUTS:
%   dXdt    - time derivative of state vector, 14×1
%
% STATE VECTOR LAYOUT (14 states, Valenti Eq. 2.7):
%   1  V_RA    Right atrial volume           [mL]
%   2  V_RV    Right ventricular volume      [mL]
%   3  V_LA    Left atrial volume            [mL]
%   4  V_LV    Left ventricular volume       [mL]
%   5  V_SAR   Systemic arterial volume      [mL]
%   6  Q_SAR   Systemic arterial flow        [mL/s]  (inductor state)
%   7  V_SC    Systemic capillary volume     [mL]
%   8  V_SVEN  Systemic venous volume        [mL]
%   9  Q_SVEN  Systemic venous flow          [mL/s]  (inductor state)
%  10  V_PAR   Pulmonary arterial volume     [mL]
%  11  Q_PAR   Pulmonary arterial flow       [mL/s]  (inductor state)
%  12  P_PC    Pulmonary capillary pressure  [mmHg]  (Valenti Eq. 2.7)
%  13  V_PVEN  Pulmonary venous volume       [mL]
%  14  Q_PVEN  Pulmonary venous flow         [mL/s]  (inductor state)
%
% GOVERNING EQUATIONS:
%   Chambers   : Valenti Eq. (2.1–2.2)   P = E(t)*(V – V0)
%   RLC compts : Valenti Eq. (2.5)       L·dQ/dt + R·Q = P_in – P_out
%                                         dV/dt = Q_in – Q_out
%   Valves     : Valenti Eq. (2.6)       non-ideal diode (open/closed R)
%   Pulm. cap  : Valenti Eq. (2.7)       dP_PC/dt = (Q_PAR – Q_COX – Q_CNO)/C_PC
%   VSD shunt  :                          Q_VSD = (P_LV – P_RV) / R_vsd
%                                         R.vsd >> large  ↔  post-surgery (closed)
%
% SIGN CONVENTIONS:
%   All flows are positive in the physiologically forward direction.
%   Q_VSD positive = LV→RV (left-to-right shunt).
%
% ASSUMPTIONS:
%   - Incompressible, Newtonian flow in all compartments.
%   - Valves are ideal diodes with two discrete resistance states.
%   - P_PC is a lumped capillary pressure state (single ODE, not two
%     separate volume ODEs) as per Valenti Eq. (2.7).
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eqs. 2.1–2.7, Table 3.3.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.1
% -----------------------------------------------------------------------

% --- Unpack state index struct (Guardrail §7.1: never hardcode indices) ---
idx    = params.idx;

V_RA   = X(idx.V_RA);    % [mL]
V_RV   = X(idx.V_RV);    % [mL]
V_LA   = X(idx.V_LA);    % [mL]
V_LV   = X(idx.V_LV);    % [mL]
V_SAR  = X(idx.V_SAR);   % [mL]
Q_SAR  = X(idx.Q_SAR);   % [mL/s]
V_SC   = X(idx.V_SC);    % [mL]
V_SVEN = X(idx.V_SVEN);  % [mL]
Q_SVEN = X(idx.Q_SVEN);  % [mL/s]
V_PAR  = X(idx.V_PAR);   % [mL]
Q_PAR  = X(idx.Q_PAR);   % [mL/s]
P_PC   = X(idx.P_PC);    % [mmHg]  pressure state (Eq. 2.7)
V_PVEN = X(idx.V_PVEN);  % [mL]
Q_PVEN = X(idx.Q_PVEN);  % [mL/s]

% --- Chamber pressures  (Eq. 2.1–2.2) ----------------------------------
[E_LV, E_RV, E_LA, E_RA] = elastance_model(t, params);

P_RA  = E_RA * (V_RA  - params.V0.RA);
P_RV  = E_RV * (V_RV  - params.V0.RV);
P_LA  = E_LA * (V_LA  - params.V0.LA);
P_LV  = E_LV * (V_LV  - params.V0.LV);

% Numerical stability guards (atrial pressures can dip briefly negative)
P_RA = max(P_RA, -5);
P_LA = max(P_LA, -5);

% --- Vascular pressures  P = (V – V0)/C  --------------------------------
P_SAR  = (V_SAR  - params.V0.SAR)  / params.C.SAR;
P_SC   = (V_SC   - params.V0.SC)   / params.C.SC;
P_SVEN = (V_SVEN - params.V0.SVEN) / params.C.SVEN;
P_PAR  = (V_PAR  - params.V0.PAR)  / params.C.PAR;
P_PVEN = (V_PVEN - params.V0.PVEN) / params.C.PVEN;

P_SVEN = max(P_SVEN, -5);
P_PVEN = max(P_PVEN, -5);

% --- Valve flows  (Eq. 2.6) — non-ideal diode ---------------------------
Q_TV  = valve_model(P_RA,  P_RV,  params);   % Tricuspid:  RA  → RV
Q_PVv = valve_model(P_RV,  P_PAR, params);   % Pulmonic:   RV  → PAR
Q_MV  = valve_model(P_LA,  P_LV,  params);   % Mitral:     LA  → LV
Q_AV  = valve_model(P_LV,  P_SAR, params);   % Aortic:     LV  → SAR

% --- VSD shunt: LV ↔ RV  ------------------------------------------------
% Q_VSD > 0  ≡  left-to-right shunt (L→R, physiologically typical for VSD)
% When params.R.vsd is very large (post-surgery), Q_VSD ≈ 0.
Q_VSD = (P_LV - P_RV) / params.R.vsd;

% =====================================================================
%  SYSTEMIC CIRCUIT  (Eq. 2.5)
%  SAR (RLC) → SC (RC) → SVEN (RLC) → RA
%  Capillary resistance split equally: R_SC/2 inlet, R_SC/2 outlet.
% =====================================================================
R_SC_half = params.R.SC / 2;

% SC→SVEN algebraic flow (outlet half of capillary)
Q_SC_out = (P_SC - P_SVEN) / R_SC_half;

% SAR inductor (Eq. 2.5) — junction pressure at SAR/SC node eliminated
dQ_SAR = (P_SAR - P_SC - Q_SAR*(params.R.SAR + R_SC_half)) / params.L.SAR;

dV_SAR  = Q_AV   - Q_SAR;
dV_SC   = Q_SAR  - Q_SC_out;

% SVEN inductor
dQ_SVEN = (P_SVEN - P_RA - params.R.SVEN*Q_SVEN) / params.L.SVEN;
dV_SVEN = Q_SC_out - Q_SVEN;

% =====================================================================
%  PULMONARY CIRCUIT  (Eq. 2.5 + 2.7)
%  PAR (RLC) → P_PC (combined capillary pressure state) → PVEN (RLC) → LA
%  Parallel capillary branches (COX, CNO) share the same P_PC upstream.
% =====================================================================
C_PC_total = params.C.PCOX + params.C.PCNO;   % combined capillary compliance

% PAR inductor
dQ_PAR = (P_PAR - P_PC - params.R.PAR*Q_PAR) / params.L.PAR;
dV_PAR = Q_PVv - Q_PAR;

% Parallel capillary outflows (both share P_PC as upstream pressure)
Q_COX = (P_PC - P_PVEN) / params.R.PCOX;
Q_CNO = (P_PC - P_PVEN) / params.R.PCNO;

% Valenti Eq. (2.7) — combined capillary pressure ODE
dP_PC = (Q_PAR - Q_COX - Q_CNO) / C_PC_total;

% PVEN inductor
dQ_PVEN = (P_PVEN - P_LA - params.R.PVEN*Q_PVEN) / params.L.PVEN;
dV_PVEN = Q_COX + Q_CNO - Q_PVEN;

% =====================================================================
%  CARDIAC CHAMBER VOLUME DERIVATIVES
%  VSD shunt: +Q_VSD into RV, -Q_VSD from LV
% =====================================================================
dV_RA = Q_SVEN  - Q_TV;
dV_RV = Q_TV    + Q_VSD  - Q_PVv;
dV_LA = Q_PVEN  - Q_MV;
dV_LV = Q_MV    - Q_AV   - Q_VSD;

% --- Assemble output via index struct (Guardrail §7.1) -------------------
dXdt = zeros(14, 1);
dXdt(idx.V_RA)   = dV_RA;
dXdt(idx.V_RV)   = dV_RV;
dXdt(idx.V_LA)   = dV_LA;
dXdt(idx.V_LV)   = dV_LV;
dXdt(idx.V_SAR)  = dV_SAR;
dXdt(idx.Q_SAR)  = dQ_SAR;
dXdt(idx.V_SC)   = dV_SC;
dXdt(idx.V_SVEN) = dV_SVEN;
dXdt(idx.Q_SVEN) = dQ_SVEN;
dXdt(idx.V_PAR)  = dV_PAR;
dXdt(idx.Q_PAR)  = dQ_PAR;
dXdt(idx.P_PC)   = dP_PC;
dXdt(idx.V_PVEN) = dV_PVEN;
dXdt(idx.Q_PVEN) = dQ_PVEN;
end
