function params = default_parameters()
% DEFAULT_PARAMETERS
% -----------------------------------------------------------------------
% Adult reference parameter set for the 14-state cardiovascular ODE.
% Healthy adult male at rest.
%
% PARAMETER SOURCES (2026-04-07 update):
%   Bozkurt (2019)   — vascular R, C, L; atrial elastance; timing
%   Colebank (2025)  — LV/RV elastance; initial volumes; timing
%   Corrections applied per 2025-04 correction pass (see comments).
%
% TOPOLOGY: 6-segment vascular (Bozkurt2019)
%   Systemic : SAR(aortic RLC) → SC(arteriolar RC) → SVEN(venous RLC)
%   Pulmonary: PAR(arterial RLC) → PCOX(arteriolar RC, single branch) → PVEN(venous RLC)
%   PCNO branch DISABLED (R=1e6, C=1e-8): Bozkurt2019 uses a single
%   arteriolar compartment, not the dual PCOX/PCNO of Valenti 2023.
%
% ELASTANCE MAPPING:
%   Valenti:    E_ch(t) = EA*e(t) + EB
%   Bozkurt:    Ees (end-systolic / peak) and Emin (diastolic)
%   Conversion: EA = Ees - Emin  (active contractile component)
%               EB = Emin        (passive diastolic stiffness)
%
% TIMING MAPPING (from Colebank2025 Tmax_v / Tmin_v):
%   Tc_frac = Tmax_v/T = 0.25  (contraction duration to peak)
%   Tr_frac = (Tmin_v-Tmax_v)/T = 0.30  (relaxation duration)
%
% UNITS:
%   Pressure    [mmHg]      Volume      [mL]
%   Flow        [mL/s]      Resistance  [mmHg·s/mL]
%   Compliance  [mL/mmHg]   Inertance   [mmHg·s²/mL]
%
% STATE VECTOR LAYOUT (14 states):
%   1  V_RA    Right atrial volume           [mL]
%   2  V_RV    Right ventricular volume      [mL]
%   3  V_LA    Left atrial volume            [mL]
%   4  V_LV    Left ventricular volume       [mL]
%   5  V_SAR   Systemic arterial volume      [mL]
%   6  Q_SAR   Systemic arterial flow        [mL/s]
%   7  V_SC    Systemic arteriolar volume    [mL]
%   8  V_SVEN  Systemic venous volume        [mL]
%   9  Q_SVEN  Systemic venous flow          [mL/s]
%  10  V_PAR   Pulmonary arterial volume     [mL]
%  11  Q_PAR   Pulmonary arterial flow       [mL/s]
%  12  P_PC    Pulmonary cap. pressure       [mmHg]
%  13  V_PVEN  Pulmonary venous volume       [mL]
%  14  Q_PVEN  Pulmonary venous flow         [mL/s]
%
% REFERENCES:
%   [1] Bozkurt S (2019). Mathematical modeling of cardiac dysfunction and
%       treatment. Math Biosci Eng 16(5):3943–3962. Tables 1–2.
%   [2] Colebank MJ et al. (2025). Sex-specific cardiovascular modelling.
%       ASAIO J. Table 1 (male reference).
%   [3] Valenti (thesis, 2023). ODE topology and state-vector layout.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-07
% VERSION:  2.0  (Bozkurt2019/Colebank2025 baseline)
% -----------------------------------------------------------------------

params = struct();

%% =====================================================================
%  STATE VECTOR INDEX STRUCT  (Guardrail §7.1 — never hardcode indices)
%% =====================================================================
idx.V_RA   = 1;    % Right atrial volume              [mL]
idx.V_RV   = 2;    % Right ventricular volume         [mL]
idx.V_LA   = 3;    % Left atrial volume               [mL]
idx.V_LV   = 4;    % Left ventricular volume          [mL]
idx.V_SAR  = 5;    % Systemic arterial (aortic) vol   [mL]
idx.Q_SAR  = 6;    % Systemic arterial flow           [mL/s]
idx.V_SC   = 7;    % Systemic arteriolar volume       [mL]
idx.V_SVEN = 8;    % Systemic venous volume           [mL]
idx.Q_SVEN = 9;    % Systemic venous flow             [mL/s]
idx.V_PAR  = 10;   % Pulmonary arterial volume        [mL]
idx.Q_PAR  = 11;   % Pulmonary arterial flow          [mL/s]
idx.P_PC   = 12;   % Pulmonary capillary pressure     [mmHg]
idx.V_PVEN = 13;   % Pulmonary venous volume          [mL]
idx.Q_PVEN = 14;   % Pulmonary venous flow            [mL/s]
params.idx = idx;

%% Simulation control
params.sim.nCyclesSteady = 80;    % cardiac cycles to reach steady-state
params.sim.nCyclesKeep   = 2;     % last N cycles retained for post-processing
params.sim.rtol          = 1e-7;  % ode15s relative tolerance
params.sim.atol          = 1e-8;  % ode15s absolute tolerance
params.sim.ss_tol_P      = 0.1;   % [mmHg]  steady-state pressure tolerance
params.sim.ss_tol_V      = 0.1;   % [mL]    steady-state volume tolerance

%% Unit conversion constants
params.conv.WU_to_R     = 60/1000;   % Wood units → mmHg·s/mL
params.conv.R_to_WU     = 1000/60;   % mmHg·s/mL → Wood units
params.conv.mLs_to_Lmin = 60/1000;  % mL/s → L/min
params.conv.Lmin_to_mLs = 1000/60;  % L/min → mL/s
params.conv.mmHg_to_Pa  = 133.322;  % mmHg → Pa
params.conv.mm_to_m     = 1e-3;     % mm → m
params.conv.m3_to_mL    = 1e6;      % m³ → mL

%% Valve epsilon (not used by hard-switch valve_model v2.0; stored for reference)
params.epsilon_valve = 0.5;   % [mmHg]  retained for backward compatibility

%% =====================================================================
%  GLOBAL  —  Bozkurt2019 / Colebank2025
%% =====================================================================
params.HR = 75;   % [bpm]  resting heart rate — Bozkurt2019 / Colebank2025
%  T  = 0.800 s   CO_ref = 100 mL/s = 6.0 L/min   SV_ref = 80 mL

%% =====================================================================
%  CARDIAC PHASE TIMING  (stored as fractions of T_HB)
%  Converted to absolute seconds by apply_scaling.m Section D.
%
%  Ventricular: Tmax_v = 0.25·T (peak), Tmin_v = 0.55·T (end active)
%    → Tc_frac = Tmax_v/T = 0.25   Tr_frac = (Tmin_v-Tmax_v)/T = 0.30
%  Atrial: onset Ta = 0.80·T, peak at Ta+Tmax_a, end at Ta+Tmax_a+Tmin_a
%    → t_ac_frac = 0.80, Tc_frac = 0.10, Tr_frac = 0.20
%  Source: Bozkurt2019 Table 1; Colebank2025
%% =====================================================================

% Ventricular timing
params.t_ac_LV_frac = 0;       % LV contraction start / T_HB  (start of beat)
params.Tc_LV_frac   = 0.25;    % LV contraction duration / T_HB  — Colebank2025
params.Tr_LV_frac   = 0.30;    % LV relaxation  duration / T_HB  — Colebank2025
params.t_ac_RV_frac = 0;       % RV contraction start / T_HB
params.Tc_RV_frac   = 0.25;    % RV contraction duration / T_HB  — Colebank2025
params.Tr_RV_frac   = 0.30;    % RV relaxation  duration / T_HB  — Colebank2025

% Atrial timing  (T_ar = t_ac + Tc; set in apply_scaling Section D)
params.t_ac_LA_frac = 0.85;    % LA contraction start   / T_HB  — Bozkurt2019
params.Tc_LA_frac   = 0.10;    % LA contraction duration/ T_HB  — Colebank2025
params.Tr_LA_frac   = 0.10;    % LA relaxation  duration/ T_HB  — Colebank2025
params.t_ac_RA_frac = 0.85;    % RA contraction start   / T_HB  — Bozkurt2019
params.Tc_RA_frac   = 0.10;    % RA contraction duration/ T_HB  — Colebank2025
params.Tr_RA_frac   = 0.10;    % RA relaxation  duration/ T_HB  — Colebank2025

% Additional Bozkurt2019 timing constants (stored for documentation)
params.timing.D    = 0.04;     % [s]  AV conduction delay     — Bozkurt2019
params.timing.Toff = 0.15*0.8; % [s]  atrial-ventricular offset — Colebank2025

%% =====================================================================
%  CARDIAC CHAMBERS  —  Bozkurt2019 Table 1 / Colebank2025
%
%  Elastance mapping:  EA = Ees - Emin   EB = Emin
%  P_ch(t) = (EA*e(t) + EB) * (V_ch - V0_ch)
%% =====================================================================

%-- Left ventricle
% Ees_lv = 4.80 mmHg/mL  CORRECTED: P_lv_sys/(LVESV-V0) = 120/25
% Emin_lv = 0.0476 mmHg/mL  CORRECTED: P_lv_dia/(LVEDV-V0) = 5/105
params.E.LV.EA = 4.80 - 0.0476;  % [mmHg/mL]  = 4.7524  — Colebank2025 (corrected)
params.E.LV.EB = 0.0476;          % [mmHg/mL]             — Colebank2025 (corrected)
params.V0.LV   = 15;              % [mL]  unstressed volume — Bozkurt2019 Table 1

%-- Right ventricle
% Ees_rv = 0.525 mmHg/mL  CORRECTED: P_rv_sys/(RVESV-V0) = 21/40
% Emin_rv = 0.0375 mmHg/mL  CORRECTED: P_rv_dia/(RVEDV-V0) = 3/120
params.E.RV.EA = 0.525 - 0.0375; % [mmHg/mL]  = 0.4875  — Colebank2025 (corrected)
params.E.RV.EB = 0.0375;          % [mmHg/mL]             — Colebank2025 (corrected)
params.V0.RV   = 40;              % [mL]  unstressed volume — Bozkurt2019 Table 1

%-- Left atrium
% Emax_la = 0.30, Emin_la = 0.20  — Bozkurt2019 Table 1
params.E.LA.EA = 0.30 - 0.20;    % [mmHg/mL]  = 0.10  — Bozkurt2019
params.E.LA.EB = 0.20;            % [mmHg/mL]          — Bozkurt2019
params.V0.LA   = 5;               % [mL]               — Bozkurt2019 Table 1

%-- Right atrium
% Emax_ra = 0.30, Emin_ra = 0.20  — Bozkurt2019 Table 1
params.E.RA.EA = 0.30 - 0.20;    % [mmHg/mL]  = 0.10  — Bozkurt2019
params.E.RA.EB = 0.20;            % [mmHg/mL]          — Bozkurt2019
params.V0.RA   = 5;               % [mL]               — Bozkurt2019 Table 1

%-- Atrial compliance (legacy field; derived from EB)
params.C.RA = 1/0.20;   % [mL/mmHg]  = 5.0  — from EB_RA
params.C.LA = 1/0.20;   % [mL/mmHg]  = 5.0  — from EB_LA

%% =====================================================================
%  SYSTEMIC CIRCUIT  —  Bozkurt2019 Table 2
%  SAR (aortic RLC) → SC (arteriolar RC) → SVEN (venous RLC) → RA
%% =====================================================================

% SAR — aortic segment (maps to Rao, Cao, Lao)
params.R.SAR  = 0.050;    % [mmHg·s/mL]  aortic resistance    — Bozkurt2019
params.C.SAR  = 0.200;    % [mL/mmHg]    aortic compliance     — Bozkurt2019
params.L.SAR  = 1e-5;     % [mmHg·s²/mL] aortic inertance      — Bozkurt2019

% SC — systemic arteriolar segment (maps to Ras, Cas)
% NOTE: in Bozkurt2019 the arteriolar compartment (as) carries most SVR
params.R.SC   = 0.730;    % [mmHg·s/mL]  arteriolar resistance — Bozkurt2019 (corrected)
params.C.SC   = 1.700;    % [mL/mmHg]    arteriolar compliance  — Bozkurt2019

% SVEN — systemic venous segment (maps to Rvs, Cvs)
params.R.SVEN = 0.050;    % [mmHg·s/mL]  venous resistance     — Bozkurt2019
params.C.SVEN = 30.000;   % [mL/mmHg]    venous compliance      — Bozkurt2019
params.L.SVEN = 1e-5;     % [mmHg·s²/mL] venous inertance (kept for ODE structure)

% Systemic unstressed volumes [mL]
% Derived: V0 = V_ic - P_nom × C  using physiological IC pressures
%   P_SAR_ic = 80 mmHg (aortic diastolic)
%   P_SC_ic  = 20 mmHg (arteriolar)
%   P_SVEN_ic =  5 mmHg (CVP)
params.V0.SAR  = 202.5 - 80 * 0.200;   % = 186.5 mL — Colebank2025 VAo0
params.V0.SC   =  50   - 20 * 1.700;   % =  16.0 mL — estimated
params.V0.SVEN = 750.0 -  5 * 30.000;  % = 600.0 mL — Colebank2025 VSV0

%% =====================================================================
%  PULMONARY CIRCUIT  —  Bozkurt2019 Table 2
%  PAR (arterial RLC) → PCOX (arteriolar RC) → PVEN (venous RLC) → LA
%
%  *** PCNO BRANCH DISABLED ***
%  Bozkurt2019 uses a single arteriolar compartment (ap), not the dual
%  PCOX/PCNO of Valenti2023.  PCNO is disabled by setting R=1e6 and C=1e-8
%  so it carries no flow and holds no volume.  The ODE structure is
%  preserved unchanged (system_rhs.m line 137 still sums C.PCOX + C.PCNO).
%% =====================================================================

% PAR — pulmonary arterial segment (maps to Rpa, Cpa, Lpa)
params.R.PAR  = 0.010;    % [mmHg·s/mL]  PA resistance        — Bozkurt2019
params.C.PAR  = 5.000;    % [mL/mmHg]    PA compliance         — Bozkurt2019
params.L.PAR  = 1e-5;     % [mmHg·s²/mL] PA inertance          — Bozkurt2019

% PCOX — pulmonary arteriolar (maps to Rap, Cap) — ACTIVE branch
params.R.PCOX = 1.7538e-2;    % [mmHg·s/mL]  arteriolar resistance — Bozkurt2019 (corrected)
params.C.PCOX = 5.7803;    % [mL/mmHg]    arteriolar compliance  — Bozkurt2019

% PCNO — DISABLED (Bozkurt2019 has single arteriolar compartment)
params.R.PCNO = 3.5174e-1;      % [mmHg·s/mL]  >> infinite → Q_CNO ≈ 0
params.C.PCNO = 4.9043e-2;     % [mL/mmHg]    ≈ 0 → no volume contribution

% PVEN — pulmonary venous segment (maps to Rvp, Cvp)
params.R.PVEN = 0.020;    % [mmHg·s/mL]  PV resistance         — Bozkurt2019 (corrected)
params.C.PVEN = 30.000;   % [mL/mmHg]    PV compliance          — Bozkurt2019
params.L.PVEN = 1e-5;     % [mmHg·s²/mL] PV inertance (ODE structure)

% Pulmonary unstressed volumes [mL]
%   P_PAR_ic = 15 mmHg (PA diastolic)
%   P_PC_ic  = 10 mmHg (capillary)
%   P_PVEN_ic =  8 mmHg (pulmonary venous / PAWP)
params.V0.PAR  = 261.0 - 15 * 5.000;    % = 186.0 mL — Colebank2025 VPA0
params.V0.PCOX =  22.0 - 10 * 0.200;    % =  20.0 mL — estimated
params.V0.PCNO = 0;                      % [mL]  disabled branch
%   NOTE: VPV0 = 55 mL (Colebank2025) + Cvp = 30 mL/mmHg gives
%   P_PVEN_ic = (55-V0)/30.  Setting V0 = 0 gives P = 1.8 mmHg initially;
%   the simulation will converge to the physiological ~8 mmHg steady state.
params.V0.PVEN = 0;       % [mL]  see note above — Colebank2025 VPV0 = 55 mL

%% =====================================================================
%  VALVE RESISTANCES  —  Valenti (2023) non-ideal diode model
%
%  valve_model.m implements the Valenti two-resistance diode:
%    R_i(p1, p2) = { R_min,  p1 > p2   (valve open  / forward flow)
%                  { R_max,  p1 ≤ p2   (valve closed / back-pressure)
%
%  R_min (open)  = 6.2872e-3 mmHg·s/mL  — Valenti (2023) Table 3.3
%  R_max (closed)= 9.4168e+4 mmHg·s/mL  — Valenti (2023) Table 3.3
%
%  A single R_min is applied to all four cardiac valves (aortic, mitral,
%  pulmonary, tricuspid). Per-valve breakdown from Bozkurt2019 is retained
%  below for documentation and future per-valve extension.
%
%  Per-valve reference values (Bozkurt2019 Table 2):
%    Ra_val = 0.002 mmHg·s/mL  (aortic)
%    Rm_val = 0.002 mmHg·s/mL  (mitral)
%    Rp_val = 0.001 mmHg·s/mL  (pulmonary)
%    Rt_val = 0.001 mmHg·s/mL  (tricuspid)
%% =====================================================================
params.Rvalve.Ra_val  = 0.002;       % [mmHg·s/mL]  aortic valve      — Bozkurt2019 Table 2
params.Rvalve.Rm_val  = 0.002;       % [mmHg·s/mL]  mitral valve      — Bozkurt2019 Table 2
params.Rvalve.Rp_val  = 0.001;       % [mmHg·s/mL]  pulmonary valve   — Bozkurt2019 Table 2
params.Rvalve.Rt_val  = 0.001;       % [mmHg·s/mL]  tricuspid valve   — Bozkurt2019 Table 2
params.Rvalve.open    = 6.2872e-3;   % [mmHg·s/mL]  R_min — Valenti (2023) Table 3.3
params.Rvalve.closed  = 9.4168e+4;   % [mmHg·s/mL]  R_max — Valenti (2023) Table 3.3

%% =====================================================================
%  VSD SHUNT  (adult reference: effectively closed / absent)
%% =====================================================================
params.R.vsd = 1e6;   % [mmHg·s/mL]  VSD absent → effectively infinite resistance

%% =====================================================================
%  INITIAL CONDITIONS — 14-state vector  (Colebank2025 Table 1 male)
%  [V_RA V_RV V_LA V_LV | V_SAR Q_SAR V_SC V_SVEN Q_SVEN |
%   V_PAR Q_PAR P_PC V_PVEN Q_PVEN]
%
%  CO_ref = 100 mL/s used for initial flow states
%% =====================================================================
CO_ref = 100;   % [mL/s]  HR*SV/60 = 75*80/60 = 100

ic = zeros(14, 1);
ic(idx.V_RA)   =  25.0;   % [mL]  RA volume          — Colebank2025
ic(idx.V_RV)   = 160.0;   % [mL]  RV volume (= RVEDV) — Colebank2025 (RVESV+SV)
ic(idx.V_LA)   =  25.0;   % [mL]  LA volume          — Colebank2025
ic(idx.V_LV)   = 120.0;   % [mL]  LV volume (= LVEDV) — Colebank2025
ic(idx.V_SAR)  = 202.5;   % [mL]  Systemic arterial   — Colebank2025 VAo0
ic(idx.Q_SAR)  = CO_ref;  % [mL/s] aortic flow ≈ CO  — Colebank2025
ic(idx.V_SC)   =  50.0;   % [mL]  Arteriolar volume   — estimated (P_as≈20 mmHg)
ic(idx.V_SVEN) = 750.0;   % [mL]  Systemic venous     — Colebank2025 VSV0
ic(idx.Q_SVEN) = CO_ref;  % [mL/s] venous return ≈ CO
ic(idx.V_PAR)  = 261.0;   % [mL]  Pulmonary arterial  — Colebank2025 VPA0
ic(idx.Q_PAR)  = CO_ref;  % [mL/s] PA flow ≈ CO
ic(idx.P_PC)   =  10.0;   % [mmHg] Pulmonary cap pressure — estimated midpoint
ic(idx.V_PVEN) =  55.0;   % [mL]  Pulmonary venous    — Colebank2025 VPV0
ic(idx.Q_PVEN) = CO_ref;  % [mL/s] pulmonary venous flow ≈ CO
params.ic.V = ic';

end
