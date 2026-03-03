function params = default_parameters()
% DEFAULT_PARAMETERS
% -----------------------------------------------------------------------
% Adult reference parameter set for the 14-state Valenti RLC cardiovascular
% model.  All values are healthy adult male at rest (BSA ~ 1.73 m²).
% This struct is the 'params_ref' argument for apply_scaling.m.
%
% INPUTS:   none
%
% OUTPUTS:
%   params  - complete parameter struct ready for apply_scaling.m
%
% UNITS:
%   Pressure    [mmHg]
%   Volume      [mL]
%   Flow        [mL/s]
%   Resistance  [mmHg·s/mL]
%   Compliance  [mL/mmHg]
%   Inertance   [mmHg·s²/mL]
%
% STATE VECTOR LAYOUT (14 states, Valenti Eq. 2.7):
%   1  V_RA    Right atrial volume           [mL]
%   2  V_RV    Right ventricular volume      [mL]
%   3  V_LA    Left atrial volume            [mL]
%   4  V_LV    Left ventricular volume       [mL]
%   5  V_SAR   Systemic arterial volume      [mL]
%   6  Q_SAR   Systemic arterial flow        [mL/s]
%   7  V_SC    Systemic capillary volume     [mL]
%   8  V_SVEN  Systemic venous volume        [mL]
%   9  Q_SVEN  Systemic venous flow          [mL/s]
%  10  V_PAR   Pulmonary arterial volume     [mL]
%  11  Q_PAR   Pulmonary arterial flow       [mL/s]
%  12  P_PC    Pulmonary capillary pressure  [mmHg]  (Valenti Eq. 2.7)
%  13  V_PVEN  Pulmonary venous volume       [mL]
%  14  Q_PVEN  Pulmonary venous flow         [mL/s]
%
% REFERENCES:
%   [1] Valenti (thesis, 2023). Table 3.3 — adult reference parameters.
%   [2] Lundquist et al. (2025). Allometric Scaling in Pediatric Cardiology.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

params = struct();

%% =====================================================================
%  STATE VECTOR INDEX STRUCT
%  Define once here; pass inside params to all functions.
%  Never hardcode numeric indices elsewhere in the codebase (Guardrail §7.1).
%
%  Units: volumes [mL], flows [mL/s], pressures [mmHg]
%% =====================================================================
idx.V_RA   = 1;    % Right atrial volume              [mL]
idx.V_RV   = 2;    % Right ventricular volume         [mL]
idx.V_LA   = 3;    % Left atrial volume               [mL]
idx.V_LV   = 4;    % Left ventricular volume          [mL]
idx.V_SAR  = 5;    % Systemic arterial volume         [mL]
idx.Q_SAR  = 6;    % Systemic arterial flow           [mL/s]
idx.V_SC   = 7;    % Systemic capillary volume        [mL]
idx.V_SVEN = 8;    % Systemic venous volume           [mL]
idx.Q_SVEN = 9;    % Systemic venous flow             [mL/s]
idx.V_PAR  = 10;   % Pulmonary arterial volume        [mL]
idx.Q_PAR  = 11;   % Pulmonary arterial flow          [mL/s]
idx.P_PC   = 12;   % Pulmonary capillary pressure     [mmHg]
idx.V_PVEN = 13;   % Pulmonary venous volume          [mL]
idx.Q_PVEN = 14;   % Pulmonary venous flow            [mL/s]
params.idx = idx;

%% Simulation control
params.sim.nCyclesSteady = 80;     % cardiac cycles integrated to reach steady-state
                                   % Set to 80: high-flow VSD (Qp/Qs~3.44, HR=150 bpm)
                                   % requires ~60-80 cycles to converge pulmonary volumes.
                                   % At 0.4s/cycle, 80 cycles = 32s wall time per simulation.
                                   % 40 required for small paediatric s<<1 scale factors
params.sim.nCyclesKeep   = 2;      % last N cycles retained for post-processing
params.sim.rtol          = 1e-7;   % ode15s relative tolerance
params.sim.atol          = 1e-8;   % ode15s absolute tolerance
% Convergence criterion: max peak-value change between last two cycles
%   < 0.1 mmHg (pressures)  and  < 0.1 mL (volumes)  → steady state reached
params.sim.ss_tol_P      = 0.1;    % [mmHg]  steady-state pressure tolerance
params.sim.ss_tol_V      = 0.1;    % [mL]    steady-state volume tolerance

%% Patient demographics (overwritten by apply_scaling / params_from_clinical)
params.HR = 75;   % [bpm]  adult resting heart rate — Source: Valenti Table 3.3

%% Unit conversion constants  (Section 6 — never scatter raw numbers in physics code)
params.conv.WU_to_R    = 60/1000;   % Wood units -> mmHg·s/mL  (1 WU = 1 mmHg/(L/min))
params.conv.R_to_WU    = 1000/60;   % mmHg·s/mL -> Wood units
params.conv.mLs_to_Lmin = 60/1000;  % mL/s -> L/min
params.conv.Lmin_to_mLs = 1000/60;  % L/min -> mL/s
params.conv.mmHg_to_Pa = 133.322;   % mmHg -> Pa  (used in Gorlin orifice eq.)
params.conv.mm_to_m    = 1e-3;      % mm   -> m   (length)
params.conv.m3_to_mL   = 1e6;       % m³   -> mL  (volume, from SI flow to mL/s)

%% Valve soft-switching parameter  (Guardrail §8.4)
%   epsilon_valve [mmHg]: width of the smooth transition zone.
%   Larger = smoother but less physiologically sharp.
%   Source: tuned to prevent solver step-size collapse; see valve_model.m.
params.epsilon_valve = 0.5;         % [mmHg]

%% =====================================================================
%  CARDIAC CHAMBERS  —  Valenti Table 3.3
%% =====================================================================

%-- Unstressed volumes [mL]  — Source: Valenti Table 3.3
params.V0.RA = 3.5385;    % [mL]  Right atrial unstressed volume    — Source: Valenti Table 3.3
params.V0.RV = 8.4067;    % [mL]  Right ventricular unstressed V    — Source: Valenti Table 3.3
params.V0.LA = 2.3085;    % [mL]  Left atrial unstressed volume     — Source: Valenti Table 3.3
params.V0.LV = 3.5385;    % [mL]  Left ventricular unstressed V     — Source: Valenti Table 3.3

%-- Atrial compliance placeholders  [mL/mmHg]
%   (not used in the elastance model; retained for legacy / postprocess)
params.C.RA = 1/0.195;    % [mL/mmHg]  — Source: Valenti Table 3.3 (EB_RA = 0.195 mmHg/mL)
params.C.LA = 1/0.27;     % [mL/mmHg]  — Source: Valenti Table 3.3 (EB_LA = 0.27 mmHg/mL)

%-- Elastance — active (EA) and passive (EB)  [mmHg/mL]
%   ODE: E_ch(t) = EA * e(t) + EB   (Valenti Eq. 2.2)

% Left ventricle
params.E.LV.EA = 2.7;       % Source: Valenti Table 3.3
params.E.LV.EB = 0.069;

% Right ventricle
params.E.RV.EA = 0.43;      % Source: Valenti Table 3.3
params.E.RV.EB = 0.041264;

% Left atrium
params.E.LA.EA = 0.38;      % Source: Valenti Table 3.3
params.E.LA.EB = 0.27;

% Right atrium
params.E.RA.EA = 0.126;     % Source: Valenti Table 3.3
params.E.RA.EB = 0.195;

%-- Cardiac phase timings (stored as fractions of T_HB)
%   Converted to absolute seconds by apply_scaling and params_from_clinical.
%   Reference: Valenti Table 2.1

% Ventricular
params.Tc_LV_frac  = 0.265;   % LV contraction duration / T_HB
params.Tr_LV_frac  = 0.40;    % LV relaxation  duration / T_HB
params.Tc_RV_frac  = 0.30;    % RV contraction duration / T_HB
params.Tr_RV_frac  = 0.40;    % RV relaxation  duration / T_HB

% Atrial
params.t_ac_LA_frac = 0.75;   % LA contraction start   / T_HB
params.Tc_LA_frac   = 0.10;   % LA contraction duration/ T_HB
params.Tr_LA_frac   = 0.80;   % LA relaxation  duration/ T_HB
params.t_ac_RA_frac = 0.80;   % RA contraction start   / T_HB
params.Tc_RA_frac   = 0.10;   % RA contraction duration/ T_HB
params.Tr_RA_frac   = 0.70;   % RA relaxation  duration/ T_HB

%% =====================================================================
%  SYSTEMIC CIRCUIT  —  Valenti Table 3.3
%  Order: SAR (arterial) → SC (capillary) → SVEN (venous) → RA
%% =====================================================================

% Systemic arterial  —  RLC segment
params.R.SAR  = 0.5911;           % [mmHg·s/mL]   Source: Valenti T3.3
params.C.SAR  = 1.3315;           % [mL/mmHg]     Source: Valenti T3.3
params.L.SAR  = 2.9643e-4;        % [mmHg·s²/mL]  Source: Valenti T3.3

% Systemic capillary  —  RC segment
params.R.SC   = 2.17e-2;          % [mmHg·s/mL]   Source: Valenti T3.3
params.C.SC   = 2.7981e-1;        % [mL/mmHg]     Source: Valenti T3.3

% Systemic venous  —  RLC segment
params.R.SVEN = 3.596e-1;         % [mmHg·s/mL]   Source: Valenti T3.3
params.C.SVEN = 75.0;             % [mL/mmHg]     Source: Valenti T3.3
params.L.SVEN = 2.0643e-5;        % [mmHg·s²/mL]  Source: Valenti T3.3

% Unstressed volumes  [mL]  —  adult estimates (Valenti T3.3)
params.V0.SAR  = 250;             % [mL]  Source: Valenti T3.3
params.V0.SC   = 50;              % [mL]  Source: Valenti T3.3
params.V0.SVEN = 3200;            % [mL]  Source: Valenti T3.3

%% =====================================================================
%  PULMONARY CIRCUIT  —  Valenti Table 3.3
%  Order: PAR (arterial) → PC (capillary, single P state) → PVEN (venous) → LA
%% =====================================================================

% Pulmonary arterial  —  RLC segment
params.R.PAR  = 7.14e-2;          % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PAR  = 6.0043;           % [mL/mmHg]     Source: Valenti T3.3
params.L.PAR  = 2.0643e-5;        % [mmHg·s²/mL]  Source: Valenti T3.3

% Oxygenated pulmonary capillary branch
params.R.PCOX = 1.7538e-2;        % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PCOX = 5.7803;           % [mL/mmHg]     Source: Valenti T3.3

% Non-oxygenated pulmonary capillary branch
params.R.PCNO = 3.5174e-1;        % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PCNO = 4.9043e-2;        % [mL/mmHg]     Source: Valenti T3.3

% Pulmonary venous  —  RLC segment
params.R.PVEN = 3.75e-2;          % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PVEN = 11.381;           % [mL/mmHg]     Source: Valenti T3.3
params.L.PVEN = 2.0643e-5;        % [mmHg·s²/mL]  Source: Valenti T3.3

% Pulmonary unstressed volumes  [mL]  — Source: Valenti T3.3
params.V0.PAR  = 150;             % [mL]  Source: Valenti T3.3
params.V0.PCOX = 80;              % [mL]  retained for postprocess reference only — Source: Valenti T3.3
params.V0.PCNO = 5;               % [mL]  retained for postprocess reference only — Source: Valenti T3.3
params.V0.PVEN = 400;             % [mL]  Source: Valenti T3.3

%% =====================================================================
%  VALVE resistances  —  Valenti Table 3.3, Eq. (2.6)
%% =====================================================================
params.Rvalve.open   = 6.2872e-3;   % R_min  [mmHg·s/mL]  — Source: Valenti Table 3.3
params.Rvalve.closed = 9.4168e+4;   % R_max  [mmHg·s/mL]  — Source: Valenti Table 3.3 (numerical guard for closed valve)

%% =====================================================================
%  VSD shunt resistance  [mmHg·s/mL]
%  Adult reference: effectively closed (not present)
%  Pre-surgery:  computed from clinical data by params_from_clinical.m
%  Post-surgery: forced to R_VSD_CLOSED in params_from_clinical.m
%% =====================================================================
params.R.vsd = 1e6;   % [mmHg·s/mL]  healthy adult: VSD not present → effectively infinite
                      % Pre-surgery:  overwritten by params_from_clinical.m from clinical data
                      % Post-surgery: forced to 1e6 mmHg·s/mL by params_from_clinical.m
                      % Source: surgical closure equivalent (infinite lumped resistance)

%% =====================================================================
%  INITIAL CONDITIONS — 14-state vector  (Valenti Eq. 2.7)
%  [V_RA V_RV V_LA V_LV | V_SAR Q_SAR V_SC V_SVEN Q_SVEN |
%   V_PAR Q_PAR P_PC V_PVEN Q_PVEN]
%  P_PC initial ≈ (V_PCOX_0 - V0.PCOX)/C_PCOX ~ 0.9 mmHg
%  Source: Valenti Table 3.3 / physiological adult estimates
%% =====================================================================
%  Each entry is physiologically motivated (not zeros) per Guardrail §8.5:
ic = zeros(14, 1);
ic(idx.V_RA)   = 150;   % [mL]    RA volume ≈ RVEDV ballpark — Source: Valenti T3.3
ic(idx.V_RV)   = 150;   % [mL]    RV volume initial          — Source: Valenti T3.3
ic(idx.V_LA)   = 150;   % [mL]    LA volume initial          — Source: Valenti T3.3
ic(idx.V_LV)   = 150;   % [mL]    LV volume ≈ LVEDV ballpark — Source: Valenti T3.3
ic(idx.V_SAR)  = 850;   % [mL]    Systemic arterial volume   — Source: Valenti T3.3
ic(idx.Q_SAR)  = 60;    % [mL/s]  SAR flow ≈ resting CO/60s  — Source: estimated
ic(idx.V_SC)   = 55;    % [mL]    Systemic capillary volume  — Source: Valenti T3.3
ic(idx.V_SVEN) = 3300;  % [mL]    Systemic venous volume     — Source: Valenti T3.3
ic(idx.Q_SVEN) = 60;    % [mL/s]  SVEN flow                  — Source: estimated
ic(idx.V_PAR)  = 200;   % [mL]    Pulmonary arterial volume  — Source: Valenti T3.3
ic(idx.Q_PAR)  = 60;    % [mL/s]  PAR flow                   — Source: estimated
ic(idx.P_PC)   = 0.9;   % [mmHg]  Pulmonary cap. pressure    — Source: Valenti Eq. 2.7
ic(idx.V_PVEN) = 450;   % [mL]    Pulmonary venous volume    — Source: Valenti T3.3
ic(idx.Q_PVEN) = 60;    % [mL/s]  PVEN flow                  — Source: estimated
params.ic.V = ic';

end
