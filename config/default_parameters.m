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

%% Simulation control
params.sim.nCyclesSteady = 20;     % cardiac cycles integrated to reach steady-state
params.sim.nCyclesKeep   = 2;      % last N cycles retained for post-processing
params.sim.rtol          = 1e-7;   % ode15s relative tolerance
params.sim.atol          = 1e-8;   % ode15s absolute tolerance

%% Patient demographics (overwritten by apply_scaling / params_from_clinical)
params.HR = 75;   % [bpm]  adult resting heart rate

%% Unit conversion constants
params.conv.WU_to_R    = 60/1000;   % Wood units -> mmHg·s/mL
params.conv.mmHg_to_Pa = 133.322;   % mmHg -> Pa  (informational)

%% =====================================================================
%  CARDIAC CHAMBERS  —  Valenti Table 3.3
%% =====================================================================

%-- Unstressed volumes (mL)
params.V0.RA = 3.5385;
params.V0.RV = 8.4067;
params.V0.LA = 2.3085;
params.V0.LV = 3.5385;

%-- Atrial compliance placeholders  [mL/mmHg]
%   (not used in the elastance model; retained for legacy / postprocess)
params.C.RA = 1/0.195;
params.C.LA = 1/0.27;

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
params.C.SAR  = 1.3315;           % [mL/mmHg]
params.L.SAR  = 2.9643e-4;        % [mmHg·s²/mL]

% Systemic capillary  —  RC segment
params.R.SC   = 2.17e-2;          % [mmHg·s/mL]   Source: Valenti T3.3
params.C.SC   = 2.7981e-1;        % [mL/mmHg]

% Systemic venous  —  RLC segment
params.R.SVEN = 3.596e-1;         % [mmHg·s/mL]   Source: Valenti T3.3
params.C.SVEN = 75.0;             % [mL/mmHg]
params.L.SVEN = 2.0643e-5;        % [mmHg·s²/mL]

% Unstressed volumes  [mL]  —  adult estimates
params.V0.SAR  = 250;
params.V0.SC   = 50;
params.V0.SVEN = 3200;

%% =====================================================================
%  PULMONARY CIRCUIT  —  Valenti Table 3.3
%  Order: PAR (arterial) → PC (capillary, single P state) → PVEN (venous) → LA
%% =====================================================================

% Pulmonary arterial  —  RLC segment
params.R.PAR  = 7.14e-2;          % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PAR  = 6.0043;           % [mL/mmHg]
params.L.PAR  = 2.0643e-5;        % [mmHg·s²/mL]

% Oxygenated pulmonary capillary branch
params.R.PCOX = 1.7538e-2;        % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PCOX = 5.7803;           % [mL/mmHg]

% Non-oxygenated pulmonary capillary branch
params.R.PCNO = 3.5174e-1;        % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PCNO = 4.9043e-2;        % [mL/mmHg]

% Pulmonary venous  —  RLC segment
params.R.PVEN = 3.75e-2;          % [mmHg·s/mL]   Source: Valenti T3.3
params.C.PVEN = 11.381;           % [mL/mmHg]
params.L.PVEN = 2.0643e-5;        % [mmHg·s²/mL]

% Pulmonary unstressed volumes  [mL]
params.V0.PAR  = 150;
params.V0.PCOX = 80;    % retained for postprocess reference only
params.V0.PCNO = 5;     % retained for postprocess reference only
params.V0.PVEN = 400;

%% =====================================================================
%  VALVE resistances  —  Valenti Table 3.3, Eq. (2.6)
%% =====================================================================
params.Rvalve.open   = 6.2872e-3;   % R_min  [mmHg·s/mL]
params.Rvalve.closed = 9.4168e+4;   % R_max  [mmHg·s/mL]  (numerical guard)

%% =====================================================================
%  VSD shunt resistance  [mmHg·s/mL]
%  Adult reference: effectively closed (not present)
%  Pre-surgery:  computed from clinical data by params_from_clinical.m
%  Post-surgery: forced to R_VSD_CLOSED in params_from_clinical.m
%% =====================================================================
params.R.vsd = 200;   % nominally ~closed; overwritten per scenario

%% =====================================================================
%  INITIAL CONDITIONS — 14-state vector  (Valenti Eq. 2.7)
%  [V_RA V_RV V_LA V_LV | V_SAR Q_SAR V_SC V_SVEN Q_SVEN |
%   V_PAR Q_PAR P_PC V_PVEN Q_PVEN]
%  P_PC initial ≈ (V_PCOX_0 - V0.PCOX)/C_PCOX ~ 0.9 mmHg
%% =====================================================================
params.ic.V = [150;  150;  150;  150; ...
               850;  60;   55;   3300; 60; ...
               200;  60;   0.9;  450;  60]';

end
