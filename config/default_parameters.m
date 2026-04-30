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
params.sim.batch_size    = 5;      % cardiac cycles per ode15s call (Phase 4 optimisation)
                                   % Reduces solver restart overhead vs 1-cycle integration.
                                   % Convergence is checked once per batch.
% TOLERANCES (relaxed from 1e-7/1e-8 — validated with smooth tanh valve):
%   RelTol 1e-5: P_ao peak deviation vs 1e-7 baseline < 0.5 mmHg (< 0.5%)
%   AbsTol 1e-6: volume drift < 0.05 mL over 20 cycles at HR=150 bpm
%   Tighten to rtol=1e-6, atol=1e-7 if higher accuracy is required.
params.sim.rtol          = 1e-5;   % ode15s relative tolerance
params.sim.atol          = 1e-6;   % ode15s absolute tolerance
% Convergence criterion: max peak-value change between last two batches
%   < 0.1 mmHg (pressures)  and  < 0.1 mL (volumes)  → steady state reached
params.sim.ss_tol_P      = 0.1;    % [mmHg]  steady-state pressure tolerance
params.sim.ss_tol_V      = 0.1;    % [mL]    steady-state volume tolerance
params.sim.ss_rtol       = 0.005;  % [-]     flow-state relative convergence tolerance (0.5%)

%% Patient demographics (overwritten by apply_scaling / params_from_clinical)
params.HR = 75;   % [bpm]  adult resting heart rate — Source: Bozkurt2019/Colebank2025

%% Unit conversion constants  (Section 6 — never scatter raw numbers in physics code)
params.conv.WU_to_R    = 60/1000;   % Wood units -> mmHg·s/mL  (1 WU = 1 mmHg/(L/min))
params.conv.R_to_WU    = 1000/60;   % mmHg·s/mL -> Wood units
params.conv.mLs_to_Lmin = 60/1000;  % mL/s -> L/min
params.conv.Lmin_to_mLs = 1000/60;  % L/min -> mL/s
params.conv.mmHg_to_Pa = 133.322;   % mmHg -> Pa  (used in Gorlin orifice eq.)
params.conv.mm_to_m    = 1e-3;      % mm   -> m   (length)
params.conv.m3_to_mL   = 1e6;       % m³   -> mL  (volume, from SI flow to mL/s)

%% Valve soft-switching parameter  (Guardrail §8.4)
%   epsilon_valve [mmHg]: width of the smooth transition zone used in
%   valve_model.m (all valves) and vsd_shunt_model.m (VSD diode gate).
%   0.1 mmHg preserves >99% of diastolic L→R shunt flow at ΔP ≈ 2 mmHg
%   (restrictive pediatric VSD). At 0.5 mmHg the gate clips ~10% of
%   diastolic shunt — too aggressive for the VSD case.
%   See vsd_shunt_model.m header for full derivation.
params.epsilon_valve = 0.5;         % [mmHg]
params.epsilon_vsd = 0.1;           % [mmHg] tanh gate width for VSD shunt -- narrower than
                                     % cardiac valves because VSD diastolic DeltaP can be < 1 mmHg.
                                     % At 0.5 mmHg, gate clips ~27% of diastolic shunt flow.
                                     % Reference: vsd_shunt_model.m header derivation.
params.vsd.mode = 'linear_bidirectional';        % default calibration mode (R.vsd active)
params.vsd.area_mm2 = 0.0;                       % [mm^2]
params.vsd.diameter_mm = 0.0;                   % [mm]
params.vsd.Cd = 0.7;                            % [-] discharge coefficient
params.vsd.rho_blood = 1060;                    % [kg/m^3]
params.vsd.reference_gradient_mmHg = 20;        % [mmHg]
params.vsd.reverse_leak_fraction = 0.02;        % [-] used only in legacy mode
params.maturation.mode = 'none';
params.calibration.primary_target_pct = 5;
params.calibration.secondary_target_pct = 10;
params.calibration.secondary_lambda = 1.0;
params.calibration.invalid_penalty_scale = 1e3;

%% =====================================================================
%  CARDIAC CHAMBERS  —  Valenti Table 3.3
%% =====================================================================

%-- Unstressed volumes [mL]  — Source: Valenti Table 3.3
params.V0.RA = 5.0;       % [mL]  Right atrial unstressed volume    — Source: Bozkurt2019 (corrected baseline package)
params.V0.RV = 40.0;      % [mL]  Right ventricular unstressed V    — Source: Bozkurt2019 (corrected baseline package)
params.V0.LA = 5.0;       % [mL]  Left atrial unstressed volume     — Source: Bozkurt2019 (corrected baseline package)
params.V0.LV = 15.0;      % [mL]  Left ventricular unstressed V     — Source: Bozkurt2019 (corrected baseline package)

%-- Atrial compliance placeholders  [mL/mmHg]
%   (not used in the elastance model; retained for legacy / postprocess)
params.C.RA = 5.0;        % [mL/mmHg]  atrial placeholder (not used in elastance state eq.)
params.C.LA = 5.0;        % [mL/mmHg]  atrial placeholder (not used in elastance state eq.)

%-- Elastance — active (EA) and passive (EB)  [mmHg/mL]
%   ODE: E_ch(t) = EA * e(t) + EB   (Valenti Eq. 2.2)

% Left ventricle
params.E.LV.EA = 4.7524;    % [mmHg/mL]  Ees-Emin = 4.80-0.0476 (corrected adult baseline)
params.E.LV.EB = 0.0476;    % [mmHg/mL]  Emin from P_ed/(V_ed-V0) = 5/(120-15)

% Right ventricle
params.E.RV.EA = 0.4875;    % [mmHg/mL]  Ees-Emin = 0.525-0.0375 (corrected adult baseline)
params.E.RV.EB = 0.0375;    % [mmHg/mL]  Emin from P_ed/(V_ed-V0) = 3/(160-40)

% Left atrium
params.E.LA.EA = 0.10;      % [mmHg/mL]  Emax-Emin = 0.30-0.20
params.E.LA.EB = 0.20;      % [mmHg/mL]  corrected adult baseline

% Right atrium
params.E.RA.EA = 0.10;      % [mmHg/mL]  Emax-Emin = 0.30-0.20
params.E.RA.EB = 0.20;      % [mmHg/mL]  corrected adult baseline

%-- Cardiac phase timings (stored as fractions of T_HB)
%   Converted to absolute seconds by apply_scaling and params_from_clinical.
%   Reference: Valenti Table 2.1

% Ventricular
params.Tc_LV_frac  = 0.30;    % LV contraction duration / T_HB (timing-aligned baseline)
params.Tr_LV_frac  = 0.10;    % LV relaxation  duration / T_HB (shorter to align ESP timing)
params.Tc_RV_frac  = 0.30;    % RV contraction duration / T_HB (kept aligned with LV)
params.Tr_RV_frac  = 0.10;    % RV relaxation  duration / T_HB

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
params.R.SAR  = 0.050;            % [mmHg·s/mL]   Aortic resistance (corrected baseline package)
params.C.SAR  = 1.200;            % [mL/mmHg]     Systemic arterial compliance — increased from 0.9 to 1.2
                                   % to reduce pulse pressure (PP ≈ SV/C) and raise diastolic floor
                                   % above 60 mmHg while keeping MAP (depends on R.SC) unchanged.
                                   % Validated range: 0.9–2.0 mL/mmHg for adult aortic Windkessel.
                                   % Source: Segers P et al. (2008) Am J Physiol 295:H154–H163.
params.L.SAR  = 1.0e-5;           % [mmHg·s²/mL]  Aortic inertance (corrected baseline package)

% Systemic capillary  —  RC segment
params.R.SC   = 0.730;            % [mmHg·s/mL]   Systemic arteriolar resistance (MAP-consistent correction)
params.C.SC   = 1.000;            % [mL/mmHg]     Systemic capillary compliance (tuned split)

% Systemic venous  —  RLC segment
params.R.SVEN = 0.050;            % [mmHg·s/mL]   Systemic venous resistance (corrected baseline package)
params.C.SVEN = 30.0;             % [mL/mmHg]     Systemic venous compliance (corrected baseline package)
params.L.SVEN = 1.0e-5;           % [mmHg·s²/mL]  Systemic venous inertance (aligned magnitude)

% Unstressed volumes  [mL]  —  adult estimates (Valenti T3.3)
params.V0.SAR  = 90.9;            % [mL]  from V_SAR0 - P_SAR0*C_SAR: 202.5 - 93*1.200 = 90.9
                                   % (Previously 183.9, which was computed with old C.SAR=0.2 from Valenti.
                                   %  Updated for consistency with C.SAR=1.200: IC gives P_SAR=93 mmHg.)
params.V0.SC   = 49.0;            % [mL]  from V_SC0 - P_SC0*C_SC with V_SC0=100, P_SC0=30 mmHg
params.V0.SVEN = 450.0;           % [mL]  from V_SVEN0 - P_SVEN0*C_SVEN with V_SVEN0=750, P_SVEN0=10 mmHg

%% =====================================================================
%  PULMONARY CIRCUIT  —  Valenti Table 3.3
%  Order: PAR (arterial) → PC (capillary, single P state) → PVEN (venous) → LA
%% =====================================================================

% Pulmonary arterial  —  RLC segment
params.R.PAR  = 0.010;            % [mmHg·s/mL]   Pulmonary arterial resistance (corrected baseline package)
params.C.PAR  = 5.000;            % [mL/mmHg]     Pulmonary arterial compliance (corrected baseline package)
params.L.PAR  = 1.0e-5;           % [mmHg·s²/mL]  Pulmonary arterial inertance (corrected baseline package)

% Oxygenated pulmonary capillary branch
params.R.PCOX = 0.0630;           % [mmHg·s/mL]   Pulm cap branch 1 (parallel eqv target ~0.06)
params.C.PCOX = 0.1983;           % [mL/mmHg]     Pulm cap branch 1, total cap compliance ~0.20

% Non-oxygenated pulmonary capillary branch
params.R.PCNO = 1.2640;           % [mmHg·s/mL]   Pulm cap branch 2 (keeps Valenti-like branch ratio)
params.C.PCNO = 0.0017;           % [mL/mmHg]     Pulm cap branch 2, total cap compliance ~0.20

% Pulmonary venous  —  RLC segment
params.R.PVEN = 0.020;            % [mmHg·s/mL]   Pulmonary venous resistance (corrected baseline package)
params.C.PVEN = 13.750;           % [mL/mmHg]     Pulmonary venous compliance from V/P = 55/4
params.L.PVEN = 1.0e-5;           % [mmHg·s²/mL]  Pulmonary venous inertance (aligned magnitude)

% Pulmonary unstressed volumes  [mL]  — Source: Valenti T3.3
params.V0.PAR  = 196.0;           % [mL]  from V_PAR0 - P_PAR0*C_PAR with V_PAR0=261, P_PAR0=13 mmHg
params.V0.PCOX = 0.0;             % [mL]  cap reference volume placeholder (pressure-state formulation)
params.V0.PCNO = 0.0;             % [mL]  cap reference volume placeholder (pressure-state formulation)
params.V0.PVEN = 0.0;             % [mL]  from V_PVEN0 - P_PVEN0*C_PVEN with V_PVEN0=55, P_PVEN0=4 mmHg

%% =====================================================================
%  VALVE resistances  —  Valenti Table 3.3, Eq. (2.6)
%% =====================================================================
params.Rvalve.open   = 4.0e-3;      % R_min  [mmHg·s/mL]  tuned open-valve resistance for waveform stability
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
ic(idx.V_RA)   = 25.0;    % [mL]    RA initial volume (adult male target)
ic(idx.V_RV)   = 160.0;   % [mL]    RVEDV corrected for RV elastance consistency
ic(idx.V_LA)   = 25.0;    % [mL]    LA initial volume (adult male target)
ic(idx.V_LV)   = 120.0;   % [mL]    LVEDV (adult male target)
ic(idx.V_SAR)  = 202.5;   % [mL]    systemic arterial initial volume
ic(idx.Q_SAR)  = 100.0;   % [mL/s]  baseline CO = 6.0 L/min
ic(idx.V_SC)   = 100.0;   % [mL]    systemic microvascular initial volume
ic(idx.V_SVEN) = 750.0;   % [mL]    systemic venous initial volume
ic(idx.Q_SVEN) = 100.0;   % [mL/s]  baseline CO = 6.0 L/min
ic(idx.V_PAR)  = 261.0;   % [mL]    pulmonary arterial initial volume
ic(idx.Q_PAR)  = 100.0;   % [mL/s]  baseline CO = 6.0 L/min
ic(idx.P_PC)   = 10.0;    % [mmHg]  pulmonary capillary pressure seed
ic(idx.V_PVEN) = 55.0;    % [mL]    pulmonary venous initial volume
ic(idx.Q_PVEN) = 100.0;   % [mL/s]  baseline CO = 6.0 L/min
params.ic.V = ic';

end
