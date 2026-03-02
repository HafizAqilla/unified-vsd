function params = apply_scaling(params_ref, patient)
% APPLY_SCALING
% -----------------------------------------------------------------------
% Patient-specific allometric scaling of cardiovascular model parameters.
%
% Scales the 14-compartment adult RLC reference (Valenti Table 3.3) to a
% paediatric patient using BSA-based power laws, then corrects initial
% conditions using the patient blood-volume formula of Lundquist (2025).
%
% INPUTS:
%   params_ref  - adult reference struct from default_parameters()
%   patient     - struct with fields:
%                   .age_years   patient age          [years]
%                   .weight_kg   body weight          [kg]
%                   .height_cm   standing height      [cm]
%                   .sex         'M' or 'F'  (reserved)
%                   .BSA         optional; computed if absent [m²]
%
% OUTPUTS:
%   params      - scaled parameter struct, ready for system_rhs.m
%
% SCALING FRAMEWORK — Layer A: Geometry / Size Prior  (Lundquist 2025)
%   Scale factor:  s = BSA_patient / BSA_ref
%
%   Resistance   R  ~ s^(-1)   (larger body → lower R; Poiseuille)
%   Compliance   C  ~ s^(+1)   (larger body → more elastic volume)
%   Inertance    L  ~ s^(-1)   (same physical basis as R)
%   Unstressed V V0 ~ s^(+1)   (proportional to body size)
%   Heart rate   HR ~ s^(-0.33)(metabolic / basal scaling)
%   Elastance    E_LV/LA ~ s^(-1)
%                E_RV/RA ~ s^(-1.5)  (right heart at lower pressures)
%
%   NOT SCALED:
%     R.vsd          — pathological defect; assigned by params_from_clinical.m
%     Rvalve.closed  — numerical guard; must remain large
%
% INITIAL CONDITION CORRECTION (Hafiz / Lundquist):
%   Total circulating blood volume (BV) is scaled as:
%     BV = weight_kg * BV_per_kg   [mL]
%   where BV_per_kg differs between infants and older patients:
%     age < 1 year : 82 mL/kg   (Ref: Lundquist 2025, citing paediatric data)
%     age ≥ 1 year : 70 mL/kg   (Ref: adult / older child reference)
%   The venous compartments (SVEN, PVEN), which hold ~85% of BV, are
%   adjusted so that the total blood volume of all volume states equals
%   BV_patient instead of BV_adult.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Table 3.3.
%   [2] Lundquist et al. (2025). Allometric Scaling in Pediatric Cardiology.
%         Eqs. 3.2, 4.1, 4.3.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

%% 1. Patient BSA  (Mosteller formula unless supplied)
BSA_ref = 1.73;   % m²  — standard adult reference (Mosteller)

if isfield(patient, 'BSA') && ~isnan(patient.BSA)
    BSA_patient = patient.BSA;
else
    % Mosteller: BSA = sqrt(height_cm * weight_kg / 3600)
    BSA_patient = sqrt(patient.height_cm * patient.weight_kg / 3600);
end

s = BSA_patient / BSA_ref;   % dimensionless relative scale factor

%% 2. Copy reference params (no field is scaled twice)
params = params_ref;
params.scaling.BSA_ref     = BSA_ref;
params.scaling.BSA_patient = BSA_patient;
params.scaling.s           = s;
params.scaling.patient     = patient;

fprintf('[apply_scaling] BSA = %.3f m²  |  s = %.3f\n', BSA_patient, s);

%% Exponent definitions  (geometry-based prior, Lundquist 2025)
eR  = -1.0;   % Resistance   R ~ s^(-1)
eC  = +1.0;   % Compliance   C ~ s^(+1)
eL  = -1.0;   % Inertance    L ~ s^(-1)
eV0 = +1.0;   % Unstressed V V0 ~ s^(+1)

% =====================================================================
%  A. CARDIAC SCALING
% =====================================================================

%-- A1. Heart rate  HR ~ s^(-0.33)  [Lundquist Eq. 3.2]
params.HR = params_ref.HR * s^(-0.33);

%-- A2–A3. Elastances  (active EA and passive EB)
% LV/LA operate at systemic pressure → scale ~ s^(-1)
% RV/RA operate at pulmonary pressure → steeper scale ~ s^(-1.5)
params.E.LV.EA = params_ref.E.LV.EA * s^(-1.0);
params.E.LV.EB = params_ref.E.LV.EB * s^(-1.0);
params.E.RV.EA = params_ref.E.RV.EA * s^(-1.5);
params.E.RV.EB = params_ref.E.RV.EB * s^(-1.5);
params.E.LA.EA = params_ref.E.LA.EA * s^(-1.0);
params.E.LA.EB = params_ref.E.LA.EB * s^(-1.0);
params.E.RA.EA = params_ref.E.RA.EA * s^(-1.5);
params.E.RA.EB = params_ref.E.RA.EB * s^(-1.5);

%-- A4. Cardiac chamber unstressed volumes  V0 ~ s^(+1)
params.V0.LV = params_ref.V0.LV * s^eV0;
params.V0.RV = params_ref.V0.RV * s^eV0;
params.V0.LA = params_ref.V0.LA * s^eV0;
params.V0.RA = params_ref.V0.RA * s^eV0;

% =====================================================================
%  B. VASCULAR SCALING
% =====================================================================

%-- B1. Resistances  R ~ s^(-1)
params.R.SAR  = params_ref.R.SAR  * s^eR;
params.R.SC   = params_ref.R.SC   * s^eR;
params.R.SVEN = params_ref.R.SVEN * s^eR;
params.R.PAR  = params_ref.R.PAR  * s^eR;
params.R.PCOX = params_ref.R.PCOX * s^eR;
params.R.PCNO = params_ref.R.PCNO * s^eR;
params.R.PVEN = params_ref.R.PVEN * s^eR;
% R.vsd is NOT scaled — assigned by params_from_clinical.m per scenario

%-- B2. Compliances  C ~ s^(+1)
params.C.SAR  = params_ref.C.SAR  * s^eC;
params.C.SC   = params_ref.C.SC   * s^eC;
params.C.SVEN = params_ref.C.SVEN * s^eC;
params.C.PAR  = params_ref.C.PAR  * s^eC;
params.C.PCOX = params_ref.C.PCOX * s^eC;
params.C.PCNO = params_ref.C.PCNO * s^eC;
params.C.PVEN = params_ref.C.PVEN * s^eC;
if isfield(params_ref.C, 'RA'), params.C.RA = params_ref.C.RA * s^eC; end
if isfield(params_ref.C, 'LA'), params.C.LA = params_ref.C.LA * s^eC; end

%-- B3. Inertances  L ~ s^(-1)
params.L.SAR  = params_ref.L.SAR  * s^eL;
params.L.SVEN = params_ref.L.SVEN * s^eL;
params.L.PAR  = params_ref.L.PAR  * s^eL;
params.L.PVEN = params_ref.L.PVEN * s^eL;

%-- B4. Vascular unstressed volumes  V0 ~ s^(+1)
params.V0.SAR  = params_ref.V0.SAR  * s^eV0;
params.V0.SC   = params_ref.V0.SC   * s^eV0;
params.V0.SVEN = params_ref.V0.SVEN * s^eV0;
params.V0.PAR  = params_ref.V0.PAR  * s^eV0;
params.V0.PVEN = params_ref.V0.PVEN * s^eV0;
if isfield(params_ref.V0,'PCOX'), params.V0.PCOX = params_ref.V0.PCOX * s^eV0; end
if isfield(params_ref.V0,'PCNO'), params.V0.PCNO = params_ref.V0.PCNO * s^eV0; end

%-- B5. Valve open resistance  area ~ s  →  R_open ~ s^(-1)
params.Rvalve.open   = params_ref.Rvalve.open * s^(-1.0);
params.Rvalve.closed = params_ref.Rvalve.closed;   % unchanged (numerical guard)

% =====================================================================
%  C. INITIAL CONDITION SCALING  (14-state vector)
%
%  State layout:
%    [V_RA V_RV V_LA V_LV | V_SAR Q_SAR V_SC V_SVEN Q_SVEN |
%     V_PAR Q_PAR P_PC V_PVEN Q_PVEN]
%
%  Primary scaling rules:
%    Volume states (mL)   : ~ s^(+1)   — body-proportional compartments
%    Flow   states (mL/s) : ~ s^(+1)   — cardiac output ~ body size
%    P_PC   state  (mmHg) : ~ s^( 0)   — physiological BP ≠ body-size
%
%  Blood-volume correction (Hafiz / Lundquist):
%    After the primary s^(+1) scaling, adjust SVEN and PVEN so that the
%    total blood volume equals weight_kg * BV_per_kg.
% =====================================================================
ic = params_ref.ic.V(:);   % 14×1

ic_scaled               = ic;
vol_idx  = [1 2 3 4 5 7 8 10 13];   % volume states
flow_idx = [6 9 11 14];              % flow / inductor states
% state 12 (P_PC) is a pressure — no body-size scaling

ic_scaled(vol_idx)  = ic(vol_idx)  * s^eV0;
ic_scaled(flow_idx) = ic(flow_idx) * s^eV0;

%-- Blood-volume correction  (Lundquist 2025; Hafiz load_profile_a.m)
% BV_per_kg: infant < 1 yr uses 82 mL/kg; older / adult uses 70 mL/kg.
age_years = patient.age_years;
if age_years < 1
    BV_per_kg = 82;   % mL/kg  Source: Lundquist (2025) / paediatric reference
else
    BV_per_kg = 70;   % mL/kg  Source: standard adult
end

BV_target = patient.weight_kg * BV_per_kg;   % intended circulating volume [mL]

% Current total blood volume = sum of all volume states after s-scaling
BV_current = sum(ic_scaled(vol_idx));          % [mL]

% Distribute the correction proportionally to SVEN (idx 8) and PVEN (idx 13),
% which together carry ≈85% of total blood volume.
if BV_current > 0
    delta_BV  = BV_target - BV_current;
    ic_scaled(8)  = ic_scaled(8)  + 0.70 * delta_BV;   % SVEN  ~70% correction
    ic_scaled(13) = ic_scaled(13) + 0.15 * delta_BV;   % PVEN  ~15% correction
end

params.ic.V = ic_scaled(:)';   % store as row vector (ode15s accepts row or col)

% =====================================================================
%  D. ABSOLUTE TIMING FIELDS  (required by elastance_model.m)
% =====================================================================
T_HB = 60 / params.HR;

params.Tc_LV   = params.Tc_LV_frac   * T_HB;
params.Tr_LV   = params.Tr_LV_frac   * T_HB;
params.Tc_RV   = params.Tc_RV_frac   * T_HB;
params.Tr_RV   = params.Tr_RV_frac   * T_HB;
params.t_ac_LA = params.t_ac_LA_frac * T_HB;
params.Tc_LA   = params.Tc_LA_frac   * T_HB;
params.t_ar_LA = params.t_ac_LA + params.Tc_LA;
params.Tr_LA   = params.Tr_LA_frac   * T_HB;
params.t_ac_RA = params.t_ac_RA_frac * T_HB;
params.Tc_RA   = params.Tc_RA_frac   * T_HB;
params.t_ar_RA = params.t_ac_RA + params.Tc_RA;
params.Tr_RA   = params.Tr_RA_frac   * T_HB;

fprintf('[apply_scaling] HR_scaled = %.1f bpm | BV_target = %.0f mL\n', ...
    params.HR, BV_target);

end  % apply_scaling
