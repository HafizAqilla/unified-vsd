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
%   All exponents verified against Lundquist (2025) Tables 2 and 3:
%     Resistance   R  ~ s^(-1)   [T1/T2] (larger body → lower R; Poiseuille)
%     Compliance   C  ~ s^(+1)   [T2]    (larger body → more elastic volume)
%     Inertance    L  ~ s^(-1)            (geometry, same basis as R)
%     Unstressed V V0 ~ s^(+1)   [T2]    (proportional to body size)
%     Heart rate   HR ~ s^(-0.33)[T2]    (metabolic / basal scaling)
%     Elastance    E_LV/LA ~ s^(-1)  [T2] "Contractility/stiffness LV/LA"
%                  E_RV/RA ~ s^(-1.5)[T2] "Contractility/stiffness RV/RA"
%     Valve open area ~ s^(+1)    [T3] → R_valve_open ~ s^(-1) (orifice flow)
%     VSD shunt area  ~ s^(+1)    [T3] → handled in params_from_clinical.m
%
%   NOT SCALED — with justification:
%     R.vsd         — pathological; assigned per patient by params_from_clinical.m.
%                    (Lundquist T3: shunt area ~ BSA is already implicit in the
%                     Gorlin calculation using the measured echo defect diameter.)
%     Rvalve.closed — numerical guard; scaled value would break the diode logic
%     epsilon_valve — Lundquist T3 gives valve opening constants ~ s^(-0.5),
%                    but those are physiological pressure thresholds in the
%                    Lundquist model; epsilon_valve in our tanh switch is a
%                    numerical continuity parameter, not equivalent.
%   BSA_ref NOTE:
%     BSA_ref = 1.73 m² matches the Valenti (2023) reference subject
%     (~70 kg / 175 cm adult). The Lundquist (2025) reference subject is
%     30yo 85 kg 180 cm → BSA ≈ 2.06 m² — only the scaling EXPONENTS
%     are borrowed from Lundquist; the base parameter values come from Valenti.
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
%   [1] Valenti (2023). Thesis: Full-order 0D cardiovascular model. Table 3.3.
%   [2] Lundquist et al. (2025). Patient-Specific Pediatric Cardiovascular LPM —
%         Scaling Across Ages and Sizes. ASAIO J.
%         DOI: 10.1097/MAT.0000000000002528.
%         Table 2 (cardiac: HR, E, V0, blood volume) and
%         Table 3 (valve open/closed area ~ BSA; shunt area ~ BSA).
%   [3] Hafiz Aqilla S (2026). Seminar Thesis, UI. Sec 3.8.4 Eq 3.7
%         (R~BSA^-1, C~BSA, V0~BSA for pediatric VSD LPM).
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

%% 1. Patient BSA  (Mosteller formula unless supplied)
BSA_ref = 1.73;   % m²  — Valenti (2023) reference subject (~70 kg/175 cm adult)
                  %       Note: Lundquist (2025) uses 85 kg/180 cm → BSA ≈ 2.06 m²
                  %       but only the scaling EXPONENTS come from Lundquist;
                  %       the base parameter values (Valenti Table 3.3) fix BSA_ref.

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

%% Exponent definitions  (verified vs Lundquist 2025, Tables 2 & 3; Hafiz Eq 3.7)
eR  = -1.0;   % Resistance   R ~ s^(-1)   [Lundquist T2, T3; Hafiz Eq 3.7]
eC  = +1.0;   % Compliance   C ~ s^(+1)   [Lundquist T2;    Hafiz Eq 3.7]
eL  = -1.0;   % Inertance    L ~ s^(-1)   [geometry, same basis as eR]
eV0 = +1.0;   % Unstressed V V0 ~ s^(+1)  [Lundquist T2 'V0'; Hafiz Eq 3.7]

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

%-- A4. Cardiac chamber unstressed volumes
%   NOTE: These s^(+1) values are OVERRIDDEN in Section C by BV_scale
%   (blood-volume ratio). Left here for reference only.
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

%-- B4. Vascular unstressed volumes
%   NOTE: These s^(+1) values are OVERRIDDEN in Section C by BV_scale.
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
%  C. INITIAL CONDITION COMPUTATION  (pressure-initialized, blood-volume-conserving)
%
%  ROOT CAUSE of paediatric failure with naive s^(+1) IC scaling:
%    V0.SVEN × s = 3200 × 0.127 = 406 mL  >  BV_infant = 303 mL
%    → venous unstressed volume exceeds total patient blood volume
%    → P_SVEN = (V_SVEN - V0.SVEN)/C_SVEN is always negative
%    → cascade pressure failure on cycle 1
%
%  FIX (two-part):
%  (1) Re-scale ALL V0 values by BV_scale = BV_patient / BV_adult_ref
%      (blood-volume ratio) instead of s^(+1) (BSA ratio).
%      Physical basis: unstressed volume is a fixed fraction of total
%      circulating blood volume across species (Guyton 1991, §15).
%      BSA^1 scaling is appropriate for R, C, L (geometry) but not for
%      V0 (reservoir filling level).
%  (2) Compute ICs from nominal filling pressures: V = V0 + P_nom × C
%      V0.SVEN is then adjusted to enforce blood conservation exactly.
%
%  Nominal pressures (size-independent physiological reference):
%    P_nom_SAR  = 65 mmHg  (diastolic aortic, conservative nominal)
%    P_nom_SC   = 15 mmHg  (systemic capillary)
%    P_nom_SVEN =  2 mmHg  (CVP; V0.SVEN adjusted to conserve BV)
%    P_nom_PAR  = 15 mmHg  (PA diastolic)
%    P_nom_PVEN =  6 mmHg  (pulmonary venous ≈ LAP)
%    P_nom_PC   =  8 mmHg  (pulmonary capillary, pressure state)
%    P_nom_LV_ED =  8 mmHg  (LVEDP nominal)
%    P_nom_RV_ED =  4 mmHg  (RVEDP nominal)
%    P_nom_LA_ED =  6 mmHg  (LAP nominal)
%    P_nom_RA_ED =  4 mmHg  (RAP nominal)
% =====================================================================

% --- Blood volume budget ----------------------------------------------
age_years = patient.age_years;
if age_years < 1
    BV_per_kg = 82;    % mL/kg  (Lundquist 2025 / paediatric reference)
else
    BV_per_kg = 70;    % mL/kg  (standard adult/older child reference)
end
BV_patient = patient.weight_kg * BV_per_kg;   % [mL] target circulating volume

% Index struct (Guardrail §7.1 — never hardcode indices)
sidx = params.idx;

% --- Re-scale ALL V0 values by blood-volume ratio (not BSA ratio) -----
%   BV_adult_ref = 70 kg × 70 mL/kg = 4900 mL  (Valenti reference person)
BV_adult_ref = 4900;                         % [mL]  Source: 70 kg × 70 mL/kg
BV_scale     = BV_patient / BV_adult_ref;    % blood-volume scaling ratio

params.V0.LV   = params_ref.V0.LV   * BV_scale;
params.V0.RV   = params_ref.V0.RV   * BV_scale;
params.V0.LA   = params_ref.V0.LA   * BV_scale;
params.V0.RA   = params_ref.V0.RA   * BV_scale;
params.V0.SAR  = params_ref.V0.SAR  * BV_scale;
params.V0.SC   = params_ref.V0.SC   * BV_scale;
params.V0.SVEN = params_ref.V0.SVEN * BV_scale;   % will be adjusted below
params.V0.PAR  = params_ref.V0.PAR  * BV_scale;
params.V0.PVEN = params_ref.V0.PVEN * BV_scale;
if isfield(params.V0,'PCOX'), params.V0.PCOX = params_ref.V0.PCOX * BV_scale; end
if isfield(params.V0,'PCNO'), params.V0.PCNO = params_ref.V0.PCNO * BV_scale; end

% --- Nominal filling pressures (size-independent) ----------------------
% P_nom_SAR uses mean arterial pressure (MAP≈93 mmHg), NOT diastolic (65 mmHg).
% Using 65 mmHg (diastolic) caused V_SAR_ic to be ~40 mL too low; blood
% conservation then raised V0.SVEN by the same amount, creating a single
% low-CO equilibrium. At MAP=93 mmHg, V_SAR starts at its true mean
% operating point and the correct adult steady state is reachable.
P_nom_SAR   = 93;   % [mmHg]  mean arterial pressure (MAP) — Source: standard adult physiology
P_nom_SC    = 15;   % [mmHg]  systemic capillary
P_nom_SVEN  =  2;   % [mmHg]  CVP — V0.SVEN will be adjusted to conserve BV
P_nom_PAR   = 15;   % [mmHg]  PA diastolic
P_nom_PVEN  =  6;   % [mmHg]  pulmonary venous
P_nom_PC    =  8;   % [mmHg]  pulmonary capillary (pressure state)
P_nom_LV_ED =  8;   % [mmHg]  LVEDP
P_nom_RV_ED =  4;   % [mmHg]  RVEDP
P_nom_LA_ED =  6;   % [mmHg]  LAP
P_nom_RA_ED =  4;   % [mmHg]  RAP

% --- Build IC vector from nominal pressures ---------------------------
ic_p = zeros(14, 1);

% Vascular compartments: V = V0 + P × C
ic_p(sidx.V_SAR)  = params.V0.SAR  + P_nom_SAR  * params.C.SAR;
ic_p(sidx.V_SC)   = params.V0.SC   + P_nom_SC   * params.C.SC;
ic_p(sidx.V_SVEN) = params.V0.SVEN + P_nom_SVEN * params.C.SVEN;   % provisional
ic_p(sidx.V_PAR)  = params.V0.PAR  + P_nom_PAR  * params.C.PAR;
ic_p(sidx.V_PVEN) = params.V0.PVEN + P_nom_PVEN * params.C.PVEN;
ic_p(sidx.P_PC)   = P_nom_PC;   % pressure state — set directly [mmHg]

% Cardiac chambers: P = E_EB × (V - V0)  →  V = V0 + P / E_EB
ic_p(sidx.V_LV) = params.V0.LV + P_nom_LV_ED / params.E.LV.EB;
ic_p(sidx.V_RV) = params.V0.RV + P_nom_RV_ED / params.E.RV.EB;
ic_p(sidx.V_LA) = params.V0.LA + P_nom_LA_ED / params.E.LA.EB;
ic_p(sidx.V_RA) = params.V0.RA + P_nom_RA_ED / params.E.RA.EB;

% Flow states: initialise to cardiac output estimate (converges within 1–2 cycles)
%   CO_adult_ref ≈ 5 L/min = 83 mL/s for a 70 kg adult.
%   Scale by BV_scale (blood-volume ratio) so that:
%     Adult  (BV_scale = 1.0): Q_init = 83    mL/s  (~5 L/min)   ✓
%     Infant (BV_scale = 0.06): Q_init =  5.1 mL/s  (~0.3 L/min) ✓
%   Using SV_est=5 mL (old code) gave Q_init=6.1 mL/s for adults — 14× too low —
%   trapping the ODE in a collapsed, low-CO equilibrium that never recovers.
CO_adult_ref_mLs = 83;                     % [mL/s]  5 L/min adult reference CO
Q_init = CO_adult_ref_mLs * BV_scale;     % [mL/s]  patient-scaled CO estimate
ic_p(sidx.Q_SAR)  = Q_init;
ic_p(sidx.Q_SVEN) = Q_init;
ic_p(sidx.Q_PAR)  = Q_init;
ic_p(sidx.Q_PVEN) = Q_init;
fprintf('[apply_scaling] Q_init=%.1f mL/s  (CO_est=%.2f L/min)  P_nom_SAR=%g mmHg\n', ...
    Q_init, Q_init*60/1000, P_nom_SAR);

% --- Blood conservation: adjust V0.SVEN so sum(V_ic) = BV_patient -----
%   SVEN is the dominant venous reservoir. Shift V0.SVEN (keeping
%   P_nom_SVEN fixed) until the total volume budget exactly matches.
vol_idx_c = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV sidx.V_SAR ...
             sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];

BV_no_SVEN    = sum(ic_p(vol_idx_c)) - ic_p(sidx.V_SVEN);   % all except SVEN
V_SVEN_target = BV_patient - BV_no_SVEN;

if V_SVEN_target < params.V0.SVEN * 0.3
    warning('apply_scaling:BVbudget', ...
        ['BV_patient (%.0f mL) is very small relative to other compartments ' ...
         '(%.0f mL). Clamping SVEN to 30%% of V0.SVEN. ' ...
         'Calibration is strongly recommended.'], BV_patient, BV_no_SVEN);
    V_SVEN_target = max(V_SVEN_target, params.V0.SVEN * 0.3);
end

% Adjust V0.SVEN to maintain P_nom_SVEN = 2 mmHg at the adjusted volume
params.V0.SVEN    = V_SVEN_target - P_nom_SVEN * params.C.SVEN;
ic_p(sidx.V_SVEN) = V_SVEN_target;

params.ic.V    = ic_p(:)';   % store as row vector (ode15s accepts row or col)
params.scaling.BV_patient = BV_patient;
params.scaling.BV_scale   = BV_scale;

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

fprintf('[apply_scaling] BSA=%.3f m² | s=%.3f | BV_patient=%.0f mL | BV_scale=%.4f\n', ...
    BSA_patient, s, BV_patient, BV_scale);
fprintf('[apply_scaling] IC: V_LV=%.1f  V_RV=%.1f  V_SVEN=%.1f mL  P_SVEN_ic=%.1f mmHg\n', ...
    ic_p(sidx.V_LV), ic_p(sidx.V_RV), ic_p(sidx.V_SVEN), P_nom_SVEN);

end  % apply_scaling
