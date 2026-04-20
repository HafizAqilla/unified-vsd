function params = apply_scaling(params_ref, patient)
% APPLY_SCALING
% -----------------------------------------------------------------------
% Patient-specific allometric scaling of cardiovascular model parameters.
%
% Scales the 14-compartment adult RLC reference (Valenti Table 3.3) to a
% paediatric patient using weight-based power laws from Zhang et al. (2019).
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
% SCALING FRAMEWORK — Zhang et al. (2019) Allometric Scaling Law (ASL)
%   Scale factor:  w = patient.weight_kg / W_ref   (W_ref = 70 kg)
%
%   All exponents from Zhang et al. (2019) Table 1:
%     Heart Rate      HR      ~ w^(-0.30)   [Zhang Table 1: β = -0.3]
%     Elastance LV/LA E_lv,la ~ w^(-0.50)   [Zhang Table 1: β = -0.5  for LV/LA]
%     Elastance RV/RA E_rv,ra ~ w^(-0.75)   [Zhang Table 1: β = -0.75 for RV/RA]
%     Unstressed Vol  V0      ~ w^(+0.80)   [Zhang Table 1: β = +0.8  blood volumes]
%     Compliance      C       ~ w^(+1.00)   [Zhang Table 1: β = +1.0  systemic+pulm]
%     Resist Systemic R_sys   ~ w^(-0.475)  [Zhang Table 1: β = -0.475 systemic vasc]
%     Resist Pulmon.  R_pul   ~ w^(-0.70)   [Zhang Table 1: β = -0.70  pulmonary vasc]
%     Valve open area ~ w^(+1)→ R_open     ~ w^(-0.90)   [Zhang Table 1 + orifice flow]
%
%   NOT SCALED — with justification:
%     R.vsd         — pathological; assigned per patient by params_from_clinical.m
%     Rvalve.closed — numerical guard; scaled value would break the diode logic
%     L.*           — inertances excluded per Zhang (2019); fixed at reference values
%
%   W_ref NOTE:
%     W_ref = 70 kg matches the Valenti (2023) reference subject
%     (~70 kg / 175 cm adult). Only the scaling EXPONENTS come from Zhang (2019);
%     the base parameter values originate from Valenti Table 3.3.
%
% INITIAL CONDITION CORRECTION (blood-volume-conserving):
%   Total circulating blood volume (BV) is scaled as:
%     BV = weight_kg * BV_per_kg   [mL]
%   where BV_per_kg differs between infants and older patients:
%     age < 1 year : 82 mL/kg   (paediatric reference)
%     age ≥ 1 year : 70 mL/kg   (standard adult/older child reference)
%   The venous compartment (SVEN) is adjusted so that the total blood
%   volume of all volume states equals BV_patient.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis: Full-order 0D cardiovascular model. Table 3.3.
%   [2] Zhang X, Haneishi H, Liu H. (2019). Multiscale modeling of the
%         cardiovascular system for infants, children, and adolescents:
%         Age-related alterations in cardiovascular parameters and hemodynamics.
%         Computers in Biology and Medicine, 108, 200–212.
%         Table 1: scaling exponents β for age-related cardiovascular parameters.
%   [3] Hafiz Aqilla S (2026). Seminar Thesis, UI.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-13
% VERSION:  2.0  (Zhang 2019 weight-based ASL replacing Lundquist BSA-based)
% -----------------------------------------------------------------------

%% 1. Weight ratio (Zhang 2019 primary scale factor)
W_ref = 70;   % kg — Valenti (2023) reference subject (~70 kg / 175 cm adult)
w = patient.weight_kg / W_ref;   % dimensionless weight ratio

%% 1b. Patient BSA (Mosteller formula) — retained for reference/logging only
BSA_ref = 1.73;   % m²  — Valenti (2023) reference subject
if isfield(patient, 'BSA') && ~isnan(patient.BSA)
    BSA_patient = patient.BSA;
else
    BSA_patient = sqrt(patient.height_cm * patient.weight_kg / 3600);
end

%% 2. Copy reference params (no field is scaled twice)
params = params_ref;
params.scaling.W_ref       = W_ref;
params.scaling.w           = w;
params.scaling.BSA_ref     = BSA_ref;
params.scaling.BSA_patient = BSA_patient;
params.scaling.patient     = patient;

fprintf('[apply_scaling] weight=%.1f kg | W_ref=%.0f kg | w=%.4f | BSA=%.3f m²\n', ...
    patient.weight_kg, W_ref, w, BSA_patient);

%% Zhang (2019) exponent definitions — Table 1
% Cardiac
eHR     = -0.30;   % Heart rate               [Zhang Table 1]
eE_lv   = -0.50;   % LV/LA elastance          [Zhang Table 1: β=-0.5]
eE_rv   = -0.75;   % RV/RA elastance          [Zhang Table 1: β=-0.75]
eV0     = +0.80;   % Unstressed volumes       [Zhang Table 1: β=+0.8]

% Vascular
eR_sys  = -0.475;  % Systemic resistance      [Zhang Table 1: β=-0.475]
eR_pul  = -0.70;   % Pulmonary resistance     [Zhang Table 1: β=-0.70]
eC      = +1.00;   % Compliance (all)         [Zhang Table 1: β=+1.0]
eRv_op  = -0.90;   % Valve open resistance    [Zhang Table 1 orifice-area → R]

% =====================================================================
%  A. CARDIAC SCALING
% =====================================================================

%-- A1. Heart rate  HR ~ w^(-0.30)  [Zhang 2019 Table 1]
params.HR = params_ref.HR * w^eHR;

%-- A2–A3. Elastances
% LV/LA: systemic pressure → β = -0.50 [Zhang Table 1]
% RV/RA: pulmonary pressure → β = -0.75 [Zhang Table 1]
params.E.LV.EA = params_ref.E.LV.EA * w^eE_lv;
params.E.LV.EB = params_ref.E.LV.EB * w^eE_lv;
params.E.LA.EA = params_ref.E.LA.EA * w^eE_lv;
params.E.LA.EB = params_ref.E.LA.EB * w^eE_lv;
params.E.RV.EA = params_ref.E.RV.EA * w^eE_rv;
params.E.RV.EB = params_ref.E.RV.EB * w^eE_rv;
params.E.RA.EA = params_ref.E.RA.EA * w^eE_rv;
params.E.RA.EB = params_ref.E.RA.EB * w^eE_rv;

%-- A4. Cardiac chamber unstressed volumes  V0 ~ w^(+0.80)
%   NOTE: These values are OVERRIDDEN in Section C by BV_scale
%   (blood-volume ratio). Left here as Zhang-scaling geometric prior.
params.V0.LV = params_ref.V0.LV * w^eV0;
params.V0.RV = params_ref.V0.RV * w^eV0;
params.V0.LA = params_ref.V0.LA * w^eV0;
params.V0.RA = params_ref.V0.RA * w^eV0;

% =====================================================================
%  B. VASCULAR SCALING
% =====================================================================

%-- B1. Resistances (Zhang 2019 Table 1)
%   Systemic: R_sys ~ w^(-0.475)  applies to SAR, SC, SVEN
%   Pulmonary: R_pul ~ w^(-0.70)  applies to PAR, PVEN, PCOX, PCNO
params.R.SAR  = params_ref.R.SAR  * w^eR_sys;
params.R.SC   = params_ref.R.SC   * w^eR_sys;
params.R.SVEN = params_ref.R.SVEN * w^eR_sys;
params.R.PAR  = params_ref.R.PAR  * w^eR_pul;
params.R.PCOX = params_ref.R.PCOX * w^eR_pul;
params.R.PCNO = params_ref.R.PCNO * w^eR_pul;
params.R.PVEN = params_ref.R.PVEN * w^eR_pul;
% R.vsd is NOT scaled — assigned by params_from_clinical.m per scenario

%-- B2. Compliances  C ~ w^(+1.00)  [Zhang Table 1]
params.C.SAR  = params_ref.C.SAR  * w^eC;
params.C.SC   = params_ref.C.SC   * w^eC;
params.C.SVEN = params_ref.C.SVEN * w^eC;
params.C.PAR  = params_ref.C.PAR  * w^eC;
params.C.PCOX = params_ref.C.PCOX * w^eC;
params.C.PCNO = params_ref.C.PCNO * w^eC;
params.C.PVEN = params_ref.C.PVEN * w^eC;
if isfield(params_ref.C, 'RA'), params.C.RA = params_ref.C.RA * w^eC; end
if isfield(params_ref.C, 'LA'), params.C.LA = params_ref.C.LA * w^eC; end

%-- B3. Inertances — NOT scaled (excluded per Zhang 2019)
%   L.SAR, L.SVEN, L.PAR, L.PVEN kept at reference values.

%-- B4. Vascular unstressed volumes  V0 ~ w^(+0.80)
%   NOTE: These values are OVERRIDDEN in Section C by BV_scale.
params.V0.SAR  = params_ref.V0.SAR  * w^eV0;
params.V0.SC   = params_ref.V0.SC   * w^eV0;
params.V0.SVEN = params_ref.V0.SVEN * w^eV0;   % will be adjusted in Section C
params.V0.PAR  = params_ref.V0.PAR  * w^eV0;
params.V0.PVEN = params_ref.V0.PVEN * w^eV0;
if isfield(params_ref.V0,'PCOX'), params.V0.PCOX = params_ref.V0.PCOX * w^eV0; end
if isfield(params_ref.V0,'PCNO'), params.V0.PCNO = params_ref.V0.PCNO * w^eV0; end

%-- B5. Valve open resistance  ~ w^(-0.90)  [orifice area ~ w^+1, Zhang Table 1]
params.Rvalve.open   = params_ref.Rvalve.open * w^eRv_op;
params.Rvalve.closed = params_ref.Rvalve.closed;   % unchanged (numerical guard)

% =====================================================================
%  C. INITIAL CONDITION COMPUTATION  (pressure-initialized, blood-volume-conserving)
%
%  Blood volume budget: BV_patient = weight_kg * BV_per_kg
%  All V0 values are re-scaled by BV_scale = BV_patient / BV_adult_ref
%  (blood-volume ratio) to ensure V0.SVEN does not exceed total BV.
%
%  Nominal filling pressures (size-independent physiological reference):
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
    BV_per_kg = 82;    % mL/kg  (paediatric reference, <1 yr)
else
    BV_per_kg = 70;    % mL/kg  (standard adult/older child reference)
end
BV_patient = patient.weight_kg * BV_per_kg;   % [mL] target circulating volume

% Index struct
sidx = params.idx;

% --- Re-scale ALL V0 values by blood-volume ratio ---------------------
%   BV_adult_ref = 70 kg × 70 mL/kg = 4900 mL  (Valenti reference person)
BV_adult_ref = 4900;                         % [mL]
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
P_nom_SAR   = 65;   % [mmHg]
P_nom_SC    = 15;   % [mmHg]
P_nom_SVEN  =  2;   % [mmHg]  CVP — V0.SVEN adjusted to conserve BV
P_nom_PAR   = 15;   % [mmHg]
P_nom_PVEN  =  6;   % [mmHg]
P_nom_PC    =  8;   % [mmHg]
P_nom_LV_ED =  8;   % [mmHg]
P_nom_RV_ED =  4;   % [mmHg]
P_nom_LA_ED =  6;   % [mmHg]
P_nom_RA_ED =  4;   % [mmHg]

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

% Flow states: initialise to cardiac output estimate
SV_assumed_mL = 5;                        % [mL] assumed stroke volume proxy for flow IC -- not patient-specific;
                                          % overwritten at steady state. Source: adult order-of-magnitude estimate.
Q_init = (params.HR / 60) * SV_assumed_mL;   % [mL/s]  initial flow IC seed
ic_p(sidx.Q_SAR)  = Q_init;
ic_p(sidx.Q_SVEN) = Q_init;
ic_p(sidx.Q_PAR)  = Q_init;
ic_p(sidx.Q_PVEN) = Q_init;

% --- Blood conservation: adjust V0.SVEN so sum(V_ic) = BV_patient -----
vol_idx_c = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV sidx.V_SAR ...
             sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];

BV_no_SVEN    = sum(ic_p(vol_idx_c)) - ic_p(sidx.V_SVEN);
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

params.ic.V    = ic_p(:)';
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

fprintf('[apply_scaling] w=%.4f | BV_patient=%.0f mL | BV_scale=%.4f\n', ...
    w, BV_patient, BV_scale);
fprintf('[apply_scaling] IC: V_LV=%.1f  V_RV=%.1f  V_SVEN=%.1f mL  P_SVEN_ic=%.1f mmHg\n', ...
    ic_p(sidx.V_LV), ic_p(sidx.V_RV), ic_p(sidx.V_SVEN), P_nom_SVEN);

end  % apply_scaling
