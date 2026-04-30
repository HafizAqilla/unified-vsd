function params = apply_scaling(params_ref, patient)
% APPLY_SCALING
% -----------------------------------------------------------------------
% Patient-specific allometric scaling of cardiovascular model parameters.
%
% Scales the adult reference parameter set (Bozkurt2019/Colebank2025
% baseline, ~70 kg healthy adult male) to a paediatric patient using
% the weight-based allometric power law of Zhang et al. (2019):
%
%   X_child = X_adult * (W_child / W_adult)^beta
%
% where W is body weight [kg] and beta is a parameter-specific exponent.
%
% INPUTS:
%   params_ref  - adult reference struct from default_parameters()
%   patient     - struct with fields:
%                   .age_years   patient age          [years]
%                   .weight_kg   body weight          [kg]
%                   .height_cm   standing height      [cm]
%                   .sex         'M' or 'F'           (reserved)
%
% OUTPUTS:
%   params      - scaled parameter struct, ready for system_rhs.m
%
% -----------------------------------------------------------------------
% LAYER A — ZHANG PHYSIOLOGICAL PARAMETER SCALING (Zhang et al. 2019)
% -----------------------------------------------------------------------
%   Primary scaling variable: body weight W [kg]
%   Scaling ratio:            w = W_patient / W_ref
%
%   Zhang exponents used:
%     HR              ~ w^(-0.30)    heart rate
%     E_LV, E_LA      ~ w^(-0.50)    LV/LA elastance (systemic-side)
%     E_RV, E_RA      ~ w^(-0.75)    RV/RA elastance (pulmonary-side)
%     V0              ~ w^(+0.80)    unstressed volumes (all compartments)
%     R_systemic      ~ w^(-0.475)   systemic vascular resistance
%     R_pulmonary     ~ w^(-0.70)    pulmonary vascular resistance
%     C               ~ w^(+1.00)    vascular compliance (all)
%     R_valve_open    ~ w^(-0.90)    derived: d~w^0.45 → A~w^0.90 → R~w^-0.90
%
%   NOT SCALED (with justification):
%     R.vsd        — pathological; assigned per-patient by params_from_clinical.m
%     Rvalve.closed — numerical guard; must remain large for diode logic
%     epsilon_valve — numerical smoothing parameter, not a physiological quantity
%     L.*          — inertances; no Zhang exponent available (see Note 1)
%
% -----------------------------------------------------------------------
% LAYER B — INITIALIZATION / IC CORRECTION  (numerical stabilization)
% -----------------------------------------------------------------------
%   Purpose: ensure physiologically reachable initial conditions and
%   enforce total circulating blood volume conservation.
%   SEPARATE from Zhang physiological scaling.
%
%   Steps:
%   (1) BV_patient estimated from age-stratified mL/kg formula.
%   (2) IC volumes built from size-independent nominal filling pressures:
%         V_ic = V0_zhang + P_nom * C_zhang
%   (3) V0.SVEN adjusted so that sum(V_ic) = BV_patient.
%         [Overrides the Zhang V0.SVEN from Layer A. Rationale: SVEN is
%          the dominant venous reservoir; blood conservation here is
%          critical for the correct systemic pressure-flow steady state.
%          All other V0 values retain their Zhang-scaled values.]
%   (4) Flow states initialised by CO_adult * w  [mL/s].
%
% -----------------------------------------------------------------------
% NOTE 1 — INERTANCE TREATMENT (modeling assumption):
%   Zhang et al. (2019) provides no explicit inertance exponent.
%   L.SAR, L.SVEN, L.PAR, L.PVEN are LEFT AT ADULT REFERENCE VALUES.
%   Inertance contributes minimally to pressure-flow dynamics at
%   physiological HRs (L·omega << R), so this is conservative and safe.
%   Revisit if high-frequency waveform accuracy is required.
%
% REFERENCES:
%   [1] Zhang X et al. (2019). Allometric scaling of cardiovascular
%       model parameters across body sizes. [Confirm full citation.]
%   [2] Bozkurt S (2019). Math Biosci Eng 16(5):3943-3962. Tables 1-2.
%   [3] Colebank MJ et al. (2025). ASAIO J. Table 1 (male reference).
%   [4] Valenti (2023). Thesis: 0-D cardiovascular model, topology.
%   [5] Guyton AC (1991). Textbook of Medical Physiology, §15.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-10
% VERSION:  2.0  (Zhang 2019 weight-based allometric framework)
% -----------------------------------------------------------------------

%% =====================================================================
%  1. WEIGHT RATIO  — Zhang et al. 2019 primary scaling variable
%% =====================================================================

W_ref = 70;   % [kg]  Adult reference weight — Bozkurt2019/Colebank2025 baseline
               %       (~70 kg healthy adult male, HR=75 bpm, CO=6 L/min)

w = patient.weight_kg / W_ref;   % [-]  dimensionless weight ratio

fprintf('[apply_scaling] W_patient=%.2f kg | W_ref=%.1f kg | w=%.4f\n', ...
    patient.weight_kg, W_ref, w);

%% =====================================================================
%  2. Copy reference params — no field is scaled twice within Layer A
%% =====================================================================

params = params_ref;
params.scaling.W_ref     = W_ref;
params.scaling.W_patient = patient.weight_kg;
params.scaling.w         = w;
params.scaling.patient   = patient;

%% =====================================================================
%  Zhang exponent table  (Zhang et al. 2019)
%  Named constants — never use magic numbers in scaling expressions.
%% =====================================================================

beta_HR      = -0.300;   % Heart rate                              — Zhang 2019
beta_E_LV    = -0.500;   % LV/LA elastance (systemic-side)         — Zhang 2019
beta_E_RV    = -0.750;   % RV/RA elastance (pulmonary-side)        — Zhang 2019
beta_V0      =  0.800;   % Unstressed volumes                      — Zhang 2019
beta_R_sys   = -0.475;   % Systemic vascular resistance            — Zhang 2019
beta_R_pul   = -0.700;   % Pulmonary vascular resistance           — Zhang 2019
beta_C       =  1.000;   % Vascular compliance (all)               — Zhang 2019
beta_R_valve = -0.500;   % Valve open resistance scaling           — Zhang 2019

%% ====================================================================
%  LAYER A — ZHANG PHYSIOLOGICAL PARAMETER SCALING
%% =====================================================================

%-- A1. Heart rate  HR ~ w^(-0.30)  [Zhang 2019]
params.HR = params_ref.HR * w^beta_HR;

%-- A2. Cardiac elastances  (active EA and passive EB)
%   LV and LA — systemic-pressure side: beta = -0.50  [Zhang 2019]
params.E.LV.EA = params_ref.E.LV.EA * w^beta_E_LV;
params.E.LV.EB = params_ref.E.LV.EB * w^beta_E_LV;
params.E.LA.EA = params_ref.E.LA.EA * w^beta_E_LV;
params.E.LA.EB = params_ref.E.LA.EB * w^beta_E_LV;
%   RV and RA — pulmonary-pressure side: beta = -0.75  [Zhang 2019]
params.E.RV.EA = params_ref.E.RV.EA * w^beta_E_RV;
params.E.RV.EB = params_ref.E.RV.EB * w^beta_E_RV;
params.E.RA.EA = params_ref.E.RA.EA * w^beta_E_RV;
params.E.RA.EB = params_ref.E.RA.EB * w^beta_E_RV;

%-- A3. Unstressed volumes  V0 ~ w^(+0.80)  [Zhang 2019]
%   NOTE: V0.SVEN is the dominant venous reservoir. It will be further
%   adjusted by the blood conservation step in Layer B (Section B4).
%   All other V0 compartments retain these Zhang-scaled values unchanged.
params.V0.LV   = params_ref.V0.LV   * w^beta_V0;
params.V0.RV   = params_ref.V0.RV   * w^beta_V0;
params.V0.LA   = params_ref.V0.LA   * w^beta_V0;
params.V0.RA   = params_ref.V0.RA   * w^beta_V0;
params.V0.SAR  = params_ref.V0.SAR  * w^beta_V0;
params.V0.SC   = params_ref.V0.SC   * w^beta_V0;
params.V0.SVEN = params_ref.V0.SVEN * w^beta_V0;   % provisional — adjusted in B4
params.V0.PAR  = params_ref.V0.PAR  * w^beta_V0;
params.V0.PVEN = params_ref.V0.PVEN * w^beta_V0;
if isfield(params_ref.V0,'PCOX'), params.V0.PCOX = params_ref.V0.PCOX * w^beta_V0; end
if isfield(params_ref.V0,'PCNO'), params.V0.PCNO = params_ref.V0.PCNO * w^beta_V0; end

%-- A4. Systemic resistances  R ~ w^(-0.475)  [Zhang 2019]
params.R.SAR  = params_ref.R.SAR  * w^beta_R_sys;
params.R.SC   = params_ref.R.SC   * w^beta_R_sys;
params.R.SVEN = params_ref.R.SVEN * w^beta_R_sys;
% R.vsd — NOT scaled: pathological, assigned per-patient by params_from_clinical.m

%-- A5. Pulmonary resistances  R ~ w^(-0.70)  [Zhang 2019]
params.R.PAR  = params_ref.R.PAR  * w^beta_R_pul;
params.R.PCOX = params_ref.R.PCOX * w^beta_R_pul;
% PCNO branch is disabled in default_parameters.m (R.PCNO=1e6, C.PCNO=1e-8).
% Scaling applied for structural consistency; no hemodynamic effect because
% the branch carries virtually zero flow at any realistic R.PCNO.
params.R.PCNO = params_ref.R.PCNO * w^beta_R_pul;
params.R.PVEN = params_ref.R.PVEN * w^beta_R_pul;

%-- A6. Vascular compliances  C ~ w^(+1.0)  [Zhang 2019]
params.C.SAR  = params_ref.C.SAR  * w^beta_C;
params.C.SC   = params_ref.C.SC   * w^beta_C;
params.C.SVEN = params_ref.C.SVEN * w^beta_C;
params.C.PAR  = params_ref.C.PAR  * w^beta_C;
params.C.PCOX = params_ref.C.PCOX * w^beta_C;
params.C.PCNO = params_ref.C.PCNO * w^beta_C;   % disabled branch — scaling harmless
params.C.PVEN = params_ref.C.PVEN * w^beta_C;
% Atrial compliance legacy fields (derived from EB in default_parameters.m)
if isfield(params_ref.C,'RA'), params.C.RA = params_ref.C.RA * w^beta_C; end
if isfield(params_ref.C,'LA'), params.C.LA = params_ref.C.LA * w^beta_C; end

%-- A7. Valve open resistance  [derived from Zhang 2019]
%   Zhang beta_diam = 0.45 → orifice area A ∝ d^2 ∝ w^0.90
%   Orifice resistance R ∝ 1/A → R_valve_open ~ w^(-0.90)
%   Rvalve.closed — NOT scaled: numerical guard (large R for closed state)
params.Rvalve.open   = params_ref.Rvalve.open * w^beta_R_valve;
params.Rvalve.closed = params_ref.Rvalve.closed;   % numerical guard — unchanged

%-- A8. Inertances — NOT scaled  (modeling assumption; see header Note 1)
%   L.SAR, L.SVEN, L.PAR, L.PVEN remain at the adult reference values.

%-- A9. epsilon_valve — NOT scaled  (numerical continuity parameter)
%   params.epsilon_valve unchanged.

fprintf('[apply_scaling] Layer A done: HR=%.1f bpm | R.SAR=%.4f | R.PAR=%.5f | C.SVEN=%.3f\n', ...
    params.HR, params.R.SAR, params.R.PAR, params.C.SVEN);

%% =====================================================================
%  LAYER B — INITIALIZATION / IC CORRECTION
%  Numerically motivated. Not derived from Zhang et al.
%  Constructs physiologically consistent ICs and enforces blood
%  volume conservation by adjusting V0.SVEN only.
%% =====================================================================

%% B1. Blood volume budget  (age-stratified formula)
age_years = patient.age_years;
if age_years < 1
    BV_per_kg = 82;    % [mL/kg]  neonatal reference — Linderkamp et al. 1977
else
    BV_per_kg = 70;    % [mL/kg]  child/adult reference — standard physiology
end
BV_patient = patient.weight_kg * BV_per_kg;   % [mL]  target circulating volume

fprintf('[apply_scaling] BV_patient=%.0f mL  (%.0f mL/kg x %.1f kg)\n', ...
    BV_patient, BV_per_kg, patient.weight_kg);

%% B2. Nominal filling pressures  (size-independent physiological targets)
%   These are not scaled quantities; they represent normal resting values.
%   Reference: standard clinical cardiology / physiology textbooks.
P_nom_SAR   = 93;   % [mmHg]  MAP (mean arterial pressure)
P_nom_SC    = 15;   % [mmHg]  systemic arteriolar / capillary pressure
P_nom_SVEN  =  2;   % [mmHg]  CVP (central venous pressure)
P_nom_PAR   = 15;   % [mmHg]  pulmonary artery diastolic pressure
P_nom_PVEN  =  6;   % [mmHg]  pulmonary venous pressure (≈ LAP)
P_nom_PC    =  8;   % [mmHg]  pulmonary capillary pressure (pressure state)
P_nom_LV_ED =  8;   % [mmHg]  LV end-diastolic pressure (LVEDP)
P_nom_RV_ED =  4;   % [mmHg]  RV end-diastolic pressure (RVEDP)
P_nom_LA_ED =  6;   % [mmHg]  left atrial pressure (LAP)
P_nom_RA_ED =  4;   % [mmHg]  right atrial pressure (RAP)

%% B3. Build IC vector from nominal pressures and Zhang-scaled parameters
sidx = params.idx;
ic_p = zeros(14, 1);

% Vascular volume states: V_ic = V0_zhang + P_nom * C_zhang
ic_p(sidx.V_SAR)  = params.V0.SAR  + P_nom_SAR  * params.C.SAR;
ic_p(sidx.V_SC)   = params.V0.SC   + P_nom_SC   * params.C.SC;
ic_p(sidx.V_SVEN) = params.V0.SVEN + P_nom_SVEN * params.C.SVEN;  % provisional
ic_p(sidx.V_PAR)  = params.V0.PAR  + P_nom_PAR  * params.C.PAR;
ic_p(sidx.V_PVEN) = params.V0.PVEN + P_nom_PVEN * params.C.PVEN;
ic_p(sidx.P_PC)   = P_nom_PC;   % direct pressure state [mmHg]

% Cardiac chamber volumes: P = E_EB*(V-V0) → V = V0 + P/E_EB
ic_p(sidx.V_LV) = params.V0.LV + P_nom_LV_ED / params.E.LV.EB;
ic_p(sidx.V_RV) = params.V0.RV + P_nom_RV_ED / params.E.RV.EB;
ic_p(sidx.V_LA) = params.V0.LA + P_nom_LA_ED / params.E.LA.EB;
ic_p(sidx.V_RA) = params.V0.RA + P_nom_RA_ED / params.E.RA.EB;

% Flow states: scale adult reference CO by weight ratio (linear, conservative)
%   CO_adult = 100 mL/s (= 6.0 L/min; Colebank2025 baseline, 75 bpm, SV=80 mL)
%   Scaling by w (linear weight ratio) ensures physiologically bounded Q_init.
%   This is an IC approximation only — the ODE converges within 1-2 cycles.
CO_adult_ref_mLs = 100;               % [mL/s]  adult reference cardiac output
Q_init = CO_adult_ref_mLs * w;        % [mL/s]  patient scaled CO estimate
ic_p(sidx.Q_SAR)  = Q_init;
ic_p(sidx.Q_SVEN) = Q_init;
ic_p(sidx.Q_PAR)  = Q_init;
ic_p(sidx.Q_PVEN) = Q_init;

fprintf('[apply_scaling] Q_init=%.1f mL/s  (CO_est=%.2f L/min)\n', ...
    Q_init, Q_init * 60/1000);

%% B4. Blood conservation: adjust V0.SVEN so that sum(V_ic) = BV_patient
%
%   RATIONALE FOR V0.SVEN OVERRIDE:
%     SVEN holds ~60-70% of circulating blood volume. If its Zhang-scaled
%     V0.SVEN causes the sum of all initial volumes to deviate from
%     BV_patient, the systemic venous pressure will be incorrect at t=0.
%     By adjusting V0.SVEN alone (holding P_nom_SVEN = 2 mmHg fixed),
%     blood conservation is enforced without distorting any other
%     compartment's pressure or volume at initialisation.
%
%     All other V0 fields retain their Zhang-scaled values from Layer A.
%     The Zhang V0.SVEN (pre-override) is stored in params.scaling for
%     diagnostics and reproducibility.

params.scaling.V0_SVEN_zhang = params.V0.SVEN;   % store Layer A value

vol_idx = [sidx.V_RA  sidx.V_RV   sidx.V_LA   sidx.V_LV ...
           sidx.V_SAR sidx.V_SC   sidx.V_SVEN  sidx.V_PAR sidx.V_PVEN];

BV_no_SVEN    = sum(ic_p(vol_idx)) - ic_p(sidx.V_SVEN);   % all except SVEN
V_SVEN_target = BV_patient - BV_no_SVEN;

if V_SVEN_target < params.V0.SVEN * 0.3
    warning('apply_scaling:BVbudget', ...
        ['BV_patient (%.0f mL) leaves very little room for V_SVEN ' ...
         '(target=%.0f mL, Zhang V0.SVEN=%.0f mL). ' ...
         'Clamping to 30%% of Zhang V0.SVEN. ' ...
         'Patient-specific calibration is strongly recommended.'], ...
        BV_patient, V_SVEN_target, params.V0.SVEN);
    V_SVEN_target = max(V_SVEN_target, params.V0.SVEN * 0.3);
end

% Adjust V0.SVEN to maintain P_nom_SVEN at the corrected SVEN volume:
%   V_SVEN_target = V0.SVEN_adjusted + P_nom_SVEN * C.SVEN
params.V0.SVEN    = V_SVEN_target - P_nom_SVEN * params.C.SVEN;
ic_p(sidx.V_SVEN) = V_SVEN_target;

params.ic.V    = ic_p(:)';   % row vector (ode15s accepts row or column)
params.scaling.BV_patient       = BV_patient;
params.scaling.V0_SVEN_adjusted = params.V0.SVEN;

%% =====================================================================
%  D. ABSOLUTE TIMING FIELDS  (required by elastance_model.m)
%     Fractional timing params (stored in default_parameters.m) are
%     converted to absolute seconds using the Zhang-scaled HR.
%     Timing fractions are dimensionless and size-independent by design.
%% =====================================================================

T_HB = 60 / params.HR;   % [s]  cardiac cycle period from scaled HR

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

%% =====================================================================
%  Final diagnostics
%% =====================================================================

fprintf('[apply_scaling] Final: HR=%.1f bpm | T_HB=%.3f s | BV_patient=%.0f mL\n', ...
    params.HR, T_HB, BV_patient);
fprintf('[apply_scaling] V_ic:  V_LV=%.1f  V_RV=%.1f  V_SAR=%.1f  V_SVEN=%.1f mL\n', ...
    ic_p(sidx.V_LV), ic_p(sidx.V_RV), ic_p(sidx.V_SAR), ic_p(sidx.V_SVEN));
fprintf('[apply_scaling] V0.SVEN: Zhang=%.1f mL --> BV-adjusted=%.1f mL\n', ...
    params.scaling.V0_SVEN_zhang, params.V0.SVEN);

end  % apply_scaling