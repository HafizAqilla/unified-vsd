function clinical = patient_profile_B()
% PATIENT_PROFILE_B
% -----------------------------------------------------------------------
% Virtual Patient Profile B — "The High-Pressure Infant" (Critical / Pre-S1P)
%
% Source: Pediatric_VSD_Clinical_Parameters.md, Table 7.2
%   Group I (Critical Infants) cohort.
%   Large non-restrictive VSD with near-systemic PAP, low Qp/Qs (1.79),
%   and systemic vasodilation (low SVR). Represents pulmonary hypertension
%   precursor — Eisenmenger physiology precursor.
%
% DEMOGRAPHICS:
%   Age: 2.0 months   Weight: 4.5 kg   BSA: 0.25 m²
%
% KEY HAEMODYNAMICS (pre-surgery):
%   HR: 150 bpm    VSD gradient: ~5 mmHg (near-equalised pressures)
%   PAP mean: 43 mmHg    MAP: 48 mmHg    Qp/Qs: 1.79
%   SVR absolute = doc SVRi(10.76 WU·m²)/BSA(0.25 m²) = 43 WU
%   PVR absolute = doc PVRi(3.11 WU·m²)/BSA(0.25 m²) = 12.44 WU
%
%   ** MODEL TRAP: Do NOT substitute Profile A SVR (98.5 WU) into Profile B.
%      Profile B has systemic vasodilation — SVR must be 43 WU to correctly
%      reproduce MAP = 48 mmHg under cardiogenic stress. **
%
% USAGE:
%   clinical = patient_profile_B();
%   main_run('pre_surgery', clinical)
%
% REFERENCES:
%   Pediatric_VSD_Clinical_Parameters.md — Table 7.2
%   Sources within that document: NCDR IMPACT, AHA cohorts (nitroprusside/
%   hydralazine studies), NeoCardio Lab catheter data.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  1.0
% -----------------------------------------------------------------------

clinical = patient_template();   % initialise all fields to NaN

%% =====================================================================
%  COMMON  — demographics
%% =====================================================================
clinical.common.age_years  = 2.0 / 12;   % [years]  2.0 months
clinical.common.weight_kg  = 4.5;        % [kg]     4.54 ± 2.36 kg
% Height estimated from BSA: h = 3600 × BSA² / w = 3600 × 0.0625 / 4.5 = 50 cm
clinical.common.height_cm  = 50;         % [cm]     estimated from Mosteller given BSA & weight
clinical.common.sex        = 1;          % 0 = female, 1 = male — assumed; not specified in source
clinical.common.BSA        = 0.25;       % [m²]     0.25 ± 0.03 m² (stated, overrides Mosteller)
clinical.common.HR         = 150;        % [bpm]    150 ± 3 bpm in large VSD infants (sympathetic activation)

%% =====================================================================
%  PRE-SURGERY  — catheterisation + echocardiography benchmarks
%% =====================================================================
pre = clinical.pre_surgery;

% ---- Shunt geometry / severity ----------------------------------------
% Gradient < 10 mmHg (near-equalised pressures, non-restrictive VSD)
pre.VSD_gradient_mmHg = 5;      % [mmHg]   representative of < 10 mmHg (source document)
% Qs derived: (MAP - RAP)/SVR = (48-5)/43 = 1.0 L/min
% Qp = QpQs × Qs = 1.79 × 1.0 = 1.79 L/min
% Q_shunt = Qp - Qs = 0.79 L/min
pre.Q_shunt_Lmin      = 0.79;   % [L/min]  derived from QpQs, MAP, SVR
pre.QpQs              = 1.79;   % [-]      1.79 ± 0.73 (Group I)
pre.VSD_diameter_mm   = NaN;    % [mm]     not specified in source document

% ---- Pulmonary circulation (near-systemic PAP) ------------------------
pre.PAP_sys_mmHg      = 63;     % [mmHg]   62.57 ± 21.0
% PAP_dia from PAP_mean: PAP_mean = (PAP_sys + 2*PAP_dia)/3
% → PAP_dia = (3*43 - 63)/2 = 33 mmHg
pre.PAP_dia_mmHg      = 33;     % [mmHg]   derived from mean and systolic
pre.PAP_mean_mmHg     = 43;     % [mmHg]   43.20 ± 9.56

% PVR = PVRi / BSA = 3.11 / 0.25 = 12.44 WU
pre.PVR_WU            = 12.44;  % [WU]     index 3.11 WU·m²; abs = index/BSA

% ---- Systemic circulation (vasodilated) --------------------------------
% SAP estimated: SAP_sys ≈ PAP_sys (non-restrictive) ≈ 65–70 mmHg
% SAP_sys=70, SAP_dia=38 → MAP = (70 + 2*38)/3 = 48.7 ≈ 48 mmHg ✓
pre.SAP_sys_mmHg      = 70;     % [mmHg]   estimated; ≈ PAP_sys for non-restrictive VSD
pre.SAP_dia_mmHg      = 38;     % [mmHg]   estimated from MAP constraint
pre.SAP_mean_mmHg     = 48;     % [mmHg]   MAP 48.40 ± 8.11 (directly stated)

% SVR = SVRi / BSA = 10.76 / 0.25 = 43.0 WU
% ** Low SVR reflects systemic vasodilation / early shock. **
pre.SVR_WU            = 43.0;   % [WU]     index 10.76 WU·m²; abs = index/BSA

% ---- Atrial pressures -------------------------------------------------
pre.RAP_mean_mmHg     = 5;      % [mmHg]   5–6 mmHg typical for simple VSD
pre.LAP_mean_mmHg     = 12;     % [mmHg]   elevated due to near-systemic PAP + high Qp/Qs;
                                %           LVEDP (10) + mitral inflow gradient (≈2 mmHg)

% ---- Ventricular volumes (echo) — estimated from physiology -----------
% Normal LVEDV for BSA 0.25 m² ≈ 20-25 mL
% Moderate volume overload (Qp/Qs 1.79) → estimated 30-35 mL
pre.LVEDV_mL          = 32;     % [mL]     ~130% of normal; moderate LV dilation
pre.LVESV_mL          = 11;     % [mL]     derived; EF ≈ (32-11)/32 = 0.66
pre.RVEDV_mL          = 32;     % [mL]     pressure-loaded by systemic-level PAP; Bernheim effect
pre.RVESV_mL          = 16;     % [mL]     RV EF ≈ 50% (impaired by pressure overload)
pre.LVEF              = 0.66;   % [-]      fraction; within hyperdynamic range for volume overload

% ---- Cardiac output ---------------------------------------------------
% Qs = (MAP - RAP) / SVR = (48-5) / 43 = 1.0 L/min
pre.CO_Lmin           = 1.0;    % [L/min]  systemic cardiac output (Qs)

% ---- LV end-diastolic pressure (assumed) ------------------------------
pre.LVEDP_mmHg        = 10;     % [mmHg]   elevated due to near-systemic PAP + Bernheim;
                                %           assumed (direct LVEDP not reported)

% ---- IC and elastance override ----------------------------------------
pre.override_IC       = true;   % [-]      activate clinical-echo IC override
pre.BV_total_mL       = 383;    % [mL]     85 mL/kg × 4.5 kg (infant <1yr; Lundquist 2025)

clinical.pre_surgery  = pre;

%% =====================================================================
%  POST-SURGERY  — NaN (not available for this profile)
%% =====================================================================
% All fields remain NaN from patient_template()

end
