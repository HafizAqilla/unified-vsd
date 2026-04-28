function clinical = patient_profile_A()
% PATIENT_PROFILE_A
% -----------------------------------------------------------------------
% Virtual Patient Profile A — "The High-Flow Infant" (Typical Surgical Candidate)
%
% Source: Pediatric_VSD_Clinical_Parameters.md, Table 7.1
%   Group II (Stable Infants) cohort + case report data.
%   Large perimembranous VSD with high Qp/Qs (3.44) and moderate PAH.
%   This is the PRIMARY benchmark case per the document recommendation.
%
% DEMOGRAPHICS:
%   Age: 1.6 months   Weight: 3.7 kg   Height: 54 cm   BSA: 0.22 m²
%
% KEY HAEMODYNAMICS (pre-surgery):
%   HR: 150 bpm    VSD diameter: 5.7 mm    Gradient: 40 mmHg
%   PAP mean: 28 mmHg    MAP: 52 mmHg    Qp/Qs: 3.44
%   Profile A SVR absolute = doc SVRi(21.66 WU·m²)/BSA(0.22 m²) = 98.5 WU
%   Note: physiological check MAP/Qs = 52/0.77 = 67.5 WU; discrepancy is
%         in the source document and arises from different CO assumptions.
%         Calibration (DO_CALIBRATION=true) will reconcile this.
%
% USAGE:
%   clinical = patient_profile_A();
%   main_run('pre_surgery', clinical)
%
% REFERENCES:
%   Pediatric_VSD_Clinical_Parameters.md — Table 7.1
%   Sources within that document: NCDR IMPACT, NeoCardio Lab, AHA cohorts.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  1.0
% -----------------------------------------------------------------------

clinical = patient_template();   % initialise all fields to NaN

%% =====================================================================
%  COMMON  — demographics
%% =====================================================================
clinical.common.age_years  = 1.6 / 12;   % [years]  1.6 months
clinical.common.weight_kg  = 3.7;        % [kg]     3.70 ± 0.70 kg
clinical.common.height_cm  = 54;         % [cm]     50th percentile for 1.6-month male
clinical.common.sex        = 1;          % 0 = female, 1 = male — AGENTS.md §3.10
clinical.common.BSA        = 0.22;       % [m²]     0.22 ± 0.02 m² (stated, overrides Mosteller)
clinical.common.HR         = 150;        % [bpm]    Sympathetic compensation; 150 ± 3 bpm in large VSD infants

%% =====================================================================
%  PRE-SURGERY  — catheterisation + echocardiography benchmarks
%% =====================================================================
pre = clinical.pre_surgery;

% ---- Shunt geometry / severity ----------------------------------------
pre.VSD_diameter_mm   = 5.7;    % [mm]     mean perimembranous defect diameter
pre.VSD_gradient_mmHg = 40;     % [mmHg]   peak instantaneous LV-RV gradient (echo Doppler)
%  Q_shunt = Qp - Qs = 2.6 - 0.77 = 1.83 L/min
pre.Q_shunt_Lmin      = 1.83;   % [L/min]  net shunt flow (Qp - Qs)
pre.QpQs              = 3.44;   % [-]      3.44 ± 0.20 (Group II)

% ---- Pulmonary circulation --------------------------------------------
pre.PAP_sys_mmHg      = 35;     % [mmHg]   34.56 ± 8.7
pre.PAP_dia_mmHg      = 15;     % [mmHg]   derived
pre.PAP_mean_mmHg     = 28;     % [mmHg]   28.2 ± 7.25

% PVR = PVRi / BSA = 2.53 / 0.22 = 11.5 WU
pre.PVR_WU            = 11.5;   % [WU]     index 2.53 WU·m²; abs = index/BSA

% ---- Systemic circulation ---------------------------------------------
pre.SAP_sys_mmHg      = 80;     % [mmHg]
pre.SAP_dia_mmHg      = 45;     % [mmHg]
pre.SAP_mean_mmHg     = 52;     % [mmHg]   MAP 52.2 ± 3.71

% SVR = SVRi / BSA = 21.66 / 0.22 = 98.5 WU
% (physiological check: MAP/Qs = 52/0.77 = 67.5 WU — see header note)
pre.SVR_WU            = 98.5;   % [WU]     index 21.66 WU·m²; abs = index/BSA

% ---- Atrial pressures -------------------------------------------------
pre.RAP_mean_mmHg     = 5;      % [mmHg]   consistently 5–6 mmHg in pediatric VSD cohorts
pre.LAP_mean_mmHg     = 10;     % [mmHg]   estimated from LVEDP + mitral gradient;
                                %           LAP ≈ LVEDP + small gradient (≈8+2=10)
                                %           Consistent with LAP=10 in MATLABterbaru profile A

% ---- Ventricular volumes (echo) ---------------------------------------
pre.LVEDV_mL          = 30;     % [mL]     ~150–210% of normal for BSA 0.22 m²; target 30 mL
pre.LVESV_mL          = 10;     % [mL]     derived; EF ≈ (30-10)/30 = 0.67
pre.RVEDV_mL          = 25;     % [mL]     ~125% of normal; Bernheim + septal shift
pre.RVESV_mL          = 12;     % [mL]
pre.LVEF              = 0.67;   % [-]      67.15 ± 4.85% (fraction, not %)

% ---- Cardiac output ---------------------------------------------------
%  Qs based on CI ~ 3.5 L/min/m² × BSA 0.22 = 0.77 L/min
pre.CO_Lmin           = 0.77;   % [L/min]  systemic cardiac output (Qs)

% ---- LV end-diastolic pressure (assumed) ------------------------------
pre.LVEDP_mmHg        = 8;      % [mmHg]   typical infant with large VSD; assumed
                                %           (direct catheter LVEDP not reported)

% ---- IC and elastance override ----------------------------------------
%  Activates override block in params_from_clinical to:
%    (a) compute V0.LV, V0.RV from echo LVESV, RVESV
%    (b) compute E.LV.EB/EA, E.RV.EB/EA from echo volumes + pressures
%    (c) distribute BV_total_mL across all 14 IC states
%  Required because BSA-allometric scaling (s=0.127) alone under-predicts
%  ventricular volumes for this < 1yr infant (BV 314 mL vs adult 4900 mL).
pre.override_IC       = true;   % [-]      activate clinical-echo IC override
pre.BV_total_mL       = 314;    % [mL]     85 mL/kg × 3.7 kg (infant <1yr; Lundquist 2025)

clinical.pre_surgery  = pre;

%% =====================================================================
%  POST-SURGERY  — NaN (not available for this profile)
%  Populate if post-closure follow-up data becomes available.
%% =====================================================================
% All fields remain NaN from patient_template()

end
