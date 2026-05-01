function clinical = patient_template()
% PATIENT_TEMPLATE
% -----------------------------------------------------------------------
% Unified clinical data input struct for the VSD lumped-parameter model.
%
% Fill in all available measurements for your patient.  Leave fields that
% are unavailable as NaN — the validation report will display NaN for those
% rows and exclude them from the RMSE calculation.
%
% The struct has three sub-sections:
%   clinical.common      — applies to both scenarios (demographics, HR)
%   clinical.pre_surgery — measurements before surgical VSD closure
%   clinical.post_surgery— measurements after surgical VSD closure
%
% Pass this struct together with a scenario flag to main_run:
%   main_run('pre_surgery',  clinical)
%   main_run('post_surgery', clinical)
%
% UNITS:  all pressures in mmHg; flows in L/min; resistances in Wood units;
%         volumes in mL; weight in kg; height in cm; BSA in m².
%
% REFERENCES:
%   [1] Clinical data dictionary: docs/clinical_data_dictionary.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

clinical = struct();

%% =====================================================================
%  COMMON  — patient demographics, measured for any scenario
%% =====================================================================
clinical.common.age_years  = NaN;   % [years]  patient age at measurement
clinical.common.weight_kg  = NaN;   % [kg]
clinical.common.height_cm  = NaN;   % [cm]
clinical.common.sex        = 'M';   % 'M' or 'F'
clinical.common.BSA        = NaN;   % [m²]  leave NaN → computed by Mosteller
clinical.common.HR         = NaN;   % [bpm]  resting heart rate
clinical.common.patient_name = '';  % [char] patient label for run folders and reports
clinical.common.patient_id   = '';  % [char] optional identifier / MRN-safe alias
clinical.common.maturation_mode = 'normal';   % 'normal' | 'pvr_fixed_day3' | 'pvr_fixed_day30' | 'none'

%% =====================================================================
%  PRE-SURGERY  — haemodynamics in the presence of the open VSD
%  Source: right-heart catheterisation and echocardiography
%% =====================================================================
pre = struct();

% ---- Shunt geometry / severity ----------------------------------------
pre.VSD_mode          = 'linear_bidirectional';   % 'linear_bidirectional' | 'orifice_bidirectional' | 'linear_left_to_right_only'
pre.VSD_diameter_mm   = NaN;   % [mm]     defect diameter (echo)
pre.VSD_gradient_mmHg = NaN;   % [mmHg]   peak instantaneous gradient (echo Doppler)
pre.Q_shunt_Lmin      = NaN;   % [L/min]  net shunt flow (cath)
pre.QpQs              = NaN;   % [-]      pulmonary-to-systemic flow ratio (cath)

% ---- Pulmonary circulation (elevated in pre-surgery) ------------------
pre.PAP_sys_mmHg      = NaN;   % [mmHg]   systolic PA pressure
pre.PAP_dia_mmHg      = NaN;   % [mmHg]   diastolic PA pressure
pre.PAP_mean_mmHg     = NaN;   % [mmHg]   mean PA pressure (cath)
pre.PVR_WU            = NaN;   % [WU]     pulmonary vascular resistance

% ---- Systemic circulation ---------------------------------------------
pre.SAP_sys_mmHg      = NaN;   % [mmHg]   systolic arterial pressure
pre.SAP_dia_mmHg      = NaN;   % [mmHg]   diastolic arterial pressure
pre.SAP_mean_mmHg     = NaN;   % [mmHg]   mean arterial pressure
pre.SVR_WU            = NaN;   % [WU]     systemic vascular resistance

% ---- Atrial pressures -------------------------------------------------
pre.RAP_mean_mmHg     = NaN;   % [mmHg]   right atrial mean pressure (cath)
pre.LAP_mean_mmHg     = NaN;   % [mmHg]   left atrial mean (or PCWP)

% ---- Ventricular volumes and function (echo / MRI) -------------------
pre.LVEDV_mL          = NaN;   % [mL]     LV end-diastolic volume
pre.LVESV_mL          = NaN;   % [mL]     LV end-systolic volume
pre.RVEDV_mL          = NaN;   % [mL]     RV end-diastolic volume
pre.RVESV_mL          = NaN;   % [mL]     RV end-systolic volume
pre.LVEF              = NaN;   % [-]      LV ejection fraction (fraction, not %)
pre.override_IC       = false; % [-]      enable echo-informed chamber tuning (V0 + E)

% ---- Cardiac output ---------------------------------------------------
pre.CO_Lmin           = NaN;   % [L/min]  systemic cardiac output (Fick or TD)

clinical.pre_surgery = pre;

%% =====================================================================
%  POST-SURGERY  — haemodynamics after VSD closure
%  Source: post-operative catheterisation / echocardiography
%% =====================================================================
post = struct();

% ---- Shunt (absent after closure) ------------------------------------
% QpQs should be ~1.0; residual shunt is modelled by a small, finite R_VSD
post.QpQs             = NaN;   % [-]      ≈ 1.0 expected; set NaN if not measured

% ---- Pulmonary circulation (normalised post-surgery) -----------------
post.PAP_sys_mmHg     = NaN;   % [mmHg]
post.PAP_dia_mmHg     = NaN;   % [mmHg]
post.PAP_mean_mmHg    = NaN;   % [mmHg]
post.PVR_WU           = NaN;   % [WU]     should be lower than pre-surgery

% ---- Systemic circulation --------------------------------------------
post.SAP_sys_mmHg     = NaN;   % [mmHg]
post.SAP_dia_mmHg     = NaN;   % [mmHg]
post.MAP_mmHg         = NaN;   % [mmHg]   mean arterial pressure
post.SVR_WU           = NaN;   % [WU]

% ---- Atrial pressures ------------------------------------------------
post.RAP_mean_mmHg    = NaN;
post.LAP_mean_mmHg    = NaN;

% ---- Ventricular volumes and function (normalised post-surgery) ------
post.LVEDV_mL         = NaN;
post.LVESV_mL         = NaN;
post.RVEDV_mL         = NaN;
post.RVESV_mL         = NaN;
post.LVEF             = NaN;
post.RVEF             = NaN;

% ---- Cardiac output --------------------------------------------------
post.CO_Lmin          = NaN;

clinical.post_surgery = post;

end
