function clinical = patient_reyna()
% PATIENT_REYNA
% -----------------------------------------------------------------------
% Reyna clinical profile for the VSD lumped-parameter model.
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
clinical.common.patient_name = 'reyna'; % [char] patient label for run folders

%% =====================================================================
%  COMMON — patient demographics, measured for any scenario
%% =====================================================================
clinical.common.age_years  = 3.17;    % [years] 3 years 2 months [cite: 80]
clinical.common.weight_kg  = 13.4;    % [kg] protocol row 4
clinical.common.height_cm  = 95.0;    % [cm] protocol row 3
clinical.common.sex        = 0;       % 0 = female, 1 = male — AGENTS.md §3.10
clinical.common.BSA        = 0.588;   % [m^2] protocol row 5
clinical.common.HR         = 119;     % [bpm] [cite: 81]

%% =====================================================================
%  PRE-SURGERY — haemodynamics in the presence of the open VSD
%  Source: right-heart catheterisation and echocardiography
%% =====================================================================
pre = struct();

% ---- Shunt geometry / severity ----------------------------------------
% VSD diameter: protocol row 7 — LV side 4.46–6.43 mm, RV side 3.63–3.70 mm.
% Use RV-side mean as the effective orifice diameter (smaller = more restrictive).
pre.VSD_diameter_mm   = 3.665;   % [mm] mean RV-side range: (3.63+3.70)/2
pre.VSD_gradient_mmHg = 69;      % [mmHg] peak systolic gradient: LV 94 - RV 25 (row 9)
pre.Q_shunt_Lmin      = 0.664;   % [L/min] Qp - Qs = 4.087 - 3.423 (rows 21–22)
pre.QpQs              = 1.194;   % [-]  protocol row 23
pre.VSD_mode          = 'orifice_bidirectional'; % [-] geometry/gradient-based shunt mode

% ---- Pulmonary circulation — catheterization (rows 16–18, 3 repeated measures)
pre.PAP_sys_mmHg      = 20;      % [mmHg] protocol row 16, repeated 20/20/20
pre.PAP_dia_mmHg      = 10;      % [mmHg] protocol row 17, repeated 10/10/10
pre.PAP_mean_mmHg     = 15;      % [mmHg] mean of 14/15/15 mmHg → 15
pre.PVR_WU            = NaN;     % [WU] protocol row 24 blank; not used as clinical target

% ---- Systemic circulation — catheterization (Right Femoral Artery, rows 12–13)
% MAP is recomputed from catheter sys/dia using MAP = dia + (sys-dia)/3.
% The protocol MAP of 95 mmHg used NIBP cuff (119/83) — different method, not used here.
% NIBP reference only: sys=119, dia=83, MAP_nibp=95 mmHg (row 10–11, 20)
pre.SAP_sys_mmHg      = 100;     % [mmHg] RFA catheter systolic  (row 12)
pre.SAP_dia_mmHg      = 57;      % [mmHg] RFA catheter diastolic (row 13)
pre.SAP_mean_mmHg     = 71.3;    % [mmHg] recomputed: 57 + (100-57)/3 = 71.3
                                  %        Source: catheter RFA, consistent with sys/dia above

pre.SVR_WU            = NaN;     % [WU] protocol row 25 blank; not used as clinical target

% ---- Atrial and ventricular filling pressures -------------------------
pre.RAP_mean_mmHg     = 5;       % [mmHg] catheter, mean of 5/5/5 mmHg (row 19)
pre.LAP_mean_mmHg     = NaN;      
pre.LVEDP_mmHg        = NaN;

% ---- Ventricular volumes and ejection fraction -----------------------
% The available LV/RV volume and EF block was confirmed to be H+1 after
% surgery, so it is not a valid pre-surgery calibration target.
pre.LVEDV_mL          = NaN;     % [mL] unavailable pre-surgery
pre.LVESV_mL          = NaN;     % [mL] unavailable pre-surgery
pre.RVEDV_mL          = NaN;     % [mL] unavailable pre-surgery
pre.RVESV_mL          = NaN;     % [mL] unavailable pre-surgery
pre.EF                = NaN;     % [-] unavailable pre-surgery

% ---- IC override flag -------------------------------------------------
% Do not tune chamber elastance/V0 from H+1 post-operative echo volumes.
pre.override_IC       = false;
pre.CO_comparator     = 'Qs_Lmin'; % [-] compare model systemic flow with protocol-derived Qs
pre.CO_uncertainty_Lmin = 0.50;    % [L/min] Fick/derived Qs uncertainty allowance

% ---- Cardiac output ---------------------------------------------------
% Qs is back-calculated from protocol Qp and Qp/Qs:
%   Qs = 4.087 / 1.194 = 3.423 L/min.
% The model reports CO_Lmin as Qs. LVCO_Lmin is reported separately as
% LVSV * HR / 1000 and should approximate Qp in VSD.
%
% We calibrate to Qs (3.423) as the CO target because:
%   - Qp and Qp/Qs are catheter/Fick entries, and Qs follows directly from them
%   - The H+1 post-operative echo volume block is excluded from pre-operative
%     fitting, so the pre-surgery objective is hemodynamic-only.
pre.CO_Lmin           = 3.423;   % [L/min] Qs = Qp/QpQs = 4.087/1.194 (rows 21 and 23)

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
post.EF               = NaN;
post.RVEF             = NaN;

% ---- Cardiac output --------------------------------------------------
post.CO_Lmin          = NaN;

clinical.post_surgery = post;

end
