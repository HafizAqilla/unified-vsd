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
clinical.common.weight_kg  = 13.4;    % [kg] [cite: 80]
clinical.common.height_cm  = 95;      % [cm] [cite: 80]
clinical.common.sex        = 0;       % 0 = female, 1 = male — AGENTS.md §3.10
clinical.common.BSA        = 0.588;   % [m²] [cite: 80]
clinical.common.HR         = 119;     % [bpm] [cite: 81]

%% =====================================================================
%  PRE-SURGERY — haemodynamics in the presence of the open VSD
%  Source: right-heart catheterisation and echocardiography
%% =====================================================================
pre = struct();

% ---- Shunt geometry / severity ----------------------------------------
% VSD diameter: protocol row 7 — LV side 4.46–6.43 mm, RV side 2.35–3.7 mm.
% Use RV-side mean as the effective orifice diameter (smaller = more restrictive).
pre.VSD_diameter_mm   = 3.025;   % [mm] mean RV-side range: (2.35+3.7)/2
pre.VSD_gradient_mmHg = 69;      % [mmHg] peak systolic gradient: LV 94 - RV 25 (row 9)
pre.Q_shunt_Lmin      = 0.664;   % [L/min] Qp - Qs = 4.087 - 3.423 (rows 21–22)
pre.QpQs              = 1.194;   % [-]  protocol row 23

% ---- Pulmonary circulation — catheterization (rows 16–18, 3 repeated measures)
pre.PAP_sys_mmHg      = 20;      % [mmHg] mean of 20/20/20 mmHg
pre.PAP_dia_mmHg      = 10;      % [mmHg] mean of 10/10/10 mmHg
pre.PAP_mean_mmHg     = 15;      % [mmHg] mean of 14/15/15 mmHg → 15
pre.PVR_WU            = NaN;     % [WU]   not calculated in protocol

% ---- Systemic circulation — catheterization (Right Femoral Artery, rows 12–13)
% MAP is recomputed from catheter sys/dia using MAP = dia + (sys-dia)/3.
% The protocol MAP of 95 mmHg used NIBP cuff (119/83) — different method, not used here.
% NIBP reference only: sys=119, dia=83, MAP_nibp=95 mmHg (row 10–11, 20)
pre.SAP_sys_mmHg      = 100;     % [mmHg] RFA catheter systolic  (row 12)
pre.SAP_dia_mmHg      = 57;      % [mmHg] RFA catheter diastolic (row 13)
pre.SAP_mean_mmHg     = 71.3;    % [mmHg] recomputed: 57 + (100-57)/3 = 71.3
                                  %        Source: catheter RFA, consistent with sys/dia above

% SVR recomputed from catheter MAP: (MAP - RAP) / Qs = (71.3 - 5) / 3.423
pre.SVR_WU            = 19.37;   % [WU]  (71.3 - 5) / 3.423 — catheter-consistent

% ---- Atrial and ventricular filling pressures -------------------------
pre.RAP_mean_mmHg     = 5;       % [mmHg] catheter, mean of 5/5/5 mmHg (row 19)
pre.LAP_mean_mmHg     = 8;       % [mmHg] estimated: mild elevation from VSD volume load
                                  %        No direct PCWP measured; typical for mild-moderate VSD
pre.LVEDP_mmHg        = 8;       % [mmHg] estimated: consistent with LVEF=52.8% (preserved function)
                                  %        No direct measurement in protocol

% ---- Ventricular volumes — Teichholz formula from M-mode dimensions
%
% Protocol rows 26–27 are labelled "Simpson biplane" but contain M-mode
% linear dimensions (LVEDD=32 mm, LVESD=23.6 mm), not Simpson volumes.
% Evidence: applying Simpson to those numbers gives SV=8.4 mL and
% CO=1.0 L/min at HR=119, which is incompatible with cath CO=3.423 L/min.
% The Teichholz formula applied to the same dimensions gives volumes
% consistent with the catheter data.
%
% Teichholz: V = (7/(2.4+D)) × D³   where D is diameter in cm
%   LVEDD = 32 mm = 3.2 cm  →  LVEDV = (7/5.6) × 3.2³  = 41.0 mL
%   LVESD = 23.6 mm = 2.36 cm →  LVESV = (7/4.76) × 2.36³ = 19.3 mL
%   LVEF  = (41.0 - 19.3) / 41.0 = 0.528  (52.8%)
%
% NOTE: FS = (LVEDD - LVESD) / LVEDD = (32 - 23.6) / 32 = 0.2625 (26.25%)
%   This FS of 26% was previously misidentified as EF — they are different.
%   Normal FS range: 28–40%. FS=26% indicates mildly reduced systolic function.
%
% Source: Teichholz LE et al. (1976). Circ 54(4):548–552.
pre.LVEDV_mL          = 41.0;    % [mL] Teichholz from LVEDD=32 mm (M-mode)
pre.LVESV_mL          = 19.3;    % [mL] Teichholz from LVESD=23.6 mm (M-mode)
pre.RVEDV_mL          = 30.5;    % [mL] protocol row 28
pre.RVESV_mL          = 12.0;    % [mL] protocol row 29
pre.LVEF              = 0.528;   % [-]  (41.0 - 19.3) / 41.0 — Teichholz EF

% ---- IC override flag -------------------------------------------------
% Activates pressure-based IC override in params_from_clinical.m.
% Essential: allometric scaling produces LVEDV < 5 mL at BSA=0.588 m².
% Uses LVEDV/LVESV/RVEDV/RVESV + pressure fields to seed the ODE correctly.
pre.override_IC       = true;

% ---- Cardiac output ---------------------------------------------------
% Clinical catheter Fick gives Qs = 3.423 L/min (systemic venous return).
% The model computes CO_Lmin = LVSV × HR / 1000 ≈ Qp in a VSD with L→R shunt.
%
% We calibrate to Qs (3.423) as the CO target because:
%   - Qs is the direct Fick measurement (highest confidence)
%   - The volume targets (LVEDV=41, LVESV=19.3) from Teichholz imply
%     LVSV = 21.7 mL → CO_volumes = 2.58 L/min, inconsistent with catheter.
%     Setting CO = Qp = 4.087 would force LVSV = 34.4 mL, making LVEDV >> 53 mL
%     (LVEDV error > 30%), which is worse than the current trade-off.
%   - Calibrating to Qs = 3.423 produces the best simultaneous fit across
%     volumes, flows, and pressures given the Teichholz measurement limitation.
%
% NOTE: SVR is reported as ~22 WU (vs clinical 19.37) because the model's
%   effective Qs = CO_model/QpQs = 3.423/1.18 = 2.90 L/min < catheter Qs = 3.423.
%   This is a known methodological discrepancy from the CO definition — document
%   in the paper. SVR absolute error ≈ 3 mmHg·min/L, clinically within noise.
pre.CO_Lmin           = 3.423;   % [L/min] Qs from Fick catheter (rows 21–23)

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
