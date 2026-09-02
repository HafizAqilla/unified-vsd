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
% Source of record for demographics and HR: the study "reyna" procedure log
% (dated 06/04/2026) — the catheterisation session that produced every pre-
% and post-closure pressure below. Full provenance (facility, MRN, case ID)
% is kept out of this tracked file per AGENTS.md Section 9.2; see
% config/private/patient_provenance.local.m (gitignored) or its checked-in
% template config/private/patient_provenance.local.m.example.
%
% These superseded an earlier "Keisya 2026-05-11 revision" (14.0 kg, 98.0 cm,
% BSA 0.6173 by Mosteller, HR 119). That revision is dated five weeks AFTER
% this catheterisation: the child had grown, so pairing May anthropometry
% with April haemodynamics mis-scaled every demographically scaled parameter.
% The measurement-day values are the correct ones for this fit.
clinical.common.age_years  = 3.17;    % [years] 3 years 2 months [cite: 80]
clinical.common.weight_kg  = 13.4;    % [kg] procedure log 06/04/2026 07.53.01
clinical.common.height_cm  = 95.0;    % [cm] procedure log 06/04/2026 07.53.09
clinical.common.sex        = 0;       % 0 = female, 1 = male — AGENTS.md §3.10
% BSA as stamped by the hospital system (DuBois: 0.007184*95^0.725*13.4^0.425
% = 0.5879). Note the prior config value used Mosteller instead; the stamped
% value is retained here so the model matches the source record exactly.
clinical.common.BSA        = 0.588;   % [m^2] procedure log 06/04/2026 07.53.20
% HR from the procedure log's own pulse row. The prior value of 119 coincides
% exactly with the NIBP SYSTOLIC on the adjacent log line (NIBP 119/83 (95)),
% which is the likely origin of the error. HR sets cycle length, so this is a
% model input, not just a reporting field: 60/119 = 0.504 s vs 60/136 = 0.441 s.
clinical.common.HR         = 136;     % [bpm] procedure log 06/04/2026 09.24.57 "Nadi 136 bpm"

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
pre.SAP_sys_mmHg      = 100;     % [mmHg] RFA catheter systolic  (row 12; log 10.39.12)
pre.SAP_dia_mmHg      = 57;      % [mmHg] RFA catheter diastolic (row 13; log 10.39.12)
% MEAN: use the catheter's OWN stamped mean, not a form-factor reconstruction.
% The procedure log records this reading as "RFA 100/57 (77)" — the transducer
% reports 77 mmHg directly. The previous value of 71.3 came from applying
% MAP = dia + (sys-dia)/3, which assumes a form factor this patient's waveform
% does not have; the same ~5-6 mmHg offset recurs post-closure (formula 75 vs
% stamped 79), so it is systematic, not noise.
% Three MAP candidates existed: NIBP cuff 95 (different method, rejected),
% form-factor 71.3 (reconstructed, rejected), catheter-stamped 77 (used).
pre.SAP_mean_mmHg     = 77;      % [mmHg] procedure log 06/04/2026 10.39.12 "RFA 100/57 (77)"

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

% ======================================================================
% SOURCE: the study "reyna" procedure log, 06/04/2026 (see
% config/private/patient_provenance.local.m for full provenance).
% Same catheterisation session as the pre-surgery block above.
% The VSD closure device was placed at 11.50.19 and released at 12.06.23;
% every value below is stamped AFTER that release (12.15-12.32), with the
% patient under the same anaesthesia and ventilator settings as the
% pre-closure readings. This is a genuine paired pre/post dataset rather
% than two separate studies, which is why the two states are directly
% comparable.
%
% Repeated measures are recorded as the hospital reported them; where three
% consecutive readings exist the modal/mean value is taken, matching the
% convention already used for the pre-surgery rows.
%
% DIRECTION CHECK (expected physiology after closure): PA pressure falls
% (20/10 mean 15 -> 17/9 mean 13), RAP unchanged (5 -> 5). Both consistent
% with removal of the left-to-right shunt.
% ======================================================================

% ---- Pulmonary circulation (normalised post-surgery) -----------------
% log 12.26.48 / 12.27.00 / 12.27.11: PA 17/8 (13), 17/9 (13), 17/9 (13)
post.PAP_sys_mmHg     = 17;    % [mmHg] post-closure PA systolic
post.PAP_dia_mmHg     = 9;     % [mmHg] post-closure PA diastolic
post.PAP_mean_mmHg    = 13;    % [mmHg] post-closure PA mean
post.PVR_WU           = NaN;   % [WU]   not measured; no post-closure CO recorded

% ---- Systemic circulation --------------------------------------------
% log 12.31.35 / 12.32.10 / 12.32.39: RFA 91/68 (79), 89/68 (78), 89/68 (79)
% RFA is used for consistency with the pre-surgery systemic rows, which also
% come from the right femoral artery rather than the descending aorta.
post.SAP_sys_mmHg     = 89;    % [mmHg] post-closure RFA systolic
post.SAP_dia_mmHg     = 68;    % [mmHg] post-closure RFA diastolic
% Catheter-stamped mean, same policy as pre.SAP_mean_mmHg above. The
% form-factor reconstruction would give 68 + (89-68)/3 = 75, again ~4 mmHg
% below the transducer's own figure - the same systematic offset seen
% pre-closure (71.3 vs 77).
post.MAP_mmHg         = 79;    % [mmHg] mean arterial pressure (maps to SAP_mean target)
post.SVR_WU           = NaN;   % [WU]   not measured; no post-closure CO recorded

% ---- Atrial pressures ------------------------------------------------
% log 12.28.14 / 12.28.24 / 12.28.35: RA 8/5 (5), 8/5 (5), 7/5 (5)
post.RAP_mean_mmHg    = 5;     % [mmHg] post-closure RA mean
post.LAP_mean_mmHg    = NaN;   % [mmHg] not measured

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
