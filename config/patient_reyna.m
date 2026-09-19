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
clinical.common.patient_id   = '';      % [char] real MRN kept in config/private/ only, see below
clinical.common.maturation_mode = 'normal'; % 'normal' | 'pvr_fixed_day3' | 'pvr_fixed_day30' | 'none'

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
% HR: the IRB-approved protocol form (Filled Protokol_VSD Pre_Reyna (2).pdf,
% IRB/65/08/ETIK/2025, row 6, "average heart rate") records 119 bpm as the
% authoritative averaged rate for this session. A prior revision of this file
% used 136 bpm ("Nadi 136 bpm", procedure log 06/04/2026 09.24.57), reasoning
% it was a single spot pulse reading close to the pressure measurements. That
% reasoning is superseded: the protocol form is the IRB-governed record and
% explicitly labels its value as an average, not a spot reading, so it takes
% precedence. HR sets cycle length, so this is a model input, not just a
% reporting field: 60/119 = 0.504 s vs 60/136 = 0.441 s — reverting this
% value invalidates every previously published fit and is expected to change
% downstream RMSE; see docs/CHANGES_SINCE_PR22.md for the reconciliation.
% Heart rate is recorded separately on each protocol form and is therefore
% stored per scenario (pre.HR, post.HR below), never shared, so the two
% states cannot inherit each other's cycle length.

% PROVENANCE ONLY — not modelled as struct fields, no established target in
% this model for outflow tract diameters or non-invasive vitals:
%   LVOT diameter 1.2 cm, RVOT diameter 1.4 cm (protocol rows 33, 35).
%   Doppler VTIs for these tracts are blank in the form, so echo Qp/Qs is
%   not computable from them.
%   Header vitals at the start of the session: SpO2 88%, RR 28/min.

%% =====================================================================
%  PRE-SURGERY — haemodynamics in the presence of the open VSD
%  Source: right-heart catheterisation and echocardiography
%% =====================================================================
pre = struct();

pre.HR                = 119;     % [bpm] pre-closure protocol form row 6

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
% MEAN: the protocol's own definition (row 20, MAP = 1/3 sys + 2/3 dia)
% applied to the catheter RFA readings of rows 12-13, so every value is
% traceable to the IRB protocol form. Applied identically post-closure.
% Not used: the catheter-stamped mean of 77 mmHg (procedure log, not part of
% the protocol record) and the form's MAP of 95 mmHg (row 20 computed from
% the NIBP cuff, a different measurement method).
pre.SAP_mean_mmHg     = (100 + 2*57)/3;   % [mmHg] = 71.33, protocol row 20 formula on rows 12-13
%
% PROVENANCE ONLY — descending aorta (DAO), not used as a fitted target:
% protocol form rows 14-15 record DAO systolic/diastolic 91/57 mmHg, a
% second catheter site from the same session. RFA 100/57 vs DAO 91/57 is a
% ~9 mmHg systolic spread between two simultaneous sites in the same
% patient, and it doubles as an empirical sigma check against the
% objective's assumed pressure uncertainty of 10 mmHg (see AGENTS.md
% sigma-weighting section): the two sites disagree by less than that
% assumed measurement noise.

pre.SVR_WU            = NaN;     % [WU] protocol row 25 blank; not used as clinical target

% ---- Atrial and ventricular filling pressures -------------------------
pre.RAP_mean_mmHg     = 5;       % [mmHg] catheter, mean of 5/5/5 mmHg (row 19)
pre.LAP_mean_mmHg     = NaN;
pre.LVEDP_mmHg        = NaN;
pre.RVEDP_mmHg        = NaN;    % [mmHg] not captured in the protocol form

% ---- Ventricular volumes and ejection fraction -----------------------
% CORRECTED TWICE (publication-readiness reconciliation, 2026-09-05 then
% 2026-09-06). The first correction moved this block from "excluded
% entirely, believed to be H+1 post-operative echo" to "pre-surgery
% consistency-only", reading the protocol form's rows 26-29 header
% ("PARAMETER VOLUME UNTUK VALIDASI MODEL - PRE RELEASE OCCLUDER") as
% same-session pre-surgery evidence. That reading was itself wrong: per
% the study owner, "pre-release occluder" means the closure device is
% already deployed and occluding the defect, simply not yet mechanically
% detached from its delivery cable — see the post-surgery block below,
% where the device is placed at 11.50.19 and released (detached) at
% 12.06.23. These volumes were measured in that placed-but-undetached
% window, i.e. the VSD is already functionally CLOSED. They belong in
% clinical.post_surgery (below), not here. This also resolves the earlier
% "internally implausible LV pair" flag (SV_LV = 8.4 mL looked far too
% small against the PRE-op Qp-implied SV of ~34 mL): a small post-closure
% stroke volume is exactly what is physiologically expected once the
% left-to-right shunt's volume load is removed, so that comparison never
% applied in the first place.
pre.LVEDV_mL          = NaN;    % [mL] not a pre-surgery measurement; see clinical.post_surgery
pre.LVESV_mL          = NaN;    % [mL] not a pre-surgery measurement; see clinical.post_surgery
pre.RVEDV_mL          = NaN;    % [mL] not a pre-surgery measurement; see clinical.post_surgery
pre.RVESV_mL          = NaN;    % [mL] not a pre-surgery measurement; see clinical.post_surgery
pre.EF                = NaN;    % [-] not a pre-surgery measurement; see clinical.post_surgery

% ---- IC override flag -------------------------------------------------
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
%   - The chamber volume block belongs to post-surgery (see above), not
%     pre-surgery, so the pre-surgery objective remains hemodynamic-only.
pre.CO_Lmin           = 3.423;   % [L/min] Qs = Qp/QpQs = 4.087/1.194 (rows 21 and 23)

clinical.pre_surgery = pre;

%% =====================================================================
%  POST-SURGERY  — haemodynamics after VSD closure
%  Source: post-operative catheterisation / echocardiography
%% =====================================================================
post = struct();

post.HR               = 108;     % [bpm] post-closure protocol form row 6 (108, 108)

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
% Same definition as pre.SAP_mean_mmHg: protocol row 17 formula
% (MAP = 1/3 sys + 2/3 dia) applied to the catheter RFA readings of rows 9-10.
post.MAP_mmHg         = (89 + 2*68)/3;    % [mmHg] = 75.0 (maps to SAP_mean target)
post.SVR_WU           = NaN;   % [WU]   not measured; no post-closure CO recorded

% ---- Atrial pressures ------------------------------------------------
% log 12.28.14 / 12.28.24 / 12.28.35: RA 8/5 (5), 8/5 (5), 7/5 (5)
post.RAP_mean_mmHg    = 5;     % [mmHg] post-closure RA mean
post.LAP_mean_mmHg    = NaN;   % [mmHg] not measured
post.LVEDP_mmHg       = NaN;   % [mmHg] not captured in the procedure log
post.RVEDP_mmHg       = NaN;   % [mmHg] not captured in the procedure log

% ---- Ventricular volumes and function ---------------------------------
% Protocol form rows 26-29, section "PARAMETER VOLUME UNTUK VALIDASI MODEL
% - PRE RELEASE OCCLUDER": measured with the closure device already
% deployed and occluding the defect, placed at 11.50.19, before mechanical
% detachment ("release") at 12.06.23. See the note in the pre-surgery
% block above for why this timing places these values here, not there.
% Reported and predicted-against as consistency-only evidence (no
% dedicated post-surgery recipe/tier config exists yet for Reyna, so
% get_calibration_targets/build_target_tiers will apply the DEFAULT tier
% policy for any future post_surgery calibration run — that default
% treats LVEDV/LVESV/LVEF as hard and RVEDV/RVESV/RVEF as soft FITTED
% targets, not consistency-only. Anyone adding a post_surgery calibration
% recipe for Reyna should decide deliberately whether that default is
% appropriate, rather than inheriting it silently.
post.LVEDV_mL         = 32.0;    % [mL] protocol row 26, pre-release occluder (= post-closure)
post.LVESV_mL         = 23.6;    % [mL] protocol row 27, pre-release occluder (= post-closure)
post.RVEDV_mL         = 30.5;    % [mL] protocol row 28, pre-release occluder (= post-closure)
post.RVESV_mL         = 12.0;    % [mL] protocol row 29, pre-release occluder (= post-closure)
post.EF               = 0.2625;  % [-] LV EF = (32-23.6)/32
post.RVEF             = 0.6066;  % [-] RV EF = (30.5-12)/30.5

% ---- Cardiac output --------------------------------------------------
post.CO_Lmin          = NaN;

clinical.post_surgery = post;

clinical.protocol = reyna_protocol_record();

end

% =======================================================================
function rec = reyna_protocol_record()
% REYNA_PROTOCOL_RECORD - every entry of the two IRB protocol forms, verbatim.
% Kept outside clinical.pre_surgery/post_surgery on purpose: nothing here is
% read by calibration. It exists so every protocol value, used or not, can be
% reported against the model with its row number (export_run_supplement).
% Sources: "Filled Protokol_VSD Pre_Reyna (2).pdf" and
% "Half-Filled Protokol_VSD Post_Reyna_09-04-2026.pdf", protocol V.03,
% IRB/65/08/ETIK/2025. Repeated readings are listed in recorded order;
% blank cells are empty ([]).
f = @(row, name, readings, unit, method, note) struct('row', row, ...
    'name', name, 'readings', readings, 'unit', unit, 'method', method, 'note', note);

rec.pre = [ ...
    f(1,  'Age',                         3 + 2/12,  'yr',    'medical record', '3 years 2 months')
    f(3,  'Height',                      95,        'cm',    'direct measurement', '')
    f(4,  'Weight',                      13.4,      'kg',    'medical record', '')
    f(5,  'Body surface area',           0.588,     'm2',    'calculated from height and weight', '')
    f(6,  'Heart rate',                  119,       'bpm',   'ECG monitor', '')
    f(7,  'VSD diameter, LV side',       [4.46 6.43], 'mm',  'echocardiography Doppler', 'range')
    f(7,  'VSD diameter, RV side',       [3.63 3.70], 'mm',  'echocardiography Doppler', 'range')
    f(9,  'LV systolic pressure',        94,        'mmHg',  'catheterization', 'recorded as gradient row')
    f(9,  'RV systolic pressure',        25,        'mmHg',  'catheterization', 'recorded as gradient row')
    f(10, 'Systolic pressure, NIBP',     119,       'mmHg',  'NIBP cuff', '')
    f(11, 'Diastolic pressure, NIBP',    83,        'mmHg',  'NIBP cuff', '')
    f(12, 'Systolic pressure, RFA',      100,       'mmHg',  'catheterization', '')
    f(13, 'Diastolic pressure, RFA',     57,        'mmHg',  'catheterization', '')
    f(14, 'Systolic pressure, DAo',      91,        'mmHg',  'catheterization', '')
    f(15, 'Diastolic pressure, DAo',     57,        'mmHg',  'catheterization', '')
    f(16, 'PA systolic pressure',        [20 20 20], 'mmHg', 'catheterization', '')
    f(17, 'PA diastolic pressure',       [10 10 10], 'mmHg', 'catheterization', '')
    f(18, 'PA mean pressure',            [14 15 15], 'mmHg', 'catheterization', '')
    f(19, 'RA mean pressure',            [5 5 5],   'mmHg',  'catheterization', '')
    f(20, 'Mean arterial pressure, form', 95,       'mmHg',  'formula on NIBP', 'form value computed from rows 10-11')
    f(21, 'Pulmonary flow Qp',           4.087,     'L/min', 'Fick (catheterization)', 'an earlier entry is crossed out')
    f(22, 'Systemic flow Qs',            [],        'L/min', 'Fick (catheterization)', 'blank')
    f(23, 'Qp/Qs',                       1.194,     '-',     'Fick (catheterization)', '')
    f(24, 'PVR',                         [],        'WU',    'calculation', 'blank')
    f(25, 'SVR',                         [],        'WU',    'calculation', 'blank')
    f(26, 'LVEDV',                       32,        'mL',    'echocardiography biplane Simpson', 'post-closure state per study decision')
    f(27, 'LVESV',                       23.6,      'mL',    'echocardiography biplane Simpson', 'post-closure state per study decision')
    f(28, 'RVEDV',                       30.5,      'mL',    'echocardiography 3D', 'post-closure state per study decision')
    f(29, 'RVESV',                       12,        'mL',    'echocardiography 3D', 'post-closure state per study decision')
    f(30, 'Stroke volume',               [],        'mL',    'echocardiography', 'blank')
    f(31, 'Cardiac output, echo',        [],        'L/min', 'echocardiography', 'blank')
    f(32, 'Ejection fraction',           [],        '%',     'echocardiography', 'blank')
    f(33, 'LVOT diameter',               1.2,       'cm',    'echocardiography', '')
    f(34, 'LVOT VTI',                    [],        'cm',    'echocardiography', 'blank')
    f(35, 'RVOT diameter',               1.4,       'cm',    'echocardiography', '')
    f(36, 'RVOT VTI',                    [],        'cm',    'echocardiography', 'blank')];

rec.post = [ ...
    f(6,  'Heart rate',                  [108 108], 'bpm',   'ECG monitor', '')
    f(7,  'Systolic pressure, NIBP',     [87 85 88], 'mmHg', 'NIBP cuff', '')
    f(8,  'Diastolic pressure, NIBP',    [56 55 55], 'mmHg', 'NIBP cuff', '')
    f(9,  'Systolic pressure, RFA',      [91 89 89], 'mmHg', 'catheterization', '')
    f(10, 'Diastolic pressure, RFA',     [68 68 68], 'mmHg', 'catheterization', '')
    f(11, 'Systolic pressure, DAo',      [104 106 103], 'mmHg', 'catheterization', 'handwritten values overwritten on form')
    f(12, 'Diastolic pressure, DAo',     [75 72 68], 'mmHg', 'catheterization', 'handwritten values overwritten on form')
    f(13, 'PA systolic pressure',        [17 17 17], 'mmHg', 'catheterization', '')
    f(14, 'PA diastolic pressure',       [8 9 9],   'mmHg',  'catheterization', '')
    f(15, 'PA mean pressure',            [13 13 13], 'mmHg', 'catheterization', '')
    f(16, 'RA mean pressure',            [5 5 5],   'mmHg',  'catheterization', '')
    f(17, 'Mean arterial pressure, form', [66 65 66], 'mmHg', 'formula on NIBP', 'each cell holds two numbers (67/66, 64/65, 66/66); lower values listed')
    f(18, 'Pulmonary flow Qp',           [],        'L/min', 'echocardiography estimate', 'blank')
    f(19, 'Systemic flow Qs',            [],        'L/min', 'echocardiography estimate', 'blank')
    f(20, 'Qp/Qs',                       [],        '-',     'echocardiography estimate', 'blank')
    f(21, 'PVR',                         [],        'WU',    'calculation', 'blank')
    f(22, 'SVR',                         [],        'WU',    'calculation', 'blank')
    f(23, 'LVEDV',                       [],        'mL',    'echocardiography', 'blank; evaluation uses pre-form row 26')
    f(24, 'LVESV',                       [],        'mL',    'echocardiography', 'blank; evaluation uses pre-form row 27')
    f(25, 'RVEDV',                       [],        'mL',    'echocardiography', 'blank; evaluation uses pre-form row 28')
    f(26, 'RVESV',                       [],        'mL',    'echocardiography', 'blank; evaluation uses pre-form row 29')
    f(27, 'Stroke volume',               [],        'mL',    'echocardiography', 'blank')
    f(28, 'Cardiac output',              [],        'L/min', 'echocardiography', 'blank')
    f(29, 'Ejection fraction',           [],        '%',     'echocardiography', 'blank')];
end
