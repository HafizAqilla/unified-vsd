function clinical = patient_profile_Razka()
% PATIENT_PROFILE_RAZKA
% -----------------------------------------------------------------------
% Patient:        RAZKA ALFA RIZKI, AN   |  MRN: 00948048
% Procedure date: 21/04/2026             |  DOB: 25/09/2021
% Facility:       RSAB Harapan Kita
% Source:         Procedure Log, 21/04/2026 (scanned, 2 pages)
%
% PROCEDURE TYPE: Kateterisasi Diagnostik (Diagnostic Cardiac Catheterization)
%   - Right and left heart catheterization via right femoral approach
%   - Catheter successfully crossed VSD ("Crossing VSD" confirmed in log)
%   - Intraoperative TEE performed — overriding aorta identified
%   - Device closure NOT performed; operator decision after TEE finding
%   - Cardiac output (Fick / thermodilution) explicitly NOT measured
%     ("no CO" noted on multiple entries)
%
% DEMOGRAPHICS:
%   DOB:    25/09/2021
%   Age at procedure (21/04/2026): 4 yr 238 days ≈ 4.64 years
%   Weight: 14.3 kg    Height: 100.4 cm    BSA: 0.629 m²
%   Sex:    NOT DOCUMENTED in procedure log — placeholder 'M'; verify
%
% KEY HAEMODYNAMICS (catheter, air rest):
%   HR:  111 bpm  (monitor 12:16:10)
%   LV:   91 / -1 / 9  mmHg  (systolic / min-diastolic / LVEDP)
%   RV:   26 /  1 / 8  mmHg  (systolic / min-diastolic / RVEDP)
%   PAp:  22 / 10 / 15 mmHg  (systolic / diastolic / mean)
%   DAO:  82 / 59 / 70 mmHg  (systolic / diastolic / mean)
%   RA:   mean ≈ 5 mmHg
%   FR  (Qp/Qs) = 1.21
%   PARI          = 0.71  — see note below
%
% NOTE — PARI:
%   Recorded verbatim from procedure log as "Hasil PARI: 0.71".
%   Exact definition (PVRi [WU·m²], PVR/SVR ratio, or other index)
%   is not stated. Stored as pre.PARI_raw only; do NOT use in simulation
%   until the definition is confirmed from the catheterisation report.
%
% NOTE — MISSING DATA:
%   CO/Fick not measured → Q_shunt_Lmin, PVR_WU, SVR_WU all NaN.
%   Echo volumes not in this log → LVEDV/LVESV/RVEDV/RVESV all NaN.
%   VSD diameter not reported → VSD_diameter_mm = NaN.
%   LAP/PCWP not measured → LAP_mean_mmHg = NaN.
%
% NOTE — AORTIC PRESSURE:
%   Monitor NIBP: 92/59 (70) mmHg matches RFA (femoral) cath 91/56 (70).
%   DAO (descending aorta) cath: 82/59 (70) — lower systolic consistent
%   with catheter position distal to aortic arch.  Simulation uses DAO.
%   LV systolic = 91 mmHg → no significant aortic valve gradient detected.
%
% USAGE:
%   clinical = patient_profile_Razka();
%   main_run('pre_surgery', clinical)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-25
% VERSION:  1.0
% -----------------------------------------------------------------------

clinical = patient_template();   % initialise all fields to NaN
clinical.common.patient_name = 'razka';    % [char] patient label for run folders
clinical.common.patient_id   = '00948048'; % [char] MRN from procedure log

%% =====================================================================
%  COMMON — demographics
%% =====================================================================
% Age: DOB 25/09/2021 → procedure 21/04/2026 = 4 yr 238 days ≈ 4.64 yr
clinical.common.age_years  = 4.64;   % [years]  derived from DOB 25/09/2021
clinical.common.weight_kg  = 14.3;   % [kg]     measured; procedure log 21/04/2026
clinical.common.height_cm  = 100.4;  % [cm]     measured; procedure log 21/04/2026
clinical.common.BSA        = 0.629;  % [m²]     stated in procedure log; overrides Mosteller
clinical.common.sex        = 1;      % NOT documented — placeholder only; verify before use
clinical.common.HR         = 111;    % [bpm]    monitor reading 21/04/2026 12:16:10

%% =====================================================================
%  PRE-SURGERY — catheterisation data, 21/04/2026, air rest
%  All pressure entries labelled "AIR REST" in procedure log.
%  Pressures reported as: systolic / diastolic (mean) or
%                          systolic / min-diastolic / EDP
%% =====================================================================
pre = clinical.pre_surgery;

% ---- Shunt severity --------------------------------------------------
% VSD catheter crossing confirmed; TEE showed overriding aorta
pre.VSD_diameter_mm   = 2.5;    % [mm]     range 2–3 mm reported; midpoint used
                                %           Source: verbal communication, 2026-04-27
pre.VSD_gradient_mmHg = NaN;    % [mmHg]   not reported (LV-RV systolic difference
                                %           exists in cath data but not stated as gradient)
pre.Q_shunt_Lmin      = NaN;    % [L/min]  CO not measured; cannot be computed
pre.QpQs              = 1.21;   % [-]      "Hasil FR: 1.21" — procedure log

% ---- Pulmonary circulation -------------------------------------------
% PAp = pulmonary artery pressure (definitive reading, final pullback)
% Supporting RPA readings: 22/8(15), 20/9(15), 20/9(14), 20/9(15)
pre.PAP_sys_mmHg      = 22;     % [mmHg]   PAp 22/10 (15) — procedure log
pre.PAP_dia_mmHg      = 10;     % [mmHg]   PAp 22/10 (15) — procedure log
pre.PAP_mean_mmHg     = 15;     % [mmHg]   PAp 22/10 (15); consistent across all RPA entries
pre.PVR_WU            = NaN;    % [WU]     CO not measured; cannot compute

% ---- Systemic circulation (DAO = Descending Aorta, cath) -------------
pre.SAP_sys_mmHg      = 82;     % [mmHg]   DAO 82/59 (70) — procedure log
pre.SAP_dia_mmHg      = 59;     % [mmHg]   DAO 82/59 (70) — procedure log
pre.SAP_mean_mmHg     = 70;     % [mmHg]   DAO 82/59 (70) — procedure log
pre.SVR_WU            = NaN;    % [WU]     CO not measured; cannot compute

% ---- Atrial pressures ------------------------------------------------
% RA readings: 6/7(5), 7/7(6), 7/7(5) — mean consistently 5–6 mmHg
pre.RAP_mean_mmHg     = 5;      % [mmHg]   RA a/v(mean): mean = 5 mmHg
pre.LAP_mean_mmHg     = NaN;    % [mmHg]   PCWP not measured in this procedure

% ---- Ventricular volumes (echo) — not reported in this document ------
pre.LVEDV_mL          = NaN;    % [mL]     not reported
pre.LVESV_mL          = NaN;    % [mL]     not reported
pre.RVEDV_mL          = NaN;    % [mL]     not reported
pre.RVESV_mL          = NaN;    % [mL]     not reported
pre.LVEF              = NaN;    % [-]      not reported

% ---- Cardiac output --------------------------------------------------
pre.CO_Lmin           = NaN;    % [L/min]  explicitly not measured ("no CO" — procedure log)

% ---- Ventricular end-diastolic pressures (direct cath) ---------------
% LV readings: 91/-1/9, 90/-1/9, 90/0/10, 91/-1/10, 90/-1/10
% Format: systolic / min-diastolic / end-diastolic pressure (EDP)
pre.LVEDP_mmHg        = 9;      % [mmHg]   LV cath; median of EDP readings (9,9,10,10,10)

% RV readings: 26/1/8, 26/1/8, 26/1/7
pre.RVEDP_mmHg        = 8;      % [mmHg]   RV cath; majority reading = 8 (two of three entries)

% ---- LV systolic pressure (direct cath) ------------------------------
pre.LVP_sys_mmHg      = 91;     % [mmHg]   LV cath; consistent across all LV entries
pre.RVP_sys_mmHg      = 26;     % [mmHg]   RV cath; consistent across all RV entries

% ---- Raw PARI (verbatim from procedure log — definition unconfirmed) -
pre.PARI_raw          = 0.71;   % [?]      "Hasil PARI: 0.71" — unit/definition not stated
                                %           Do NOT use in simulation until definition confirmed

clinical.pre_surgery  = pre;

%% =====================================================================
%  POST-SURGERY — not applicable
%  This was a diagnostic catheterisation only; VSD closure was not
%  performed.  All post_surgery fields remain NaN from patient_template().
%% =====================================================================

end
