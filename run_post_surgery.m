% RUN_POST_SURGERY
% -----------------------------------------------------------------------
% Post-surgery simulation for the patient defined in run_patient_case.m.
%
% PURPOSE:
%   Continues from the calibrated pre-surgery model.  The calibration step
%   in main_run adjusted the free parameters (R.SAR, R.PAR, E.LV.EA, etc.)
%   to match pre-operative haemodynamics.  Those adjusted parameters carry
%   over to the post-surgery model as the starting point, representing the
%   patient's intrinsic cardiac physiology independent of the shunt.
%
% WORKFLOW:
%   1. Load calibrated pre-surgery params_cal from the most recent .mat file
%      saved by main_run (results/tables/params_calibrated_pre_surgery_*.mat).
%   2. Build the post-surgery clinical struct with expected post-op targets.
%   3. Call main_run('post_surgery', clinical) — this internally:
%        a. Applies allometric scaling (unchanged patient demographics)
%        b. Calls params_from_clinical (sets R.vsd = 1e6, maps post-op SVR/PVR)
%        c. Runs baseline simulation
%        d. Runs fmincon calibration on post-surgery free parameters
%           (see calibration_param_sets.m → case 'post_surgery')
%        e. Generates plots and validation report
%
% PATIENT:
%   3-year-old female | 13.4 kg | 95 cm
%   Same demographics as run_patient_case.m, consistent with the calibrated
%   pre-surgery result.
%
% POST-SURGERY HAEMODYNAMIC TARGETS:
%   After successful VSD closure in a child with pre-op Qp/Qs ≈ 1.19 and
%   MAP ≈ 95 mmHg, the expected early post-operative state is:
%     - Qp/Qs ≈ 1.0   (shunt abolished; residual leak < 0.1 L/min expected)
%     - SVR:  normalises toward physiological range for age/weight
%     - PAP:  expected to normalise rapidly in mildly elevated pre-op cases
%     - LVEF: maintained (pre-op LVEF not impaired)
%
%   All targets set to NaN where no post-operative measurement is available.
%   main_run will simulate without calibration targets for those rows
%   (validation table will mark them as NaN, excluded from RMSE).
%
% KEY DESIGN DECISION — param warm-start from calibrated pre-surgery file:
%   The calibrated pre-surgery params_cal already encodes the patient's
%   personalised R/C/E values.  Using these as the initial point for the
%   post-surgery optimisation (rather than allometric-scaled params0) places
%   the post-surgery fmincon at a physiologically informed starting location.
%
%   Implementation: the calibrated params_cal are INJECTED into main_run by
%   adding a 'CalibStart' entry to the clinical struct.  main_run detects
%   this and uses it in place of the allometric-scaled params0.
%
%   ALTERNATIVE — if this injection mechanism is not yet implemented, the
%   simplest route is to run main_run normally; main_run will re-scale from
%   demographics and re-calibrate from scratch.  Results will still be
%   physiologically consistent because the patient demographics are identical.
%
% HOW TO RUN:
%   1. Confirm that results/tables/ contains a
%      params_calibrated_pre_surgery_*.mat file from a completed pre-surgery run.
%   2. Optionally update the POST-SURGERY CLINICAL TARGETS block below with
%      actual post-operative measurements if available.
%   3. Run this script (F5 or "Run") in MATLAB.
%
% REFERENCES:
%   [1] run_patient_case.m          — original pre-surgery patient definition
%   [2] calibration_param_sets.m   — post-surgery free-parameter list
%   [3] params_from_clinical.m     — VSD closure: R.vsd = 1e6 mmHg·s/mL
%   [4] docs/theory_notes.md       — governing equations
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-14
% VERSION:  1.0
% -----------------------------------------------------------------------

clear; clc;
fprintf('=====================================================\n');
fprintf('  POST-SURGERY SIMULATION — run_post_surgery.m\n');
fprintf('=====================================================\n\n');

%% =====================================================================
%  SECTION 0 — Locate most-recent calibrated pre-surgery .mat file
%  The pre-surgery calibration run saves:
%    results/tables/params_calibrated_pre_surgery_<YYYYMMDD_HHMMSS>.mat
%  We sort by name (which embeds timestamp) and take the last (newest).
%% =====================================================================
root        = fileparts(mfilename('fullpath'));
addpath(genpath(root));

tables_dir  = fullfile(root, 'results', 'tables');
latest_mat  = fullfile(tables_dir, 'params_calibrated_pre_surgery_20260426_170809.mat');

if ~isfile(latest_mat)
    error(['run_post_surgery:noCalibFile\n' ...
           'The specified calibration file was not found:\n  %s\n\n'], latest_mat);
end

fprintf('[run_post_surgery] Loading calibrated parameters from:\n');
fprintf('  %s\n\n', latest_mat);

loaded       = load(latest_mat);   % fields: params_cal, calib_out, report
params_pre   = loaded.params_cal;  % calibrated pre-surgery parameter struct

%% =====================================================================
%  SECTION 1 — Build post-surgery clinical struct
%
%  Patient demographics are IDENTICAL to run_patient_case.m.
%  The post_surgery sub-struct captures expected / measured post-op targets.
%  Set any field to NaN where no measurement is available — it will be
%  excluded from the RMSE but will still appear in the validation table.
%% =====================================================================
clinical = patient_template();

% ---- Common demographics (same patient, different time point) ----------
clinical.common.age_years  = 3 + (2 / 12);   % [years]  3 years 2 months
clinical.common.weight_kg  = 13.4;            % [kg]
clinical.common.height_cm  = 95;              % [cm]
clinical.common.sex        = 'F';             % Female
% Heart rate: post-op HR commonly decreases slightly once shunt load removed.
% Set to NaN to use the allometric-scaled value, or supply if measured.
clinical.common.HR         = NaN;             % [bpm]  NaN → allometric default

% ---- Pre-surgery data (required by params_from_clinical for allometric
%      reference; values are those used in the original calibration run) ----
clinical.pre_surgery.VSD_diameter_mm   = 6.0;    % [mm]
clinical.pre_surgery.VSD_gradient_mmHg = 94;     % [mmHg]
clinical.pre_surgery.PAP_sys_mmHg      = 20;     % [mmHg]
clinical.pre_surgery.PAP_dia_mmHg      = 10;     % [mmHg]
clinical.pre_surgery.PAP_mean_mmHg     = 15;     % [mmHg]
clinical.pre_surgery.SAP_sys_mmHg      = 119;    % [mmHg]
clinical.pre_surgery.SAP_dia_mmHg      = 83;     % [mmHg]
clinical.pre_surgery.SAP_mean_mmHg     = 95;     % [mmHg]
clinical.pre_surgery.RAP_mean_mmHg     = 5;      % [mmHg]
clinical.pre_surgery.QpQs              = 1.194;  % [-]

% =====================================================================
%  POST-SURGERY CLINICAL TARGETS
%  *** EDIT THIS BLOCK with actual post-operative measurements. ***
%
%  Physiological expectations for this patient (pre-op Qp/Qs ≈ 1.19):
%  - MAP:   expected to remain near pre-op (~90–100 mmHg) immediately
%           post-closure.  Mild drop possible due to reduced preload.
%  - Qp/Qs: should return to ≈ 1.0 (shunt abolished).
%  - PAP:   low pre-op PAP → expected to normalise rapidly.
%  - SVR:   may increase slightly as LV afterload redistributes.
%  - LVEF:  expected to be preserved (pre-op function normal).
%
%  Source: expected values derived from:
%    Saxena A et al. (2007). J Thorac Cardiovasc Surg 134:955–960.
%    Hoffman JIE (2016). Congenital Heart Disease in Adults. Ch. 6.
% =====================================================================
post = clinical.post_surgery;

% ---- Shunt (abolished post-closure) ----------------------------------
post.QpQs              = 0.9;     % [-]      target: no residual shunt
                                  %          set NaN if not measured

% ---- Pulmonary circulation (expected normalisation) ------------------
% Low pre-op PAP (mean 15 mmHg) → expected to remain or decrease slightly.
post.PAP_sys_mmHg      = 17;    % [mmHg]  set when (echo/RHC) available
post.PAP_dia_mmHg      = 8.67;    % [mmHg]
post.PAP_mean_mmHg     = 13;    % [mmHg]
post.PVR_WU            = NaN;    % [WU]    expected ↓; set if measured

% ---- Systemic circulation --------------------------------------------
% With shunt closed, systemic CO = pulmonary CO → SAP similar or slightly ↑.
post.SAP_sys_mmHg      = 104.3;    % [mmHg]  set when available
post.SAP_dia_mmHg      = 71.67;    % [mmHg]
post.MAP_mmHg          = 66;    % [mmHg]  mean arterial pressure post-op
post.SVR_WU            = NaN;    % [WU]    set if SVR measured

% ---- Atrial pressures ------------------------------------------------
post.RAP_mean_mmHg     = 5;    % [mmHg]
post.RAP_max_mmHg      = 7.67;    % [mmHg]
post.RAP_min_mmHg      = 5;    % [mmHg]

post.LAP_mean_mmHg     = NaN;    % [mmHg]

% ---- Ventricular volumes (expected remodelling) ----------------------
% With volume overload removed, LVEDV/RVEDV typically decrease over weeks.
% Immediate post-op echo volumes may still be elevated.
post.LVEDV_mL          = 32;    % [mL]  set from post-op echo if available
post.LVESV_mL          = 23.6;    % [mL]
post.RVEDV_mL          = 30.5;    % [mL]
post.RVESV_mL          = 12;    % [mL]
post.LVEF              = 0.618;  % [-]   fraction (61.8%)
post.RVEF              = NaN;    % [-]   fraction

post.LVSV_mL           = 26.4;   % [mL]  left ventricular stroke volume

% ---- Cardiac output --------------------------------------------------
post.CO_Lmin           = 3.0;    % [L/min]  set from Fick/TD if measured

clinical.post_surgery = post;

%% =====================================================================
%  SECTION 2 — Warm-start: inject calibrated pre-surgery params
%
%  The fields that carry over from the pre-surgery calibration:
%    Ventricular elastances:  E.LV.EA, E.LV.EB, E.RV.EA, E.RV.EB
%    Unstressed volumes:      V0.LV, V0.RV
%    Vascular resistances:    R.SAR, R.SC, R.SVEN, R.PAR, R.PCOX, R.PVEN
%    Compliances:             C.SAR, C.PAR, C.SVEN, C.PVEN
%
%  These represent the patient's personalised cardiovascular physiology
%  which does not change instantaneously with surgical closure.
%
%  The clinical.pre_surgery.CalibParams field is used to signal to
%  main_run that a warm-start struct is available.
%% =====================================================================
clinical.pre_surgery.CalibParams = params_pre;  % warm-start payload

fprintf('[run_post_surgery] Warm-start parameters loaded from pre-surgery calibration.\n');
fprintf('  Pre-surgery calibration J* = %.6f\n', loaded.calib_out.fbest);
fprintf('  Key calibrated parameters:\n');
fprintf('    E.LV.EA = %.4f mmHg/mL  |  E.LV.EB = %.4f mmHg/mL\n', ...
        params_pre.E.LV.EA, params_pre.E.LV.EB);
fprintf('    E.RV.EA = %.4f mmHg/mL  |  E.RV.EB = %.4f mmHg/mL\n', ...
        params_pre.E.RV.EA, params_pre.E.RV.EB);
fprintf('    R.SAR   = %.4f mmHg·s/mL | R.PAR  = %.4f mmHg·s/mL\n', ...
        params_pre.R.SAR, params_pre.R.PAR);
fprintf('    R.vsd(pre) = %.2e (will be set to 1e6 for post_surgery)\n\n', ...
        params_pre.R.vsd);

%% =====================================================================
%  SECTION 3 — Run the full post-surgery pipeline via main_run
%
%  main_run('post_surgery', clinical) internally executes:
%    Step 1: apply_scaling        — allometric scaling for patient
%    Step 2: params_from_clinical — R.vsd = 1e6; map post-op SVR/PVR targets
%    Step 3: baseline simulation  — integrate_system with closed VSD
%    Step 4: calibration          — fmincon on post-surgery param set
%    Step 5: validation report    — comparison vs clinical.post_surgery targets
%    Step 6: plots                — haemodynamic waveforms
%
%  To skip calibration and run a quick forward simulation only,
%  set DO_CALIBRATION = false in main_run.m before running.
%% =====================================================================
fprintf('[run_post_surgery] Launching main_run(''post_surgery'', clinical)...\n\n');
main_run('post_surgery', clinical);

fprintf('\n[run_post_surgery] Post-surgery simulation complete.\n');
fprintf('Results saved to: %s\n', fullfile(root, 'results', 'tables'));
