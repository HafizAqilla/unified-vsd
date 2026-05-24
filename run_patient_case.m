% RUN_PATIENT_CASE
% -----------------------------------------------------------------------
% Patient entry point — Reyna
% Simulation and calibration setup for Reyna VSD patient.
%
% Follows the same simple style as MATLAB/VSD/run_patient_case.m:
%   1. Clear workspace
%   2. Fill in patient clinical data inline
%   3. Call main_run
%
% USAGE:
%   Run from the MATLAB command window:
%        >> run run_patient_case
%
% OUTPUTS:
%   All outputs are saved automatically by main_run under results/.
%   Console log is written to results/console_reyna_pre_surgery.txt.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-22
% VERSION:  2.1  (inline patient data — matches MATLAB/VSD/run_patient_case.m style)
% -----------------------------------------------------------------------

%% ---- Housekeeping ------------------------------------------------------
clear; clc;

% Build a clean MATLAB path: project root only, strip shadow/worktree dirs.
root = fileparts(mfilename('fullpath'));
restoredefaultpath();
root_paths = strsplit(genpath(root), pathsep);
is_shadow = contains(root_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.git'   filesep], 'IgnoreCase', true);
addpath(strjoin(root_paths(~is_shadow), pathsep));

% Ensure results directory exists
results_dir = fullfile(root, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

% =========================================================================
%  PATIENT — REYNA
% =========================================================================
clinical = patient_template();   % initialise all fields to NaN

% ---- Demographics -------------------------------------------------------
clinical.common.patient_name = 'reyna';
clinical.common.age_years    = 3 + (2 / 12); % [years]  3 years 2 months
clinical.common.weight_kg    = 14.0;          % [kg]     Keisya revision 2026-05-11
clinical.common.height_cm    = 98.0;          % [cm]     Keisya revision 2026-05-11
clinical.common.sex          = 0;             % 0 = female  — AGENTS.md §3.10
clinical.common.BSA          = 0.6173;        % [m²]  Mosteller: sqrt(14.0*98.0/3600)
clinical.common.HR           = 119;           % [bpm]

% ---- VSD defect ---------------------------------------------------------
clinical.pre_surgery.VSD_diameter_mm   = 3.025;  % [mm]    mean RV-side (2.35+3.7)/2
clinical.pre_surgery.VSD_gradient_mmHg = 69;     % [mmHg]  peak: LV 94 - RV 25
clinical.pre_surgery.QpQs              = 1.194;  % [-]
clinical.pre_surgery.Q_shunt_Lmin      = 0.664;  % [L/min] Qp - Qs = 4.087 - 3.423

% ---- Pulmonary pressures (catheter) -------------------------------------
clinical.pre_surgery.PAP_sys_mmHg   = 20;    % [mmHg]
clinical.pre_surgery.PAP_dia_mmHg   = 10;    % [mmHg]
clinical.pre_surgery.PAP_mean_mmHg  = 15;    % [mmHg]
clinical.pre_surgery.PVR_WU         = NaN;   % [WU]    not calculated in protocol

% ---- Systemic pressures (RFA catheter) ----------------------------------
clinical.pre_surgery.SAP_sys_mmHg   = 100;   % [mmHg]  RFA catheter systolic
clinical.pre_surgery.SAP_dia_mmHg   = 57;    % [mmHg]  RFA catheter diastolic
clinical.pre_surgery.SAP_mean_mmHg  = 71.3;  % [mmHg]  57 + (100-57)/3
clinical.pre_surgery.SVR_WU         = 19.37; % [WU]    (71.3-5) / 3.423

% ---- Atrial pressures ---------------------------------------------------
clinical.pre_surgery.RAP_mean_mmHg  = 5;     % [mmHg]
clinical.pre_surgery.LAP_mean_mmHg  = 8;     % [mmHg]  estimated (no PCWP)

% ---- Ventricular pressures & filling ------------------------------------
clinical.pre_surgery.LVEDP_mmHg     = 8;     % [mmHg]  estimated; no direct measurement
clinical.pre_surgery.LVP_sys_mmHg   = NaN;   % [mmHg]
clinical.pre_surgery.RVP_sys_mmHg   = NaN;   % [mmHg]

% ---- Echo volumes (Teichholz from M-mode, LVEDD=32 mm, LVESD=23.6 mm) --
clinical.pre_surgery.LVEDV_mL   = 41.0;   % [mL]   (7/5.6) * 3.2^3
clinical.pre_surgery.LVESV_mL   = 19.3;   % [mL]   (7/4.76) * 2.36^3
clinical.pre_surgery.RVEDV_mL   = 30.5;   % [mL]   protocol row 28
clinical.pre_surgery.RVESV_mL   = 12.0;   % [mL]   protocol row 29
clinical.pre_surgery.LVEF       = 0.528;  % [-]    (41.0-19.3)/41.0

% ---- Cardiac output (Qs, Fick catheter) ---------------------------------
clinical.pre_surgery.CO_Lmin    = 3.423;  % [L/min]  systemic Fick (rows 21-23)

% ---- Run ----------------------------------------------------------------
fprintf('Starting simulation for patient: %s (%.1f kg, %.2f yr)\n', ...
    clinical.common.patient_name, clinical.common.weight_kg, clinical.common.age_years);

diary(fullfile(results_dir, sprintf('console_%s_pre_surgery.txt', clinical.common.patient_name)));
main_run('pre_surgery', clinical);
diary off

fprintf('\nDone. Results saved under: %s\n', results_dir);
