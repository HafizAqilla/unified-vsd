% RUN_PATIENT_CASE
% -----------------------------------------------------------------------
% Operational entry point for real patient datasets.
%
% USAGE:
%   Run from the MATLAB command window:
%        >> run run_patient_case
%   Or:
%        >> delete(gcp('nocreate')); clear functions; run run_patient_case
%
% OUTPUTS:
%   - Console validation table: simulated vs. clinical
%   - results/tables/calib_diagnostics_pre_surgery_<timestamp>.mat
%   - results/tables/params_calibrated_pre_surgery.mat
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-13
% VERSION:  1.2
% -----------------------------------------------------------------------

%% ---- housekeeping ------------------------------------------------------
% Kill stale parallel pool so workers reload fresh .m files after any edits.
delete(gcp('nocreate'));
clear functions                  % flush main-process function cache
clearvars -except;

% Build clean path: project root only, strip all shadow/worktree dirs.
% This prevents MATLAB from resolving .m files from stale
% .claude/worktrees copies that may be earlier on the path.
root = fileparts(mfilename('fullpath'));
restoredefaultpath();            % wipe entire path back to MATLAB builtins
root_paths = strsplit(genpath(root), pathsep);
is_shadow = contains(root_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.git'   filesep], 'IgnoreCase', true);
addpath(strjoin(root_paths(~is_shadow), pathsep));

% Note: parallel pool is created on-demand by fmincon if needed.
% addAttachedFiles() in run_calibration.m ensures workers always load the
% correct objective_calibration.m regardless of path resolution order.
% Calibration is currently set to serial (UseParallel=false) to avoid
% worker path/cache mismatch issues.

%% =====================================================================
%  PATIENT DEMOGRAPHICS  — Load from patient profile
%% =====================================================================

% Enable GSA so the pipeline runs both pre- and post-calibration Sobol.
setenv('UNIFIED_VSD_DO_GSA', '1');

% Choose which patient to run by calling their profile function here:
clinical = patient_reyna();

%% =====================================================================
%  RUN PIPELINE
%% =====================================================================
fprintf('\n===== run_patient_case: starting pre-surgery pipeline =====\n');
fprintf('  Patient: %.1f kg | age %.1f mo | %.0f cm | sex=%d\n', ...
    clinical.common.weight_kg, clinical.common.age_years * 12, ...
    clinical.common.height_cm, clinical.common.sex);

main_run('pre_surgery', clinical);

fprintf('\n===== run_patient_case: done =====\n');
