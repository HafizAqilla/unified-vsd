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

% Patient roster — add or remove entries here to control which patients run.
patient_fns = {
    % @patient_reyna
    @patient_profile_Razka
};
patient_labels = {
    
    % 'reyna'
    'razka'
};

%% =====================================================================
%  RUN PIPELINE — iterates over all patients in roster above
%% =====================================================================
for p = 1:numel(patient_fns)
    label   = patient_labels{p};
    clinical = patient_fns{p}();

    fprintf('\n===== run_patient_case: patient=%s =====\n', label);
    fprintf('  Patient: %.1f kg | age %.1f mo | %.0f cm | sex=%d\n', ...
        clinical.common.weight_kg, clinical.common.age_years * 12, ...
        clinical.common.height_cm, clinical.common.sex);

    results_dir = fullfile(root, 'results');
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    diary(fullfile(results_dir, sprintf('console_%s_pre_surgery.txt', label)));
    main_run('pre_surgery', clinical);
    diary off

    % Rename calibrated params so next patient does not overwrite.
    src = fullfile(root, 'results', 'tables', 'params_calibrated_pre_surgery.mat');
    dst = fullfile(root, 'results', 'tables', sprintf('params_calibrated_pre_surgery_%s.mat', label));
    if isfile(src)
        movefile(src, dst);
        fprintf('[run_patient_case] Saved calibrated params to: %s\n', dst);
    end

    fprintf('\n===== run_patient_case: done patient=%s =====\n', label);
end

fprintf('\n===== run_patient_case: all patients complete =====\n');
