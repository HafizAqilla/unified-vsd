%% compare_scaling_modes.m
% =========================================================================
% Compare Zhang weight allometry against Lundquist BSA allometry.
%
% PURPOSE:
%   Runs the acceptance trio with each scaling law under identical runtime
%   switches. This is intentionally portfolio-based: the better scaling
%   choice should preserve clinical fit and fitted-parameter plausibility
%   across full-data, sparse-cath, and synthetic benchmark cases.
%
% OUTPUTS:
%   Standard main_run output folders under results/runs/.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% =========================================================================

clear; clc;
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

previous_scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
previous_plots = getenv('UNIFIED_VSD_DO_PLOTS');
previous_overlay = getenv('UNIFIED_VSD_DO_OVERLAY');
previous_gsa = getenv('UNIFIED_VSD_DO_GSA');

setenv('UNIFIED_VSD_DO_PLOTS', '0');
setenv('UNIFIED_VSD_DO_OVERLAY', '0');
setenv('UNIFIED_VSD_DO_GSA', '0');

scaling_modes = {'zhang', 'lundquist_bsa'};
case_builders = {@patient_reyna, @patient_profile_Razka, @patient_profile_A};
case_labels = {'Reyna full-data', 'Razka sparse-cath', 'Patient A synthetic'};

cleanup_obj = onCleanup(@() restore_environment(previous_scaling_mode, ...
    previous_plots, previous_overlay, previous_gsa));

fprintf('============================================================\n');
fprintf('  UNIFIED VSD MODEL - Scaling Mode Comparison\n');
fprintf('============================================================\n');
fprintf('Runs: %d scaling modes x %d cases = %d main_run calls.\n', ...
    numel(scaling_modes), numel(case_builders), numel(scaling_modes) * numel(case_builders));
fprintf('GSA and plotting are disabled for comparability and speed.\n\n');

for i_mode = 1:numel(scaling_modes)
    scaling_mode = scaling_modes{i_mode};
    setenv('UNIFIED_VSD_SCALING_MODE', scaling_mode);
    fprintf('\n============================================================\n');
    fprintf('  SCALING MODE: %s\n', scaling_mode);
    fprintf('============================================================\n');

    for i_case = 1:numel(case_builders)
        fprintf('\n--- Case %d/%d: %s | scaling=%s ---\n', ...
            i_case, numel(case_builders), case_labels{i_case}, scaling_mode);
        clinical = case_builders{i_case}();
        main_run('pre_surgery', clinical);
    end
end

fprintf('\n============================================================\n');
fprintf('  Scaling comparison finished.\n');
fprintf('  Inspect each run''s validation report and plausibility summary.\n');
fprintf('  Prefer the mode that is robust across Reyna, Razka, and Patient A.\n');
fprintf('============================================================\n');

% RESTORE_ENVIRONMENT - restore user MATLAB session environment variables.
function restore_environment(previous_scaling_mode, previous_plots, previous_overlay, previous_gsa)
restore_one('UNIFIED_VSD_SCALING_MODE', previous_scaling_mode);
restore_one('UNIFIED_VSD_DO_PLOTS', previous_plots);
restore_one('UNIFIED_VSD_DO_OVERLAY', previous_overlay);
restore_one('UNIFIED_VSD_DO_GSA', previous_gsa);
end

% RESTORE_ONE - restore or clear one environment variable.
function restore_one(name, value)
if isempty(value)
    setenv(name, '');
else
    setenv(name, value);
end
end
