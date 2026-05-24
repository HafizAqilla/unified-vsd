%% run_virtual_patients.m
% =========================================================================
% Quick-start script: runs both virtual patient profiles from the
% Pediatric_VSD_Clinical_Parameters.md benchmark document.
%
% USAGE — in MATLAB Command Window:
%   run('scripts/run_virtual_patients.m')    % run both profiles
%   PROFILE = 'A'; run('scripts/run_virtual_patients.m')   % Profile A only
%   PROFILE = 'B'; run('scripts/run_virtual_patients.m')   % Profile B only
%
% PROFILES:
%   A — "The High-Flow Infant"    : 3.7 kg, 1.6 mo, Qp/Qs 3.44, PAP mean 28 mmHg
%   B — "The High-Pressure Infant": 4.5 kg, 2.0 mo, Qp/Qs 1.79, PAP mean 43 mmHg
%
% SOURCE:
%   Pediatric_VSD_Clinical_Parameters.md, Tables 7.1 and 7.2
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% =========================================================================

%% ---- setup -------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

if ~exist('PROFILE', 'var')
    PROFILE = 'both';   % 'A', 'B', or 'both'
end

%% ---- Profile A — High-Flow Infant --------------------------------------
if any(strcmpi(PROFILE, {'A', 'both'}))
    fprintf('\n========================================\n');
    fprintf('  PROFILE A — High-Flow Infant\n');
    fprintf('  3.7 kg | 1.6 mo | Qp/Qs 3.44 | PAP mean 28 mmHg\n');
    fprintf('========================================\n');

    clinical_A = patient_profile_A();
    main_run('pre_surgery', clinical_A);
end

%% ---- Profile B — High-Pressure Infant (removed 2026-05-24: orphaned config)
% Use patient_template() to create new synthetic benchmark profiles.

fprintf('\n[run_virtual_patients] Done.\n');
