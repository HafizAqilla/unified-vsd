function run_quick_regression_suite()
% RUN_QUICK_REGRESSION_SUITE
% -----------------------------------------------------------------------
% Run a lightweight no-plot, no-GSA regression sweep across benchmark
% patient profiles.
%
% This helper is intended for overnight or checkpoint verification after
% structural calibration changes.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-02
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(build_clean_project_path(root));

setenv('UNIFIED_VSD_DO_PLOTS', '0');
setenv('UNIFIED_VSD_DO_OVERLAY', '0');
setenv('UNIFIED_VSD_DO_GSA', '0');

cases = {
    'pre_surgery', patient_reyna()
    'pre_surgery', patient_profile_Razka()
    'pre_surgery', patient_profile_A()
    };

for i = 1:size(cases, 1)
    scenario = cases{i, 1};
    clinical = cases{i, 2};
    fprintf('\n[run_quick_regression_suite] Case %d/%d: %s (%s)\n', ...
        i, size(cases, 1), resolve_case_label(clinical), scenario);
    main_run(scenario, clinical);
end

function project_path = build_clean_project_path(root)
% BUILD_CLEAN_PROJECT_PATH — exclude shadow worktrees from regression runs.
root_paths = strsplit(genpath(root), pathsep);
root_paths = root_paths(~cellfun('isempty', root_paths));
is_shadow = contains(root_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, root_paths);
project_path = strjoin(root_paths(~is_shadow & is_existing), pathsep);
end
end

function label = resolve_case_label(clinical)
% RESOLVE_CASE_LABEL — robust patient label for benchmark logging.
label = 'unknown_case';
if ~isstruct(clinical) || ~isfield(clinical, 'common')
    return;
end

if isfield(clinical.common, 'patient_name') && ~isempty(clinical.common.patient_name)
    label = clinical.common.patient_name;
elseif isfield(clinical.common, 'patient_id') && ~isempty(clinical.common.patient_id)
    label = clinical.common.patient_id;
elseif isfield(clinical.common, 'weight_kg') && isfield(clinical.common, 'height_cm')
    label = sprintf('patient_%0.1fkg_%0.0fcm', ...
        clinical.common.weight_kg, clinical.common.height_cm);
end
end
