% RUN_POST_SURGERY
% -----------------------------------------------------------------------
% Entry-point script for running the post-surgery VSD closure scenario.
%
% This script mirrors Keisya's post-surgery workflow while using the current
% run-folder architecture:
%   1. Load the patient profile.
%   2. Apply editable post-operative targets.
%   3. Load the latest pre-to-post seed from a completed pre-surgery run.
%   4. Inject the calibrated pre-op parameters into main_run('post_surgery').
%
% HOW TO RUN:
%   >> run run_post_surgery
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-18
% VERSION:  1.0
% -----------------------------------------------------------------------

clear; clc;

root = fileparts(mfilename('fullpath'));
restoredefaultpath();
root_paths = strsplit(genpath(root), pathsep);
is_shadow = contains(root_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(root_paths, [filesep '.git' filesep], 'IgnoreCase', true);
addpath(strjoin(root_paths(~is_shadow), pathsep));

fprintf('=====================================================\n');
fprintf('  POST-SURGERY SIMULATION - run_post_surgery.m\n');
fprintf('=====================================================\n\n');

% Keep the post-op handoff quick by default. Remove this line if a full
% post-op PCE sensitivity pass is required.
setenv('UNIFIED_VSD_DO_GSA', '0');

clinical = patient_reyna();
clinical = apply_post_surgery_targets(clinical);

patient_label = resolve_patient_label(clinical);
seed_file = find_latest_pre_to_post_seed(root, patient_label);
params_pre = load_pre_surgery_params(seed_file);
clinical.pre_surgery.CalibParams = params_pre;

fprintf('[run_post_surgery] Patient: %s\n', patient_label);
fprintf('[run_post_surgery] Loaded pre-op seed:\n  %s\n', seed_file);
fprintf('[run_post_surgery] Injecting calibrated pre-op physiology; main_run will close R.vsd.\n\n');

main_run('post_surgery', clinical);

fprintf('\n[run_post_surgery] Post-surgery run complete.\n');

function clinical = apply_post_surgery_targets(clinical)
% APPLY_POST_SURGERY_TARGETS - set editable post-op targets [clinical units].
post = clinical.post_surgery;

post.QpQs          = 1.0;      % [-] closed VSD target, no residual shunt
post.PAP_sys_mmHg  = 17.0;     % [mmHg]
post.PAP_dia_mmHg  = 8.67;     % [mmHg]
post.PAP_mean_mmHg = 13.0;     % [mmHg]
post.PVR_WU        = NaN;      % [WU]

post.SAP_sys_mmHg  = 104.3;    % [mmHg]
post.SAP_dia_mmHg  = 71.67;    % [mmHg]
post.MAP_mmHg      = post.SAP_dia_mmHg + ...
    (post.SAP_sys_mmHg - post.SAP_dia_mmHg) / 3; % [mmHg]
post.SVR_WU        = NaN;      % [WU]

post.RAP_mean_mmHg = 5.0;      % [mmHg]
post.RAP_sys_mmHg  = 7.67;     % [mmHg]
post.RAP_dia_mmHg  = 5.0;      % [mmHg]
post.LAP_mean_mmHg = NaN;      % [mmHg]

post.LVEDV_mL      = 32.0;     % [mL]
post.EF            = 0.618;    % [-] fraction, not percent
post.LVESV_mL      = post.LVEDV_mL * (1 - post.EF); % [mL]
post.RVEDV_mL      = 30.5;     % [mL]
post.RVESV_mL      = 12.0;     % [mL]
post.RVEF          = (post.RVEDV_mL - post.RVESV_mL) / post.RVEDV_mL; % [-]
post.CO_Lmin       = 3.0;      % [L/min]

clinical.post_surgery = post;
end

function label = resolve_patient_label(clinical)
% RESOLVE_PATIENT_LABEL - return a run-folder-safe patient label.
label = 'patient';
if isfield(clinical, 'common') && isfield(clinical.common, 'patient_name') && ...
        ~isempty(clinical.common.patient_name)
    label = lower(regexprep(char(clinical.common.patient_name), '[^a-zA-Z0-9]+', '_'));
end
end

function seed_file = find_latest_pre_to_post_seed(root, patient_label)
% FIND_LATEST_PRE_TO_POST_SEED - locate newest pre-op handoff MAT file.
patient_patterns = {
    fullfile(root, 'results', 'runs', sprintf('*_%s_pre_surgery', patient_label), 'mat', 'pre_to_post_seed_latest.mat')
    fullfile(root, 'results', 'runs', sprintf('*_%s_pre_surgery', patient_label), 'mat', 'pre_to_post_seed_*.mat')
    };
fallback_patterns = {
    fullfile(root, 'results', 'runs', '*_pre_surgery', 'mat', 'pre_to_post_seed_latest.mat')
    fullfile(root, 'results', 'runs', '*_pre_surgery', 'mat', 'pre_to_post_seed_*.mat')
    fullfile(root, 'results', 'tables', 'params_calibrated_pre_surgery*.mat')
    };

seed_file = newest_matching_file(patient_patterns);
if isempty(seed_file)
    seed_file = newest_matching_file(fallback_patterns);
end
if isempty(seed_file)
    error('run_post_surgery:noPreOpSeed', ...
        ['No pre-surgery seed was found. Run main_run(''pre_surgery'', clinical) ', ...
         'first, then rerun run_post_surgery.m.']);
end
end

function params_pre = load_pre_surgery_params(seed_file)
% LOAD_PRE_SURGERY_PARAMS - read calibrated pre-op params from modern or legacy MAT files.
loaded = load(seed_file);

if isfield(loaded, 'pre_to_post_seed')
    seed = loaded.pre_to_post_seed;
    if isfield(seed, 'accepted_candidate') && isfield(seed.accepted_candidate, 'params')
        params_pre = seed.accepted_candidate.params;
    elseif isfield(seed, 'params_calibrated')
        params_pre = seed.params_calibrated;
    elseif isfield(seed, 'params_post_seed_closed_vsd')
        params_pre = seed.params_post_seed_closed_vsd;
    else
        error('run_post_surgery:badSeedPayload', ...
            'pre_to_post_seed does not contain calibrated parameter fields.');
    end
elseif isfield(loaded, 'params_cal')
    params_pre = loaded.params_cal;
else
    error('run_post_surgery:badLegacyPayload', ...
        'MAT file does not contain pre_to_post_seed or params_cal.');
end
end

function file_path = newest_matching_file(patterns)
% NEWEST_MATCHING_FILE - return newest file across a list of wildcard patterns.
listing_cells = cell(numel(patterns), 1);
for pattern_idx = 1:numel(patterns)
    listing_cells{pattern_idx} = dir(patterns{pattern_idx});
end
candidates = vertcat(listing_cells{:});
file_path = '';
if isempty(candidates)
    return;
end
[~, order] = sort([candidates.datenum], 'descend');
latest = candidates(order(1));
file_path = fullfile(latest.folder, latest.name);
end
