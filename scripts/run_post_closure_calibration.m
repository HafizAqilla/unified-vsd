function batch_dir = run_post_closure_calibration(varargin)
% RUN_POST_CLOSURE_CALIBRATION
% -----------------------------------------------------------------------
% Post-closure calibration for the paper (Section "Post-closure
% calibration"). Runs AFTER the zero-shot prediction and does not change
% it: the pre-closure procedure is applied to the post-closure record so
% that the parameters of the two states can be compared.
%
% Same settings as run_all (Zhang scaling, clinical seeding, 6 starts,
% multi-start seed 101, 300 evaluations / 40 iterations per stage). The
% GSA screen is not repeated: the free parameters are the pre-closure
% GSA-selected set without vsd.Cd, fixed in
% config/calibration_recipes/reyna_post_surgery.m together with the
% targets (seven catheter pressures and Qs; Qp/Qs and volumes compared only).
%
% USAGE:
%   >> run_post_closure_calibration
%
% OUTPUTS (results/batches/postcal_<timestamp>/):
%   run_post_closure_calibration.log, manifest.json
%   <timestamp>_reyna_post_surgery_postcal/   main_run output folder
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-10-07
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath(genpath(root));

p = inputParser;
addParameter(p, 'MultistartSeed', 101);
addParameter(p, 'NumStarts', 6);
addParameter(p, 'MaxFunEvals', 300);
addParameter(p, 'MaxIterations', 40);
parse(p, varargin{:});
opt = p.Results;

batch_dir = fullfile(root, 'results', 'batches', ...
    ['postcal_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
mkdir(batch_dir);
diary(fullfile(batch_dir, 'run_post_closure_calibration.log'));
diary on;
cleanup = onCleanup(@() diary('off'));

env = struct( ...
    'UNIFIED_VSD_RUN_MODE', 'exploratory', ...
    'UNIFIED_VSD_SCALING_MODE', 'zhang', ...
    'UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS', '1', ...
    'UNIFIED_VSD_FREEZE_CLINICAL_PROFILE', '1', ...
    'UNIFIED_VSD_DO_GSA', '0', ...
    'UNIFIED_VSD_NUM_STARTS', num2str(opt.NumStarts), ...
    'UNIFIED_VSD_MULTISTART_SEED', num2str(opt.MultistartSeed), ...
    'UNIFIED_VSD_MAX_FUN_EVALS', num2str(opt.MaxFunEvals), ...
    'UNIFIED_VSD_MAX_ITERATIONS', num2str(opt.MaxIterations), ...
    'UNIFIED_VSD_FMINCON_PARALLEL', '0', ...
    'UNIFIED_VSD_DO_PLOTS', '0', ...
    'UNIFIED_VSD_DO_OVERLAY', '0', ...
    'UNIFIED_VSD_RUNS_PARENT', batch_dir, ...
    'UNIFIED_VSD_RUN_TAG', 'postcal');
for f = fieldnames(env)'
    setenv(f{1}, env.(f{1}));
end
set(groot, 'DefaultFigureVisible', 'off');

manifest = struct();
manifest.analysis = 'post-closure calibration (after the zero-shot prediction)';
manifest.started = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
manifest.recipe = 'config/calibration_recipes/reyna_post_surgery.m';
manifest.settings = env;
[~, manifest.git_commit] = system(sprintf('git -C "%s" rev-parse HEAD', root));
[~, manifest.git_branch] = system(sprintf('git -C "%s" rev-parse --abbrev-ref HEAD', root));
[~, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', root));
manifest.git_commit = strtrim(manifest.git_commit);
manifest.git_branch = strtrim(manifest.git_branch);
manifest.git_dirty_files = strtrim(dirty);
manifest.matlab_version = version;
fid = fopen(fullfile(batch_dir, 'manifest.json'), 'w');
fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
fclose(fid);

fprintf('=== run_post_closure_calibration %s ===\n', batch_dir);
main_run('post_surgery', patient_reyna());
diary(fullfile(batch_dir, 'run_post_closure_calibration.log'));
diary on;
fprintf('=== run_post_closure_calibration finished %s ===\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end
