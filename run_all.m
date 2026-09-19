function batch_dir = run_all(varargin)
% RUN_ALL
% -----------------------------------------------------------------------
% Reproduces the paper results from the clinical record: for each seed,
% one full pre-closure pipeline (scaling -> GSA screening -> multi-start
% calibration of the GSA-selected parameters -> post-calibration GSA ->
% identifiability) followed by the zero-shot post-closure prediction and a
% raw-data export. Runs are sequential; a failed seed is logged and the
% next seed still runs. Every run is exported as soon as it finishes.
%
% OUTPUTS (results/batches/<batch_id>/):
%   run_all.log            full console output of the whole batch
%   manifest.json          code version, software versions, settings, seeds
%   runs.csv               one row per seed: folder, status, times (updated
%                          after every seed, so it survives an interruption)
%   figure_data/           combined CSVs across seeds (aggregate_supplement)
%   figures/               publication figures (R), when R is available
%
% USAGE:
%   >> run_all                       % seeds 101, 202, 303
%   >> run_all('Seeds', 101)         % one seed
%   >> run_all('FiguresOnly', 'results/batches/<batch_id>')
% -----------------------------------------------------------------------

root = fileparts(mfilename('fullpath'));
cd(root);
addpath(genpath(root));

p = inputParser;
addParameter(p, 'Seeds', [101 202 303]);
addParameter(p, 'NumStarts', 6);
addParameter(p, 'GsaN', 128);
addParameter(p, 'MaxFunEvals', 300);
addParameter(p, 'MaxIterations', 40);
addParameter(p, 'FiguresOnly', '');
parse(p, varargin{:});
opt = p.Results;

if ~isempty(opt.FiguresOnly)
    batch_dir = opt.FiguresOnly;
    make_figures(batch_dir);
    return;
end

batch_id = ['paper_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
batch_dir = fullfile(root, 'results', 'batches', batch_id);
mkdir(batch_dir);
diary(fullfile(batch_dir, 'run_all.log'));
diary on;
cleanup = onCleanup(@() diary('off'));
fprintf('=== run_all %s started %s ===\n', batch_id, timestamp_now());

settings = struct( ...
    'UNIFIED_VSD_RUN_MODE', 'exploratory', ...
    'UNIFIED_VSD_SCALING_MODE', 'zhang', ...
    'UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS', '1', ...
    'UNIFIED_VSD_FREEZE_CLINICAL_PROFILE', '1', ...
    'UNIFIED_VSD_DO_GSA', '1', ...
    'UNIFIED_VSD_GSA_PCE_N', num2str(opt.GsaN), ...
    'UNIFIED_VSD_NUM_STARTS', num2str(opt.NumStarts), ...
    'UNIFIED_VSD_MAX_FUN_EVALS', num2str(opt.MaxFunEvals), ...
    'UNIFIED_VSD_MAX_ITERATIONS', num2str(opt.MaxIterations), ...
    'UNIFIED_VSD_FMINCON_PARALLEL', '0', ...
    'UNIFIED_VSD_DO_PLOTS', '0', ...
    'UNIFIED_VSD_DO_OVERLAY', '0');
for f = fieldnames(settings)'
    setenv(f{1}, settings.(f{1}));
end
set(groot, 'DefaultFigureVisible', 'off');

manifest = build_manifest(root, batch_id, opt, settings);
write_json(fullfile(batch_dir, 'manifest.json'), manifest);
fprintf('Code: %s @ %s (uncommitted changes: %s)\n', manifest.git_branch, ...
    manifest.git_commit, string(manifest.git_dirty));

clinical = patient_reyna();
runs = table();
for k = 1:numel(opt.Seeds)
    seed = opt.Seeds(k);
    setenv('UNIFIED_VSD_MULTISTART_SEED', num2str(seed));
    t_start = datetime('now');
    fprintf('\n=== Seed %d (%d of %d) started %s ===\n', seed, k, numel(opt.Seeds), timestamp_now());
    before = list_run_dirs(root);
    run_dir = ''; status = 'failed'; message = ''; export_status = 'not attempted';
    try
        main_run('pre_surgery', clinical);
        run_dir = newest_new_run_dir(root, before);
        status = 'completed';
        try
            export_run_supplement(run_dir, clinical);
            export_status = 'exported';
        catch err
            export_status = ['export failed: ' err.message];
            fprintf(2, '%s\n', getReport(err, 'extended'));
        end
    catch err
        message = err.message;
        fprintf(2, '%s\n', getReport(err, 'extended'));
    end
    t_end = datetime('now');
    row = table(seed, {run_dir}, {status}, {export_status}, {message}, ...
        {char(t_start, 'yyyy-MM-dd HH:mm:ss')}, {char(t_end, 'yyyy-MM-dd HH:mm:ss')}, ...
        seconds(t_end - t_start) / 60, ...
        'VariableNames', {'Seed','RunFolder','Status','Export','Error','Start','End','WallMinutes'});
    runs = [runs; row]; %#ok<AGROW>
    writetable(runs, fullfile(batch_dir, 'runs.csv'));
    fprintf('=== Seed %d %s (%s) in %.1f min ===\n', seed, status, export_status, row.WallMinutes);
end

ok = strcmp(runs.Export, 'exported');
if any(ok)
    aggregate_supplement(runs.RunFolder(ok), runs.Seed(ok), fullfile(batch_dir, 'figure_data'));
    make_figures(batch_dir);
end
fprintf('\n=== run_all %s finished %s: %d of %d seeds exported ===\n', ...
    batch_id, timestamp_now(), nnz(ok), height(runs));
end

% =======================================================================
function make_figures(batch_dir)
rscript = find_rscript();
script = fullfile(fileparts(mfilename('fullpath')), 'scripts', 'R', 'make_figures.R');
if isempty(rscript)
    fprintf(2, 'R not found; figures skipped. Run later with run_all(''FiguresOnly'', ...).\n');
    return;
end
cmd = sprintf('"%s" "%s" "%s" "%s"', rscript, script, ...
    fullfile(batch_dir, 'figure_data'), fullfile(batch_dir, 'figures'));
[status, out] = system(cmd);
fprintf('%s\n', out);
fprintf('[run_all] R figures exit status %d\n', status);
end

function r = find_rscript()
r = '';
cands = [{'Rscript'}, arrayfun(@(d) fullfile(d.folder, d.name, 'bin', 'Rscript.exe'), ...
    dir('C:\Program Files\R\R-*'), 'UniformOutput', false)'];
for c = cands
    if strcmp(c{1}, 'Rscript')
        [s, ~] = system('Rscript --version');
        if s == 0, r = 'Rscript'; return; end
    elseif exist(c{1}, 'file')
        r = c{1}; return;
    end
end
end

function m = build_manifest(root, batch_id, opt, settings)
m = struct();
m.batch_id = batch_id;
m.started = timestamp_now();
m.seeds = opt.Seeds;
m.num_starts = opt.NumStarts;
m.gsa_design_size = opt.GsaN;
m.max_function_evaluations = opt.MaxFunEvals;
m.max_iterations = opt.MaxIterations;
m.gsa_threshold = 0.10;
m.environment = settings;
[~, m.git_commit] = system(sprintf('git -C "%s" rev-parse HEAD', root));
[~, m.git_branch] = system(sprintf('git -C "%s" rev-parse --abbrev-ref HEAD', root));
[~, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', root));
m.git_commit = strtrim(m.git_commit);
m.git_branch = strtrim(m.git_branch);
m.git_dirty = ~isempty(strtrim(dirty));
m.git_dirty_files = strtrim(dirty);
m.matlab_version = version;
v = ver;
m.toolboxes = arrayfun(@(t) sprintf('%s %s', t.Name, t.Version), v, 'UniformOutput', false);
m.computer = computer;
m.uqlab_path = '';
d = dir(fullfile(root, 'toolbox', 'UQLab*'));
if ~isempty(d), m.uqlab_path = d(1).name; end
rs = find_rscript();
m.r_version = '';
if ~isempty(rs)
    [~, out] = system(sprintf('"%s" --version', rs));
    m.r_version = strtrim(out);
end
end

function d = list_run_dirs(root)
x = dir(fullfile(root, 'results', 'runs', '*_pre_surgery'));
d = fullfile({x.folder}, {x.name});
end

function run_dir = newest_new_run_dir(root, before)
after = list_run_dirs(root);
new = setdiff(after, before);
assert(~isempty(new), 'run_all:noRunFolder', 'main_run did not create a run folder.');
info = cellfun(@dir, new, 'UniformOutput', false);
[~, i] = max(cellfun(@(s) max([s.datenum]), info));
run_dir = new{i};
end

function t = timestamp_now()
t = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end

function write_json(path, s)
fid = fopen(path, 'w');
fprintf(fid, '%s', jsonencode(s, 'PrettyPrint', true));
fclose(fid);
end
