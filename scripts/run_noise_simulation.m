function batch_dir = run_noise_simulation(varargin)
% RUN_NOISE_SIMULATION
% -----------------------------------------------------------------------
% Measurement-noise (Monte Carlo) simulation for the paper (Section 2.6).
%
% For each replicate k = 1..N:
%   1. the pre-closure clinical record is perturbed within its declared
%      measurement uncertainty (draw_measurement_noise);
%   2. the model is recalibrated to that record with the SAME procedure as
%      the main result: same scaling, same clinical seeding, same staged
%      multi-start calibration, same number of starts and the same fixed
%      multi-start seed. The fitted parameter set is the one selected by the
%      reference run's GSA screen (loaded from its mask) and is not
%      re-screened, so only the measurements change between replicates;
%   3. the recalibrated parameters are applied with the defect closed to
%      predict the post-closure state, which is compared with the
%      UNPERTURBED post-closure record (still held out, never fitted).
%
% The spread of the calibrated parameters and of the post-closure
% prediction across replicates is the measurement-noise uncertainty.
%
% Each replicate is saved as soon as it finishes; rerunning with the same
% 'BatchDir' skips replicates already done, so an interrupted batch resumes.
%
% USAGE:
%   >> run_noise_simulation('NumReplicates', 2, 'Workers', 0)   % pilot
%   >> run_noise_simulation('NumReplicates', 50, 'Workers', 6)  % paper
%   >> run_noise_simulation('BatchDir', 'results/batches/noise_...')  % resume
%   >> run_noise_simulation('BatchDir', ..., 'NumReplicates', 20)    % resume, first 20 only
%   >> aggregate_noise_simulation('results/batches/noise_...')  % tables only
%
% NAME-VALUE OPTIONS:
%   NumReplicates  number of perturbed records                   [50]
%   NoiseSeed      seed for the measurement-noise draws           [20260926]
%   MultistartSeed seed for the scrambled starts, fixed for all   [101]
%   NumStarts      multi-start count (6 = main result)            [6]
%   ReferenceRun   run folder whose GSA mask is reused
%                  [results/runs/20260920_005419_reyna_pre_surgery]
%   Workers        0 = serial; >0 = parfor over that many workers [0]
%   BatchDir       existing batch folder to resume                ['']
%   OutputRoot     parent folder of new batches         [results/batches]
%
% OUTPUTS (results/batches/noise_<timestamp>/):
%   manifest.json, run_noise_simulation.log
%   noise_spec.csv            perturbed targets, measured values, sigma
%   noise_draws.csv           the perturbed record of every replicate
%   mask_reference.mat        the fixed parameter mask
%   replicates/<run>_mcNNN/   main_run output folder of each replicate
%   replicates/mcNNN.mat      compact result of each replicate
%   summary/                  written by aggregate_noise_simulation
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-26
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
cd(root);
addpath(genpath(root));

p = inputParser;
addParameter(p, 'NumReplicates', 50);
addParameter(p, 'NoiseSeed', 20260926);
addParameter(p, 'MultistartSeed', 101);
addParameter(p, 'NumStarts', 6);
addParameter(p, 'MaxFunEvals', 300);
addParameter(p, 'MaxIterations', 40);
addParameter(p, 'ReferenceRun', fullfile(root, 'results', 'runs', ...
    '20260920_005419_reyna_pre_surgery'));
addParameter(p, 'Workers', 0);
addParameter(p, 'BatchDir', '');
addParameter(p, 'OutputRoot', fullfile(root, 'results', 'batches'));
parse(p, varargin{:});
opt = p.Results;

clinical = patient_reyna();

% ---- batch folder: new, or resume an existing one -------------------------
if isempty(opt.BatchDir)
    batch_id = ['noise_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
    batch_dir = fullfile(opt.OutputRoot, batch_id);
    mkdir(batch_dir);
    is_resume = false;
else
    batch_dir = char(opt.BatchDir);
    assert(isfolder(batch_dir), 'run_noise_simulation:noBatch', 'No batch folder %s', batch_dir);
    saved = jsondecode(fileread(fullfile(batch_dir, 'manifest.json')));
    % A resumed batch must keep the settings it was started with. The one
    % exception is a smaller NumReplicates: the draws are sequential from a
    % fixed stream, so replicates 1..M are the same for any total >= M and
    % truncating to the first M keeps an unbiased, pre-drawn sample. The
    % reduction is recorded in the manifest.
    requested = opt.NumReplicates;
    opt.NumReplicates = saved.num_replicates;
    if ~ismember('NumReplicates', p.UsingDefaults) && requested ~= saved.num_replicates
        assert(requested < saved.num_replicates, 'run_noise_simulation:grow', ...
            'A resumed batch can only be reduced (%d planned, %d requested).', ...
            saved.num_replicates, requested);
        if ~isfield(saved, 'num_replicates_initial')
            saved.num_replicates_initial = saved.num_replicates;
        end
        saved.num_replicates = requested;
        saved.num_replicates_changed = timestamp_now();
        saved.num_replicates_reason = ['reduced to the first ' num2str(requested) ...
            ' pre-drawn replicates because of compute cost (about 84 min per recalibration)'];
        write_json(fullfile(batch_dir, 'manifest.json'), saved);
        opt.NumReplicates = requested;
    end
    opt.NoiseSeed = saved.noise_seed;
    opt.MultistartSeed = saved.multistart_seed;
    opt.NumStarts = saved.num_starts;
    opt.MaxFunEvals = saved.max_function_evaluations;
    opt.MaxIterations = saved.max_iterations;
    opt.ReferenceRun = saved.reference_run;
    is_resume = true;
end
rep_dir = fullfile(batch_dir, 'replicates');
if ~exist(rep_dir, 'dir'), mkdir(rep_dir); end
diary(fullfile(batch_dir, 'run_noise_simulation.log'));
diary on;
cleanup = onCleanup(@() diary('off'));
fprintf('=== run_noise_simulation %s %s ===\n', batch_dir, timestamp_now());

% ---- fixed parameter mask from the reference run --------------------------
mask_file = fullfile(batch_dir, 'mask_reference.mat');
if ~is_resume
    ref = load(fullfile(opt.ReferenceRun, 'mat', 'pre_to_post_seed_latest.mat'));
    ref = ref.pre_to_post_seed;
    optMask = ref.optMask; %#ok<NASGU>
    sobol_ST_threshold = ref.sobol_ST_threshold; %#ok<NASGU>
    reference_active_names = ref.calib_out.names(:);
    reference_names_all = ref.calib_out.names_all(:); %#ok<NASGU>
    reference_primary_metrics = ref.primary_metrics(:); %#ok<NASGU>
    reference_x_active = ref.parameter_vector_calibrated_active(:); %#ok<NASGU>
    save(mask_file, 'optMask', 'sobol_ST_threshold', 'reference_active_names', ...
        'reference_names_all', 'reference_primary_metrics', 'reference_x_active');
    fprintf('Reference run: %s\nFixed active set (%d): %s\n', opt.ReferenceRun, ...
        numel(reference_active_names), strjoin(reference_active_names', ', '));
end
M = load(mask_file);

% ---- perturbed records (all drawn up front: independent of scheduling) -----
[records, draws, spec] = draw_measurement_noise(clinical, opt.NumReplicates, opt.NoiseSeed);
if ~is_resume
    writetable(spec, fullfile(batch_dir, 'noise_spec.csv'));
    writetable(draws, fullfile(batch_dir, 'noise_draws.csv'));
    manifest = build_manifest(root, batch_dir, opt);
    write_json(fullfile(batch_dir, 'manifest.json'), manifest);
end
disp(spec);

% ---- replicates -------------------------------------------------------------
pending = [];
for k = 1:opt.NumReplicates
    if ~isfile(replicate_file(rep_dir, k))
        pending(end + 1) = k; %#ok<AGROW>
    end
end
fprintf('%d of %d replicates pending.\n', numel(pending), opt.NumReplicates);

env = struct( ...
    'UNIFIED_VSD_RUN_MODE', 'exploratory', ...
    'UNIFIED_VSD_SCALING_MODE', 'zhang', ...
    'UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS', '1', ...
    'UNIFIED_VSD_FREEZE_CLINICAL_PROFILE', '1', ...
    'UNIFIED_VSD_DO_GSA', '1', ...
    'UNIFIED_VSD_MASK_FILE', mask_file, ...
    'UNIFIED_VSD_SKIP_FINAL_GSA', '1', ...
    'UNIFIED_VSD_NUM_STARTS', num2str(opt.NumStarts), ...
    'UNIFIED_VSD_MULTISTART_SEED', num2str(opt.MultistartSeed), ...
    'UNIFIED_VSD_MAX_FUN_EVALS', num2str(opt.MaxFunEvals), ...
    'UNIFIED_VSD_MAX_ITERATIONS', num2str(opt.MaxIterations), ...
    'UNIFIED_VSD_FMINCON_PARALLEL', '0', ...
    'UNIFIED_VSD_DO_PLOTS', '0', ...
    'UNIFIED_VSD_DO_OVERLAY', '0', ...
    'UNIFIED_VSD_RUNS_PARENT', rep_dir);

pending_records = records(pending);
if opt.Workers > 0 && numel(pending) > 1
    pool = gcp('nocreate');
    if isempty(pool) || pool.NumWorkers ~= opt.Workers
        delete(pool);
        parpool('Processes', opt.Workers);
    end
    parfor j = 1:numel(pending)
        run_replicate(pending(j), pending_records{j}, clinical, env, rep_dir, M);
    end
else
    for j = 1:numel(pending)
        run_replicate(pending(j), pending_records{j}, clinical, env, rep_dir, M);
        % main_run switches its own diary on and off; resume the batch log.
        diary(fullfile(batch_dir, 'run_noise_simulation.log'));
        diary on;
    end
end

aggregate_noise_simulation(batch_dir);
fprintf('=== run_noise_simulation finished %s ===\n', timestamp_now());
end

% =======================================================================
function run_replicate(k, clinical_k, clinical_original, env, rep_dir, M)
% One replicate: recalibrate to the perturbed record, predict post-closure.
% Failures are saved (status 'failed') rather than thrown, so one bad
% replicate never stops the batch and is still reported.
fields = fieldnames(env);
for i = 1:numel(fields)
    setenv(fields{i}, env.(fields{i}));
end
tag = sprintf('mc%03d', k);
setenv('UNIFIED_VSD_RUN_TAG', tag);
set(groot, 'DefaultFigureVisible', 'off');

r = struct('replicate', k, 'status', 'failed', 'message', '', 'run_dir', '', ...
    'wall_minutes', NaN);
t0 = tic;
try
    main_run('pre_surgery', clinical_k);
    hits = dir(fullfile(rep_dir, ['*_' tag]));
    assert(numel(hits) == 1, 'run_noise_simulation:runFolder', ...
        'Expected one run folder for %s, found %d.', tag, numel(hits));
    run_dir = fullfile(hits(1).folder, hits(1).name);
    r.run_dir = run_dir;

    S = load(fullfile(run_dir, 'mat', 'pre_to_post_seed_latest.mat'));
    S = S.pre_to_post_seed;
    active = S.calib_out.names(:);
    assert(isequal(active, M.reference_active_names(:)), ...
        'run_noise_simulation:activeSetChanged', ...
        'Replicate %d fitted a different parameter set: %s', k, strjoin(active', ', '));
    assert(isequal(S.primary_metrics(:), M.reference_primary_metrics(:)), ...
        'run_noise_simulation:primaryChanged', ...
        'Replicate %d used different primary metrics: %s', k, strjoin(S.primary_metrics, ', '));

    ac = S.accepted_candidate;
    r.parameter_names = active;
    r.parameter_values = S.parameter_vector_calibrated_active(:);
    r.parameter_initial = S.parameter_vector_initial_active(:);
    r.calibration_status = ac.calibration_status.label;
    r.rollback_applied = isfield(ac, 'rollback_applied') && ac.rollback_applied;
    r.multistart_best_start = S.calib_out.multistart.best_start_index;
    r.multistart_table = S.calib_out.multistart.table;

    % Pre-closure fit: against the perturbed record it was fitted to, and
    % against the original measured record.
    gate = ac.report.full_metric_gate.table;
    r.pre_gate_table = gate;
    r.pre_vs_original = compare_metrics(ac.metrics, get_calibration_targets('pre_surgery', clinical_original));
    r.pre_vs_perturbed = compare_metrics(ac.metrics, get_calibration_targets('pre_surgery', clinical_k));
    r.chi2_per_obs = ac.report.chi_squared.chi2_per_obs;

    % Post-closure prediction against the UNPERTURBED held-out record.
    pred = evaluate_post_closure_prediction( ...
        fullfile(run_dir, 'mat', 'params_accepted_candidate_pre_surgery.mat'), ...
        clinical_original, false);
    r.post_table = pred.table;
    r.post_within10 = pred.within10;
    r.post_n = pred.n;
    r.post_chi2_per_n = pred.chi2 / max(pred.n, 1);
    r.post_metrics = pred.metrics;
    r.status = 'completed';
catch err
    r.message = getReport(err, 'extended', 'hyperlinks', 'off');
    fprintf(2, '[run_noise_simulation] replicate %d failed:\n%s\n', k, r.message);
end
r.wall_minutes = toc(t0) / 60;
save_replicate(replicate_file(rep_dir, k), r);
fprintf('[run_noise_simulation] replicate %d %s in %.1f min\n', k, r.status, r.wall_minutes);
end

function save_replicate(path, r) %#ok<INUSD>
save(path, 'r');
end

function f = replicate_file(rep_dir, k)
f = fullfile(rep_dir, sprintf('mc%03d.mat', k));
end

function T = compare_metrics(metrics, targets)
names = {}; measured = []; model = [];
for i = 1:numel(targets)
    t = targets(i);
    if ~isfinite(t.ClinicalValue) || ~isfield(metrics, t.Metric), continue; end
    names{end + 1, 1} = t.Metric; %#ok<AGROW>
    measured(end + 1, 1) = t.ClinicalValue; %#ok<AGROW>
    model(end + 1, 1) = metrics.(t.Metric); %#ok<AGROW>
end
T = table(names, measured, model, 100 * (model - measured) ./ abs(measured), ...
    'VariableNames', {'Metric', 'Measured', 'Model', 'Error_pct'});
end

function m = build_manifest(root, batch_dir, opt)
m = struct();
m.analysis = 'measurement-noise Monte Carlo simulation';
m.batch_dir = batch_dir;
m.started = timestamp_now();
m.num_replicates = opt.NumReplicates;
m.noise_seed = opt.NoiseSeed;
m.multistart_seed = opt.MultistartSeed;
m.num_starts = opt.NumStarts;
m.max_function_evaluations = opt.MaxFunEvals;
m.max_iterations = opt.MaxIterations;
m.reference_run = opt.ReferenceRun;
m.workers = opt.Workers;
m.noise_model = ['independent Gaussian, sigma = declared 1-sigma uncertainty of each ' ...
    'measured pre-closure target; MAP and shunt flow recomputed from the perturbed ' ...
    'values; inputs and post-closure record unperturbed; physically impossible draws redrawn'];
[~, m.git_commit] = system(sprintf('git -C "%s" rev-parse HEAD', root));
[~, m.git_branch] = system(sprintf('git -C "%s" rev-parse --abbrev-ref HEAD', root));
[~, dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=no', root));
m.git_commit = strtrim(m.git_commit);
m.git_branch = strtrim(m.git_branch);
m.git_dirty = ~isempty(strtrim(dirty));
m.matlab_version = version;
end

function t = timestamp_now()
t = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end

function write_json(path, s)
fid = fopen(path, 'w');
fprintf(fid, '%s', jsonencode(s, 'PrettyPrint', true));
fclose(fid);
end
