%% test_protocol_seed_resolution.m
% Unit tests for protocol loading and pre-to-post seed resolution.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'config'), '-begin');
addpath(fullfile(project_root, 'src', 'utils'), '-begin');

fprintf('================================================\n');
fprintf('  UNIFIED VSD MODEL - Protocol Seed Test\n');
fprintf('================================================\n\n');

n_pass = 0;
n_fail = 0;

%% Test 1: Publication protocol defaults to Zhang and disables fallback.
protocol_path = fullfile(project_root, 'config', 'run_protocols', ...
    'reyna_pre_post_publishable.m');
protocol = load_run_protocol(protocol_path, 'RootDir', project_root);
if strcmp(protocol.Mode, 'publication') && strcmp(protocol.ScalingMode, 'zhang') && ...
        strcmp(protocol.ScalingRole, 'primary_prior') && ~protocol.AllowSeedFallback
    fprintf('  [PASS] Publication protocol loads with Zhang primary and no fallback.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Publication protocol defaults mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 2: Publication protocol without seed path errors.
did_error = false;
try
    resolve_pre_to_post_seed(protocol, project_root);
catch ME
    did_error = strcmp(ME.identifier, ...
        'resolve_pre_to_post_seed:publicationSeedRequired');
end
if did_error
    fprintf('  [PASS] Publication mode requires pinned seed.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Publication mode did not require pinned seed.\n');
    n_fail = n_fail + 1;
end

%% Test 3: Explicit pinned seed resolves.
tmp_root = fullfile(tempdir(), 'unified_vsd_seed_resolution_test');
if exist(tmp_root, 'dir')
    rmdir(tmp_root, 's');
end
mkdir(tmp_root);
seed_file = fullfile(tmp_root, 'seed.mat');
pre_to_post_seed = struct();
pre_to_post_seed.source_scenario = 'pre_surgery';
pre_to_post_seed.timestamp = '20260620_000000';
pre_to_post_seed.calibration_status = struct('label', 'ACCEPT');
pre_to_post_seed.scaling_policy = resolve_scaling_policy('', 'publication');
save(seed_file, 'pre_to_post_seed');

protocol.PreSeedPath = seed_file;
seed_resolution = resolve_pre_to_post_seed(protocol, project_root);
if strcmp(seed_resolution.Status, 'resolved') && ...
        strcmp(seed_resolution.Source, 'pinned') && ...
        strcmp(seed_resolution.Metadata.CalibrationStatus, 'ACCEPT')
    fprintf('  [PASS] Pinned seed resolves with metadata.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Pinned seed resolution mismatch.\n');
    n_fail = n_fail + 1;
end

%% Test 4: Exploratory fallback finds latest seed.
fallback_dir = fullfile(tmp_root, 'results', 'runs', ...
    '20260620_000000_reyna_pre_surgery', 'mat');
mkdir(fallback_dir);
fallback_seed = fullfile(fallback_dir, 'pre_to_post_seed_latest.mat');
save(fallback_seed, 'pre_to_post_seed');
protocol_fallback = protocol;
protocol_fallback.Mode = 'exploratory';
protocol_fallback.PreSeedPath = '';
protocol_fallback.AllowSeedFallback = true;
seed_fallback = resolve_pre_to_post_seed(protocol_fallback, tmp_root);
if strcmp(seed_fallback.Status, 'resolved') && ...
        strcmp(seed_fallback.Source, 'latest_fallback') && ...
        contains(seed_fallback.ResolvedPath, 'pre_to_post_seed_latest.mat')
    fprintf('  [PASS] Exploratory fallback seed discovery works.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Exploratory fallback seed discovery mismatch.\n');
    n_fail = n_fail + 1;
end

fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
fprintf('==========================================\n');

assert(n_fail == 0, 'test_protocol_seed_resolution:failed', ...
    '%d protocol seed test(s) failed.', n_fail);

