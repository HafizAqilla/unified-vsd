%% test_protocol_manifest.m
% Unit tests for protocol manifest writer.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'config'), '-begin');
addpath(fullfile(project_root, 'src', 'utils'), '-begin');

fprintf('===============================================\n');
fprintf('  UNIFIED VSD MODEL - Protocol Manifest Test\n');
fprintf('===============================================\n\n');

protocol = load_run_protocol(fullfile(project_root, 'config', 'run_protocols', ...
    'reyna_pre_post_publishable.m'), 'RootDir', project_root);
seed_resolution = struct('Status', 'not_requested', 'ResolvedPath', '');
out_dir = fullfile(tempdir(), 'unified_vsd_protocol_manifest_test');
if exist(out_dir, 'dir')
    rmdir(out_dir, 's');
end
manifest_file = write_protocol_manifest(out_dir, protocol, seed_resolution);
raw = fileread(manifest_file);

assert(isfile(manifest_file), 'Protocol manifest was not written.');
assert(contains(raw, '"ProtocolID"'), 'Manifest missing ProtocolID.');
assert(contains(raw, 'reyna_pre_post_publishable_v1'), ...
    'Manifest missing protocol identifier.');
assert(contains(raw, '"ScalingMode"'), 'Manifest missing scaling mode.');
assert(contains(raw, 'zhang'), 'Manifest missing Zhang policy.');

fprintf('  [PASS] protocol manifest JSON includes protocol and scaling policy.\n');

