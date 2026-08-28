%% test_scaling_method_provenance.m
% Unit tests for scaling provenance artifact writer.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
addpath(fullfile(project_root, 'config'), '-begin');
addpath(fullfile(project_root, 'src', 'utils'), '-begin');

fprintf('=================================================\n');
fprintf('  UNIFIED VSD MODEL - Scaling Provenance Test\n');
fprintf('=================================================\n\n');

out_dir = fullfile(tempdir(), 'unified_vsd_scaling_provenance_test');
if exist(out_dir, 'dir')
    rmdir(out_dir, 's');
end
mkdir(out_dir);

policy = resolve_scaling_policy('', 'publication');
provenance_file = write_scaling_method_provenance(out_dir, policy, 'pass');
tbl = readtable(provenance_file, 'TextType', 'string');

assert(isfile(provenance_file), 'Provenance CSV was not written.');
assert(height(tbl) == 1, 'Provenance CSV should contain one row.');
assert(tbl.ScalingMode(1) == "zhang", 'Expected Zhang scaling mode.');
assert(tbl.ScalingRole(1) == "primary_prior", 'Expected primary prior role.');
assert(contains(tbl.DOI(1), "10.1016"), 'Expected Zhang DOI.');
assert(contains(tbl.ImplementationVariant(1), "zhang_weight"), ...
    'Expected Zhang implementation variant.');
assert(tbl.BaselineGateStatus(1) == "pass", 'Expected baseline gate status.');

fprintf('  [PASS] scaling provenance CSV includes citation and implementation caveats.\n');

