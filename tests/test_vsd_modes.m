%% test_vsd_modes.m
% Regression test for VSD mode switching.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

params = default_parameters();

params.vsd.mode = 'linear_bidirectional';
params.R.vsd = 10;
Q_lin = vsd_shunt_model(80, 20, params);
assert(Q_lin > 0, 'Linear bidirectional mode should give positive L->R flow.');

params.vsd.mode = 'orifice_bidirectional';
params.vsd.area_mm2 = pi * (4 / 2)^2;
Q_orifice_small = vsd_shunt_model(80, 20, params);
params.vsd.area_mm2 = pi * (6 / 2)^2;
Q_orifice_large = vsd_shunt_model(80, 20, params);
assert(Q_orifice_large > Q_orifice_small, 'Larger VSD area should increase shunt flow.');

Q_reverse = vsd_shunt_model(20, 80, params);
assert(Q_reverse < 0, 'Bidirectional orifice mode should allow reverse shunt.');

params.vsd.mode = 'linear_left_to_right_only';
Q_legacy_reverse = vsd_shunt_model(20, 80, params);
assert(Q_legacy_reverse > -1e-3, 'Legacy mode should suppress strong reverse shunt.');

fprintf('[PASS] VSD mode switching behaves as expected.\n');
