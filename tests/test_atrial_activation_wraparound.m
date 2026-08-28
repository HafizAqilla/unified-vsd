%% test_atrial_activation_wraparound.m
% Regression test for atrial activation across heartbeat wraparound.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Atrial Wraparound Test\n');
fprintf('==========================================\n\n');

params = default_parameters();
params.HR = 60;
T_HB = 60 / params.HR;

params.Tc_LV = params.Tc_LV_frac * T_HB;
params.Tr_LV = params.Tr_LV_frac * T_HB;
params.Tc_RV = params.Tc_RV_frac * T_HB;
params.Tr_RV = params.Tr_RV_frac * T_HB;

params.t_ac_LA = 0.75 * T_HB;
params.Tc_LA = 0.10 * T_HB;
params.t_ar_LA = params.t_ac_LA + params.Tc_LA;
params.Tr_LA = 0.80 * T_HB;

params.t_ac_RA = 0.75 * T_HB;
params.Tc_RA = 0.10 * T_HB;
params.t_ar_RA = params.t_ac_RA + params.Tc_RA;
params.Tr_RA = 0.80 * T_HB;

t_active_after_wrap = 1.05 * T_HB;
t_rest_after_tail = 1.70 * T_HB;

[~, ~, E_LA, E_RA] = elastance_model([t_active_after_wrap, t_rest_after_tail], params);

active_tol = 1e-8;
rest_tol = 1e-10;

if E_LA(1) > params.E.LA.EB + active_tol && E_RA(1) > params.E.RA.EB + active_tol
    fprintf('  [PASS] Atrial relaxation remains active after heartbeat wraparound.\n');
else
    error('Atrial activation dropped to baseline immediately after heartbeat wraparound.');
end

if abs(E_LA(2) - params.E.LA.EB) < rest_tol && abs(E_RA(2) - params.E.RA.EB) < rest_tol
    fprintf('  [PASS] Atrial activation returns to passive baseline after relaxation tail.\n');
else
    error('Atrial activation did not return to passive baseline after the relaxation tail.');
end

fprintf('\n  ALL ATRIAL WRAPAROUND TESTS PASSED\n');
