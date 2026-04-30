%% test_validity_gates.m
% Smoke test for the hard validity gate helper.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

params = default_parameters();
params = recompute_test_timing(params);

sim = struct();
sim.t = linspace(0, 0.8, 8)';
sim.V = repmat(params.ic.V, numel(sim.t), 1);
sim.V(:, params.idx.V_LV) = -1;
sim.ss_reached = false;

metrics = struct('PVR', 50, 'SVR', 0.01, 'LVEF', 1.2, 'RVEF', 0.0, ...
    'QpQs', 1.5, 'Qs_mean_mLs', 0);
validity = evaluate_simulation_validity(sim, params, metrics, 'pre_surgery');

assert(~validity.is_valid, 'Synthetic invalid simulation should fail validity gates.');
assert(validity.penalty > 0, 'Invalid simulation should incur positive penalty.');
assert(any(strcmp(validity.failed_flags, 'svr_out_of_bounds')), ...
    'Synthetic low SVR should trip the SVR bound.');

clinical = patient_profile_A();
sim_ok = sim;
sim_ok.V(:, params.idx.V_LV) = params.ic.V(params.idx.V_LV);
sim_ok.ss_reached = true;
metrics_ok = struct( ...
    'PVR', clinical.pre_surgery.PVR_WU, ...
    'SVR', clinical.pre_surgery.SVR_WU, ...
    'LVEF', clinical.pre_surgery.LVEF, ...
    'RVEF', 0.55, ...
    'QpQs', clinical.pre_surgery.QpQs, ...
    'Qs_mean_mLs', 5);
validity_ok = evaluate_simulation_validity(sim_ok, params, metrics_ok, 'pre_surgery', clinical);
assert(~validity_ok.flags.svr_out_of_bounds, ...
    'Clinical high-SVR pediatric target should not be rejected as invalid.');
assert(validity_ok.metric_bounds.SVR(2) >= clinical.pre_surgery.SVR_WU, ...
    'Scenario-aware SVR upper bound should cover the clinical target.');

fprintf('[PASS] Validity gates reject non-physiological simulation states.\n');

function params = recompute_test_timing(params)
T_HB = 60 / params.HR;
params.Tc_LV   = params.Tc_LV_frac   * T_HB;
params.Tr_LV   = params.Tr_LV_frac   * T_HB;
params.Tc_RV   = params.Tc_RV_frac   * T_HB;
params.Tr_RV   = params.Tr_RV_frac   * T_HB;
params.t_ac_LA = params.t_ac_LA_frac * T_HB;
params.Tc_LA   = params.Tc_LA_frac   * T_HB;
params.t_ar_LA = params.t_ac_LA + params.Tc_LA;
params.Tr_LA   = params.Tr_LA_frac   * T_HB;
params.t_ac_RA = params.t_ac_RA_frac * T_HB;
params.Tc_RA   = params.Tc_RA_frac   * T_HB;
params.t_ar_RA = params.t_ac_RA + params.Tc_RA;
params.Tr_RA   = params.Tr_RA_frac   * T_HB;
end
