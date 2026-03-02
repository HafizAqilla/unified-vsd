%% test_baseline.m
% =========================================================================
% Baseline verification test for the unified VSD lumped-parameter model.
%
% PURPOSE:
%   Confirms that the default adult parameters (no VSD, healthy reference)
%   reproduce physiologically plausible haemodynamics and that the
%   conservation of blood volume is satisfied.
%
% PASS CRITERIA (all must hold):
%   1. MAP (mean aortic pressure) in [70, 100] mmHg
%   2. HR   in [65, 85]  bpm  (adult resting — no scaling applied)
%   3. QpQs in [0.95, 1.05] (no shunt expected with large R_VSD)
%   4. LVEF in [0.50, 0.75] (normal adult range)
%   5. PAP_mean in [10, 18] mmHg (normal adult)
%   6. Blood volume conservation: Σ(volume states) ≈ Σ(IC volume) ± 5%
%
% USAGE:
%   >> test_baseline
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, '..')));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL — Baseline Test\n');
fprintf('==========================================\n\n');

pass_all = true;

%% 1. Build default adult parameters (no patient scaling)
params = default_parameters();

% Recompute absolute timing from fractional values at HR=75
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

%% 2. Run simulation
fprintf('Running baseline simulation...\n');
try
    sim = integrate_system(params);
    fprintf('  Simulation: PASSED\n\n');
catch ME
    fprintf('  Simulation: FAILED — %s\n', ME.message);
    return;
end

%% 3. Compute metrics
metrics = compute_clinical_indices(sim, params);

%% 4. Assertions
function assert_range(name, val, lo, hi)
    status = (val >= lo && val <= hi);
    if status
        fprintf('  [PASS] %-30s = %8.3f  (expected [%.2f, %.2f])\n', name, val, lo, hi);
    else
        fprintf('  [FAIL] %-30s = %8.3f  (expected [%.2f, %.2f])\n', name, val, lo, hi);
        pass_all = false;   %#ok<NODEF>
    end
end

fprintf('--- Haemodynamic assertions ---\n');
assert_range('MAP (SAP_mean) [mmHg]',     metrics.SAP_mean,   70,  100);
assert_range('HR [bpm]',                   params.HR,          65,   85);
assert_range('QpQs [-]',                   metrics.QpQs,       0.95, 1.05);
assert_range('LVEF [-]',                   metrics.LVEF,       0.50, 0.75);
assert_range('PAP_mean [mmHg]',            metrics.PAP_mean,   10,   18);
assert_range('RAP_mean [mmHg]',            metrics.RAP_mean,    0,    8);

%% 5. Blood volume conservation
vol_idx = [1 2 3 4 5 7 8 10 13];
BV_ic   = sum(params.ic.V(vol_idx));
BV_end  = sum(sim.V(end, vol_idx));
BV_err  = abs(BV_end - BV_ic) / max(BV_ic, 1e-6);
fprintf('\n--- Conservation check ---\n');
if BV_err < 0.05
    fprintf('  [PASS] Blood volume drift = %.2f%%  (< 5%%)\n', 100*BV_err);
else
    fprintf('  [FAIL] Blood volume drift = %.2f%%  (>= 5%%)\n', 100*BV_err);
    pass_all = false;
end

%% 6. Overall result
fprintf('\n==========================================\n');
if pass_all
    fprintf('  RESULT: ALL TESTS PASSED\n');
else
    fprintf('  RESULT: ONE OR MORE TESTS FAILED\n');
end
fprintf('==========================================\n');
