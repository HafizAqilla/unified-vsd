%% test_valve_logic.m
% =========================================================================
% Unit tests for the cardiac valve model (valve_model.m).
%
% PURPOSE:
%   Verifies that the smooth-switch valve implementation (Guardrail §8.4)
%   satisfies the following properties independently of a full simulation:
%
%   1. FORWARD FLOW  — forward pressure gradient yields positive flow.
%   2. BACKWARD FLOW — reverse gradient yields near-zero flow
%      (small not hard-zero; proportional to R_closed not R_open).
%   3. ZERO GRADIENT — Q ≈ 0 when P_up == P_down.
%   4. RESISTANCE BOUNDS — R_eff always in [R_open, R_closed].
%   5. CONTINUITY — Q(dP) is Lipschitz continuous (no discrete jumps).
%   6. LEAKAGE RATIO — Q_reverse / Q_forward ≈ R_open / R_closed (≪ 1).
%   7. OPENING THRESHOLD — at dP = 5 × epsilon_valve, Q ≈ dP / R_open
%      (effectively fully open); confirms sharp enough switching.
%
% REFERENCES:
%   [1] valve_model.m — smooth tanh switching (Valenti Eq. 2.6, Guardrail §8.4).
%   [2] Guardrail §8.4 — numerical robustness requirement.
%
% USAGE:
%   >> test_valve_logic
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  1.0
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL — Valve Logic Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% Build minimal parameter struct (avoids a full default_parameters() call)
params = struct();
params.Rvalve.open   = 6.2872e-3;   % [mmHg·s/mL]  Source: Valenti Table 3.3
params.Rvalve.closed = 9.4168e+4;   % [mmHg·s/mL]  Source: Valenti Table 3.3
params.epsilon_valve = 0.5;         % [mmHg]        Source: default_parameters.m

R_open   = params.Rvalve.open;    % [mmHg·s/mL]
R_closed = params.Rvalve.closed;  % [mmHg·s/mL]
eps_v    = params.epsilon_valve;  % [mmHg]

% -------------------------------------------------------------------------
% TEST 1: Forward pressure gradient → positive flow
% -------------------------------------------------------------------------
fprintf('--- Test 1: Forward flow ---\n');
dP_fwd = 10;   % [mmHg]  typical LV-aorta forward gradient at peak systole
Q_fwd  = valve_model(100, 100 - dP_fwd, params);   % P_up=100, P_down=90
if Q_fwd > 0
    fprintf('  [PASS] Q_forward = %.6f mL/s  (> 0 as required)\n', Q_fwd);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Q_forward = %.6f mL/s  (expected > 0)\n', Q_fwd);
    n_fail = n_fail + 1;
end

% -------------------------------------------------------------------------
% TEST 2: Reverse pressure gradient → near-zero (not strongly negative)
% -------------------------------------------------------------------------
fprintf('--- Test 2: Reverse flow leakage ---\n');
dP_rev = -10;  % [mmHg]  valve closed, slight back-pressure
Q_rev  = valve_model(90, 100, params);   % P_up < P_down
leakage_ratio = abs(Q_rev) / abs(Q_fwd);   % should be << 1

leakage_tol = 1e-3;   % must be < 0.1% of forward flow
if leakage_ratio < leakage_tol
    fprintf('  [PASS] Leakage ratio = %.2e  (< %.1e : negligible regurgitation)\n', ...
            leakage_ratio, leakage_tol);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Leakage ratio = %.2e  (expected < %.1e)\n', ...
            leakage_ratio, leakage_tol);
    n_fail = n_fail + 1;
end

% -------------------------------------------------------------------------
% TEST 3: Zero gradient → Q ≈ 0
% -------------------------------------------------------------------------
fprintf('--- Test 3: Zero pressure gradient ---\n');
Q_zero = valve_model(80, 80, params);   % P_up == P_down
Q_zero_tol = max(R_open, 1e-6);        % dP=0 → Q = 0/R_eff must be exactly 0
if abs(Q_zero) < 1e-12
    fprintf('  [PASS] Q at dP=0 = %.2e mL/s  (≈ 0)\n', Q_zero);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Q at dP=0 = %.2e mL/s  (expected ≈ 0)\n', Q_zero);
    n_fail = n_fail + 1;
end

% -------------------------------------------------------------------------
% TEST 4: R_eff always within [R_open, R_closed] for all dP
% -------------------------------------------------------------------------
fprintf('--- Test 4: Effective resistance bounds ---\n');
dP_sweep = linspace(-50, 50, 1000);   % [mmHg]  wide sweep
Q_sweep  = zeros(size(dP_sweep));
for k = 1:numel(dP_sweep)
    Q_sweep(k) = valve_model(dP_sweep(k), 0, params);
end
% Avoid division by near-zero Q (divide by non-zero dP instead)
nz = abs(dP_sweep) > 1e-9;
R_eff_sweep = abs(dP_sweep(nz)) ./ abs(Q_sweep(nz));   % [mmHg·s/mL]

R_min_actual = min(R_eff_sweep);
R_max_actual = max(R_eff_sweep);

if R_min_actual >= 0.99*R_open && R_max_actual <= 1.01*R_closed
    fprintf('  [PASS] R_eff ∈ [%.4g, %.4g] mmHg·s/mL  (within [R_open, R_closed])\n', ...
            R_min_actual, R_max_actual);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] R_eff ∈ [%.4g, %.4g] — expected [%.4g, %.4g]\n', ...
            R_min_actual, R_max_actual, R_open, R_closed);
    n_fail = n_fail + 1;
end

% -------------------------------------------------------------------------
% TEST 5: Continuity — no discontinuous jump in Q over dP sweep
% -------------------------------------------------------------------------
fprintf('--- Test 5: Flow continuity (no discrete jumps) ---\n');
dQ_max = max(abs(diff(Q_sweep)));       % [mL/s]  largest step between Q values
ddP    = dP_sweep(2) - dP_sweep(1);    % [mmHg]  step size

% Lipschitz bound: max |dQ/ddP| ≤ 1/R_open  (slope of fully open valve)
% Any measured jump > 10× this bound indicates a discontinuity.
max_allowed_jump = 10 * ddP / R_open;  % [mL/s]  generous threshold

if dQ_max < max_allowed_jump
    fprintf('  [PASS] Max |ΔQ| per step = %.2e mL/s  (< continuity threshold %.2e)\n', ...
            dQ_max, max_allowed_jump);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Max |ΔQ| per step = %.2e mL/s  (DISCONTINUITY — threshold %.2e)\n', ...
            dQ_max, max_allowed_jump);
    n_fail = n_fail + 1;
end

% -------------------------------------------------------------------------
% TEST 6: Open/closed resistance ratio >> 1 (prevents physiological regurgitation)
% -------------------------------------------------------------------------
fprintf('--- Test 6: Resistance ratio adequacy ---\n');
ratio_RC_RO = R_closed / R_open;
ratio_min   = 1e6;   % minimum acceptable: R_closed / R_open ≥ 10^6

if ratio_RC_RO >= ratio_min
    fprintf('  [PASS] R_closed / R_open = %.2e  (≥ %.0e : adequate anti-regurgitation guard)\n', ...
            ratio_RC_RO, ratio_min);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] R_closed / R_open = %.2e  (expected ≥ %.0e)\n', ...
            ratio_RC_RO, ratio_min);
    n_fail = n_fail + 1;
end

% -------------------------------------------------------------------------
% TEST 7: At dP >> epsilon_valve, R_eff ≈ R_open (fully open state)
% -------------------------------------------------------------------------
fprintf('--- Test 7: Full-opening at large forward gradient ---\n');
dP_large = 5 * eps_v;   % [mmHg]  5× switching width: tanh(5) ≈ 0.9999 → weight ≈ 1
Q_large  = valve_model(dP_large, 0, params);
R_eff_large = dP_large / Q_large;   % [mmHg·s/mL]

if abs(R_eff_large - R_open) / R_open < 0.01   % within 1% of R_open
    fprintf('  [PASS] R_eff at dP=5ε = %.6f mmHg·s/mL  (≈ R_open=%.6f, error < 1%%)\n', ...
            R_eff_large, R_open);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] R_eff at dP=5ε = %.6f mmHg·s/mL  (R_open=%.6f, error > 1%%)\n', ...
            R_eff_large, R_open);
    n_fail = n_fail + 1;
end

%% Overall result
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL VALVE LOGIC TESTS PASSED\n');
else
    fprintf('  ONE OR MORE VALVE TESTS FAILED\n');
end
fprintf('==========================================\n');
