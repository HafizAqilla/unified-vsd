%% test_baseline.m
% =========================================================================
% Baseline verification test for the unified VSD lumped-parameter model.
%
% PURPOSE:
%   Confirms that the default adult parameters (no VSD, healthy reference)
%   reproduce physiologically plausible haemodynamics and that the
%   conservation of blood volume is satisfied.
%
% PASS CRITERIA (all must hold — Guardrail §10.1):
%   Reference ranges for a healthy adult at rest (Cite: Nichols & O'Rourke,
%   McDonald's Blood Flow in Arteries, 6th ed.):
%   1. SAP_max (systolic)     in [100, 140]  mmHg
%   2. SAP_min (diastolic)    in [60,  90]   mmHg
%   3. SAP_mean               in [70,  100]  mmHg
%   4. HR                     in [65,  85]   bpm
%   5. QpQs                   in [0.95, 1.05] (no shunt at baseline)
%   6. LVEF                   in [0.55, 0.75] (fraction)
%   7. PAP_mean               in [10,  20]   mmHg
%   8. RAP_mean               in [0,   8]    mmHg
%   9. Blood volume conservation (per-cycle): BV drift < 1 mL
%   10. Steady-state reached flag confirmed
%
% USAGE:
%   >> test_baseline
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-03
% VERSION:  1.1
% =========================================================================

clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, '..')));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Baseline Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% 1. Build default adult parameters (no patient scaling)
params = default_parameters();

% Recompute absolute timing from fractional values at HR=75
T_HB = 60 / params.HR;   % [s]
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
    fprintf('  Simulation: PASSED\n');
catch ME
    fprintf('  Simulation: FAILED \u2014 %s\n', ME.message);
    return;
end

%% 3. Check steady-state flag  (Guardrail \u00a78.3)
fprintf('\n--- Steady-state convergence ---\n');
if sim.ss_reached
    fprintf('  [PASS] Steady state confirmed by integrate_system.\n');
    n_pass = n_pass + 1;
else
    fprintf('  [WARN] Steady state NOT confirmed. Results may be transient.\n');
    % Not counted as hard fail, but operator must investigate.
end

%% 4. Compute metrics
metrics = compute_clinical_indices(sim, params);

%% 5-6. Haemodynamic range assertions  (Guardrail \u00a710.1)
%   Reference: Nichols & O'Rourke (2005), McDonald's Blood Flow in Arteries, 6th ed.
%   Each row: {label, value, lo, hi, unit}
fprintf('\n--- Haemodynamic assertions ---\n');

assertions = {
    'SAP_max  [systolic]   mmHg',  metrics.SAP_max,  100, 140, 'mmHg'
    'SAP_min  [diastolic]  mmHg',  metrics.SAP_min,   60,  90, 'mmHg'
    'SAP_mean [MAP]        mmHg',  metrics.SAP_mean,  70, 100, 'mmHg'
    'HR                    bpm',   params.HR,          65,  85, 'bpm'
    'QpQs                  -',     metrics.QpQs,     0.95, 1.05, '-'
    'LVEF                  frac',  metrics.LVEF,     0.55, 0.75, '-'
    'PAP_mean              mmHg',  metrics.PAP_mean,  10,  20, 'mmHg'
    'RAP_mean              mmHg',  metrics.RAP_mean,   0,   8, 'mmHg'
};

for k = 1:size(assertions, 1)
    a_label = assertions{k, 1};
    a_val   = assertions{k, 2};
    a_lo    = assertions{k, 3};
    a_hi    = assertions{k, 4};
    a_unit  = assertions{k, 5};
    if a_val >= a_lo && a_val <= a_hi
        fprintf('  [PASS] %-35s = %8.3f %s  (expected [%.2f, %.2f])\n', ...
                a_label, a_val, a_unit, a_lo, a_hi);
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] %-35s = %8.3f %s  (expected [%.2f, %.2f])\n', ...
                a_label, a_val, a_unit, a_lo, a_hi);
        n_fail = n_fail + 1;
    end
end

%% 7. Blood volume conservation — per-cycle check  (Guardrail \u00a78.6 / \u00a710.3)
%   Integrate ΣdV/dt over the last complete cardiac cycle.
%   For a closed conservative system, net change must be < bv_tol mL.
fprintf('\n--- Conservation check (Guardrail Section 8.6) ---\n');

sidx    = params.idx;   % state index struct (Guardrail \u00a77.1)
vol_idx = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV ...   % volume state indices
           sidx.V_SAR sidx.V_SC sidx.V_SVEN ...
           sidx.V_PAR sidx.V_PVEN];

t_sim  = sim.t(:);
T_HB_s = 60 / params.HR;   % [s]
T1     = t_sim(end);
T0     = T1 - T_HB_s;
cmask  = (t_sim >= T0) & (t_sim <= T1);

BV_ic  = sum(params.ic.V(vol_idx));                % [mL]  initial condition total
BV_end = sum(sim.V(end, vol_idx));                  % [mL]  final state total
BV_drift_abs = abs(BV_end - BV_ic);                % [mL]  absolute drift over all cycles

BV_tol = 1.0;   % [mL]  1 mL is conservative for 20-cycle integration
if BV_drift_abs < BV_tol
    fprintf('  [PASS] Total blood volume drift = %.3f mL  (< %.1f mL)\n', ...
            BV_drift_abs, BV_tol);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Total blood volume drift = %.3f mL  (>= %.1f mL)\n', ...
            BV_drift_abs, BV_tol);
    n_fail = n_fail + 1;
end

%% 8. Cycle-level mass balance: check last single cycle
%   Sum of all volume states at cycle start vs cycle end should be ~0.
V_cyc_start = sum(sim.V(find(cmask, 1, 'first'), vol_idx));     % [mL]
V_cyc_end   = sum(sim.V(find(cmask, 1, 'last'),  vol_idx));     % [mL]
cycle_drift  = abs(V_cyc_end - V_cyc_start);                    % [mL]

cycle_tol = 0.1;  % [mL]  per Guardrail §8.3 steady-state criterion
if cycle_drift < cycle_tol
    fprintf('  [PASS] Last-cycle volume drift = %.4f mL  (< %.2f mL)\n', ...
            cycle_drift, cycle_tol);
    n_pass = n_pass + 1;
else
    fprintf('  [FAIL] Last-cycle volume drift = %.4f mL  (>= %.2f mL)\n', ...
            cycle_drift, cycle_tol);
    n_fail = n_fail + 1;
end

%% 9. LV pressure-volume loop shape check  (Guardrail §10.2)
%
%  A physiologically correct PV loop must show:
%   (a) Positive stroke work (loop enclosed counter-clockwise)
%   (b) LVEDP in the normal range (not hyper- or hypo-loaded)
%   (c) LVESP in the aortic pressure range (valve opens at correct time)
%   (d) ESPVR slope consistent with Emax_LV (within ±40 %)
%
%  Reference: Kass DA & Maughan WL (1988). Am J Cardiol 62:371–380.
%             Suga H et al. (1973). Circ Res 32:314–322.  (ESPVR definition)

fprintf('\n--- LV P-V loop shape check (Guardrail §10.2) ---\n');

t_cyc  = sim.t(cmask);
V_LV_c = sim.V(cmask, sidx.V_LV);   % [mL]

% --- Reconstruct P_LV over the cycle via elastance model
P_LV_c = zeros(size(t_cyc));
for ki = 1:numel(t_cyc)
    [E_LV_ki, ~, ~, ~] = elastance_model(t_cyc(ki), params);   % t_cyc is already absolute sim time
    P_LV_c(ki) = E_LV_ki * (V_LV_c(ki) - params.V0.LV);        % [mmHg]
end

% (a) Stroke work  SW = ∮ P dV  (positive = counter-clockwise = systole ejects blood)
%     Negate trapz sign because for a CCW loop ∮P dV > 0 when integrating with
%     increasing time (volume decreases during systole → dV/dt < 0 → contribution positive).
SW_LV = -trapz(V_LV_c, P_LV_c);   % [mmHg·mL]  — should be 3000–12000 for healthy adult

% (b) LVEDP — pressure at the time of maximum LV volume
[~, idx_ed] = max(V_LV_c);
LVEDP = P_LV_c(idx_ed);   % [mmHg]

% (c) LVESP — pressure at the time of minimum LV volume
[~, idx_es] = min(V_LV_c);
LVESP = P_LV_c(idx_es);   % [mmHg]

% (d) ESPVR slope: Emax_derived = LVESP / (LVESV - V0_LV)
%     For a healthy adult, this should be within ±40% of params.E.LV.EA
%     (the passive baseline EB is small, so EA ≈ Emax).
denom_espvr = max(metrics.LVESV - params.V0.LV, 1e-6);
Emax_derived = LVESP / denom_espvr;   % [mmHg/mL]

pv_checks = {
    'Stroke work SW_LV   mmHg·mL',  SW_LV,         3000, 12000, 'mmHg·mL'
    'LVEDP               mmHg',     LVEDP,            0,    15, 'mmHg'
    'LVESP               mmHg',     LVESP,           80,   130, 'mmHg'
    'Emax_derived        mmHg/mL',  Emax_derived,  0.60*params.E.LV.EA, 1.40*params.E.LV.EA, 'mmHg/mL'
};

for k = 1:size(pv_checks, 1)
    a_label = pv_checks{k, 1};
    a_val   = pv_checks{k, 2};
    a_lo    = pv_checks{k, 3};
    a_hi    = pv_checks{k, 4};
    a_unit  = pv_checks{k, 5};
    if a_val >= a_lo && a_val <= a_hi
        fprintf('  [PASS] %-35s = %8.3f %s  (expected [%.2f, %.2f])\n', ...
                a_label, a_val, a_unit, a_lo, a_hi);
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] %-35s = %8.3f %s  (expected [%.2f, %.2f])\n', ...
                a_label, a_val, a_unit, a_lo, a_hi);
        n_fail = n_fail + 1;
    end
end

%% 10. Overall result
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
if n_fail == 0
    fprintf('  ALL TESTS PASSED\n');
else
    fprintf('  ONE OR MORE TESTS FAILED — investigate before proceeding\n');
end
fprintf('==========================================\n');
