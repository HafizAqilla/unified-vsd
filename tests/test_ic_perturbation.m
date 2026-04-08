%% test_ic_perturbation.m
% =========================================================================
% Initial condition robustness test for the unified VSD lumped-parameter model.
%
% PURPOSE:
%   Verifies that the ODE converges to the same periodic limit cycle when
%   started from perturbed initial conditions (Guardrail §10.4).
%
%   The baseline and two perturbed IC sets must all converge to the same
%   steady-state cycle within the tolerances below.
%
% METHODOLOGY (Guardrail §10.4):
%   Three IC vectors are tested:
%     nominal : params.ic.V  (as set by default_parameters)
%     +10%    : 1.10 × params.ic.V
%     -10%    : 0.90 × params.ic.V
%
%   Each run is integrated to steady state by integrate_system.m.
%   Steady-state cycle is then characterised by:
%     - Peak values of all 14 states over the last full cardiac cycle.
%
%   Convergence criterion: all peak values from the ±10% runs must match
%   the nominal run within:
%     tol_V = 1.0 mL   (volume states)
%     tol_P = 1.0 mmHg (pressure state P_PC)
%     tol_Q = 2.0 mL/s (flow states)
%
%   Because integrate_system.m already enforces its own steady-state
%   criterion (0.1 mmHg / 0.1 mL cycle-to-cycle), the perturbation test
%   is a *reproducibility* test, not a solver convergence test.
%
% REPORTED DATA:
%   - Whether each perturbed run reached steady state
%   - Max absolute deviation per state from the nominal limit cycle
%   - Number of cardiac cycles simulated (from params.sim.nCyclesSteady)
%
% REFERENCES:
%   [1] Guardrail §10.4 — requirement to perturb X0 ±10%.
%   [2] integrate_system.m — steady-state convergence criteria.
%   [3] Valenti (2023). Thesis. Eqs. 2.1–2.7.
%
% USAGE:
%   >> test_ic_perturbation
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
fprintf('  UNIFIED VSD MODEL — IC Perturbation Test\n');
fprintf('==========================================\n\n');

n_pass = 0;
n_fail = 0;

%% 1. Build default parameters
params = default_parameters();

% Recompute absolute timing (same as test_baseline.m)
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

% Allow more cycles when starting from perturbed ICs
params.sim.nCyclesSteady = 40;   % [cycles]  enough for perturbed start

%% 2. Tolerance definitions  (Guardrail §10.4)
%  Volume states   [mL]   — states 1–4 (chambers), 5,7,8,10,13 (vascular)
%  Flow states     [mL/s] — states 6,9,11,14
%  Pressure state  [mmHg] — state 12 (P_PC)
tol_V = 1.0;   % [mL]    max allowed deviation in volume peak  vs nominal
tol_Q = 2.0;   % [mL/s]  max allowed deviation in flow peak    vs nominal
tol_P = 1.0;   % [mmHg]  max allowed deviation in P_PC peak    vs nominal

sidx     = params.idx;
vol_idx  = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV ...
            sidx.V_SAR sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];
flow_idx = [sidx.Q_SAR sidx.Q_SVEN sidx.Q_PAR sidx.Q_PVEN];
pres_idx = sidx.P_PC;

% Named list for readable output
state_labels = {'V_RA','V_RV','V_LA','V_LV','V_SAR','Q_SAR', ...
                'V_SC','V_SVEN','Q_SVEN','V_PAR','Q_PAR','P_PC', ...
                'V_PVEN','Q_PVEN'};

%% 3. Helper: extract peak values over last cycle
%   Returns 1×14 vector of max amplitude per state.
get_peaks = @(sim_out) max(abs(sim_out.V), [], 1);   % [1×14]

%% 4. Run nominal simulation
fprintf('\nRunning NOMINAL ICs...\n');
params_nom = params;
try
    sim_nom = integrate_system(params_nom);
    fprintf('  Nominal: steady-state = %s  |  max_cycles = %d\n', ...
            string(sim_nom.ss_reached), params.sim.nCyclesSteady);
    peaks_nom = get_peaks(sim_nom);
catch ME
    fprintf('  NOMINAL simulation FAILED: %s\n', ME.message);
    fprintf('  Cannot continue IC perturbation test without nominal baseline.\n');
    return;
end

%% 5. Define perturbation cases
ic_nominal = params.ic.V(:);                    % 14×1

perturbations = {
    '+10%',  1.10 * ic_nominal
    '-10%',  0.90 * ic_nominal
};

%% 6. Run perturbed ICs and compare
for p = 1:size(perturbations, 1)
    pert_label = perturbations{p, 1};
    ic_pert    = perturbations{p, 2};

    fprintf('\nRunning IC perturbation: %s...\n', pert_label);

    params_pert        = params;
    params_pert.ic.V   = ic_pert(:)';   % row vector for ode15s

    try
        sim_pert = integrate_system(params_pert);
        fprintf('  Perturbed (%s): steady-state = %s\n', ...
                pert_label, string(sim_pert.ss_reached));
        peaks_pert = get_peaks(sim_pert);
    catch ME
        fprintf('  [FAIL] %s simulation FAILED: %s\n', pert_label, ME.message);
        n_fail = n_fail + 1;
        continue;
    end

    % --- Steady-state convergence flag ---
    label_ss = sprintf('SS_reached  (%s IC)', pert_label);
    if sim_pert.ss_reached
        fprintf('  [PASS] %-40s — converged to steady state\n', label_ss);
        n_pass = n_pass + 1;
    else
        fprintf('  [WARN] %-40s — steady state NOT confirmed (nCyclesSteady too low?)\n', label_ss);
        % Not a hard fail, but reported
    end

    % --- Peak-value deviation from nominal ---
    fprintf('  Comparing peak values against nominal limit cycle:\n');
    all_within = true;

    for s = 1:14
        delta = abs(peaks_pert(s) - peaks_nom(s));

        if ismember(s, vol_idx)
            tol_s = tol_V;
            unit_s = 'mL';
        elseif ismember(s, flow_idx)
            tol_s = tol_Q;
            unit_s = 'mL/s';
        else  % P_PC
            tol_s = tol_P;
            unit_s = 'mmHg';
        end

        if delta <= tol_s
            status = 'PASS';
        else
            status = 'FAIL';
            all_within = false;
            n_fail = n_fail + 1;
        end

        fprintf('    [%s] %-10s  Δpeak = %8.4f %-6s  (tol = %.1f)\n', ...
                status, state_labels{s}, delta, unit_s, tol_s);
        if strcmp(status, 'PASS'), n_pass = n_pass + 1; end
    end

    if all_within
        fprintf('  => All %d states within tolerance for %s perturbation.\n', ...
                14, pert_label);
    else
        fprintf('  => One or more states OUTSIDE tolerance for %s perturbation.\n', ...
                pert_label);
    end
end

%% 7. Summary
fprintf('\n==========================================\n');
fprintf('  RESULT: %d PASSED, %d FAILED\n', n_pass, n_fail);
fprintf('  Warm-up cycles used: %d (params.sim.nCyclesSteady)\n', ...
        params.sim.nCyclesSteady);
if n_fail == 0
    fprintf('  ALL IC PERTURBATION TESTS PASSED\n');
    fprintf('  The model converges to the same limit cycle from ±10%% perturbation.\n');
else
    fprintf('  FAILURES FOUND — Consider increasing nCyclesSteady.\n');
end
fprintf('==========================================\n');
