%% tune_baseline_sapmin_lvesp.m
% Constrained multi-objective tuning around the current baseline.
% Focus objective: SAP_min and LVESP.
% Hard constraints: preserve already-passed baseline gates.

clear; clc;

root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

base = default_parameters();

x0 = [base.Tc_LV_frac, base.Tr_LV_frac, base.C.SAR, base.Rvalve.open];
lb = [0.26, 0.06, 0.85, 0.0030];
ub = [0.38, 0.16, 1.30, 0.0060];

n_starts = 10;
starts = zeros(n_starts, numel(x0));
starts(1, :) = x0;
rng(42);
for i = 2:n_starts
    starts(i, :) = lb + rand(1, numel(x0)) .* (ub - lb);
end

if exist('fmincon', 'file') ~= 2
    error('tune_baseline:missingFmincon', 'fmincon is required for this tuning script.');
end

opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'off', ...
    'MaxIterations', 40, ...
    'MaxFunctionEvaluations', 500);

best = struct('x', [], 'f', inf, 'exitflag', -999, 'aux', []);

for si = 1:n_starts
    x_init = starts(si, :);
    problem.objective = @(x) objective_fn(x, base);
    problem.nonlcon   = @(x) constraints_fn(x, base);
    problem.x0 = x_init;
    problem.lb = lb;
    problem.ub = ub;
    problem.options = opts;
    problem.solver = 'fmincon';

    try
        [x_sol, f_sol, exitflag] = fmincon(problem);
        [c_sol, ~, aux_sol] = constraints_fn(x_sol, base);
        feasible = all(c_sol <= 0);

        fprintf('start %2d | f = %.6f | exitflag = %d | feasible = %d\n', ...
            si, f_sol, exitflag, feasible);

        if feasible && f_sol < best.f
            best.x = x_sol;
            best.f = f_sol;
            best.exitflag = exitflag;
            best.aux = aux_sol;
        end
    catch ME
        fprintf('start %2d | failed: %s\n', si, ME.message);
    end
end

if isempty(best.x)
    error('tune_baseline:noFeasibleSolution', 'No feasible tuned solution found in this run.');
end

fprintf('\n=== BEST FEASIBLE SOLUTION ===\n');
fprintf('Tc_frac       = %.6f\n', best.x(1));
fprintf('Tr_frac       = %.6f\n', best.x(2));
fprintf('C.SAR         = %.6f\n', best.x(3));
fprintf('C.SC (derived)= %.6f\n', max(1.9 - best.x(3), 0.2));
fprintf('Rvalve.open   = %.6f\n', best.x(4));
fprintf('Objective     = %.6f\n', best.f);

disp(best.aux.metrics);
fprintf('LVEDP = %.3f, LVESP = %.3f, SW = %.3f, Emax_derived = %.3f\n', ...
    best.aux.LVEDP, best.aux.LVESP, best.aux.SW_LV, best.aux.Emax_derived);

save(fullfile(project_root, 'results', 'tables', 'tune_baseline_sapmin_lvesp_best.mat'), 'best');

%% ------------------------------------------------------------------------
function f = objective_fn(x, base)
    [~, ~, aux] = constraints_fn(x, base);

    sap_min = aux.metrics.SAP_min;
    lvesp = aux.LVESP;

    % Multi-objective target with range penalties + soft center pull.
    sap_range_pen = max(0, 60 - sap_min).^2 + max(0, sap_min - 90).^2;
    esp_range_pen = max(0, 80 - lvesp).^2 + max(0, lvesp - 130).^2;

    sap_center_pen = ((sap_min - 70) / 10).^2;
    esp_center_pen = ((lvesp - 100) / 20).^2;

    f = 4.0 * sap_range_pen + 2.0 * esp_range_pen + 0.25 * sap_center_pen + 0.10 * esp_center_pen;
end

%% ------------------------------------------------------------------------
function [c, ceq, aux] = constraints_fn(x, base)
    p = apply_candidate(base, x);
    sim = integrate_system(p);
    m = compute_clinical_indices(sim, p);

    [LVEDP, LVESP, SW_LV, Emax_derived, BV_drift_abs, cycle_drift] = pv_and_conservation(sim, p);

    % Preserve already-passed gates as hard constraints.
    c = [
        100 - m.SAP_max
        m.SAP_max - 140
        70 - m.SAP_mean
        m.SAP_mean - 100
        0.95 - m.QpQs
        m.QpQs - 1.05
        0.55 - m.LVEF
        m.LVEF - 0.75
        10 - m.PAP_mean
        m.PAP_mean - 20
        0 - m.RAP_mean
        m.RAP_mean - 8
        LVEDP - 15
        3000 - SW_LV
        SW_LV - 12000
        0.60 * p.E.LV.EA - Emax_derived
        Emax_derived - 1.40 * p.E.LV.EA
        BV_drift_abs - 1.0
        cycle_drift - 0.1
    ];

    ceq = [];

    aux = struct();
    aux.metrics = m;
    aux.LVEDP = LVEDP;
    aux.LVESP = LVESP;
    aux.SW_LV = SW_LV;
    aux.Emax_derived = Emax_derived;
    aux.BV_drift_abs = BV_drift_abs;
    aux.cycle_drift = cycle_drift;
end

%% ------------------------------------------------------------------------
function p = apply_candidate(base, x)
    p = base;

    p.Tc_LV_frac = x(1);
    p.Tr_LV_frac = x(2);
    p.Tc_RV_frac = x(1);
    p.Tr_RV_frac = x(2);

    p.C.SAR = x(3);
    p.C.SC  = max(1.9 - x(3), 0.2);

    p.Rvalve.open = x(4);

    T = 60 / p.HR;
    p.Tc_LV   = p.Tc_LV_frac   * T;
    p.Tr_LV   = p.Tr_LV_frac   * T;
    p.Tc_RV   = p.Tc_RV_frac   * T;
    p.Tr_RV   = p.Tr_RV_frac   * T;
    p.t_ac_LA = p.t_ac_LA_frac * T;
    p.Tc_LA   = p.Tc_LA_frac   * T;
    p.t_ar_LA = p.t_ac_LA + p.Tc_LA;
    p.Tr_LA   = p.Tr_LA_frac   * T;
    p.t_ac_RA = p.t_ac_RA_frac * T;
    p.Tc_RA   = p.Tc_RA_frac   * T;
    p.t_ar_RA = p.t_ac_RA + p.Tc_RA;
    p.Tr_RA   = p.Tr_RA_frac   * T;
end

%% ------------------------------------------------------------------------
function [LVEDP, LVESP, SW_LV, Emax_derived, BV_drift_abs, cycle_drift] = pv_and_conservation(sim, p)
    sidx = p.idx;

    t = sim.t(:);
    T = 60 / p.HR;
    mask = t >= (t(end) - T);

    tc = t(mask);
    Vlv = sim.V(mask, sidx.V_LV);
    Vsar = sim.V(mask, sidx.V_SAR);

    Elv = zeros(size(tc));
    Plv = zeros(size(tc));
    for k = 1:numel(tc)
        [Elv(k), ~, ~, ~] = elastance_model(tc(k), p);
        Plv(k) = Elv(k) * (Vlv(k) - p.V0.LV);
    end

    Psar = (Vsar - p.V0.SAR) ./ p.C.SAR;
    dPav = Plv - Psar;
    wav  = 0.5 + 0.5 * tanh(dPav ./ p.epsilon_valve);
    Rav  = wav .* p.Rvalve.open + (1 - wav) .* p.Rvalve.closed;
    Qav  = dPav ./ Rav;

    [~, ied] = max(Vlv);
    LVEDP = Plv(ied);

    idx_eject = find(Qav > 1.0);
    if ~isempty(idx_eject)
        [~, k_es] = max(Elv(idx_eject));
        ies = idx_eject(k_es);
    else
        [~, ies] = min(Vlv);
    end

    LVESP = Plv(ies);
    LVESV_es = Vlv(ies);

    SW_LV = -trapz(Vlv, Plv);

    denom = max(LVESV_es - p.V0.LV, 1e-6);
    Emax_derived = LVESP / denom;

    vol_idx = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV sidx.V_SAR sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];
    BV_ic = sum(p.ic.V(vol_idx));
    BV_end = sum(sim.V(end, vol_idx));
    BV_drift_abs = abs(BV_end - BV_ic);

    cyc_start = find(mask, 1, 'first');
    cyc_end   = find(mask, 1, 'last');
    cycle_drift = abs(sum(sim.V(cyc_end, vol_idx)) - sum(sim.V(cyc_start, vol_idx)));
end
