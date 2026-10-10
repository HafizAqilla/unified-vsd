function T = verify_solver_convergence(run_dir)
% VERIFY_SOLVER_CONVERGENCE
% -----------------------------------------------------------------------
% Numerical verification of the reported states (pre-closure calibrated
% and post-closure prediction): (1) ODE tolerance and maximum step,
% (2) number of cycles after the steady-state flag, (3) whether the
% -5 mmHg pressure floor is active, (4) volume conservation over a cycle.
% Writes <run_dir>/supplement/solver_verification.csv.
% -----------------------------------------------------------------------
if nargin < 1
    run_dir = 'results/runs/20260920_005419_reyna_pre_surgery';
end
S = load(fullfile(run_dir, 'mat', 'params_accepted_candidate_pre_surgery.mat'));
f = fieldnames(S); P = [];
for i = 1:numel(f)
    v = S.(f{i});
    if isstruct(v) && isfield(v, 'R') && isfield(v, 'vsd'), P = v; break; end
    if isstruct(v) && isfield(v, 'params') && isstruct(v.params), P = v.params; break; end
end
assert(~isempty(P), 'verify_solver_convergence:noParams', 'No parameter struct found.');
clinical = patient_reyna();
Ppost = P; Ppost.R.vsd = 1e6; Ppost.vsd.area_mm2 = 0;
Ppost.HR = clinical_heart_rate(clinical, 'post_surgery');
Ppost = recompute_cycle_timing(Ppost);

states = {'pre_calibrated', P; 'post_prediction', Ppost};
settings = { ...
    'reference',        1e-5, 1e-6, 20, 0; ...
    'tol x10 tighter',  1e-6, 1e-7, 20, 0; ...
    'tol x100 tighter', 1e-7, 1e-8, 20, 0; ...
    'MaxStep T/50',     1e-5, 1e-6, 50, 0; ...
    '+20 cycles',       1e-5, 1e-6, 20, 20};
metrics = {'SAP_max','SAP_min','SAP_mean','PAP_max','PAP_min','PAP_mean', ...
    'RAP_mean','LAP_mean','CO_Lmin','QpQs','LVEDV','LVESV','RVEDV','RVESV','LVEF','RVEF'};
rows = {};
for s = 1:size(states, 1)
    ref = [];
    for k = 1:size(settings, 1)
        p = states{s, 2};
        p.sim.rtol = settings{k, 2}; p.sim.atol = settings{k, 3};
        p.sim.maxstep_div = settings{k, 4};
        p.sim.nCyclesSteady = 80 + settings{k, 5};
        if settings{k, 5} > 0
            p.sim.ss_tol_V = 0; p.sim.ss_tol_P = 0;   % force the full cycle count
        end
        sim = integrate_system(p);
        m = compute_clinical_indices(sim, p);
        y = cellfun(@(n) m.(n), metrics);
        % pressure floor and conservation over the last cycle
        X = sim.V; idx = p.idx;
        PRA = local_pressure(sim, p, 'RA'); PLA = local_pressure(sim, p, 'LA');
        PSV = (X(:, idx.V_SVEN) - p.V0.SVEN) / p.C.SVEN;
        PPV = (X(:, idx.V_PVEN) - p.V0.PVEN) / p.C.PVEN;
        floor_active = any([PRA; PLA; PSV; PPV] < -5);
        vol_idx = [idx.V_RA idx.V_RV idx.V_LA idx.V_LV idx.V_SAR idx.V_SC idx.V_SVEN idx.V_PAR idx.V_PVEN];
        Vtot = sum(X(:, vol_idx), 2) + (p.C.PCOX + p.C.PCNO) * X(:, idx.P_PC);
        drift = max(Vtot) - min(Vtot);
        if k == 1, ref = y; end
        rel = max(abs(y - ref) ./ max(abs(ref), 1e-9)) * 100;
        rows(end + 1, :) = {states{s, 1}, settings{k, 1}, logical(sim.ss_reached), rel, ...
            min([PRA; PLA]), min([PSV; PPV]), floor_active, drift}; %#ok<AGROW>
        fprintf('%-16s %-17s ss=%d maxRelChange=%.4f%% minPatrial=%.2f minPven=%.2f floor=%d Vdrift=%.4f mL\n', rows{end, :});
    end
end
T = cell2table(rows, 'VariableNames', {'State','Setting','SteadyState', ...
    'MaxRelChange_pct','MinAtrialP_mmHg','MinVenousP_mmHg','FloorActive','TotalVolumeRange_mL'});
writetable(T, fullfile(run_dir, 'supplement', 'solver_verification.csv'));
end

function P = local_pressure(sim, p, ch)
[E_LV, E_RV, E_LA, E_RA] = elastance_model(sim.t, p); %#ok<ASGLU>
switch ch
    case 'RA', P = E_RA .* (sim.V(:, p.idx.V_RA) - p.V0.RA);
    case 'LA', P = E_LA .* (sim.V(:, p.idx.V_LA) - p.V0.LA);
end
end
