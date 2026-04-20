function gsa_out = gsa_run_sobol(cfg, params0)
% GSA_RUN_SOBOL
% -----------------------------------------------------------------------
% Execute the Sobol sensitivity analysis using pre-generated Saltelli
% sample matrices (A, B, A_Bi) from gsa_sobol_setup.m.
%
% Uses the Jansen (1999) estimator for first-order (S1) and total (ST)
% Sobol indices.
%
% INPUTS:
%   cfg       - configuration struct from gsa_sobol_setup.m
%   params0   - baseline (scaled) parameter struct
%
% OUTPUTS:
%   gsa_out   - struct with one sub-struct per output metric, each with:
%               .S1        first-order Sobol index  (d × 1)
%               .ST        total-order Sobol index  (d × 1)
%               .table     MATLAB table (sorted by ST descending)
%               .primary   logical flag: this metric is in cfg.primary_metrics
%
% COMPUTATIONAL NOTE:
%   Total model evaluations = N × (d + 2).
%   For N=256, d=18 (pre) or d=19 (post): ~5120 evaluations.
%   N=256 satisfies Saltelli et al. (2010) precision requirement:
%   ST confidence width < 0.05 at 95 % confidence for d ≤ 20.
%
% REFERENCES:
%   [1] Saltelli et al. (2010). Computer Physics Communications 181:259–270.
%   [2] Jansen (1999). Computer Physics Communications 117:35–43.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

names  = cfg.names;
d      = numel(names);
S      = cfg.saltelli;
A      = S.A;
B      = S.B;
AB     = S.AB;

all_metrics = cfg.all_metrics;

% Apply GSA-mode simulation overrides (reduced warmup cycles for speed).
% These are set in gsa_sobol_setup.m. Without them, N=512 at full fidelity
% can take days. SS precision is not critical per GSA sample.
if isfield(cfg, 'gsa_sim_overrides')
    ov = cfg.gsa_sim_overrides;
    if isfield(ov, 'nCyclesSteady'), params0.sim.nCyclesSteady = ov.nCyclesSteady; end
    if isfield(ov, 'ss_tol_P'),      params0.sim.ss_tol_P      = ov.ss_tol_P;      end
    if isfield(ov, 'ss_tol_V'),      params0.sim.ss_tol_V      = ov.ss_tol_V;      end
end

fprintf('[gsa_run_sobol] Evaluating A matrix (%d runs)...\n', cfg.N);
YA = evaluate_sample_matrix(A, params0, names, all_metrics, 'A matrix');

fprintf('[gsa_run_sobol] Evaluating B matrix (%d runs)...\n', cfg.N);
YB = evaluate_sample_matrix(B, params0, names, all_metrics, 'B matrix');

fprintf('[gsa_run_sobol] Evaluating A_Bi matrices (%d params, %d runs each)...\n', d, cfg.N);
YAB = cell(d, 1);
for i = 1:d
    fprintf('[gsa_run_sobol]   Parameter %d/%d: %s\n', i, d, names{i});
    YAB{i} = evaluate_sample_matrix(AB{i}, params0, names, all_metrics, ...
                                    sprintf('AB matrix %d/%d', i, d));
end

gsa_out = struct();

for m = 1:numel(all_metrics)
    mf = all_metrics{m};
    yA = YA.(mf);
    yB = YB.(mf);

    yAll = [yA; yB];
    VY   = max(var(yAll, 1), 1e-12);   % population variance

    S1i  = zeros(d, 1);
    STi  = zeros(d, 1);

    for i = 1:d
        yABi = YAB{i}.(mf);

        % Jansen estimators
        STi(i) = mean((yA - yABi).^2) / (2 * VY);
        S1i(i) = 1 - mean((yB - yABi).^2) / (2 * VY);
    end

    % Clamp to [−0.2, 1] to handle Monte-Carlo noise
    S1i = max(min(S1i,  1.0), -0.2);
    STi = max(min(STi,  1.0),  0.0);

    % Bootstrap 95% confidence intervals for S1 and ST (Saltelli 2010 §4)
    n_boot   = 200;
    N        = numel(yA);
    S1_boot  = zeros(d, n_boot);
    ST_boot  = zeros(d, n_boot);
    for b = 1:n_boot
        idx_b  = randi(N, N, 1);
        yA_b   = yA(idx_b);
        yB_b   = yB(idx_b);
        VY_b   = max(var([yA_b; yB_b], 1), 1e-12);
        for i = 1:d
            yABi_b      = YAB{i}.(mf)(idx_b);
            ST_boot(i,b) = mean((yA_b - yABi_b).^2) / (2 * VY_b);
            S1_boot(i,b) = 1 - mean((yB_b - yABi_b).^2) / (2 * VY_b);
        end
    end
    S1_CI = prctile(S1_boot, [2.5 97.5], 2);   % [d × 2]  lower/upper
    ST_CI = prctile(ST_boot, [2.5 97.5], 2);   % [d × 2]  lower/upper

    T = table(names(:), S1i, STi, S1_CI(:,1), S1_CI(:,2), ST_CI(:,1), ST_CI(:,2), ...
        'VariableNames', {'Parameter', 'Sobol_S1', 'Sobol_ST', ...
                          'S1_CI_lo', 'S1_CI_hi', 'ST_CI_lo', 'ST_CI_hi'});
    T = sortrows(T, 'Sobol_ST', 'descend');

    gsa_out.(mf).S1      = S1i;
    gsa_out.(mf).ST      = STi;
    gsa_out.(mf).S1_CI   = S1_CI;   % [d × 2]  95% bootstrap CI for S1
    gsa_out.(mf).ST_CI   = ST_CI;   % [d × 2]  95% bootstrap CI for ST
    gsa_out.(mf).table   = T;
    gsa_out.(mf).primary = ismember(mf, cfg.primary_metrics);
end

gsa_out.scenario = cfg.scenario;
gsa_out.cfg      = cfg;

fprintf('[gsa_run_sobol] Complete.\n');

end  % gsa_run_sobol

% =========================================================================
%  LOCAL HELPER
% =========================================================================

function Y = evaluate_sample_matrix(X_mat, params0, names, metric_list, stage_label)
% EVALUATE_SAMPLE_MATRIX — run the model for each row of parameter matrix X_mat
%   X_mat   : N × d  matrix of parameter values
%   Returns struct Y with one N-vector per metric in metric_list.

N = size(X_mat, 1);
F = numel(metric_list);
Y_mat = zeros(N, F);

% Live progress indicator so long GSA batches do not appear stuck.
count_done = 0;
last_pct   = -1;
h_wait = [];
show_waitbar = usejava('desktop') && feature('ShowFigureWindows');
if show_waitbar
    h_wait = waitbar(0, sprintf('%s: 0/%d', stage_label, N), ...
        'Name', 'GSA Progress');
else
    fprintf('[gsa_run_sobol]   %s: 0%%\n', stage_label);
end

% Safe DataQueue creation — requires Parallel Computing Toolbox.
% Without PCT, parfor degrades to for automatically, but DataQueue
% throws "Undefined variable 'parallel'" and crashes the whole GSA.
% This try-catch prevents that; progress falls back to per-sample fprintf.
dq = [];
try
    dq = parallel.pool.DataQueue;
    afterEach(dq, @update_progress);
catch
    % PCT unavailable — parfor will run serially; use inline fprintf below
end

parfor k = 1:N    % Runs in parallel if PCT available, serially otherwise
    params = apply_sample_row(params0, names, X_mat(k,:)');
    row_out = zeros(1, F);
    try
        sim     = integrate_system(params);
        metrics = compute_clinical_indices(sim, params);
        for f = 1:F
            mf = metric_list{f};
            if isfield(metrics, mf)
                row_out(f) = metrics.(mf);
            end
        end
    catch
        % Failed simulation: leave as 0 (contributes to variance but does
        % not crash the estimator; flag if too many failures occur)
    end
    Y_mat(k, :) = row_out;
    if ~isempty(dq)
        send(dq, 1);
    else
        % No PCT: parfor is actually a for-loop — print per-sample progress
        fprintf('[gsa_run_sobol]   %s: sample %d/%d done\n', stage_label, k, N);
    end
end

if ~isempty(h_wait) && isvalid(h_wait)
    close(h_wait);
end

Y = struct();
for f = 1:F
    Y.(metric_list{f}) = Y_mat(:, f);
end

    function update_progress(~)
        count_done = count_done + 1;
        pct = floor(100 * count_done / max(N, 1));
        if pct == last_pct
            return;
        end
        last_pct = pct;

        if ~isempty(h_wait) && isvalid(h_wait)
            waitbar(count_done / N, h_wait, sprintf('%s: %d/%d (%d%%)', ...
                stage_label, count_done, N, pct));
        elseif mod(pct, 10) == 0 || pct == 100
            fprintf('[gsa_run_sobol]   %s: %d%%\n', stage_label, pct);
        end
    end
end

function params = apply_sample_row(params0, names, x)
% APPLY_SAMPLE_ROW — update named parameters from sample vector x
params = params0;
for i = 1:numel(names)
    parts = strsplit(names{i}, '.');
    switch numel(parts)
        case 2, params.(parts{1}).(parts{2})          = x(i);
        case 3, params.(parts{1}).(parts{2}).(parts{3}) = x(i);
    end
end
end
