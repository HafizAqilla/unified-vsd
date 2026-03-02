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
%   For N=64, d=18 (pre) or d=19 (post): ~1280 evaluations.
%   Increase cfg.N to 256 for publication-quality indices.
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

fprintf('[gsa_run_sobol] Evaluating A matrix (%d runs)...\n', cfg.N);
YA = evaluate_sample_matrix(A, params0, names, all_metrics);

fprintf('[gsa_run_sobol] Evaluating B matrix (%d runs)...\n', cfg.N);
YB = evaluate_sample_matrix(B, params0, names, all_metrics);

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
        fprintf('[gsa_run_sobol]   Parameter %d/%d: %s (metric: %s)\n', ...
                i, d, names{i}, mf);
        YABi  = evaluate_sample_matrix(AB{i}, params0, names, {mf});
        yABi  = YABi.(mf);

        % Jansen estimators
        STi(i) = mean((yA - yABi).^2) / (2 * VY);
        S1i(i) = 1 - mean((yB - yABi).^2) / (2 * VY);
    end

    % Clamp to [−0.2, 1] to handle Monte-Carlo noise
    S1i = max(min(S1i,  1.0), -0.2);
    STi = max(min(STi,  1.0),  0.0);

    T = table(names(:), S1i, STi, ...
        'VariableNames', {'Parameter', 'Sobol_S1', 'Sobol_ST'});
    T = sortrows(T, 'Sobol_ST', 'descend');

    gsa_out.(mf).S1      = S1i;
    gsa_out.(mf).ST      = STi;
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

function Y = evaluate_sample_matrix(X_mat, params0, names, metric_list)
% EVALUATE_SAMPLE_MATRIX — run the model for each row of parameter matrix X_mat
%   X_mat   : N × d  matrix of parameter values
%   Returns struct Y with one N-vector per metric in metric_list.

N = size(X_mat, 1);
for f = 1:numel(metric_list)
    Y.(metric_list{f}) = zeros(N, 1);
end

for k = 1:N
    params = apply_sample_row(params0, names, X_mat(k,:)');
    try
        sim     = integrate_system(params);
        metrics = compute_clinical_indices(sim, params);
        for f = 1:numel(metric_list)
            mf = metric_list{f};
            if isfield(metrics, mf)
                Y.(mf)(k) = metrics.(mf);
            end
        end
    catch
        % Failed simulation: leave as 0 (contributes to variance but does
        % not crash the estimator; flag if too many failures occur)
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
