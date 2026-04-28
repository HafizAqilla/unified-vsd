function gsa_out = gsa_run_pce(cfg, params0)
% GSA_RUN_PCE — see header for full documentation.
% VERSION 2.2: Shared experimental design — all metrics use one ODE batch.
%
% KEY CHANGE FROM v2.1:
%   Previously, UQLab re-ran 200 ODE simulations per metric (200 × n_metrics
%   total). Now, 200 Halton samples are evaluated ONCE and all metric values
%   are extracted in a single batch. Each PCE then receives pre-computed
%   (X_train, Y_train) instead of a FullModel, reducing ODE cost from
%   ~5000 to 200 per GSA pass.

names       = cfg.names;
d           = numel(names);
all_metrics = cfg.all_metrics;

gsa_out = struct();

%% =====================================================================
%  STEP 0 — Generate shared experimental design (one ODE batch for all metrics)
%
%  Previously: 200 ODE runs × n_metrics per GSA pass (~5000 total).
%  Now: 200 ODE runs once; each PCE receives pre-computed (X_train, Y_train).
%  Speedup: ~25× for pre_surgery (25 metrics), ~20× for post_surgery.
%% =====================================================================
N_train = cfg.PCEOpts.ExpDesign.NSamples;  % [dimensionless] number of training samples

%% Per-metric autosave setup (checkpoint file for crash recovery)
autosave_dir  = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'gsa');
if ~exist(autosave_dir, 'dir'), mkdir(autosave_dir); end

p_weight = round(params0.scaling.patient.weight_kg, 1);
autosave_file = fullfile(autosave_dir, ...
    sprintf('gsa_pce_%s_%.1fkg_checkpoint.mat', cfg.scenario, p_weight));
fprintf('[gsa_run_pce] Autosave checkpoint: %s\n', autosave_file);

% Checkpoint validation: discard if patient x0 changed (stale bounds).
stale = true;
if exist(autosave_file, 'file')
    tmp = load(autosave_file, 'gsa_out');
    ck = tmp.gsa_out;
    if isfield(ck, 'cfg') && isfield(ck.cfg, 'x0') && ...
            numel(ck.cfg.x0) == numel(cfg.x0)
        rel_diff = max(abs(ck.cfg.x0 - cfg.x0) ./ max(abs(cfg.x0), 1e-12));
        if rel_diff <= 1e-4
            stale = false;
            fprintf('[Crash Recovery] Loading previous partial results...\n');
            gsa_out = ck;
        else
            fprintf('[gsa_run_pce] Checkpoint x0 mismatch (max rel diff=%.2e) — discarding stale checkpoint, retraining.\n', rel_diff);
        end
    else
        fprintf('[gsa_run_pce] Checkpoint incompatible (size/field mismatch) — discarding.\n');
    end
end

% Check whether the shared training batch is already in the checkpoint.
need_batch = ~isfield(gsa_out, 'X_train') || ~isfield(gsa_out, 'Y_all');

if need_batch
    fprintf('\n[gsa_run_pce] Generating shared Halton design (%d samples, d=%d)...\n', ...
            N_train, d);
    X_train = uq_getSample(cfg.Input, N_train, 'Halton');   % [N_train × d]
    Y_all   = nan(N_train, numel(all_metrics));              % [N_train × n_metrics]

    % Apply reduced-fidelity sim overrides to params0 for the GSA batch.
    % This reduces warm-up from 80 cycles (default) to 10 cycles per sample,
    % cutting ODE cost by ~8x without materially affecting Sobol index rankings.
    if isfield(cfg, 'gsa_sim_overrides')
        ov = cfg.gsa_sim_overrides;
        if isfield(ov, 'nCyclesSteady'), params0.sim.nCyclesSteady = ov.nCyclesSteady; end
        if isfield(ov, 'ss_tol_P'),      params0.sim.ss_tol_P      = ov.ss_tol_P;      end
        if isfield(ov, 'ss_tol_V'),      params0.sim.ss_tol_V      = ov.ss_tol_V;      end
        fprintf('[gsa_run_pce] Sim overrides applied: nCyclesSteady=%d, ss_tol_P=%.1f, ss_tol_V=%.1f\n', ...
                params0.sim.nCyclesSteady, params0.sim.ss_tol_P, params0.sim.ss_tol_V);
    end

    % Single batch: run the full ODE model once per sample, extract all metrics.
    for n = 1:N_train
        if mod(n, 50) == 0 || n == 1
            fprintf('[gsa_run_pce] ODE batch: sample %d/%d...\n', n, N_train);
        end
        params_n = params0;
        for i = 1:d
            parts = strsplit(names{i}, '.');
            switch numel(parts)
                case 2; params_n.(parts{1}).(parts{2})            = X_train(n, i);
                case 3; params_n.(parts{1}).(parts{2}).(parts{3}) = X_train(n, i);
            end
        end
        try
            sim_n     = integrate_system(params_n);
            metrics_n = compute_clinical_indices(sim_n, params_n);
            for m = 1:numel(all_metrics)
                mf = all_metrics{m};
                if isfield(metrics_n, mf)
                    Y_all(n, m) = metrics_n.(mf);
                else
                    Y_all(n, m) = 0;
                end
            end
        catch
            Y_all(n, :) = 0;   % failed simulation: fill row with zeros
        end
    end

    gsa_out.X_train = X_train;
    gsa_out.Y_all   = Y_all;
    gsa_out.scenario = cfg.scenario;
    gsa_out.cfg      = cfg;
    save(autosave_file, 'gsa_out', '-v7.3');
    fprintf('[gsa_run_pce] ODE batch complete. Checkpoint saved.\n');
else
    X_train = gsa_out.X_train;
    Y_all   = gsa_out.Y_all;
    fprintf('[gsa_run_pce] Loaded pre-computed ODE batch from checkpoint (%d samples).\n', N_train);
end

%% =====================================================================
%  Loop over each output metric — fit one PCE per QoI using shared data
%  (no additional ODE evaluations; PCE fitting only)
%% =====================================================================
for m = 1:numel(all_metrics)
    mf = all_metrics{m};
    fprintf('\n[gsa_run_pce] --- Metric %d/%d: %s ---\n', m, numel(all_metrics), mf);
    if isfield(gsa_out, mf) && isfield(gsa_out.(mf), 'surrogate')
        fprintf('[gsa_run_pce] Loaded %s from checkpoint, skipping training...\n', mf);
        continue;
    end

    % ------------------------------------------------------------------
    % STEP 1: Train PCE surrogate using pre-computed (X_train, Y_train)
    %         No FullModel needed — data already evaluated in batch above.
    % ------------------------------------------------------------------
    Y_train = Y_all(:, m);   % [N_train × 1] scalar output for this metric

    PO = cfg.PCEOpts;
    PO = rmfield_safe(PO, 'FullModel');    % ensure no stale FullModel reference
    PO.ExpDesign.Sampling = 'user';        % required when providing X,Y manually
    PO.ExpDesign.X = X_train;             % [N_train × d] pre-computed inputs
    PO.ExpDesign.Y = Y_train;             % [N_train × 1] pre-computed outputs

    fprintf('[gsa_run_pce] Fitting PCE for %s (pre-computed data, method: %s)...\n', ...
            mf, PO.Method);

    myPCE = uq_createModel(PO, '-private');

    if isfield(myPCE, 'Error') && isfield(myPCE.Error, 'LOO')
        fprintf('[gsa_run_pce] LOO error for %s: %.4f\n', mf, myPCE.Error.LOO);
    end

    % ------------------------------------------------------------------
    % STEP 2: Compute PCE-Sobol indices analytically
    % ------------------------------------------------------------------
    SobolPCE              = struct();
    SobolPCE.Type         = 'Sensitivity';
    SobolPCE.Method       = 'Sobol';
    SobolPCE.Model        = myPCE;
    SobolPCE.Input        = cfg.Input;
    SobolPCE.Sobol.Order  = 1;

    mySobolAnalysis = uq_createAnalysis(SobolPCE, '-private');
    res             = mySobolAnalysis.Results;

    % ------------------------------------------------------------------
    % STEP 3: Package results (same format as gsa_run_sobol)
    % ------------------------------------------------------------------
    S1i = res.FirstOrder(:);
    STi = res.Total(:);
    S1i = max(min(S1i,  1.0), -0.2);
    STi = max(min(STi,  1.0),  0.0);

    T = table(names(:), S1i, STi, ...
        'VariableNames', {'Parameter', 'Sobol_S1', 'Sobol_ST'});
    T = sortrows(T, 'Sobol_ST', 'descend');

    gsa_out.(mf).S1        = S1i;
    gsa_out.(mf).ST        = STi;
    gsa_out.(mf).table     = T;
    gsa_out.(mf).primary   = ismember(mf, cfg.primary_metrics);
    gsa_out.(mf).surrogate = myPCE;

    % AUTOSAVE after every completed metric for crash recovery.
    gsa_out.scenario = cfg.scenario;
    gsa_out.cfg      = cfg;
    save(autosave_file, 'gsa_out', '-v7.3');
    fprintf('[gsa_run_pce] Checkpoint saved (%d/%d metrics done).\n', ...
            m, numel(all_metrics));
end

gsa_out.scenario = cfg.scenario;
gsa_out.cfg      = cfg;

fprintf('\n[gsa_run_pce] Complete. All %d metrics processed.\n', numel(all_metrics));

end  % gsa_run_pce


% =========================================================================
%  LOCAL HELPER
% =========================================================================
function s = rmfield_safe(s, f)
% RMFIELD_SAFE — remove field f from struct s only if it exists.
if isfield(s, f)
    s = rmfield(s, f);
end
end
