function gsa_out = gsa_run_pce(cfg, params0)
% GSA_RUN_PCE
% -----------------------------------------------------------------------
% Run PCE-based global sensitivity analysis with one shared ODE batch.
%
% The expensive model evaluations are performed once for a shared Halton
% training design. Each output metric then trains its own PCE surrogate from
% the pre-computed (X_train, Y_train) data, so adding metrics does not add
% another ODE batch.
% -----------------------------------------------------------------------

names = cfg.names;
d = numel(names);
all_metrics = cfg.all_metrics;

gsa_out = struct();

%% Generate shared experimental design
N_train = cfg.PCEOpts.ExpDesign.NSamples;

autosave_dir = getenv('UNIFIED_VSD_GSA_DIR');
if isempty(autosave_dir)
    autosave_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'gsa');
end
if ~exist(autosave_dir, 'dir')
    mkdir(autosave_dir);
end

autosave_file = getenv('UNIFIED_VSD_GSA_CHECKPOINT_FILE');
if isempty(autosave_file)
    p_weight = round(params0.scaling.patient.weight_kg, 1);
    autosave_file = fullfile(autosave_dir, ...
        sprintf('gsa_pce_%s_%.1fkg_checkpoint.mat', cfg.scenario, p_weight));
end
fprintf('[gsa_run_pce] Autosave checkpoint: %s\n', autosave_file);

checkpoint_signature = build_checkpoint_signature(cfg, N_train, all_metrics);

% Checkpoint validation: names, bounds, metrics, sample count, and nominal
% point must all match.  A matching x0 alone is insufficient because a
% stale checkpoint can silently reuse an obsolete GSA mask or metric set.
if exist(autosave_file, 'file')
    tmp = load(autosave_file, 'gsa_out');
    ck = tmp.gsa_out;
    [checkpoint_ok, checkpoint_reason] = checkpoint_matches(ck, checkpoint_signature);
    if checkpoint_ok
        fprintf('[Crash Recovery] Loading matching GSA checkpoint.\n');
        gsa_out = ck;
    else
        fprintf('[gsa_run_pce] Checkpoint rejected (%s); retraining.\n', checkpoint_reason);
    end
end

need_batch = ~isfield(gsa_out, 'X_train') || ~isfield(gsa_out, 'Y_all');

if need_batch
    fprintf('\n[gsa_run_pce] Generating shared Halton design (%d samples, d=%d)...\n', ...
        N_train, d);
    X_train = uq_getSample(cfg.Input, N_train, 'Halton');
    Y_all = nan(N_train, numel(all_metrics));

    if isfield(cfg, 'gsa_sim_overrides')
        ov = cfg.gsa_sim_overrides;
        if isfield(ov, 'nCyclesSteady')
            params0.sim.nCyclesSteady = ov.nCyclesSteady;
        end
        if isfield(ov, 'ss_tol_P')
            params0.sim.ss_tol_P = ov.ss_tol_P;
        end
        if isfield(ov, 'ss_tol_V')
            params0.sim.ss_tol_V = ov.ss_tol_V;
        end
        fprintf('[gsa_run_pce] Sim overrides applied: nCyclesSteady=%d, ss_tol_P=%.1f, ss_tol_V=%.1f\n', ...
            params0.sim.nCyclesSteady, params0.sim.ss_tol_P, params0.sim.ss_tol_V);
    end

    for n = 1:N_train
        if mod(n, 50) == 0 || n == 1
            fprintf('[gsa_run_pce] ODE batch: sample %d/%d...\n', n, N_train);
        end
        params_n = params0;
        for i = 1:d
            parts = strsplit(names{i}, '.');
            switch numel(parts)
                case 2
                    params_n.(parts{1}).(parts{2}) = X_train(n, i);
                case 3
                    params_n.(parts{1}).(parts{2}).(parts{3}) = X_train(n, i);
            end
        end
        try
            sim_n = integrate_system(params_n);
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
            Y_all(n, :) = 0;
        end
    end

    gsa_out.X_train = X_train;
    gsa_out.Y_all = Y_all;
    gsa_out.scenario = cfg.scenario;
    gsa_out.cfg = cfg;
    gsa_out.checkpoint_signature = checkpoint_signature;
    save(autosave_file, 'gsa_out', '-v7.3');
    fprintf('[gsa_run_pce] ODE batch complete. Checkpoint saved.\n');
else
    X_train = gsa_out.X_train;
    Y_all = gsa_out.Y_all;
    fprintf('[gsa_run_pce] Loaded pre-computed ODE batch from checkpoint (%d samples).\n', N_train);
end

%% Fit one PCE per output metric from shared data
for m = 1:numel(all_metrics)
    mf = all_metrics{m};
    fprintf('\n[gsa_run_pce] --- Metric %d/%d: %s ---\n', m, numel(all_metrics), mf);
    if isfield(gsa_out, mf) && isfield(gsa_out.(mf), 'surrogate')
        fprintf('[gsa_run_pce] Loaded %s from checkpoint, skipping training...\n', mf);
        continue;
    end

    Y_train = Y_all(:, m);

    PO = cfg.PCEOpts;
    PO = rmfield_safe(PO, 'FullModel');
    PO.ExpDesign.Sampling = 'user';
    PO.ExpDesign.X = X_train;
    PO.ExpDesign.Y = Y_train;

    fprintf('[gsa_run_pce] Fitting PCE for %s (pre-computed data, method: %s)...\n', ...
        mf, PO.Method);

    myPCE = uq_createModel(PO, '-private');

    if isfield(myPCE, 'Error') && isfield(myPCE.Error, 'LOO')
        fprintf('[gsa_run_pce] LOO error for %s: %.4f\n', mf, myPCE.Error.LOO);
    end

    SobolPCE = struct();
    SobolPCE.Type = 'Sensitivity';
    SobolPCE.Method = 'Sobol';
    SobolPCE.Model = myPCE;
    SobolPCE.Input = cfg.Input;
    SobolPCE.Sobol.Order = 1;

    mySobolAnalysis = uq_createAnalysis(SobolPCE, '-private');
    res = mySobolAnalysis.Results;

    S1i = res.FirstOrder(:);
    STi = res.Total(:);
    S1i = max(min(S1i, 1.0), -0.2);
    STi = max(min(STi, 1.0), 0.0);

    T = table(names(:), S1i, STi, ...
        'VariableNames', {'Parameter', 'Sobol_S1', 'Sobol_ST'});
    T = sortrows(T, 'Sobol_ST', 'descend');

    gsa_out.(mf).S1 = S1i;
    gsa_out.(mf).ST = STi;
    gsa_out.(mf).table = T;
    gsa_out.(mf).primary = ismember(mf, cfg.primary_metrics);
    gsa_out.(mf).surrogate = myPCE;

    gsa_out.scenario = cfg.scenario;
    gsa_out.cfg = cfg;
    save(autosave_file, 'gsa_out', '-v7.3');
    fprintf('[gsa_run_pce] Checkpoint saved (%d/%d metrics done).\n', ...
        m, numel(all_metrics));
end

gsa_out.scenario = cfg.scenario;
gsa_out.cfg = cfg;
gsa_out.checkpoint_signature = checkpoint_signature;

fprintf('\n[gsa_run_pce] Complete. All %d metrics processed.\n', numel(all_metrics));

end

function signature = build_checkpoint_signature(cfg, N_train, all_metrics)
% BUILD_CHECKPOINT_SIGNATURE - identity of the expensive GSA batch.
signature = struct();
signature.version = 1;
signature.scenario = cfg.scenario;
signature.names = cfg.names(:)';
signature.x0 = cfg.x0(:);
signature.lb = cfg.lb(:);
signature.ub = cfg.ub(:);
signature.all_metrics = all_metrics(:)';
signature.N_train = N_train;
end

function [tf, reason] = checkpoint_matches(checkpoint, expected)
tf = false;
reason = 'missing_signature';
if ~isstruct(checkpoint) || ~isfield(checkpoint, 'checkpoint_signature')
    return;
end
actual = checkpoint.checkpoint_signature;
required = fieldnames(expected);
for idx = 1:numel(required)
    if ~isfield(actual, required{idx})
        reason = ['signature_missing_' required{idx}];
        return;
    end
end
if ~strcmp(actual.scenario, expected.scenario) || actual.N_train ~= expected.N_train
    reason = 'scenario_or_sample_count_mismatch';
    return;
end
if ~isequal(actual.names, expected.names) || ~isequal(actual.all_metrics, expected.all_metrics)
    reason = 'parameter_or_metric_names_mismatch';
    return;
end
for field_name = {'x0','lb','ub'}
    field = field_name{1};
    if numel(actual.(field)) ~= numel(expected.(field)) || ...
            max(abs(actual.(field)(:) - expected.(field)(:))) > 1e-10 * ...
            max(1, max(abs(expected.(field))))
        reason = [field '_mismatch'];
        return;
    end
end
if ~isfield(checkpoint, 'X_train') || ~isfield(checkpoint, 'Y_all') || ...
        size(checkpoint.X_train, 1) ~= expected.N_train || ...
        size(checkpoint.X_train, 2) ~= numel(expected.names) || ...
        size(checkpoint.Y_all, 1) ~= expected.N_train || ...
        size(checkpoint.Y_all, 2) ~= numel(expected.all_metrics)
    reason = 'batch_dimensions_mismatch';
    return;
end
tf = true;
reason = 'match';
end

function s = rmfield_safe(s, f)
if isfield(s, f)
    s = rmfield(s, f);
end
end
