function gsa_out = gsa_run_pce(cfg, params0)
% GSA_RUN_PCE — see header for full documentation.
% VERSION 2.1: Added per-metric autosave for crash recovery.

names       = cfg.names;
d           = numel(names);
all_metrics = cfg.all_metrics;

gsa_out = struct();

%% Per-metric autosave setup (checkpoint file for crash recovery)
autosave_dir  = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'gsa');
if ~exist(autosave_dir, 'dir'), mkdir(autosave_dir); end

p_weight = round(params0.scaling.patient.weight_kg, 1);
autosave_file = fullfile(autosave_dir, ...
    sprintf('gsa_pce_%s_%.1fkg_checkpoint.mat', cfg.scenario, p_weight));
fprintf('[gsa_run_pce] Autosave checkpoint: %s\n', autosave_file);
if exist(autosave_file, 'file')
    tmp = load(autosave_file, 'gsa_out');
    ck = tmp.gsa_out;
    % Validate checkpoint: if cfg.x0 changed (patient data was corrected),
    % the surrogate bounds are stale and must be discarded.
    stale = false;
    if isfield(ck, 'cfg') && isfield(ck.cfg, 'x0') && ...
            numel(ck.cfg.x0) == numel(cfg.x0)
        rel_diff = max(abs(ck.cfg.x0 - cfg.x0) ./ max(abs(cfg.x0), 1e-12));
        if rel_diff > 1e-4
            stale = true;
            fprintf('[gsa_run_pce] Checkpoint x0 mismatch (max rel diff=%.2e) — discarding stale checkpoint, retraining.\n', rel_diff);
        end
    else
        stale = true;
        fprintf('[gsa_run_pce] Checkpoint incompatible (size/field mismatch) — discarding.\n');
    end
    if ~stale
        fprintf('[Crash Recovery] Loading previous partial results...\n');
        gsa_out = ck;
    end
end

%% =====================================================================
%  Loop over each output metric
%  (following the SoBioS pattern: one PCE surrogate per QoI)
%% =====================================================================
for m = 1:numel(all_metrics)
    mf = all_metrics{m};
    fprintf('\n[gsa_run_pce] --- Metric %d/%d: %s ---\n', m, numel(all_metrics), mf);
    if isfield(gsa_out, mf) && isfield(gsa_out.(mf), 'surrogate')
        fprintf('[gsa_run_pce] Loaded %s from checkpoint, skipping training...\n', mf);
        continue;
    end

    % ------------------------------------------------------------------
    % STEP 1: Build UQLab model wrapper for this metric
    % ------------------------------------------------------------------
    ModelOpts.mFile            = 'gsa_pce_vsd_qoi';
    ModelOpts.isVectorized     = false;
    ModelOpts.Parameters.cfg   = cfg;
    ModelOpts.Parameters.mf    = mf;
    ModelOpts.Parameters.names = names;
    ModelOpts.Parameters.p0    = params0;

    myModel = uq_createModel(ModelOpts, '-private');

    % ------------------------------------------------------------------
    % STEP 2: Train PCE surrogate for this metric
    % ------------------------------------------------------------------
    PO           = cfg.PCEOpts;
    PO.FullModel = myModel;

    fprintf('[gsa_run_pce] Training PCE for %s (%d samples, method: %s)...\n', ...
            mf, PO.ExpDesign.NSamples, PO.Method);

    myPCE = uq_createModel(PO, '-private');

    if isfield(myPCE, 'Error') && isfield(myPCE.Error, 'LOO')
        fprintf('[gsa_run_pce] LOO error for %s: %.4f\n', mf, myPCE.Error.LOO);
    end

    % ------------------------------------------------------------------
    % STEP 3: Compute PCE-Sobol indices analytically
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
    % STEP 4: Package results (same format as gsa_run_sobol)
    % ------------------------------------------------------------------
    S1i = res.FirstOrder(:);
    STi = res.Total(:);
    S1i = max(min(S1i,  1.0), -0.2);
    STi = max(min(STi,  1.0),  0.0);

    T = table(names(:), S1i, STi, ...
        'VariableNames', {'Parameter', 'Sobol_S1', 'Sobol_ST'});
    T = sortrows(T, 'Sobol_ST', 'descend');

    gsa_out.(mf).S1      = S1i;
    gsa_out.(mf).ST      = STi;
    gsa_out.(mf).table   = T;
    gsa_out.(mf).primary = ismember(mf, cfg.primary_metrics);
    gsa_out.(mf).surrogate = myPCE;

    % ------------------------------------------------------------------
    % AUTOSAVE: write checkpoint after every completed metric
    %   If MATLAB crashes, load this file to recover partial results.
    % ------------------------------------------------------------------
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
