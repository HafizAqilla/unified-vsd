function summary = aggregate_noise_simulation(batch_dir)
% AGGREGATE_NOISE_SIMULATION
% -----------------------------------------------------------------------
% Collects the replicate results of run_noise_simulation into long-format
% CSVs (one row per replicate x quantity) and summary tables. Every
% replicate is reported, including failed ones and those whose calibration
% status is not 'accepted'; nothing is filtered by result quality.
%
% OUTPUTS (<batch_dir>/summary/):
%   noise_runs.csv                 one row per replicate: status, gate,
%                                  chi2/N, post within-10% count, time
%   noise_parameters_long.csv      calibrated value of every fitted parameter
%   noise_pre_long.csv             pre-closure model vs original and
%                                  perturbed record
%   noise_post_long.csv            post-closure prediction vs held-out record
%   noise_summary_parameters.csv   mean, SD, CV, 2.5/50/97.5 percentiles,
%                                  and the unperturbed reference value
%   noise_summary_pre.csv          same for each pre-closure target
%   noise_summary_post.csv         same for each post-closure target, plus
%                                  the fraction of replicates within 10% and
%                                  whether the measurement lies inside the
%                                  2.5-97.5 percentile band of the prediction
%   noise_summary.txt              headline numbers
%
% Statistics use completed replicates only; the count of failed replicates
% is stated in noise_summary.txt.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-26
% VERSION:  1.0
% -----------------------------------------------------------------------

rep_dir = fullfile(batch_dir, 'replicates');
out_dir = fullfile(batch_dir, 'summary');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
M = load(fullfile(batch_dir, 'mask_reference.mat'));

files = dir(fullfile(rep_dir, 'mc*.mat'));
assert(~isempty(files), 'aggregate_noise_simulation:none', 'No replicate results in %s', rep_dir);

runs = table(); P = {}; PRE = {}; POST = {};
for f = 1:numel(files)
    S = load(fullfile(files(f).folder, files(f).name));
    r = S.r;
    ok = strcmp(r.status, 'completed');
    row = table(r.replicate, {r.status}, {get_or(r, 'calibration_status', '')}, ...
        get_or(r, 'rollback_applied', false), get_or(r, 'multistart_best_start', NaN), ...
        NaN, NaN, get_or(r, 'chi2_per_obs', NaN), ...
        get_or(r, 'post_within10', NaN), get_or(r, 'post_n', NaN), ...
        get_or(r, 'post_chi2_per_n', NaN), r.wall_minutes, {r.run_dir}, ...
        {first_line(r.message)}, ...
        'VariableNames', {'Replicate','Status','CalibrationStatus','RollbackApplied', ...
        'SelectedStart','GovernedRMSE_vsPerturbed','GovernedWithinGate', ...
        'Chi2PerObs','PostWithin10','PostN','PostChi2PerN','WallMinutes','RunFolder','Error'});
    if ok
        gate = r.pre_gate_table;
        g = logical(gate.InPrimaryRMSE);
        row.GovernedRMSE_vsPerturbed = sqrt(mean((gate.Error_pct(g) / 100) .^ 2));
        row.GovernedWithinGate = nnz(gate.WithinGate(g));
        P{end + 1} = table(repmat(r.replicate, numel(r.parameter_names), 1), ...
            r.parameter_names(:), r.parameter_values(:), r.parameter_initial(:), ...
            'VariableNames', {'Replicate','Parameter','Calibrated','Initial'}); %#ok<AGROW>
        orig = renamevars(r.pre_vs_original, {'Measured','Error_pct'}, ...
            {'MeasuredOriginal','ErrorVsOriginal_pct'});
        pert = renamevars(removevars(r.pre_vs_perturbed, 'Model'), ...
            {'Measured','Error_pct'}, {'MeasuredPerturbed','ErrorVsPerturbed_pct'});
        pre = outerjoin(orig, pert, 'Keys', 'Metric', 'MergeKeys', true, 'Type', 'left');
        pre.Governed = ismember(pre.Metric, gate.Metric(g));
        pre = addvars(pre, repmat(r.replicate, height(pre), 1), 'Before', 1, ...
            'NewVariableNames', 'Replicate');
        PRE{end + 1} = pre; %#ok<AGROW>
        post = r.post_table;
        post = addvars(post, repmat(r.replicate, height(post), 1), 'Before', 1, ...
            'NewVariableNames', 'Replicate');
        POST{end + 1} = post; %#ok<AGROW>
    end
    runs = [runs; row]; %#ok<AGROW>
end
runs = sortrows(runs, 'Replicate');
writetable(runs, fullfile(out_dir, 'noise_runs.csv'));

n_total = height(runs);
n_ok = nnz(strcmp(runs.Status, 'completed'));
summary = struct('n_total', n_total, 'n_completed', n_ok);
if n_ok == 0
    write_text(out_dir, runs, summary, table(), table());
    return;
end

P = vertcat(P{:}); PRE = vertcat(PRE{:}); POST = vertcat(POST{:});
writetable(P, fullfile(out_dir, 'noise_parameters_long.csv'));
writetable(PRE, fullfile(out_dir, 'noise_pre_long.csv'));
writetable(POST, fullfile(out_dir, 'noise_post_long.csv'));

% ---- parameters ------------------------------------------------------------
sp = describe(P, 'Parameter', 'Calibrated', M.reference_active_names(:));
ref = nan(height(sp), 1);
for i = 1:height(sp)
    ref(i) = M.reference_x_active(strcmp(M.reference_active_names, sp.Parameter{i}));
end
sp = addvars(sp, ref, 'After', 'Parameter', 'NewVariableNames', 'Reference');
writetable(sp, fullfile(out_dir, 'noise_summary_parameters.csv'));

% ---- pre-closure ------------------------------------------------------------
spre = describe(PRE, 'Metric', 'Model', unique(PRE.Metric, 'stable'));
[~, ia] = unique(PRE.Metric, 'stable');
spre = addvars(spre, PRE.MeasuredOriginal(ia), PRE.Governed(ia), 'After', 'Metric', ...
    'NewVariableNames', {'Measured','Governed'});
writetable(spre, fullfile(out_dir, 'noise_summary_pre.csv'));

% ---- post-closure -----------------------------------------------------------
metrics = unique(POST.Metric, 'stable');
spost = describe(POST, 'Metric', 'Predicted', metrics);
measured = nan(numel(metrics), 1); frac10 = nan(numel(metrics), 1);
covered = false(numel(metrics), 1);
for i = 1:numel(metrics)
    sel = strcmp(POST.Metric, metrics{i});
    measured(i) = POST.Measured(find(sel, 1));
    frac10(i) = mean(abs(POST.Error_pct(sel)) <= 10);
    covered(i) = measured(i) >= spost.P2_5(i) && measured(i) <= spost.P97_5(i);
end
spost = addvars(spost, measured, 'After', 'Metric', 'NewVariableNames', 'Measured');
spost.MeanError_pct = 100 * (spost.Mean - measured) ./ abs(measured);
spost.FractionWithin10 = frac10;
spost.MeasuredInside95Band = covered;
writetable(spost, fullfile(out_dir, 'noise_summary_post.csv'));

summary.parameters = sp;
summary.pre = spre;
summary.post = spost;
write_text(out_dir, runs, summary, sp, spost);
disp(sp); disp(spost);
end

% =======================================================================
function T = describe(D, key, value, order)
n = numel(order);
Mean = nan(n, 1); SD = nan(n, 1); CV_pct = nan(n, 1);
P2_5 = nan(n, 1); P50 = nan(n, 1); P97_5 = nan(n, 1); N = zeros(n, 1);
for i = 1:n
    x = D.(value)(strcmp(D.(key), order{i}));
    x = x(isfinite(x));
    N(i) = numel(x);
    if isempty(x), continue; end
    Mean(i) = mean(x); SD(i) = std(x);
    CV_pct(i) = 100 * SD(i) / abs(Mean(i));
    q = prctile(x, [2.5 50 97.5]);
    P2_5(i) = q(1); P50(i) = q(2); P97_5(i) = q(3);
end
T = table(order(:), N, Mean, SD, CV_pct, P2_5, P50, P97_5, ...
    'VariableNames', {key, 'N', 'Mean', 'SD', 'CV_pct', 'P2_5', 'P50', 'P97_5'});
end

function write_text(out_dir, runs, summary, sp, spost)
fid = fopen(fullfile(out_dir, 'noise_summary.txt'), 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, 'Measurement-noise simulation summary\n');
fprintf(fid, 'Replicates: %d total, %d completed, %d failed\n', summary.n_total, ...
    summary.n_completed, summary.n_total - summary.n_completed);
if summary.n_completed == 0, return; end
ok = strcmp(runs.Status, 'completed');
labels = unique(runs.CalibrationStatus(ok));
for i = 1:numel(labels)
    fprintf(fid, '  status %-30s %d\n', labels{i}, nnz(strcmp(runs.CalibrationStatus(ok), labels{i})));
end
fprintf(fid, 'Governed RMSE vs perturbed record: median %.4f (range %.4f-%.4f)\n', ...
    median(runs.GovernedRMSE_vsPerturbed(ok)), min(runs.GovernedRMSE_vsPerturbed(ok)), ...
    max(runs.GovernedRMSE_vsPerturbed(ok)));
fprintf(fid, 'Governed targets within 10%%: median %g of 9 (range %g-%g)\n', ...
    median(runs.GovernedWithinGate(ok)), min(runs.GovernedWithinGate(ok)), max(runs.GovernedWithinGate(ok)));
fprintf(fid, 'Post-closure targets within 10%%: median %g of %g (range %g-%g)\n', ...
    median(runs.PostWithin10(ok)), median(runs.PostN(ok)), min(runs.PostWithin10(ok)), ...
    max(runs.PostWithin10(ok)));
fprintf(fid, '\nParameters (CV %% across replicates):\n');
for i = 1:height(sp)
    fprintf(fid, '  %-22s ref %.4g  mean %.4g  CV %.1f%%  95%% [%.4g, %.4g]\n', sp.Parameter{i}, ...
        sp.Reference(i), sp.Mean(i), sp.CV_pct(i), sp.P2_5(i), sp.P97_5(i));
end
fprintf(fid, '\nPost-closure prediction:\n');
for i = 1:height(spost)
    fprintf(fid, '  %-10s measured %.4g  predicted %.4g [%.4g, %.4g]  within10 %.0f%%  covered %d\n', ...
        spost.Metric{i}, spost.Measured(i), spost.Mean(i), spost.P2_5(i), spost.P97_5(i), ...
        100 * spost.FractionWithin10(i), spost.MeasuredInside95Band(i));
end
end

function v = get_or(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function s = first_line(msg)
s = '';
if ~isempty(msg)
    lines = splitlines(string(msg));
    s = char(lines(1));
end
end
