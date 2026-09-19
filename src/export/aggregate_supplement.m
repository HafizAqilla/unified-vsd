function aggregate_supplement(run_dirs, seeds, out_dir)
% AGGREGATE_SUPPLEMENT
% -----------------------------------------------------------------------
% Combines the per-run supplement CSVs of several seeds into long tables
% with a Seed column (<name>_all.csv) and adds cross-seed summaries
% (mean, sample SD, min, max; n = number of runs). These are the only
% inputs the figure scripts read.
%
% INPUTS:
%   run_dirs - cell array of run folders that finished and were exported   [-]
%   seeds    - multi-start seed of each run, same order                     [-]
%   out_dir  - destination folder                                            [-]
% -----------------------------------------------------------------------

if ~exist(out_dir, 'dir'), mkdir(out_dir); end
files = {'protocol_comparison_pre','protocol_comparison_post','inputs', ...
    'parameters_calibration','gsa_indices_initial','gsa_indices_final', ...
    'gsa_surrogate_quality','gsa_mask','multistart','identifiability_columns', ...
    'identifiability_pairs','identifiability_reduction','timeseries_pre_baseline', ...
    'timeseries_pre_calibrated','timeseries_post_prediction'};

all_tbl = struct();
for f = 1:numel(files)
    parts = {};
    for r = 1:numel(run_dirs)
        path = fullfile(run_dirs{r}, 'supplement', [files{f} '.csv']);
        if ~exist(path, 'file'), continue; end
        T = readtable(path, 'TextType', 'char');
        T = addvars(T, repmat(seeds(r), height(T), 1), 'Before', 1, 'NewVariableNames', 'Seed');
        parts{end+1} = T; %#ok<AGROW>
    end
    if isempty(parts), continue; end
    T = vertcat(parts{:});
    writetable(T, fullfile(out_dir, [files{f} '_all.csv']));
    all_tbl.(files{f}) = T;
end

% ---- run metadata ------------------------------------------------------
meta = {};
for r = 1:numel(run_dirs)
    path = fullfile(run_dirs{r}, 'supplement', 'run_metadata.json');
    if exist(path, 'file')
        m = jsondecode(fileread(path));
        meta(end+1, :) = {seeds(r), m.run_folder, m.calibration_status, ...
            m.n_calibrated_parameters, m.heart_rate_pre_bpm, m.heart_rate_post_bpm, ...
            m.post_within_10pct, m.post_n_compared, m.post_chi2_per_n, ...
            m.steady_state_pre_calibrated, m.steady_state_post_prediction}; %#ok<AGROW>
    end
end
if ~isempty(meta)
    writetable(cell2table(meta, 'VariableNames', {'Seed','RunFolder','CalibrationStatus', ...
        'NCalibratedParameters','HeartRatePre','HeartRatePost','PostWithin10Pct', ...
        'PostNCompared','PostChi2PerN','SteadyStatePre','SteadyStatePost'}), ...
        fullfile(out_dir, 'run_metadata_all.csv'));
end

% ---- cross-seed summaries -----------------------------------------------
if isfield(all_tbl, 'protocol_comparison_pre')
    write_protocol_summary(all_tbl.protocol_comparison_pre, fullfile(out_dir, 'protocol_comparison_pre_summary.csv'));
end
if isfield(all_tbl, 'protocol_comparison_post')
    write_protocol_summary(all_tbl.protocol_comparison_post, fullfile(out_dir, 'protocol_comparison_post_summary.csv'));
end
if isfield(all_tbl, 'parameters_calibration')
    P = all_tbl.parameters_calibration;
    S = summarise(P, 'Parameter', {'ValueUsed','RatioToScaled'});
    first = P(P.Seed == P.Seed(1), {'Parameter','unit','baseline_scaled','seeded_value','lb','ub','GSASelected','Calibrated'});
    S = outerjoin(first, S, 'Keys', 'Parameter', 'MergeKeys', true);
    writetable(S, fullfile(out_dir, 'parameters_calibration_summary.csv'));
end
if isfield(all_tbl, 'gsa_indices_final')
    G = all_tbl.gsa_indices_final;
    G.Key = strcat(G.Output, '|', G.Parameter);
    S = summarise(G, 'Key', {'S1','ST'});
    parts = split(S.Key, '|');
    S.Output = parts(:, 1); S.Parameter = parts(:, 2); S.Key = [];
    writetable(movevars(S, {'Output','Parameter'}, 'Before', 1), ...
        fullfile(out_dir, 'gsa_indices_final_summary.csv'));
end
fprintf('[aggregate_supplement] %d runs combined into %s\n', numel(run_dirs), out_dir);
end

% =======================================================================
function write_protocol_summary(T, path)
T = T(~cellfun(@isempty, T.ModelMetric), :);
T.Key = strcat(T.Variable, '|', T.Role);
S = summarise(T, 'Key', {'ModelValue','ErrorPct','AbsErrorPct','ModelBaseline','ErrorPctBaseline'});
parts = split(S.Key, '|');
S.Variable = parts(:, 1); S.Role = parts(:, 2);
first = T(T.Seed == T.Seed(1), {'Key','Form','ProtocolRow','Unit','Method','ModelMetric','ComparedValue','ReadingsMean'});
S = outerjoin(first, S, 'Keys', 'Key', 'MergeKeys', true);
w = varfun(@(x) sum(x == 1), T(:, {'Key','Within10Pct'}), 'GroupingVariables', 'Key');
S = outerjoin(S, w(:, {'Key','GroupCount','Fun_Within10Pct'}), 'Keys', 'Key', 'MergeKeys', true);
S = renamevars(S, {'GroupCount','Fun_Within10Pct'}, {'NRuns','NRunsWithin10Pct'});
S.Key = [];
writetable(movevars(S, {'Variable','Role'}, 'Before', 1), path);
end

function S = summarise(T, key, vars)
keys = unique(T.(key), 'stable');
S = table(keys, 'VariableNames', {key});
for v = 1:numel(vars)
    [mu, sd, lo, hi] = deal(nan(numel(keys), 1));
    for k = 1:numel(keys)
        x = T.(vars{v})(strcmp(T.(key), keys{k}));
        x = x(isfinite(x));
        if isempty(x), continue; end
        mu(k) = mean(x); lo(k) = min(x); hi(k) = max(x);
        if numel(x) > 1, sd(k) = std(x); end
    end
    S.([vars{v} '_mean']) = mu;
    S.([vars{v} '_sd']) = sd;
    S.([vars{v} '_min']) = lo;
    S.([vars{v} '_max']) = hi;
end
end
