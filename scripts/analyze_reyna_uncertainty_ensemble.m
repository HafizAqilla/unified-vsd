function analyze_reyna_uncertainty_ensemble(csv_path)
% ANALYZE_REYNA_UNCERTAINTY_ENSEMBLE
% -----------------------------------------------------------------------
% Reviews the latest Reyna uncertainty ensemble and classifies selected
% metrics as prior-sensitive or likely structural, based on ensemble
% spread and correlation to the uncertainty knobs.
%
% INPUTS:
%   csv_path  - optional path to a reyna_uncertainty_ensemble_*.csv file  [-]
%
% OUTPUTS:
%   results/tables/reyna_uncertainty_analysis_<timestamp>.csv
%   results/tables/reyna_uncertainty_analysis_<timestamp>.txt
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-08
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));
results_dir = fullfile(root, 'results', 'tables');

if nargin < 1 || isempty(csv_path)
    listing = dir(fullfile(results_dir, 'reyna_uncertainty_ensemble_*.csv'));
    if isempty(listing)
        error('analyze_reyna_uncertainty_ensemble:noEnsemble', ...
            'No Reyna uncertainty ensemble CSV found. Run run_reyna_uncertainty_ensemble first.');
    end
    [~, order] = sort([listing.datenum], 'descend');
    csv_path = fullfile(listing(order(1)).folder, listing(order(1)).name);
end

ensemble_tbl = readtable(csv_path);
clinical = patient_reyna();
target_struct = resolve_target_struct(clinical, 'pre_surgery');

metric_names = {'CO_Lmin','SVR','PAP_min'};
n_metric = numel(metric_names);

metric_col = strings(n_metric, 1);
target_col = nan(n_metric, 1);
mean_col = nan(n_metric, 1);
min_col = nan(n_metric, 1);
max_col = nan(n_metric, 1);
spread_pct_col = nan(n_metric, 1);
median_abs_err_col = nan(n_metric, 1);
min_abs_err_col = nan(n_metric, 1);
dom_knob_col = strings(n_metric, 1);
dom_corr_col = nan(n_metric, 1);
assessment_col = strings(n_metric, 1);

knob_names = {'BVScale','SVRScale','PVRScale'};

for idx = 1:n_metric
    metric_name = metric_names{idx};
    values = ensemble_tbl.(metric_name);
    target = target_struct.(metric_name);
    abs_err_pct = 100 * abs(values - target) ./ max(abs(target), 1e-9);

    corr_values = nan(numel(knob_names), 1);
    for knob_idx = 1:numel(knob_names)
        knob_vals = ensemble_tbl.(knob_names{knob_idx});
        corr_values(knob_idx) = safe_corr(knob_vals, values);
    end
    [max_abs_corr, max_idx] = max(abs(corr_values));

    metric_col(idx) = string(metric_name);
    target_col(idx) = target;
    mean_col(idx) = mean(values, 'omitnan');
    min_col(idx) = min(values);
    max_col(idx) = max(values);
    spread_pct_col(idx) = 100 * (max(values) - min(values)) / max(abs(target), 1e-9);
    median_abs_err_col(idx) = median(abs_err_pct, 'omitnan');
    min_abs_err_col(idx) = min(abs_err_pct);
    dom_knob_col(idx) = string(knob_names{max_idx});
    dom_corr_col(idx) = corr_values(max_idx);
    assessment_col(idx) = classify_metric_behavior( ...
        min_abs_err_col(idx), spread_pct_col(idx), max_abs_corr);
end

analysis_tbl = table(metric_col, target_col, mean_col, min_col, max_col, ...
    spread_pct_col, median_abs_err_col, min_abs_err_col, dom_knob_col, ...
    dom_corr_col, assessment_col, ...
    'VariableNames', {'Metric','ClinicalTarget','EnsembleMean','EnsembleMin', ...
    'EnsembleMax','SpreadPctOfTarget','MedianAbsErrorPct','MinAbsErrorPct', ...
    'DominantKnob','DominantCorrelation','Assessment'});

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_csv = fullfile(results_dir, sprintf('reyna_uncertainty_analysis_%s.csv', timestamp));
out_txt = fullfile(results_dir, sprintf('reyna_uncertainty_analysis_%s.txt', timestamp));
writetable(analysis_tbl, out_csv);

fid = fopen(out_txt, 'w');
if fid < 0
    error('analyze_reyna_uncertainty_ensemble:openFailed', ...
        'Unable to write uncertainty analysis summary.');
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Reyna Uncertainty Ensemble Analysis\n');
fprintf(fid, '==================================\n');
fprintf(fid, 'Source CSV: %s\n\n', csv_path);
for idx = 1:height(analysis_tbl)
    fprintf(fid, '%s\n', analysis_tbl.Metric(idx));
    fprintf(fid, '  Target: %.4f\n', analysis_tbl.ClinicalTarget(idx));
    fprintf(fid, '  Range: %.4f to %.4f\n', analysis_tbl.EnsembleMin(idx), analysis_tbl.EnsembleMax(idx));
    fprintf(fid, '  SpreadPctOfTarget: %.2f\n', analysis_tbl.SpreadPctOfTarget(idx));
    fprintf(fid, '  MedianAbsErrorPct: %.2f\n', analysis_tbl.MedianAbsErrorPct(idx));
    fprintf(fid, '  MinAbsErrorPct: %.2f\n', analysis_tbl.MinAbsErrorPct(idx));
    fprintf(fid, '  DominantKnob: %s (corr=%.3f)\n', ...
        analysis_tbl.DominantKnob(idx), analysis_tbl.DominantCorrelation(idx));
    fprintf(fid, '  Assessment: %s\n\n', analysis_tbl.Assessment(idx));
end

fprintf('[analyze_reyna_uncertainty_ensemble] Analysis exported:\n  %s\n  %s\n', ...
    out_csv, out_txt);
disp(analysis_tbl);
end

function r = safe_corr(x, y)
valid = isfinite(x) & isfinite(y);
if nnz(valid) < 3
    r = NaN;
    return;
end
c = corrcoef(x(valid), y(valid));
if numel(c) < 4
    r = NaN;
else
    r = c(1, 2);
end
end

function label = classify_metric_behavior(min_abs_err_pct, spread_pct, max_abs_corr)
if min_abs_err_pct <= 5
    label = "prior_sensitive_recoverable";
elseif min_abs_err_pct >= 20 && spread_pct < 5
    label = "likely_structural_or_strongly_constrained";
elseif spread_pct >= 10 || max_abs_corr >= 0.5
    label = "prior_sensitive";
else
    label = "likely_structural";
end
end

function target_struct = resolve_target_struct(clinical, scenario)
targets = get_calibration_targets(scenario, clinical);
metric_names = {'CO_Lmin','SVR','PAP_min'};
target_struct = struct();
for idx = 1:numel(metric_names)
    metric_name = metric_names{idx};
    match_idx = find(strcmp({targets.Metric}, metric_name), 1, 'first');
    if isempty(match_idx)
        error('analyze_reyna_uncertainty_ensemble:missingTarget', ...
            'Metric target mapping for %s was not found.', metric_name);
    end
    target_struct.(metric_name) = targets(match_idx).ClinicalValue;
end
end

function project_path = build_clean_project_path(root)
project_paths = strsplit(genpath(root), pathsep);
project_paths = project_paths(~cellfun('isempty', project_paths));
is_shadow = contains(project_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, project_paths);
project_path = strjoin(project_paths(~is_shadow & is_existing), pathsep);
end
