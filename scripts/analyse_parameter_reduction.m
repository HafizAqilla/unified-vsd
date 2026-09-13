function out = analyse_parameter_reduction(run_dir, scenario, max_report)
% ANALYSE_PARAMETER_REDUCTION
% -----------------------------------------------------------------------
% Which parameters can be dropped, and how much does dropping them help?
%
% Motivation
% ----------
% The binding limitation on every claim from this model is that the fit is
% underdetermined: N = 9 governed observations against p = 12 free
% parameters gives dof = -3, so a low chi2/N is guaranteed and proves
% nothing (results doc §4, §9).
%
% There are two levers. Raising N needs the joint pre/post fit. Lowering p
% has never been tried, and the identifiability report already says which
% parameters are redundant. This script turns that report into a concrete
% reduced set.
%
% Method
% ------
% The scaled sensitivity matrix S (n_metrics x n_parameters) is built ONCE at
% the calibrated operating point, then subsets are scored by cond(S(:,idx)).
% Because the columns are fixed, every subset is evaluated by linear algebra
% alone -- no re-simulation and no recalibration. This makes an exhaustive
% search over subset sizes cheap, where a calibration-based search would cost
% hours per candidate.
%
% What the score means
% --------------------
% cond(S) is the amplification of parameter uncertainty implied by the
% metrics: low is identifiable, > 1e3 is flagged near-dependent. Column norm
% is how strongly a parameter moves the metrics at all; a near-zero column is
% a parameter the data cannot constrain in any combination.
%
% Selecting a subset by conditioning alone is not sufficient on its own -- a
% dropped parameter is FIXED at its calibrated value, which is a modelling
% commitment, not a free win. The reduced set must then be re-run and its
% gate count compared. This script narrows the candidates; it does not
% license skipping that run.
%
% INPUTS:
%   run_dir    - a completed run folder containing mat/ and tables/         [-]
%   scenario   - scenario string (default 'pre_surgery')                    [-]
%   max_report - how many best subsets to print per size (default 3)        [-]
%
% OUTPUTS:
%   out.full_cond    - cond(S) over all parameters                          [-]
%   out.by_size      - best subset found at each size, with cond and dof    [-]
%   out.column_norms - per-parameter column norm table                      [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-30
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(scenario); scenario = 'pre_surgery'; end
if nargin < 3 || isempty(max_report); max_report = 3; end

gate_csv = fullfile(run_dir, 'tables', ...
    sprintf('full_metric_gate_%s.csv', scenario));
cand_mat = fullfile(run_dir, 'mat', ...
    sprintf('params_accepted_candidate_%s.mat', scenario));
assert(exist(gate_csv, 'file') == 2, 'missing %s', gate_csv);
assert(exist(cand_mat, 'file') == 2, 'missing %s', cand_mat);

gate_table = readtable(gate_csv);
S_in = load(cand_mat);
[params_best, calib_out] = unpack(S_in, run_dir, scenario);

report = analyse_parameter_identifiability(params_best, scenario, ...
    calib_out, gate_table, '');
S = report.sensitivity_matrix;
names = report.parameter_table.Parameter;

keep = all(isfinite(S), 1);
if ~all(keep)
    fprintf('  dropping %d parameter(s) with non-finite sensitivity\n', ...
        nnz(~keep));
    S = S(:, keep);
    names = names(keep);
end

n_metrics = size(S, 1);
p_full = size(S, 2);
out = struct();
out.full_cond = cond(S);
out.n_metrics = n_metrics;
out.column_norms = table(names, vecnorm(S, 2, 1)', ...
    'VariableNames', {'Parameter', 'ColumnNorm'});
out.column_norms = sortrows(out.column_norms, 'ColumnNorm');

fprintf('\n=== PARAMETER REDUCTION ANALYSIS (%s) ===\n', scenario);
fprintf('  metrics N = %d, full parameter set p = %d\n', n_metrics, p_full);
fprintf('  cond(S) over the full set : %.3g   (dof = %d)\n', ...
    out.full_cond, n_metrics - p_full);
fprintf('\n  Column norms (weakest first -- these constrain the fit least):\n');
for i = 1:height(out.column_norms)
    fprintf('    %-22s %10.4f\n', out.column_norms.Parameter{i}, ...
        out.column_norms.ColumnNorm(i));
end

% Exhaustive over subset sizes where it is cheap, greedy above that.
rows = {};
for k = min(p_full, n_metrics):-1:2
    if nchoosek(p_full, k) <= 5000
        combos = nchoosek(1:p_full, k);
    else
        combos = greedy_subsets(S, k, 400);
    end
    best_c = inf; best_idx = [];
    scored = nan(size(combos, 1), 1);
    for r = 1:size(combos, 1)
        c = cond(S(:, combos(r, :)));
        scored(r) = c;
        if c < best_c; best_c = c; best_idx = combos(r, :); end
    end
    [~, order] = sort(scored);
    fprintf('\n  p = %d  (dof = %+d)  best cond = %.4g\n', ...
        k, n_metrics - k, best_c);
    for r = 1:min(max_report, numel(order))
        idx = combos(order(r), :);
        fprintf('     %.4g : %s\n', scored(order(r)), ...
            strjoin(names(idx)', ', '));
    end
    rows{end+1} = struct('p', k, 'dof', n_metrics - k, ...
        'cond', best_c, 'names', {names(best_idx)'}); %#ok<AGROW>
end
out.by_size = rows;
end

% =======================================================================
function [params_best, calib_out] = unpack(S_in, run_dir, scenario)
% The accepted-candidate package stores the parameter struct and the case
% profile, but not the active-parameter name list. The names are recovered
% from the identifiability CSV the same run exported, which is by
% construction exactly the set that run treated as free -- so the reduction
% analysis is scored over the same p the run reported.
params_best = [];
if isfield(S_in, 'accepted_candidate') && isfield(S_in.accepted_candidate, 'params')
    params_best = S_in.accepted_candidate.params;
end
assert(~isempty(params_best), 'no accepted_candidate.params in candidate mat');

calib_out = struct();
if isfield(S_in, 'case_profile')
    calib_out.caseProfile = S_in.case_profile;
end

ident_csv = fullfile(run_dir, 'tables', ...
    sprintf('parameter_identifiability_%s.csv', scenario));
assert(exist(ident_csv, 'file') == 2, ...
    'missing %s (needed for the active parameter names)', ident_csv);
tbl = readtable(ident_csv);
calib_out.names = tbl.Parameter(:)';
end

% =======================================================================
function combos = greedy_subsets(S, k, n_seeds)
% Backward elimination from random seeds: drop the column whose removal
% leaves the best-conditioned remainder, until k columns remain.
p = size(S, 2);
combos = zeros(n_seeds, k);
for s = 1:n_seeds
    idx = 1:p;
    if s > 1
        idx = idx(randperm(p));
        idx = sort(idx(1:max(k, p - mod(s, 3))));
    end
    while numel(idx) > k
        best_c = inf; drop_at = 1;
        for j = 1:numel(idx)
            trial = idx; trial(j) = [];
            c = cond(S(:, trial));
            if c < best_c; best_c = c; drop_at = j; end
        end
        idx(drop_at) = [];
    end
    combos(s, :) = sort(idx);
end
combos = unique(combos, 'rows');
end
