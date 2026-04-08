%% check_sobol_ranking_stability.m
% Batch 3 utility: compare Sobol ST ranking stability across sample sizes.
%
% Default workflow:
%   1) Build pre-surgery patient A parameters
%   2) Run Sobol GSA at N_small and N_large
%   3) Compare parameter ranking by mean ST over primary metrics
%   4) Report Spearman rho and top-k overlap

clear; clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

scenario = 'pre_surgery';
clinical = patient_profile_A();

% Batch 3 comparison pair. Adjust as needed.
N_small = 512;
N_large = 1024;
top_k = 8;

fprintf('[check_sobol_ranking_stability] Scenario=%s | N_small=%d | N_large=%d\n', ...
	scenario, N_small, N_large);

% Build baseline params exactly as main pipeline.
params_ref = default_parameters();
patient.age_years = clinical.common.age_years;
patient.weight_kg = clinical.common.weight_kg;
patient.height_cm = clinical.common.height_cm;
patient.sex = clinical.common.sex;
if isfield(clinical.common, 'BSA') && ~isnan(clinical.common.BSA)
	patient.BSA = clinical.common.BSA;
end

params0 = apply_scaling(params_ref, patient);
params0 = params_from_clinical(params0, clinical, scenario);

% Run GSA at two N values.
cfg_small = gsa_sobol_setup(params0, scenario, N_small);
out_small = gsa_run_sobol(cfg_small, params0);

cfg_large = gsa_sobol_setup(params0, scenario, N_large);
out_large = gsa_run_sobol(cfg_large, params0);

% Build aggregate ranking from primary metrics (mean ST across primary outputs).
[names_small, score_small] = aggregate_primary_st(cfg_small, out_small);
[names_large, score_large] = aggregate_primary_st(cfg_large, out_large);

if ~isequal(names_small, names_large)
	error('check_sobol_ranking_stability:nameMismatch', ...
		'Parameter name order mismatch between N_small and N_large outputs.');
end

all_names = names_small;

[score_small_sorted, order_small] = sort(score_small, 'descend');
[score_large_sorted, order_large] = sort(score_large, 'descend');

rank_small = zeros(size(order_small));
rank_large = zeros(size(order_large));
rank_small(order_small) = 1:numel(order_small);
rank_large(order_large) = 1:numel(order_large);

rho = corr(rank_small(:), rank_large(:), 'Type', 'Spearman', 'Rows', 'complete');

top_k = min(top_k, numel(all_names));
top_small = all_names(order_small(1:top_k));
top_large = all_names(order_large(1:top_k));
overlap_count = numel(intersect(top_small, top_large));
overlap_frac = overlap_count / max(top_k, 1);

fprintf('\n--- Batch 3 Ranking Stability Summary ---\n');
fprintf('Spearman rank correlation: %.4f\n', rho);
fprintf('Top-%d overlap: %d/%d (%.1f%%)\n', top_k, overlap_count, top_k, 100 * overlap_frac);

T_small = table(all_names(order_small), score_small_sorted, rank_small(order_small), ...
	'VariableNames', {'Parameter', 'MeanST_Nsmall', 'Rank_Nsmall'});
T_large = table(all_names(order_large), score_large_sorted, rank_large(order_large), ...
	'VariableNames', {'Parameter', 'MeanST_Nlarge', 'Rank_Nlarge'});

T_compare = table(all_names, score_small, score_large, rank_small, rank_large, abs(rank_small - rank_large), ...
	'VariableNames', {'Parameter', 'MeanST_Nsmall', 'MeanST_Nlarge', 'Rank_Nsmall', 'Rank_Nlarge', 'AbsRankDelta'});
T_compare = sortrows(T_compare, 'Rank_Nlarge', 'ascend');

disp(T_compare(1:min(15, height(T_compare)), :));

out_dir = fullfile(project_root, 'results', 'tables');
if ~exist(out_dir, 'dir')
	mkdir(out_dir);
end

save(fullfile(out_dir, sprintf('sobol_rank_stability_%s_N%d_vs_N%d.mat', scenario, N_small, N_large)), ...
	'scenario', 'N_small', 'N_large', 'top_k', 'rho', 'overlap_count', 'overlap_frac', ...
	'T_small', 'T_large', 'T_compare', 'cfg_small', 'cfg_large');

fprintf('\nSaved Batch 3 stability artifact in results/tables/.\n');

function [names, mean_st] = aggregate_primary_st(cfg, gsa_out)
% AGGREGATE_PRIMARY_ST - mean total-order index over scenario primary metrics.
names = cfg.names(:);
n_param = numel(names);
n_primary = numel(cfg.primary_metrics);
ST = zeros(n_param, n_primary);

for i = 1:n_primary
	mf = cfg.primary_metrics{i};
	if ~isfield(gsa_out, mf) || ~isfield(gsa_out.(mf), 'ST')
		error('check_sobol_ranking_stability:missingMetric', ...
			'Missing ST for primary metric: %s', mf);
	end
	st_i = gsa_out.(mf).ST(:);
	if numel(st_i) ~= n_param
		error('check_sobol_ranking_stability:lengthMismatch', ...
			'ST length mismatch for metric %s.', mf);
	end
	ST(:, i) = st_i;
end

mean_st = mean(ST, 2, 'omitnan');
end
