function polish_post_surgery_direct_fit()
% POLISH_POST_SURGERY_DIRECT_FIT
% -----------------------------------------------------------------------
% Direct local polishing pass for the Reyna post-surgery fit.
%
% This script is intentionally separate from main_run.m. It explores whether
% the current post-surgery target set can be fit to RMSE 0.08-0.10 before
% changing the production calibration workflow.
%
% INPUTS:
%   Uses the latest reyna post-surgery run package in results/runs/.
%
% OUTPUTS:
%   results/runs/<latest>/
%       mat/post_surgery_direct_polish.mat
%       tables/post_surgery_direct_polish_validation.csv
%
% ASSUMPTIONS:
%   - VSD remains surgically closed through params.R.vsd = 1e6.
%   - The polished parameters are exploratory until reviewed.
%
% SIGN CONVENTIONS:
%   - No new flow sign convention is introduced.
%
% REFERENCES:
%   [1] docs/theory_notes.md. Scenario definitions and validation targets.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-18
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));

run_dir = find_latest_post_surgery_run(root);
package_file = newest_file(fullfile(run_dir, 'mat', 'run_package_post_surgery_*.mat'));
best_file = fullfile(run_dir, 'mat', 'params_best_candidate_post_surgery.mat');
random_file = fullfile(run_dir, 'mat', 'random_search_post_tmp.mat');

loaded_package = load(package_file, 'run_package');
loaded_best = load(best_file, 'best_candidate');
clinical = loaded_package.run_package.clinical;       % [-]
scenario = 'post_surgery';                            % [-]

if isfile(random_file)
    loaded_random = load(random_file, 'best_params');
    params_start = loaded_random.best_params;          % [-]
else
    params_start = loaded_best.best_candidate.params;  % [-]
end

targets = get_calibration_targets(scenario, clinical);
target_names = {targets.Metric};
target_values = [targets.ClinicalValue];
target_mask = isfinite(target_values);
target_names = target_names(target_mask);
target_values = target_values(target_mask);

x0 = pack_params(params_start);
[lb, ub] = polishing_bounds(x0);

objective = @(x) score_candidate(unpack_params(params_start, x), ...
    target_names, target_values);
score0 = objective(x0);
fprintf('[polish_post_surgery_direct_fit] Start score = %.4f\n', score0);

opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 700, ...
    'MaxIterations', 90, ...
    'OptimalityTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'FiniteDifferenceType', 'forward');

[xbest, fbest, exitflag, output] = fmincon(objective, x0, [], [], [], [], lb, ub, [], opts);

params_polished = unpack_params(params_start, xbest);
[score_polished, metrics_polished, sim_polished, validation_table] = ...
    score_candidate(params_polished, target_names, target_values);

fprintf('[polish_post_surgery_direct_fit] Final score = %.4f\n', score_polished);
fprintf('[polish_post_surgery_direct_fit] CO=%.3f L/min | MAP=%.2f mmHg | LVEF=%.3f | LVEDV=%.2f mL\n', ...
    metrics_polished.CO_Lmin, metrics_polished.SAP_mean, ...
    metrics_polished.LVEF, metrics_polished.LVEDV);

out_file = fullfile(run_dir, 'mat', 'post_surgery_direct_polish.mat');
save(out_file, 'params_polished', 'metrics_polished', 'sim_polished', ...
    'validation_table', 'score0', 'score_polished', 'x0', 'xbest', ...
    'fbest', 'exitflag', 'output', 'target_names', 'target_values');
writetable(validation_table, fullfile(run_dir, 'tables', ...
    'post_surgery_direct_polish_validation.csv'));

fprintf('[polish_post_surgery_direct_fit] Saved: %s\n', out_file);
end

function run_dir = find_latest_post_surgery_run(root)
% FIND_LATEST_POST_SURGERY_RUN - newest Reyna post-op run folder.
listing = dir(fullfile(root, 'results', 'runs', '*_reyna_post_surgery'));
if isempty(listing)
    error('polish_post_surgery_direct_fit:noPostRun', ...
        'No reyna post-surgery run folder was found.');
end
[~, order] = sort([listing.datenum], 'descend');
run_dir = fullfile(listing(order(1)).folder, listing(order(1)).name);
end

function file_path = newest_file(pattern)
% NEWEST_FILE - newest file matching a wildcard pattern.
listing = dir(pattern);
if isempty(listing)
    error('polish_post_surgery_direct_fit:noFile', ...
        'No file matched pattern: %s', pattern);
end
[~, order] = sort([listing.datenum], 'descend');
file_path = fullfile(listing(order(1)).folder, listing(order(1)).name);
end

function x = pack_params(params)
% PACK_PARAMS - transform positive parameters to log10 space where useful.
x = [
    log10(params.R.SAR)
    log10(params.R.SVEN)
    log10(params.R.PAR)
    log10(params.C.SAR)
    log10(params.C.PAR)
    log10(params.E.LV.EA)
    log10(params.E.LV.EB)
    params.V0.LV
    log10(params.E.RV.EA)
    log10(params.E.RV.EB)
    params.V0.RV
    log10(params.E.RA.EA)
    log10(params.E.LA.EA)
    ];
end

function [lb, ub] = polishing_bounds(x0)
% POLISHING_BOUNDS - broad exploratory bounds around physiologic positive ranges.
lb = x0;
ub = x0;

lb(1) = log10(0.05);  ub(1) = log10(2.50);  % R.SAR [mmHg*s/mL]
lb(2) = log10(0.02);  ub(2) = log10(0.90);  % R.SVEN [mmHg*s/mL]
lb(3) = log10(0.02);  ub(3) = log10(1.50);  % R.PAR [mmHg*s/mL]
lb(4) = log10(0.25);  ub(4) = log10(1.80);  % C.SAR [mL/mmHg]
lb(5) = log10(0.80);  ub(5) = log10(5.00);  % C.PAR [mL/mmHg]
lb(6) = log10(8.00);  ub(6) = log10(80.0);  % E.LV.EA [mmHg/mL]
lb(7) = log10(0.03);  ub(7) = log10(1.20);  % E.LV.EB [mmHg/mL]
lb(8) = 0.0;          ub(8) = 22.0;         % V0.LV [mL]
lb(9) = log10(1.00);  ub(9) = log10(10.0);  % E.RV.EA [mmHg/mL]
lb(10) = log10(0.03); ub(10) = log10(0.90); % E.RV.EB [mmHg/mL]
lb(11) = 2.0;         ub(11) = 22.0;        % V0.RV [mL]
lb(12) = log10(0.05); ub(12) = log10(3.00); % E.RA.EA [mmHg/mL]
lb(13) = log10(0.05); ub(13) = log10(3.00); % E.LA.EA [mmHg/mL]
end

function params = unpack_params(params_template, x)
% UNPACK_PARAMS - write candidate vector back into a parameter struct.
params = params_template;

R_SAR_base = params_template.R.SAR;          % [mmHg*s/mL]
R_SC_base = params_template.R.SC;            % [mmHg*s/mL]
R_PAR_base = params_template.R.PAR;          % [mmHg*s/mL]
R_PCOX_base = params_template.R.PCOX;        % [mmHg*s/mL]
R_PCNO_base = params_template.R.PCNO;        % [mmHg*s/mL]
R_PVEN_base = params_template.R.PVEN;        % [mmHg*s/mL]

params.R.SAR = 10 ^ x(1);                    % [mmHg*s/mL]
params.R.SC = R_SC_base * params.R.SAR / max(R_SAR_base, 1e-12);
params.R.SVEN = 10 ^ x(2);                   % [mmHg*s/mL]

params.R.PAR = 10 ^ x(3);                    % [mmHg*s/mL]
pulmonary_scale = params.R.PAR / max(R_PAR_base, 1e-12);
params.R.PCOX = R_PCOX_base * pulmonary_scale;
params.R.PCNO = R_PCNO_base * pulmonary_scale;
params.R.PVEN = R_PVEN_base * pulmonary_scale;

params.C.SAR = 10 ^ x(4);                    % [mL/mmHg]
params.C.PAR = 10 ^ x(5);                    % [mL/mmHg]
params.E.LV.EA = 10 ^ x(6);                  % [mmHg/mL]
params.E.LV.EB = 10 ^ x(7);                  % [mmHg/mL]
params.V0.LV = x(8);                         % [mL]
params.E.RV.EA = 10 ^ x(9);                  % [mmHg/mL]
params.E.RV.EB = 10 ^ x(10);                 % [mmHg/mL]
params.V0.RV = x(11);                        % [mL]
params.E.RA.EA = 10 ^ x(12);                 % [mmHg/mL]
params.E.LA.EA = 10 ^ x(13);                 % [mmHg/mL]
params.R.vsd = 1e6;                          % [mmHg*s/mL]
end

function [score, metrics, sim, validation_table] = score_candidate(params, target_names, target_values)
% SCORE_CANDIDATE - RMSE over finite post-surgery validation targets.
sim = integrate_system(params);
metrics = compute_clinical_indices(sim, params);

model_values = nan(numel(target_names), 1);
error_frac = nan(numel(target_names), 1);
for target_idx = 1:numel(target_names)
    metric_name = target_names{target_idx};
    if isfield(metrics, metric_name)
        model_values(target_idx) = metrics.(metric_name);
        error_frac(target_idx) = (model_values(target_idx) - target_values(target_idx)) / ...
            max(abs(target_values(target_idx)), 1e-9);
    end
end

finite_error = error_frac(isfinite(error_frac));
if isempty(finite_error)
    score = 1e6;
else
    score = sqrt(mean(finite_error .^ 2));
end

validation_table = table(target_names(:), target_values(:), model_values(:), ...
    100 * error_frac(:), 'VariableNames', ...
    {'Metric','Clinical','Model','Error_pct'});

if ~sim.ss_reached
    score = score + 10;
end
end
