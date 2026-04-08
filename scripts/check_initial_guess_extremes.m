%% check_initial_guess_extremes.m
% Robustness check: run calibration from multiple extreme x0 perturbations.
% This script helps answer whether optimisation performance is sensitive to
% initial guess selection in calibration_param_sets.

clear; clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
project_paths = strsplit(genpath(project_root), pathsep);
is_shadow = contains(project_paths, [filesep '.claude' filesep]);
addpath(strjoin(project_paths(~is_shadow), pathsep));

% Choose benchmark profile/scenario (dummy/synthetic-compatible).
scenario = 'pre_surgery';
clinical = patient_profile_A();

% Build baseline params exactly as in main pipeline.
params_ref = default_parameters();
patient.age_years  = clinical.common.age_years;
patient.weight_kg  = clinical.common.weight_kg;
patient.height_cm  = clinical.common.height_cm;
patient.sex        = clinical.common.sex;
if isfield(clinical.common, 'BSA') && ~isnan(clinical.common.BSA)
	patient.BSA = clinical.common.BSA;
end
params0 = apply_scaling(params_ref, patient);
params0 = params_from_clinical(params0, clinical, scenario);

% Calibration config and objective.
calib = calibration_param_sets(scenario, params0);
obj = @(x) safe_objective(x, params0, clinical, calib, scenario);

if exist('fmincon', 'file') ~= 2
	error('check_initial_guess_extremes:missingFmincon', ...
		'fmincon is required for this script.');
end

opts = optimoptions('fmincon', ...
	'Algorithm', 'sqp', ...
	'UseParallel', false, ...
	'FiniteDifferenceStepSize', 1e-5, ...
	'FiniteDifferenceType', 'forward', ...
	'Display', 'none', ...
	'MaxFunctionEvaluations', 1200, ...
	'MaxIterations', 80, ...
	'OptimalityTolerance', 1e-6, ...
	'StepTolerance', 1e-8);

% Extreme x0 scales to stress-test convergence.
scale_set = [0.30, 0.50, 1.00, 1.80, 2.50];
n = numel(scale_set);

rows = strings(n, 1);
J_start = zeros(n, 1);
J_best = zeros(n, 1);
exitflag = zeros(n, 1);
improvement_pct = zeros(n, 1);
status = strings(n, 1);

for i = 1:n
	s = scale_set(i);
	rows(i) = sprintf('x0_scale=%.2f', s);
	x0_try = calib.x0 .* s;
	x0_try = min(max(x0_try, calib.lb), calib.ub);

	try
		J0 = obj(x0_try);

		if J0 >= 1e8
			xbest = x0_try;
			fbest = J0;
			ef = -998;
			status(i) = "skipped:invalid_start";
			fprintf('[%s] SKIPPED | invalid start objective (J0=%.3e)\n', rows(i), J0);
		else
			[xbest, fbest, ef] = fmincon(obj, x0_try, [], [], [], [], calib.lb, calib.ub, [], opts);
			status(i) = "ok";
		end
	catch ME
		J0 = 1e9;
		xbest = x0_try;
		fbest = 1e9;
		ef = -999;
		status(i) = "failed:" + string(ME.identifier);
		fprintf('[%s] FAILED | %s\n', rows(i), ME.message);
	end

	J_start(i) = J0;
	J_best(i) = fbest;
	exitflag(i) = ef;
	improvement_pct(i) = 100 * (J0 - fbest) / max(abs(J0), 1e-12);

	fprintf('[%s] J0=%.6f -> J*=%.6f | exitflag=%d | improve=%.2f%%\n', ...
		rows(i), J0, fbest, ef, improvement_pct(i));

	%#ok<NASGU>
	if i == 1
		xbest_ref = xbest;
	end
end

T = table(rows, J_start, J_best, exitflag, improvement_pct, status, ...
	'VariableNames', {'Case', 'J0', 'Jbest', 'Exitflag', 'Improve_pct', 'Status'});
disp(T);

out_dir = fullfile(project_root, 'results', 'tables');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_file = fullfile(out_dir, sprintf('initial_guess_extremes_%s.mat', scenario));
save(out_file, 'T', 'scale_set', 'calib');
fprintf('Saved robustness summary: %s\n', out_file);

function J = safe_objective(x, params0, clinical, calib, scenario)
% SAFE_OBJECTIVE - objective wrapper that prevents single-run crashes.
try
	J = objective_calibration(x, params0, clinical, calib, scenario);
	if ~isfinite(J)
		J = 1e9;
	end
catch
	J = 1e9;
end
end
