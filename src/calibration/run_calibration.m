function [params_best, calib_out] = run_calibration(params0, clinical, scenario, optMask, fastMode, pce_surrogate, primaryMetrics)
% RUN_CALIBRATION
% -----------------------------------------------------------------------
% Staged pediatric calibration:
%   Stage A - vascular/shunt fit
%   Stage B - chamber volume/function fit
%   Stage C - joint polish on <=5 sensitive parameters
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  5.0
% -----------------------------------------------------------------------

if nargin < 4, optMask = []; end
if nargin < 5, fastMode = false; end
if nargin < 6, pce_surrogate = []; end
if nargin < 7, primaryMetrics = {}; end

if exist('fmincon', 'file') ~= 2
    error('run_calibration:missingFmincon', ...
        'fmincon is not available. Install Optimization Toolbox to run calibration.');
end

% Optional parallel finite-difference mode for fmincon.
% Enable with: setenv('UNIFIED_VSD_FMINCON_PARALLEL','1')
do_parallel = false;
par_env = getenv('UNIFIED_VSD_FMINCON_PARALLEL');
if ~isempty(par_env)
    do_parallel = any(strcmpi(strtrim(par_env), {'1','true','yes','on'}));
end

calib = calibration_param_sets(scenario, params0, optMask, primaryMetrics);
J0 = objective_calibration(calib.x0, params0, clinical, make_stage_calib(calib, params0, calib.names, calib.metricFields), scenario, pce_surrogate);

params_stage = params0;
stage_history_cell = cell(1, 3);

[params_stage, stage_history_cell{1}] = optimize_stage(params_stage, clinical, scenario, calib, ...
    calib.stageA.names, calib.stageA.metricFields, do_parallel, pce_surrogate, fastMode, 'A');
[params_stage, stage_history_cell{2}] = optimize_stage(params_stage, clinical, scenario, calib, ...
    calib.stageB.names, calib.stageB.metricFields, do_parallel, pce_surrogate, fastMode, 'B');

stageC_names = select_stage_c_names(calib, pce_surrogate, scenario);
[params_best, stage_history_cell{3}] = optimize_stage(params_stage, clinical, scenario, calib, ...
    stageC_names, calib.metricFields, do_parallel, pce_surrogate, fastMode, 'C');
last_stage = stage_history_cell{3};

final_calib = make_stage_calib(calib, params_best, stageC_names, calib.metricFields);
final_x = pack_x(params_best, final_calib.names);
fbest = objective_calibration(final_x, params_best, clinical, final_calib, scenario, pce_surrogate);

calib_out = struct();
calib_out.names = final_calib.names;
calib_out.names_all = calib.names_all;
calib_out.mask = calib.mask;
calib_out.x0 = pack_x(params0, final_calib.names);
calib_out.xbest = final_x;
calib_out.x0_all = calib.x0_all;
calib_out.xbest_all = pack_x(params_best, calib.names_all);
calib_out.J0 = J0;
calib_out.fbest = fbest;
calib_out.lb = final_calib.lb;
calib_out.ub = final_calib.ub;
calib_out.scenario = scenario;
calib_out.primaryMetrics = calib.primaryMetrics;
calib_out.improvement = J0 - fbest;
calib_out.best_stage = 3;
calib_out.best_restart = 0;
calib_out.stage_history = stage_history_cell;
calib_out.exitflag = last_stage.exitflag;
calib_out.output = last_stage.output;
calib_out.objective_breakdown = build_objective_breakdown(params_best, clinical, scenario, calib);
calib_out.use_parallel = do_parallel;

end

function [params_out, stage_out] = optimize_stage(params_in, clinical, scenario, calib, stage_names, stage_metrics, do_parallel, pce_surrogate, fastMode, stage_label)
stage_calib = make_stage_calib(calib, params_in, stage_names, stage_metrics);
params_out = params_in;
stage_out = struct('label', stage_label, 'names', {stage_calib.names}, ...
    'x0', stage_calib.x0, 'xbest', stage_calib.x0, 'fval', NaN, ...
    'exitflag', NaN, 'output', struct(), 'skipped', false);

if isempty(stage_calib.names)
    stage_out.skipped = true;
    return;
end

obj = @(x) objective_calibration(x, params_in, clinical, stage_calib, scenario, pce_surrogate);
f0 = obj(stage_calib.x0);
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'HessianApproximation', 'lbfgs', ...
    'FiniteDifferenceStepSize', 1e-5, ...
    'FiniteDifferenceType', 'forward', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', ternary(fastMode, 1200, 3000), ...
    'MaxIterations', ternary(fastMode, 100, 250), ...
    'OptimalityTolerance', 1e-5, ...
    'StepTolerance', 1e-6, ...
    'UseParallel', do_parallel);

try
    [xbest, fval, exitflag, output] = fmincon( ...
        obj, stage_calib.x0, [], [], [], [], stage_calib.lb, stage_calib.ub, [], opts);
catch ME
    if do_parallel
        fprintf('[run_calibration] Stage %s parallel attempt failed: %s\n', ...
            stage_label, ME.message);
        fprintf('[run_calibration] Stage %s retrying in serial mode.\n', stage_label);
        opts = optimoptions(opts, 'UseParallel', false);
        try
            [xbest, fval, exitflag, output] = fmincon( ...
                obj, stage_calib.x0, [], [], [], [], stage_calib.lb, stage_calib.ub, [], opts);
        catch retryME
            stage_out.skipped = true;
            stage_out.output = struct( ...
                'message', retryME.message, ...
                'parallel_error', ME.message, ...
                'retried_serial', true);
            fprintf('[run_calibration] Stage %s skipped after serial retry: %s\n', ...
                stage_label, retryME.message);
            return;
        end
    else
        stage_out.skipped = true;
        stage_out.output = struct('message', ME.message);
        fprintf('[run_calibration] Stage %s skipped: %s\n', stage_label, ME.message);
        return;
    end
end

if ~isfinite(fval) || fval > f0
    stage_out.skipped = true;
    stage_out.fval = f0;
    stage_out.output = output;
    fprintf('[run_calibration] Stage %s rejected: objective %.6g did not improve over %.6g.\n', ...
        stage_label, fval, f0);
    return;
end

for i = 1:numel(stage_calib.names)
    params_out = set_param_by_name(params_out, stage_calib.names{i}, xbest(i));
end

stage_out.xbest = xbest;
stage_out.fval = fval;
stage_out.exitflag = exitflag;
stage_out.output = output;
end

function stage_calib = make_stage_calib(calib, params_ref, requested_names, stage_metrics)
stage_mask = calib.mask & ismember(calib.names_all, requested_names(:)');
stage_names = calib.names_all(stage_mask);

stage_calib = calib;
stage_calib.names = stage_names;
stage_calib.metricFields = stage_metrics(:)';
stage_calib.x0 = pack_x(params_ref, stage_names);
stage_calib.lb = calib.lb_all(stage_mask);
stage_calib.ub = calib.ub_all(stage_mask);
end

function names = select_stage_c_names(calib, gsa_out, scenario)
active_names = calib.names_all(calib.mask);
if isempty(gsa_out) || ~isstruct(gsa_out)
    names = fallback_stage_c_names(active_names, scenario);
    return;
end

scores = zeros(numel(active_names), 1);
for i = 1:numel(active_names)
    scores(i) = aggregate_parameter_score(active_names{i}, gsa_out, calib.primaryMetrics);
end
[~, order] = sort(scores, 'descend');
keep = min(5, numel(order));
names = active_names(order(1:keep));
if keep == 0
    names = fallback_stage_c_names(active_names, scenario);
end
end

function score = aggregate_parameter_score(param_name, gsa_out, primary_metrics)
score = 0;
for k = 1:numel(primary_metrics)
    mf = primary_metrics{k};
    if ~isfield(gsa_out, mf) || ~isfield(gsa_out.(mf), 'table')
        continue;
    end
    tbl = gsa_out.(mf).table;
    if ~ismember('Parameter', tbl.Properties.VariableNames) || ~ismember('ST', tbl.Properties.VariableNames)
        continue;
    end
    idx = find(strcmp(tbl.Parameter, param_name), 1, 'first');
    if ~isempty(idx)
        score = score + tbl.ST(idx);
    end
end
end

function names = fallback_stage_c_names(active_names, scenario)
switch scenario
    case 'pre_surgery'
        preferred = {'R.vsd','R.PAR','R.SAR','E.LV.EA','V0.LV'};
    otherwise
        preferred = {'R.SAR','R.PAR','C.SAR','E.LV.EA','V0.LV'};
end
names = preferred(ismember(preferred, active_names));
if isempty(names)
    names = active_names(1:min(5, numel(active_names)));
end
end

function breakdown = build_objective_breakdown(params, clinical, scenario, calib)
try
    metrics = compute_clinical_indices(integrate_system(params), params);
catch
    breakdown = table(cell(0,1), cell(0,1), zeros(0,1), zeros(0,1), ...
        'VariableNames', {'Metric','Tier','RelativeError','WeightedContribution'});
    return;
end
targets = get_calibration_targets(scenario, clinical);
target_names = {targets.Metric};

rows = cell(numel(calib.metricFields), 1);
weighted = nan(numel(calib.metricFields), 1);
unweighted = nan(numel(calib.metricFields), 1);
tiers = cell(numel(calib.metricFields), 1);
for k = 1:numel(calib.metricFields)
    mf = calib.metricFields{k};
    rows{k} = mf;
    idx = find(strcmp(target_names, mf), 1, 'first');
    if isempty(idx) || ~isfield(metrics, mf)
        continue;
    end
    y_clin = targets(idx).ClinicalValue;
    if isnan(y_clin)
        continue;
    end
    err_rel = abs(metrics.(mf) - y_clin) / max(abs(y_clin), 1e-6);
    unweighted(k) = err_rel;
    if ismember(mf, calib.primaryMetrics)
        weighted(k) = calib.weights.(mf) * (err_rel / calib.primaryTarget)^2;
        tiers{k} = 'primary';
    else
        weighted(k) = calib.secondaryLambda * calib.weights.(mf) * (err_rel / calib.secondaryTarget)^2;
        tiers{k} = 'secondary';
    end
end

breakdown = table(rows, tiers, unweighted, weighted, ...
    'VariableNames', {'Metric','Tier','RelativeError','WeightedContribution'});
end

function x = pack_x(params, names)
x = zeros(numel(names), 1);
for i = 1:numel(names)
    x(i) = get_param_by_name(params, names{i});
end
end

function v = get_param_by_name(params, name)
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts)
    v = v.(parts{k});
end
end

function params = set_param_by_name(params, name, value)
parts = strsplit(name, '.');
switch numel(parts)
    case 2
        params.(parts{1}).(parts{2}) = value;
    case 3
        params.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_param_by_name:unsupportedDepth', ...
            'Only 2- or 3-level dot notation supported. Got: %s', name);
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
