function [params_best, calib_out] = run_calibration(params0, clinical, scenario, optMask, fastMode, pce_surrogate, primaryMetrics, caseProfile, registryContext)
% RUN_CALIBRATION
% -----------------------------------------------------------------------
% Staged pediatric calibration:
%   Stage A - vascular/shunt fit
%   Stage B - chamber volume/function fit
%   Stage C - joint polish on <=5 sensitive parameters
%   Stage D - optional systemic-output polish for adaptive patient cases
%   Stage E - optional plausibility polish without losing validation fit
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  5.2
% -----------------------------------------------------------------------

if nargin < 4, optMask = []; end
if nargin < 5, fastMode = false; end
if nargin < 6, pce_surrogate = []; end
if nargin < 7, primaryMetrics = {}; end
if nargin < 8 || isempty(caseProfile), caseProfile = struct(); end
if nargin < 9 || isempty(registryContext), registryContext = struct(); end

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

calib = calibration_param_sets(scenario, params0, optMask, primaryMetrics, caseProfile, registryContext);
J0 = objective_calibration(calib.x0, params0, clinical, make_stage_calib(calib, params0, calib.names, calib.metricFields), scenario, pce_surrogate);

params_stage = params0;
stage_history_cell = cell(1, 5);

[params_stage, stage_history_cell{1}] = optimize_stage(params_stage, clinical, scenario, calib, ...
    calib.stageA.names, calib.stageA.metricFields, do_parallel, pce_surrogate, fastMode, 'A');
[params_stage, stage_history_cell{2}] = optimize_stage(params_stage, clinical, scenario, calib, ...
    calib.stageB.names, calib.stageB.metricFields, do_parallel, pce_surrogate, fastMode, 'B');

stageC_names = select_stage_c_names(calib, pce_surrogate, scenario);
[params_best, stage_history_cell{3}] = optimize_stage(params_stage, clinical, scenario, calib, ...
    stageC_names, calib.metricFields, do_parallel, pce_surrogate, fastMode, 'C');
last_stage = stage_history_cell{3};
final_stage_names = stageC_names;
final_stage_metrics = calib.metricFields;
final_calib_source = calib;
best_stage = 3;

if should_run_systemic_polish(calib, scenario)
    stageD_calib = apply_systemic_polish_profile(calib);
    stageD_names = select_systemic_polish_names(stageD_calib);
    stageD_metrics = select_systemic_polish_metrics(stageD_calib);
    [params_polished, stage_history_cell{4}] = optimize_stage(params_best, clinical, scenario, stageD_calib, ...
        stageD_names, stageD_metrics, do_parallel, pce_surrogate, fastMode, 'D');

    if ~stage_history_cell{4}.skipped
        [accept_polish, polish_acceptance] = accept_systemic_polish( ...
            params_best, params_polished, clinical, scenario, calib, pce_surrogate);
        stage_history_cell{4}.acceptance = polish_acceptance;
    else
        accept_polish = false;
    end

    if ~stage_history_cell{4}.skipped && accept_polish
        params_best = params_polished;
        last_stage = stage_history_cell{4};
        final_stage_names = stage_history_cell{4}.names;
        final_stage_metrics = stageD_metrics;
        final_calib_source = stageD_calib;
        best_stage = 4;
    elseif ~stage_history_cell{4}.skipped
        stage_history_cell{4}.skipped = true;
        stage_history_cell{4}.output.rejection_reason = ...
            'Rejected by RMSE/primary-gate guard after systemic polish.';
        fprintf(['[run_calibration] Stage D rejected by RMSE/primary-gate guard: ', ...
            'J %.6g -> %.6g, RMSE %.6g -> %.6g, primary_fail %d -> %d, ', ...
            'systemic_pair %.6g -> %.6g.\n'], ...
            polish_acceptance.J_before, polish_acceptance.J_after, ...
            polish_acceptance.rmse_before, polish_acceptance.rmse_after, ...
            polish_acceptance.primary_fail_before, polish_acceptance.primary_fail_after, ...
            polish_acceptance.systemic_pair_score_before, ...
            polish_acceptance.systemic_pair_score_after);
    end
end

if should_run_plausibility_polish(final_calib_source)
    stageE_calib = apply_plausibility_polish_profile(final_calib_source);
    stageE_names = select_plausibility_polish_names(stageE_calib);
    stageE_metrics = stageE_calib.metricFields;
    [params_plausible, stage_history_cell{5}] = optimize_stage(params_best, clinical, scenario, stageE_calib, ...
        stageE_names, stageE_metrics, do_parallel, pce_surrogate, fastMode, 'E');

    if ~stage_history_cell{5}.skipped
        [accept_plausibility, plausibility_acceptance] = accept_plausibility_polish( ...
            params_best, params_plausible, clinical, scenario, stageE_calib);
        stage_history_cell{5}.acceptance = plausibility_acceptance;
    else
        accept_plausibility = false;
    end

    if ~stage_history_cell{5}.skipped && accept_plausibility
        params_best = params_plausible;
        last_stage = stage_history_cell{5};
        final_stage_names = stage_history_cell{5}.names;
        final_stage_metrics = stageE_metrics;
        final_calib_source = stageE_calib;
        best_stage = 5;
    elseif ~stage_history_cell{5}.skipped
        stage_history_cell{5}.skipped = true;
        stage_history_cell{5}.output.rejection_reason = ...
            'Rejected by plausibility/RMSE guard after plausibility polish.';
        fprintf(['[run_calibration] Stage E rejected by plausibility guard: ', ...
            'warnings %d -> %d, RMSE %.6g -> %.6g, primary_fail %d -> %d, ', ...
            'systemic_pair %.6g -> %.6g.\n'], ...
            plausibility_acceptance.warning_before, plausibility_acceptance.warning_after, ...
            plausibility_acceptance.rmse_before, plausibility_acceptance.rmse_after, ...
            plausibility_acceptance.primary_fail_before, plausibility_acceptance.primary_fail_after, ...
            plausibility_acceptance.systemic_pair_score_before, ...
            plausibility_acceptance.systemic_pair_score_after);
    end
end

stage_history_cell = trim_empty_stage_history(stage_history_cell);

final_calib = make_stage_calib(final_calib_source, params_best, final_stage_names, final_stage_metrics);
final_x = pack_x(params_best, final_calib.referenceParams, final_calib.names, final_calib.caseProfile);
active_xbest = pack_x(params_best, calib.referenceParams, calib.names, calib.caseProfile);
fbest = objective_calibration(final_x, params_best, clinical, final_calib, scenario, pce_surrogate);

calib_out = struct();
calib_out.names = final_calib.names;
calib_out.names_all = calib.names_all;
calib_out.mask = calib.mask;
calib_out.x0 = pack_x(params0, final_calib.referenceParams, final_calib.names, final_calib.caseProfile);
calib_out.xbest = final_x;
calib_out.x0_all = calib.x0_all;
calib_out.x0_active = calib.x0;
calib_out.xbest_all = pack_x(params_best, calib.referenceParams, calib.names_all, calib.caseProfile);
calib_out.J0 = J0;
calib_out.fbest = fbest;
calib_out.lb = final_calib.lb;
calib_out.ub = final_calib.ub;
calib_out.scenario = scenario;
calib_out.primaryMetrics = calib.primaryMetrics;
calib_out.caseProfile = final_calib_source.caseProfile;
calib_out.parameterRegistry = calib.parameterRegistry;
calib_out.parameterRegistryActive = calib.parameterRegistryActive;
calib_out.xbest_active = active_xbest;
calib_out.parameterPlausibility = evaluate_parameter_plausibility(active_xbest, calib.parameterRegistryActive);
calib_out.improvement = J0 - fbest;
calib_out.best_stage = best_stage;
calib_out.best_restart = 0;
calib_out.stage_history = stage_history_cell;
calib_out.exitflag = last_stage.exitflag;
calib_out.output = last_stage.output;
calib_out.objective_breakdown = build_objective_breakdown(params_best, clinical, scenario, final_calib_source);
calib_out.use_parallel = do_parallel;

end

function tf = should_run_systemic_polish(calib, scenario)
tf = strcmp(scenario, 'pre_surgery') && ...
    isfield(calib, 'caseProfile') && isstruct(calib.caseProfile) && ...
    isfield(calib.caseProfile, 'systemicPolishEnabled') && ...
    calib.caseProfile.systemicPolishEnabled && ...
    isfield(calib.caseProfile, 'systemicPolishNames') && ...
    ~isempty(calib.caseProfile.systemicPolishNames);
end

function calib = apply_systemic_polish_profile(calib)
if isfield(calib.caseProfile, 'systemicPolishWeights') && ...
        isstruct(calib.caseProfile.systemicPolishWeights)
    calib.caseProfile.systemicLoadWeights = calib.caseProfile.systemicPolishWeights;
end
end

function names = select_systemic_polish_names(calib)
active_names = calib.names_all(calib.mask);
preferred = calib.caseProfile.systemicPolishNames(:)';
names = preferred(ismember(preferred, active_names));
end

function metrics = select_systemic_polish_metrics(calib)
if isfield(calib.caseProfile, 'systemicPolishMetrics') && ...
        ~isempty(calib.caseProfile.systemicPolishMetrics)
    metrics = calib.caseProfile.systemicPolishMetrics(:)';
else
    metrics = {'SAP_mean','SAP_max','CO_Lmin','SVR','QpQs','PAP_mean','LVEDV','LVEF'};
end
end

function tf = should_run_plausibility_polish(calib)
tf = isfield(calib, 'caseProfile') && isstruct(calib.caseProfile) && ...
    isfield(calib.caseProfile, 'plausibilityPolishEnabled') && ...
    calib.caseProfile.plausibilityPolishEnabled;
end

function calib = apply_plausibility_polish_profile(calib)
param_lambda = profile_scalar(calib.caseProfile, 'plausibilityPolishParamLambda', 2.25);
boundary_lambda = profile_scalar(calib.caseProfile, 'plausibilityPolishBoundaryLambda', 120.0);
calib.paramPlausibilityLambda = max(calib.paramPlausibilityLambda, param_lambda);
calib.boundaryPlausibilityLambda = max(calib.boundaryPlausibilityLambda, boundary_lambda);
end

function names = select_plausibility_polish_names(calib)
% SELECT_PLAUSIBILITY_POLISH_NAMES - all active knobs may trade off gently.
names = calib.names_all(calib.mask);
end

function [accept, diagnostics] = accept_systemic_polish(params_before, params_after, clinical, scenario, calib, pce_surrogate)
% ACCEPT_SYSTEMIC_POLISH - accept validation improvement without gate damage.
diagnostics = struct();
diagnostics.J_before = objective_for_params(params_before, clinical, scenario, calib, pce_surrogate);
diagnostics.J_after = objective_for_params(params_after, clinical, scenario, calib, pce_surrogate);
diagnostics.before = validation_summary_for_params(params_before, clinical, scenario, ...
    calib.primaryMetrics, calib.caseProfile);
diagnostics.after = validation_summary_for_params(params_after, clinical, scenario, ...
    calib.primaryMetrics, calib.caseProfile);
diagnostics.rmse_before = diagnostics.before.rmse;
diagnostics.rmse_after = diagnostics.after.rmse;
diagnostics.primary_fail_before = diagnostics.before.primary_fail_count;
diagnostics.primary_fail_after = diagnostics.after.primary_fail_count;
diagnostics.co_abs_error_before = diagnostics.before.co_abs_error_pct;
diagnostics.co_abs_error_after = diagnostics.after.co_abs_error_pct;
diagnostics.svr_abs_error_before = diagnostics.before.svr_abs_error_pct;
diagnostics.svr_abs_error_after = diagnostics.after.svr_abs_error_pct;
diagnostics.systemic_pair_score_before = diagnostics.before.systemic_pair_score;
diagnostics.systemic_pair_score_after = diagnostics.after.systemic_pair_score;

rmse_ok = isfinite(diagnostics.rmse_after) && isfinite(diagnostics.rmse_before) && ...
    diagnostics.rmse_after < diagnostics.rmse_before - 1e-6;
primary_gate_ok = diagnostics.primary_fail_after <= diagnostics.primary_fail_before;
systemic_pair_ok = diagnostics.systemic_pair_score_after <= ...
    diagnostics.systemic_pair_score_before * (1 + 1e-6);
accept = rmse_ok && primary_gate_ok && systemic_pair_ok;
diagnostics.objective_improved = isfinite(diagnostics.J_after) && isfinite(diagnostics.J_before) && ...
    diagnostics.J_after <= diagnostics.J_before * (1 + 1e-6);
diagnostics.rmse_ok = rmse_ok;
diagnostics.primary_gate_ok = primary_gate_ok;
diagnostics.systemic_pair_ok = systemic_pair_ok;
diagnostics.accept = accept;
end

function [accept, diagnostics] = accept_plausibility_polish(params_before, params_after, clinical, scenario, calib)
% ACCEPT_PLAUSIBILITY_POLISH - reduce warnings while preserving fit quality.
diagnostics = struct();
diagnostics.before = validation_summary_for_params(params_before, clinical, scenario, ...
    calib.primaryMetrics, calib.caseProfile);
diagnostics.after = validation_summary_for_params(params_after, clinical, scenario, ...
    calib.primaryMetrics, calib.caseProfile);
diagnostics.plausibility_before = plausibility_summary_for_params(params_before, calib);
diagnostics.plausibility_after = plausibility_summary_for_params(params_after, calib);

diagnostics.rmse_before = diagnostics.before.rmse;
diagnostics.rmse_after = diagnostics.after.rmse;
diagnostics.primary_fail_before = diagnostics.before.primary_fail_count;
diagnostics.primary_fail_after = diagnostics.after.primary_fail_count;
diagnostics.systemic_pair_score_before = diagnostics.before.systemic_pair_score;
diagnostics.systemic_pair_score_after = diagnostics.after.systemic_pair_score;
diagnostics.warning_before = diagnostics.plausibility_before.n_warning;
diagnostics.warning_after = diagnostics.plausibility_after.n_warning;
diagnostics.warning_fraction_before = diagnostics.plausibility_before.warning_fraction;
diagnostics.warning_fraction_after = diagnostics.plausibility_after.warning_fraction;

rmse_tolerance = profile_scalar(calib.caseProfile, 'plausibilityPolishRmseTolerance', 0.005);
systemic_tolerance = profile_scalar(calib.caseProfile, 'plausibilityPolishSystemicPairTolerance', 0.05);
warning_ok = diagnostics.warning_after < diagnostics.warning_before || ...
    diagnostics.warning_fraction_after <= 0.40;
rmse_ok = diagnostics.rmse_after <= diagnostics.rmse_before + rmse_tolerance;
primary_gate_ok = diagnostics.primary_fail_after <= diagnostics.primary_fail_before;
systemic_pair_ok = diagnostics.systemic_pair_score_after <= ...
    diagnostics.systemic_pair_score_before * (1 + systemic_tolerance);

accept = warning_ok && rmse_ok && primary_gate_ok && systemic_pair_ok;
diagnostics.warning_ok = warning_ok;
diagnostics.rmse_ok = rmse_ok;
diagnostics.primary_gate_ok = primary_gate_ok;
diagnostics.systemic_pair_ok = systemic_pair_ok;
diagnostics.accept = accept;
end

function plausibility = plausibility_summary_for_params(params, calib)
x_active = pack_x(params, calib.referenceParams, calib.names, calib.caseProfile);
plausibility = evaluate_parameter_plausibility(x_active, calib.parameterRegistryActive);
end

function J = objective_for_params(params, clinical, scenario, calib, pce_surrogate)
global_calib = make_stage_calib(calib, params, calib.names, calib.metricFields);
x = pack_x(params, global_calib.referenceParams, global_calib.names, global_calib.caseProfile);
J = objective_calibration(x, params, clinical, global_calib, scenario, pce_surrogate);
end

function summary = validation_summary_for_params(params, clinical, scenario, primary_metrics, case_profile)
summary = struct('rmse', Inf, 'primary_fail_count', Inf, ...
    'co_abs_error_pct', Inf, 'svr_abs_error_pct', Inf, ...
    'systemic_pair_score', Inf);
if nargin < 5
    case_profile = struct();
end
try
    metrics = compute_clinical_indices(integrate_system(params), params);
catch
    return;
end

targets = get_calibration_targets(scenario, clinical);
clinical_values = nan(numel(targets), 1);
model_values = nan(numel(targets), 1);
metric_names = cell(numel(targets), 1);
for idx = 1:numel(targets)
    metric_name = targets(idx).Metric;
    metric_names{idx} = metric_name;
    clinical_values(idx) = targets(idx).ClinicalValue;
    if isfield(metrics, metric_name)
        model_values(idx) = metrics.(metric_name);
    end
end

valid = isfinite(clinical_values) & isfinite(model_values);
if isfield(case_profile, 'targetTiers') && isstruct(case_profile.targetTiers) && ...
        isfield(case_profile.targetTiers, 'excluded_from_primary_rmse')
    valid = valid & ~ismember(metric_names, ...
        case_profile.targetTiers.excluded_from_primary_rmse(:)');
end
if ~any(valid)
    return;
end
relative_error = (model_values(valid) - clinical_values(valid)) ./ ...
    max(abs(clinical_values(valid)), 1e-9);
summary.rmse = sqrt(mean(relative_error.^2));

primary_mask = ismember(metric_names, primary_metrics(:)');
primary_valid = primary_mask & valid;
if any(primary_valid)
    primary_error_pct = 100 * abs(model_values(primary_valid) - clinical_values(primary_valid)) ./ ...
        max(abs(clinical_values(primary_valid)), 1e-9);
    summary.primary_fail_count = sum(primary_error_pct > 5);
else
    summary.primary_fail_count = 0;
end

summary.co_abs_error_pct = metric_abs_error_pct(metric_names, clinical_values, model_values, 'CO_Lmin');
summary.svr_abs_error_pct = metric_abs_error_pct(metric_names, clinical_values, model_values, 'SVR');
systemic_errors = [summary.co_abs_error_pct, summary.svr_abs_error_pct];
systemic_errors = systemic_errors(isfinite(systemic_errors));
if ~isempty(systemic_errors)
    summary.systemic_pair_score = sum((systemic_errors / 100).^2);
end
end

function abs_error_pct = metric_abs_error_pct(metric_names, clinical_values, model_values, metric_name)
abs_error_pct = Inf;
idx = find(strcmp(metric_names, metric_name), 1, 'first');
if isempty(idx) || ~isfinite(clinical_values(idx)) || ~isfinite(model_values(idx))
    return;
end
abs_error_pct = 100 * abs(model_values(idx) - clinical_values(idx)) / ...
    max(abs(clinical_values(idx)), 1e-9);
end

function value = profile_scalar(profile, field_name, default_value)
value = default_value;
if isstruct(profile) && isfield(profile, field_name) && ...
        isfinite(profile.(field_name))
    value = profile.(field_name);
end
end

function history = trim_empty_stage_history(history)
last_nonempty = find(~cellfun(@isempty, history), 1, 'last');
if isempty(last_nonempty)
    history = {};
else
    history = history(1:last_nonempty);
end
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
max_fun_evals = env_numeric_or_default('UNIFIED_VSD_MAX_FUN_EVALS', ...
    ternary(fastMode, 1200, 3000));
max_iterations = env_numeric_or_default('UNIFIED_VSD_MAX_ITERATIONS', ...
    ternary(fastMode, 100, 250));
opts = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'HessianApproximation', 'lbfgs', ...
    'FiniteDifferenceStepSize', 1e-5, ...
    'FiniteDifferenceType', 'forward', ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', max_fun_evals, ...
    'MaxIterations', max_iterations, ...
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
    params_out = set_calibration_param_value( ...
        params_out, calib.referenceParams, stage_calib.names{i}, xbest(i), calib.caseProfile);
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
stage_calib.x0 = pack_x(params_ref, calib.referenceParams, stage_names, calib.caseProfile);
stage_calib.lb = calib.lb_all(stage_mask);
stage_calib.ub = calib.ub_all(stage_mask);
stage_calib.parameterRegistryActive = calib.parameterRegistry(stage_mask, :);
end

function names = select_stage_c_names(calib, gsa_out, scenario)
active_names = calib.names_all(calib.mask);
if isempty(gsa_out) || ~isstruct(gsa_out)
    names = fallback_stage_c_names(active_names, scenario, calib.caseProfile);
    return;
end

scores = zeros(numel(active_names), 1);
for i = 1:numel(active_names)
    scores(i) = aggregate_parameter_score(active_names{i}, gsa_out, calib.primaryMetrics, calib.caseProfile);
end
[~, order] = sort(scores, 'descend');
keep = min(5, numel(order));
names = active_names(order(1:keep));
if keep == 0
    names = fallback_stage_c_names(active_names, scenario, calib.caseProfile);
end
end

function score = aggregate_parameter_score(param_name, gsa_out, primary_metrics, case_profile)
score = 0;
member_names = expand_calibration_parameter_names(param_name, case_profile);
for k = 1:numel(primary_metrics)
    mf = primary_metrics{k};
    if ~isfield(gsa_out, mf) || ~isfield(gsa_out.(mf), 'table')
        continue;
    end
    tbl = gsa_out.(mf).table;
    if ~ismember('Parameter', tbl.Properties.VariableNames) || ~ismember('ST', tbl.Properties.VariableNames)
        continue;
    end
    metric_score = 0;
    for member_idx = 1:numel(member_names)
        member_mask = strcmp(tbl.Parameter, member_names{member_idx});
        if any(member_mask)
            metric_score = max(metric_score, max(tbl.ST(member_mask)));
        end
    end
    score = score + metric_score;
end
end

function names = fallback_stage_c_names(active_names, scenario, case_profile)
if nargin >= 3 && isstruct(case_profile) && ...
        isfield(case_profile, 'stageCPreferredNames') && ...
        ~isempty(case_profile.stageCPreferredNames)
    preferred = case_profile.stageCPreferredNames(:)';
    names = preferred(ismember(preferred, active_names));
    if ~isempty(names)
        return;
    end
end

switch scenario
    case 'pre_surgery'
        preferred = {'R.vsd','vsd.Cd','R.PAR','R.SAR','E.LV.EA','V0.LV'};
    otherwise
        preferred = {'R.SAR','R.PAR','E.LV.EA','V0.LV','R.SVEN'};
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
systemic_bundle_metrics = {'SAP_mean','RAP_mean','CO_Lmin','SVR'};
scenario_src = struct();
if isstruct(clinical) && isfield(clinical, scenario)
    scenario_src = clinical.(scenario);
end
has_full_systemic_bundle = isfield(scenario_src, 'CO_Lmin') && ...
    isfield(scenario_src, 'SAP_mean_mmHg') && isfield(scenario_src, 'RAP_mean_mmHg') && ...
    all(isfinite([scenario_src.CO_Lmin, scenario_src.SAP_mean_mmHg, scenario_src.RAP_mean_mmHg]));

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
    if has_full_systemic_bundle && ismember(mf, systemic_bundle_metrics)
        err_rel = abs(metrics.(mf) - y_clin) / max(abs(y_clin), 1e-6);
        unweighted(k) = err_rel;
        tiers{k} = 'systemic_bundle';
        weighted(k) = NaN;
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

function x = pack_x(params, reference_params, names, case_profile)
x = zeros(numel(names), 1);
for i = 1:numel(names)
    x(i) = get_calibration_param_value(params, reference_params, names{i}, case_profile);
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end

function value = env_numeric_or_default(name, default_value)
value = default_value;
raw = getenv(name);
if isempty(raw)
    return;
end
candidate = str2double(raw);
if isfinite(candidate) && candidate > 0
    value = candidate;
end
end
