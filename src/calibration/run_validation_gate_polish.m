function [params_best, gate_out] = run_validation_gate_polish(params_start, clinical, scenario, calib, fastMode)
% RUN_VALIDATION_GATE_POLISH
% -----------------------------------------------------------------------
% Final validation-facing polish for patient calibration.
%
% This stage optimizes the same governed clinical error surface that is
% reported by validation_report.m. It is intentionally separate from the
% mechanistic staged objective: the staged objective finds a physiologic
% basin, while this polish prevents a candidate with lower internal loss but
% worse reported validation RMSE from being kept.
%
% INPUTS:
%   params_start - calibrated parameter struct before final polish       [-]
%   clinical     - unified clinical data struct                          [-]
%   scenario     - scenario string: 'pre_surgery' | 'post_surgery'        [-]
%   calib        - calibration configuration from calibration_param_sets  [-]
%   fastMode     - true for reduced evaluation budget                    [-]
%
% OUTPUTS:
%   params_best  - accepted polished parameters, or params_start         [-]
%   gate_out     - diagnostics for acceptance/rejection                  [-]
%
% ASSUMPTIONS:
%   - Registry bounds define hard plausibility limits.
%   - Derived-validation metrics such as PVR/SVR remain report-only.
%   - Warnings are reported but do not reject a physiologic improvement.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%   [2] AGENTS.md, Sections 9-10.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-21
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 5
    fastMode = false;
end

params_best = params_start;
gate_out = init_gate_output(calib);

if isempty(gate_out.names)
    gate_out.skipped = true;
    gate_out.output.message = 'No active parameters available for validation polish.';
    return;
end

before = evaluate_gate_candidate(params_start, clinical, scenario, calib);
gate_out.before = before;

if before.valid && before.primary_fail_count == 0 && before.primary_rmse <= gate_acceptance_rmse(calib)
    gate_out.skipped = true;
    gate_out.output.message = 'Existing candidate already satisfies the validation gate.';
    return;
end

stage_calib = make_gate_calib(calib, params_start, gate_out.names);
x_current = pack_gate_x(params_start, stage_calib.referenceParams, stage_calib.names, stage_calib.caseProfile);
x_current = clamp_x(x_current, stage_calib.lb, stage_calib.ub);

starts = build_gate_starts(x_current, stage_calib, calib);
max_fun_evals = env_numeric_or_default('UNIFIED_VSD_GATE_POLISH_MAX_FUN_EVALS', ...
    ternary(fastMode, 350, 900));
max_iterations = env_numeric_or_default('UNIFIED_VSD_GATE_POLISH_MAX_ITERATIONS', ...
    ternary(fastMode, 50, 120));
max_starts = env_numeric_or_default('UNIFIED_VSD_GATE_POLISH_STARTS', ...
    ternary(fastMode, 1, min(2, size(starts, 2))));
max_starts = min(max_starts, size(starts, 2));

obj = @(x) validation_gate_objective( ...
    x, params_start, clinical, scenario, stage_calib, before);

opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize', 1e-4, ...
    'Display', 'off', ...
    'MaxFunctionEvaluations', max_fun_evals, ...
    'MaxIterations', max_iterations, ...
    'OptimalityTolerance', 1e-5, ...
    'StepTolerance', 1e-6);

best_score = obj(x_current);
best_x = x_current;
best_exitflag = NaN;
best_output = struct();
restart_scores = nan(max_starts, 1);
restart_exitflags = nan(max_starts, 1);

for start_idx = 1:max_starts
    x0 = starts(:, start_idx);
    try
        [x_trial, f_trial, exitflag, output] = fmincon( ...
            obj, x0, [], [], [], [], stage_calib.lb, stage_calib.ub, [], opts);
    catch ME
        restart_scores(start_idx) = Inf;
        restart_exitflags(start_idx) = NaN;
        best_output.last_error = ME.message;
        continue;
    end

    restart_scores(start_idx) = f_trial;
    restart_exitflags(start_idx) = exitflag;
    if isfinite(f_trial) && f_trial < best_score
        best_score = f_trial;
        best_x = x_trial(:);
        best_exitflag = exitflag;
        best_output = output;
    end
end

params_trial = apply_gate_x(params_start, stage_calib.referenceParams, ...
    stage_calib.names, best_x, stage_calib.caseProfile);
after = evaluate_gate_candidate(params_trial, clinical, scenario, calib);
plausibility = evaluate_parameter_plausibility( ...
    pack_gate_x(params_trial, calib.referenceParams, calib.names, calib.caseProfile), ...
    calib.parameterRegistryActive);

[accept, acceptance] = accept_gate_candidate(before, after, plausibility, calib);

gate_out.after = after;
gate_out.acceptance = acceptance;
gate_out.plausibility = plausibility;
gate_out.x0 = x_current;
gate_out.xbest = best_x;
gate_out.fval = best_score;
gate_out.exitflag = best_exitflag;
gate_out.output = best_output;
gate_out.restart_scores = restart_scores;
gate_out.restart_exitflags = restart_exitflags;
gate_out.skipped = false;

if accept
    params_best = params_trial;
else
    gate_out.output.rejection_reason = acceptance.reason;
end
end

function gate_out = init_gate_output(calib)
gate_out = struct();
gate_out.label = 'F';
gate_out.names = calib.names_all(calib.mask);
gate_out.x0 = [];
gate_out.xbest = [];
gate_out.fval = NaN;
gate_out.exitflag = NaN;
gate_out.output = struct();
gate_out.skipped = false;
gate_out.before = struct();
gate_out.after = struct();
gate_out.acceptance = struct();
gate_out.plausibility = struct();
end

function stage_calib = make_gate_calib(calib, params_ref, requested_names)
stage_mask = calib.mask & ismember(calib.names_all, requested_names(:)');
stage_calib = calib;
stage_calib.names = calib.names_all(stage_mask);
stage_calib.x0 = pack_gate_x(params_ref, calib.referenceParams, ...
    stage_calib.names, calib.caseProfile);
stage_calib.lb = calib.lb_all(stage_mask);
stage_calib.ub = calib.ub_all(stage_mask);
stage_calib.parameterRegistryActive = calib.parameterRegistry(stage_mask, :);
end

function starts = build_gate_starts(x_current, stage_calib, calib)
x_seed = pack_gate_x_from_seed(calib, stage_calib.names);
x_seed = clamp_x(x_seed, stage_calib.lb, stage_calib.ub);
x_mid = 0.5 * (stage_calib.lb(:) + stage_calib.ub(:));
starts = [x_current(:), x_seed(:), x_mid(:)];

if numel(stage_calib.names) >= 3
    x_low_resistance = x_current(:);
    resistance_mask = startsWith(stage_calib.names(:), 'R.') | ...
        startsWith(stage_calib.names(:), 'group.R');
    x_low_resistance(resistance_mask) = max(stage_calib.lb(resistance_mask), ...
        0.85 * x_low_resistance(resistance_mask));
    starts(:, end + 1) = clamp_x(x_low_resistance, stage_calib.lb, stage_calib.ub); %#ok<AGROW>
end

[~, unique_idx] = unique(round(starts', 10), 'rows', 'stable');
starts = starts(:, unique_idx);
end

function x_seed = pack_gate_x_from_seed(calib, names)
x_seed = nan(numel(names), 1);
for idx = 1:numel(names)
    seed_idx = find(strcmp(calib.names, names{idx}), 1, 'first');
    if ~isempty(seed_idx)
        x_seed(idx) = calib.x0(seed_idx);
    else
        all_idx = find(strcmp(calib.names_all, names{idx}), 1, 'first');
        x_seed(idx) = calib.x0_all(all_idx);
    end
end
end

function score = validation_gate_objective(x, params_base, clinical, scenario, calib, before)
params = apply_gate_x(params_base, calib.referenceParams, ...
    calib.names, x(:), calib.caseProfile);
candidate = evaluate_gate_candidate(params, clinical, scenario, calib);
if ~candidate.valid
    score = calib.invalidPenaltyScale;
    return;
end

score = 0;
valid_rows = candidate.rows.ValidForScore;
if any(valid_rows)
    rel_error = candidate.rows.RelativeError(valid_rows);
    gate = candidate.rows.Gate(valid_rows);
    weight = candidate.rows.Weight(valid_rows);
    scaled_error = abs(rel_error) ./ max(gate, 1e-9);
    excess = max(0, abs(rel_error) - gate) ./ max(gate, 1e-9);
    score = score + mean(weight .* scaled_error.^2);
    score = score + 2.0 * mean(weight .* excess.^2);
end

primary_excess = max(0, candidate.max_primary_abs_error_pct - gate_acceptance_pct(calib)) / 100;
score = score + 5.0 * candidate.primary_fail_count + 20.0 * primary_excess^2;

if isfinite(before.primary_rmse) && candidate.primary_rmse > before.primary_rmse
    score = score + 5.0 * (candidate.primary_rmse - before.primary_rmse)^2;
end

score = score + 0.03 * parameter_log_drift_penalty(x, calib);
score = score + candidate.validity_penalty;
end

function candidate = evaluate_gate_candidate(params, clinical, scenario, calib)
candidate = struct();
candidate.valid = false;
candidate.physiology_valid = false;
candidate.validity_penalty = calib.invalidPenaltyScale;
candidate.primary_rmse = Inf;
candidate.full_rmse = Inf;
candidate.hard_rmse = Inf;
candidate.soft_rmse = Inf;
candidate.primary_fail_count = Inf;
candidate.max_primary_abs_error_pct = Inf;
candidate.rows = table();
candidate.metrics = struct();

try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
    validity = evaluate_simulation_validity(sim, params, metrics, scenario, clinical);
catch
    return;
end

targets = get_calibration_targets(scenario, clinical);
target_tiers = calib.targetTiers;
rows = build_gate_rows(targets, metrics, target_tiers, calib);

primary_mask = rows.ValidForRMSE & rows.IncludedInPrimaryRMSE;
full_mask = rows.ValidForRMSE;
hard_mask = rows.ValidForRMSE & strcmp(rows.Tier, 'hard');
soft_mask = rows.ValidForRMSE & strcmp(rows.Tier, 'soft');

candidate.valid = true;
candidate.physiology_valid = validity.is_valid;
candidate.validity_penalty = validity.penalty;
candidate.metrics = metrics;
candidate.validity = validity;
candidate.rows = rows;
candidate.primary_rmse = rmse_from_rows(rows, primary_mask);
candidate.full_rmse = rmse_from_rows(rows, full_mask);
candidate.hard_rmse = rmse_from_rows(rows, hard_mask);
candidate.soft_rmse = rmse_from_rows(rows, soft_mask);

primary_gate_mask = rows.ValidForRMSE & ismember(rows.Metric, calib.primaryMetrics(:));
if any(primary_gate_mask)
    abs_error_pct = abs(rows.RelativeError(primary_gate_mask)) * 100;
    candidate.primary_fail_count = sum(abs_error_pct > gate_acceptance_pct(calib));
    candidate.max_primary_abs_error_pct = max(abs_error_pct);
else
    candidate.primary_fail_count = 0;
    candidate.max_primary_abs_error_pct = 0;
end
end

function rows = build_gate_rows(targets, metrics, target_tiers, calib)
n_targets = numel(targets);
metric_col = cell(n_targets, 1);
tier_col = cell(n_targets, 1);
clinical_col = nan(n_targets, 1);
model_col = nan(n_targets, 1);
relative_error_col = nan(n_targets, 1);
gate_col = nan(n_targets, 1);
weight_col = zeros(n_targets, 1);
valid_rmse_col = false(n_targets, 1);
valid_score_col = false(n_targets, 1);
included_primary_col = false(n_targets, 1);
included_cal_col = false(n_targets, 1);

for idx = 1:n_targets
    metric_name = targets(idx).Metric;
    metric_col{idx} = metric_name;
    clinical_col(idx) = targets(idx).ClinicalValue;
    if isfield(metrics, metric_name)
        model_col(idx) = metrics.(metric_name);
    end

    [tier, included_cal, included_primary] = target_tier_metadata(target_tiers, metric_name);
    tier_col{idx} = tier;
    included_cal_col(idx) = included_cal;
    included_primary_col(idx) = included_primary;

    if isfinite(clinical_col(idx)) && isfinite(model_col(idx))
        relative_error_col(idx) = (model_col(idx) - clinical_col(idx)) / ...
            max(abs(clinical_col(idx)), 1e-9);
        valid_rmse_col(idx) = true;
    end

    [gate_col(idx), weight_col(idx), valid_score_col(idx)] = ...
        gate_policy(metric_name, tier, included_cal, included_primary, calib);
end

valid_score_col = valid_score_col & valid_rmse_col;

rows = table(metric_col, tier_col, clinical_col, model_col, relative_error_col, ...
    gate_col, weight_col, valid_rmse_col, valid_score_col, included_cal_col, ...
    included_primary_col, ...
    'VariableNames', {'Metric','Tier','Clinical','Model','RelativeError', ...
    'Gate','Weight','ValidForRMSE','ValidForScore','IncludedInCalibration', ...
    'IncludedInPrimaryRMSE'});
end

function [tier, included_cal, included_primary] = target_tier_metadata(target_tiers, metric_name)
tier = 'unavailable';
included_cal = false;
included_primary = false;
if ~isstruct(target_tiers) || ~isfield(target_tiers, 'table') || isempty(target_tiers.table)
    return;
end
tier_tbl = target_tiers.table;
idx = find(strcmp(tier_tbl.Metric, metric_name), 1, 'first');
if isempty(idx)
    return;
end
tier = tier_tbl.Tier{idx};
included_cal = tier_tbl.IncludedInCalibration(idx);
included_primary = tier_tbl.IncludedInPrimaryRMSE(idx);
end

function [gate, weight, score_metric] = gate_policy(metric_name, tier, included_cal, included_primary, calib)
gate = gate_acceptance_pct(calib) / 100;
weight = 0;
score_metric = false;

if any(strcmp(tier, {'unavailable','consistency_check_only', ...
        'derived_validation','validation_holdout'}))
    return;
end
if ~(included_cal || included_primary)
    return;
end

score_metric = true;
switch tier
    case 'hard'
        gate = gate_acceptance_pct(calib) / 100;
        weight = 1.60;
    case 'soft'
        gate = gate_secondary_pct(calib) / 100;
        weight = 0.90;
    case 'validation_only'
        gate = gate_secondary_pct(calib) / 100;
        weight = 0.40;
    otherwise
        gate = gate_secondary_pct(calib) / 100;
        weight = 0.50;
end

if ismember(metric_name, calib.primaryMetrics(:)')
    weight = weight * 1.50;
end
if strcmp(metric_name, 'CO_Lmin')
    weight = weight * 1.60;
end
if strcmp(metric_name, 'PAP_mean') || strcmp(metric_name, 'QpQs')
    weight = weight * 1.25;
end
end

function value = rmse_from_rows(rows, mask)
if ~any(mask)
    value = Inf;
    return;
end
value = sqrt(mean(rows.RelativeError(mask).^2));
end

function [accept, diagnostics] = accept_gate_candidate(before, after, plausibility, calib)
diagnostics = struct();
diagnostics.accept = false;
diagnostics.reason = '';
diagnostics.primary_rmse_before = before.primary_rmse;
diagnostics.primary_rmse_after = after.primary_rmse;
diagnostics.primary_fail_before = before.primary_fail_count;
diagnostics.primary_fail_after = after.primary_fail_count;
diagnostics.max_primary_before = before.max_primary_abs_error_pct;
diagnostics.max_primary_after = after.max_primary_abs_error_pct;
diagnostics.plausibility_fail = plausibility.n_fail;
diagnostics.physiology_valid = after.physiology_valid;

if ~after.valid || ~after.physiology_valid
    diagnostics.reason = 'invalid_or_nonphysiological';
    accept = false;
    return;
end
if plausibility.n_fail > 0
    diagnostics.reason = 'hard_parameter_plausibility_failure';
    accept = false;
    return;
end

rmse_improved = after.primary_rmse < before.primary_rmse - gate_rmse_delta(calib);
primary_fail_improved = after.primary_fail_count < before.primary_fail_count;
max_error_improved = after.max_primary_abs_error_pct < before.max_primary_abs_error_pct - 1.0;
not_worse = after.primary_rmse <= before.primary_rmse + 1e-6 && ...
    after.primary_fail_count <= before.primary_fail_count;

accept = (rmse_improved || primary_fail_improved || max_error_improved) && not_worse;
diagnostics.rmse_improved = rmse_improved;
diagnostics.primary_fail_improved = primary_fail_improved;
diagnostics.max_error_improved = max_error_improved;
diagnostics.not_worse = not_worse;
diagnostics.accept = accept;
if accept
    diagnostics.reason = 'accepted_validation_gate_improvement';
else
    diagnostics.reason = 'no_validation_gate_improvement';
end
end

function penalty = parameter_log_drift_penalty(x, calib)
penalty = 0;
if ~isfield(calib, 'parameterRegistryActive') || isempty(calib.parameterRegistryActive)
    return;
end
baseline = calib.parameterRegistryActive.baseline_scaled(:);
x = x(:);
valid = isfinite(x) & isfinite(baseline) & x > 0 & baseline > 0;
if any(valid)
    penalty = sum(log(x(valid) ./ baseline(valid)).^2);
end
end

function params = apply_gate_x(params, reference_params, names, x, case_profile)
for idx = 1:numel(names)
    params = set_calibration_param_value( ...
        params, reference_params, names{idx}, x(idx), case_profile);
end
end

function x = pack_gate_x(params, reference_params, names, case_profile)
x = zeros(numel(names), 1);
for idx = 1:numel(names)
    x(idx) = get_calibration_param_value(params, reference_params, names{idx}, case_profile);
end
end

function x = clamp_x(x, lb, ub)
x = min(max(x(:), lb(:)), ub(:));
end

function value = gate_acceptance_pct(calib)
value = 10;
if isfield(calib, 'caseProfile') && isfield(calib.caseProfile, 'acceptancePrimaryErrorPct') && ...
        isfinite(calib.caseProfile.acceptancePrimaryErrorPct)
    value = calib.caseProfile.acceptancePrimaryErrorPct;
end
end

function value = gate_secondary_pct(calib)
value = 15;
if isfield(calib, 'caseProfile') && isfield(calib.caseProfile, 'acceptanceSecondaryErrorPct') && ...
        isfinite(calib.caseProfile.acceptanceSecondaryErrorPct)
    value = calib.caseProfile.acceptanceSecondaryErrorPct;
end
end

function value = gate_acceptance_rmse(calib)
value = gate_acceptance_pct(calib) / 100;
end

function value = gate_rmse_delta(calib)
value = 0.003;
if isfield(calib, 'caseProfile') && isfield(calib.caseProfile, 'validationPolishRmseDelta') && ...
        isfinite(calib.caseProfile.validationPolishRmseDelta)
    value = calib.caseProfile.validationPolishRmseDelta;
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

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
