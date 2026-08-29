function tests = test_parameter_identifiability()
% TEST_PARAMETER_IDENTIFIABILITY
% -----------------------------------------------------------------------
% Contract tests for identifiability_tables_from_matrix (the pure-math
% reduction step) and a smoke test for analyse_parameter_identifiability
% (the forward-simulation wrapper), added by PRD
% reyna_statistical_calibration_v1 Phase 3.
%
% GSA/Sobol ranks parameters one at a time by marginal influence. It cannot
% detect that two retained parameters are collinear -- e.g. C.SAR and
% group.R_sys_scale trading off along the systemic RC time constant
% identified in docs/reyna_zhang_fullmetric_results_20260828.md section 3.5.
% These tests exercise the analytic cases the PRD specifies directly:
% duplicated columns, orthogonal columns, and a zero (inactive) column.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 3)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function test_duplicated_columns_are_perfectly_correlated(tc)
% Two exactly-duplicated columns must show |rho| = 1 and be flagged.
S = [1 1 0; 2 2 1; 3 3 -1; 4 4 2];
names = {'A', 'B', 'C'};
out = identifiability_tables_from_matrix(S, names);

row = strcmp(out.pair_table.ParameterA, 'A') & strcmp(out.pair_table.ParameterB, 'B');
verifyTrue(tc, any(row));
verifyEqual(tc, out.pair_table.Correlation(row), 1.0, 'AbsTol', 1e-9);
verifyTrue(tc, out.pair_table.Collinear(row));

param_row = strcmp(out.parameter_table.Parameter, 'A');
verifyTrue(tc, out.parameter_table.Flagged(param_row));
end

function test_orthogonal_columns_have_near_zero_correlation_and_good_condition(tc)
% Orthogonal columns (a Hadamard-like design) must show correlations near
% zero and a condition number near 1 (well-conditioned).
S = [1 1 1; 1 -1 1; 1 1 -1; 1 -1 -1];
names = {'A', 'B', 'C'};
out = identifiability_tables_from_matrix(S, names);

verifyEqual(tc, max(abs(out.pair_table.Correlation)), 0, 'AbsTol', 1e-9, ...
    'Orthogonal columns must show zero pairwise correlation.');
verifyFalse(tc, any(out.pair_table.Collinear), ...
    'No pair should be flagged collinear when columns are orthogonal.');
verifyLessThan(tc, out.condition_number, 2, ...
    'An orthogonal design should be well-conditioned (cond near 1).');
end

function test_zero_column_is_flagged_inactive_without_dividing_by_zero(tc)
S = [1 0 2; 2 0 -1; 3 0 3; -1 0 1];
names = {'A', 'ZeroCol', 'C'};
out = identifiability_tables_from_matrix(S, names);

row = strcmp(out.parameter_table.Parameter, 'ZeroCol');
verifyTrue(tc, any(row));
verifyTrue(tc, out.parameter_table.Inactive(row), ...
    'A parameter with zero sensitivity to every governed metric must be flagged inactive.');
verifyEqual(tc, out.parameter_table.ColumnNorm(row), 0, 'AbsTol', 1e-12);
verifyTrue(tc, all(isfinite(out.parameter_table.ColumnNorm)), ...
    'No NaN/Inf may result from a zero column (division-by-zero guard).');
% A zero column makes S exactly rank-deficient, so cond(S) = Inf is the
% mathematically correct answer here, not a computation failure. The
% guarantee under test is that computing it never throws.
verifyTrue(tc, isnumeric(out.condition_number) && isscalar(out.condition_number), ...
    'condition_number must be a plain numeric scalar (Inf/NaN allowed), never a thrown error.');
end

function test_all_nan_column_is_dropped_not_zero_filled(tc)
% A column that never evaluated (every perturbation failed) must be
% excluded from the reported table entirely -- it is a different finding
% from "evaluated but inactive".
S = [1 NaN 2; 2 NaN -1; 3 NaN 3];
names = {'A', 'NeverEvaluated', 'C'};
out = identifiability_tables_from_matrix(S, names);

verifyFalse(tc, any(strcmp(out.parameter_table.Parameter, 'NeverEvaluated')), ...
    'A parameter that never evaluated must not appear as a false "inactive" row.');
verifyEqual(tc, height(out.parameter_table), 2);
end

function test_output_has_one_row_per_active_parameter_plus_pair_table(tc)
S = [1 2 3; 4 5 6; 7 8 9; 1 1 1];
names = {'P1', 'P2', 'P3'};
out = identifiability_tables_from_matrix(S, names);

verifyEqual(tc, height(out.parameter_table), 3);
verifyEqual(tc, height(out.pair_table), 3, ...
    'Three parameters must produce C(3,2) = 3 pairs.');
verifyTrue(tc, all(ismember({'Parameter','ColumnNorm','Inactive', ...
    'MaxAbsCorrelation','MostCorrelatedWith','Flagged'}, ...
    out.parameter_table.Properties.VariableNames)));
verifyTrue(tc, all(ismember({'ParameterA','ParameterB','Correlation','Collinear'}, ...
    out.pair_table.Properties.VariableNames)));
end

function test_pair_table_covers_every_unordered_pair_exactly_once(tc)
S = magic(4);
names = {'W', 'X', 'Y', 'Z'};
out = identifiability_tables_from_matrix(S, names);
verifyEqual(tc, height(out.pair_table), 6, 'C(4,2) = 6 pairs expected.');

pairs = strcat(out.pair_table.ParameterA, '__', out.pair_table.ParameterB);
verifyEqual(tc, numel(unique(pairs)), 6, 'No pair may be duplicated.');
end

function test_empty_matrix_returns_empty_tables_without_error(tc)
out = identifiability_tables_from_matrix([], {});
verifyEqual(tc, height(out.parameter_table), 0);
verifyEqual(tc, height(out.pair_table), 0);
verifyTrue(tc, isnan(out.condition_number));
end

function test_single_parameter_has_no_pairs(tc)
S = [1; 2; 3];
names = {'Only'};
out = identifiability_tables_from_matrix(S, names);
verifyEqual(tc, height(out.parameter_table), 1);
verifyEqual(tc, height(out.pair_table), 0);
end

function test_analyse_parameter_identifiability_smoke_reyna_operating_point(tc)
% End-to-end smoke test at the current Reyna calibration operating point.
% Not a numeric-value assertion (the forward model is expensive and its
% exact operating point will move as the campaign iterates) -- this checks
% the wrapper actually produces a well-formed report against a real
% params_best/calib_out/gate_table triple, exercising the perturbation and
% simulation path that the analytic tests above cannot reach.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');

params_ref = default_parameters();
patient = struct('age_years', clinical.common.age_years, ...
    'weight_kg', clinical.common.weight_kg, 'height_cm', clinical.common.height_cm, ...
    'sex', clinical.common.sex, 'BSA', clinical.common.BSA, ...
    'maturation_mode', 'normal', 'scaling_mode', 'zhang');
params_scaled = apply_scaling(params_ref, patient);
params0 = params_from_clinical(params_scaled, clinical, 'pre_surgery', params_scaled, profile);

% Use a small, cheap subset (2 parameters) rather than the full active set
% so this smoke test stays fast; it is exercising the wrapper's plumbing,
% not producing a publishable identifiability report.
calib_out = struct('names', {{'R.SVEN', 'C.SAR'}}, 'caseProfile', profile);

gate_table = table({'RAP_mean'; 'SAP_mean'}, {'mmHg'; 'mmHg'}, ...
    {'hard'; 'hard'}, [5.0; 71.3], [5.0; 71.3], [0; 0], [0; 0], ...
    [0.25; 3.565], [0; 0], [0; 0], true(2,1), true(2,1), true(2,1), ...
    true(2,1), {'none'; 'none'}, ...
    'VariableNames', {'Metric','Unit','Tier','Clinical','Model', ...
    'Error_pct','AbsError_pct','Sigma','ZScore','ZScoreSquared', ...
    'InObjective','InPrimaryRMSE','WithinGate','WithinExcellent','Flag'});

report = analyse_parameter_identifiability(params0, 'pre_surgery', calib_out, gate_table, '');

verifyTrue(tc, isstruct(report));
verifyTrue(tc, ismember('parameter_table', fieldnames(report)));
verifyTrue(tc, ismember('sensitivity_matrix', fieldnames(report)));
verifyEqual(tc, size(report.sensitivity_matrix), [2, 2]);
end
