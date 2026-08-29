function tests = test_chi_squared_report()
% TEST_CHI_SQUARED_REPORT
% -----------------------------------------------------------------------
% Contract tests for compute_chi_squared_report, the discrepancy-principle
% goodness-of-fit statistic added by PRD reyna_statistical_calibration_v1
% Phase 2.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 2)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function tbl = make_table(z_scores, in_primary_rmse, names)
% MAKE_TABLE - minimal synthetic gate table with the columns
% compute_chi_squared_report requires.
n = numel(z_scores);
if nargin < 2 || isempty(in_primary_rmse)
    in_primary_rmse = true(n, 1);
end
if nargin < 3 || isempty(names)
    names = arrayfun(@(i) sprintf('M%d', i), 1:n, 'UniformOutput', false)';
end
z_scores = z_scores(:);
tbl = table(names(:), z_scores, z_scores.^2, logical(in_primary_rmse(:)), ...
    'VariableNames', {'Metric', 'ZScore', 'ZScoreSquared', 'InPrimaryRMSE'});
end

function test_all_metrics_at_one_sigma_gives_chi2_per_obs_one(tc)
tbl = make_table([1, -1, 1, -1, 1, -1, 1, -1, 1]);
report = compute_chi_squared_report(tbl, 7);
verifyEqual(tc, report.n_obs, 9);
verifyEqual(tc, report.chi2, 9, 'AbsTol', 1e-9);
verifyEqual(tc, report.chi2_per_obs, 1.0, 'AbsTol', 1e-9);
verifyEqual(tc, report.interpretation, 'consistent');
end

function test_all_metrics_at_three_sigma_is_underfit(tc)
tbl = make_table(3 * ones(9, 1));
report = compute_chi_squared_report(tbl, 7);
verifyEqual(tc, report.chi2_per_obs, 9.0, 'AbsTol', 1e-9);
verifyEqual(tc, report.interpretation, 'underfit');
end

function test_all_metrics_at_tenth_sigma_is_overfit(tc)
tbl = make_table(0.1 * ones(9, 1));
report = compute_chi_squared_report(tbl, 7);
verifyEqual(tc, report.chi2_per_obs, 0.01, 'AbsTol', 1e-9);
verifyEqual(tc, report.interpretation, 'overfit');
end

function test_dof_le_zero_gives_nan_reduced_chi2_no_error(tc)
tbl = make_table([1, 1, 1]);
% 3 observations, 7 parameters -> dof = -4.
report = compute_chi_squared_report(tbl, 7);
verifyEqual(tc, report.dof, -4);
verifyTrue(tc, isnan(report.chi2_reduced));
verifyTrue(tc, isfinite(report.chi2), 'chi2 itself must still compute.');
end

function test_negative_dof_is_reported_not_clamped_to_zero(tc)
% A clamped dof of 0 reads as "exactly determined" when the truth may be
% "more free parameters than observations" -- the single most important
% caveat on any fit reported here. Reporting 0 for N=9,p=12 also prints a
% line that is arithmetically false. Regression guard for that behaviour.
tbl = make_table(ones(9, 1));
report = compute_chi_squared_report(tbl, 12);
verifyEqual(tc, report.dof, -3, ...
    'dof must be the true N - p, including when negative.');
verifyEqual(tc, report.dof_note, 'over_parameterised');
verifyTrue(tc, isnan(report.chi2_reduced), ...
    'reduced chi2 must still refuse to divide by a non-positive dof.');
end

function test_interpretation_is_qualified_when_underdetermined(tc)
% chi2/N inside the "consistent" band means nothing when the model has at
% least as many free parameters as observations: residuals that small are
% guaranteed by construction. The label must not read as validation.
tbl = make_table(ones(9, 1));           % chi2/N = 1.0, squarely "consistent"
under = compute_chi_squared_report(tbl, 12);   % dof = -3
verifyEqual(tc, under.interpretation, 'consistent_but_underdetermined');

% With real degrees of freedom the plain band label is correct and kept.
ok = compute_chi_squared_report(tbl, 3);       % dof = 6
verifyEqual(tc, ok.interpretation, 'consistent');
end

function test_underfit_verdict_is_never_softened_by_low_dof(tc)
% Failing to match the data DESPITE having excess freedom is a genuine and
% interpretable failure -- qualifying it would weaken a real signal rather
% than prevent an overclaim.
tbl = make_table(3 * ones(9, 1));       % chi2/N = 9, well above 2.0
report = compute_chi_squared_report(tbl, 12);  % dof = -3
verifyEqual(tc, report.interpretation, 'underfit', ...
    'an underfit verdict must survive dof <= 0 unqualified.');
end

function test_dof_le_two_sets_insufficient_dof_note(tc)
tbl = make_table(ones(9, 1));
report = compute_chi_squared_report(tbl, 7); % dof = 2
verifyEqual(tc, report.dof, 2);
verifyEqual(tc, report.dof_note, 'insufficient_dof');

report_ok = compute_chi_squared_report(tbl, 3); % dof = 6
verifyEqual(tc, report_ok.dof, 6);
verifyEqual(tc, report_ok.dof_note, 'sufficient_dof', ...
    ['dof_note must never be empty-char: an empty field made struct2table ', ...
     'throw a row-count-mismatch error in export_chi_squared_summary, ', ...
     'which crashed a completed 6-start calibration run before its ', ...
     'results could be saved (2026-08-29). See compute_chi_squared_report.m.']);
end

function test_only_governed_rows_are_counted(tc)
z = [1, 1, 1, 5]; % last row is an outlier that must be excluded
in_primary = [true, true, true, false];
tbl = make_table(z, in_primary);
report = compute_chi_squared_report(tbl, 2);
verifyEqual(tc, report.n_obs, 3);
verifyEqual(tc, report.chi2, 3, 'AbsTol', 1e-9, ...
    'The non-governed outlier row must not inflate chi2.');
end

function test_nonfinite_rows_are_skipped_and_counted(tc)
tbl = make_table([1, 1, NaN, 1]);
report = compute_chi_squared_report(tbl, 2);
verifyEqual(tc, report.n_obs, 3);
verifyEqual(tc, report.n_skipped, 1);
verifyEqual(tc, report.chi2, 3, 'AbsTol', 1e-9);
end

function test_reproduces_the_measured_reyna_candidate(tc)
% Cross-check against the manually computed values in
% docs/reyna_statistical_calibration_prd.md section 1.2: chi2 = 13.22 at
% N=9, p=7, chi2/N = 1.47.
names = {'RAP_mean';'PAP_min';'PAP_max';'PAP_mean';'SAP_min';'SAP_max'; ...
    'SAP_mean';'QpQs';'CO_Lmin'};
z = [0.15; 1.99; 1.73; 0.67; 1.74; 0.93; 1.16; 0.47; 0.58];
tbl = make_table(z, true(9,1), names);
report = compute_chi_squared_report(tbl, 7);
verifyEqual(tc, report.chi2, 13.22, 'AbsTol', 0.05);
verifyEqual(tc, report.chi2_per_obs, 1.47, 'AbsTol', 0.01);
verifyEqual(tc, report.dof, 2);
verifyEqual(tc, report.worst_metric, 'PAP_min', ...
    'PAP_min at z=1.99 must be identified as the worst residual, not SAP_min at 17%% error.');
end

function test_empty_or_missing_table_does_not_error(tc)
report = compute_chi_squared_report(table(), 7);
verifyEqual(tc, report.n_obs, 0);
verifyEqual(tc, report.interpretation, 'unavailable');

report2 = compute_chi_squared_report([], 7);
verifyEqual(tc, report2.n_obs, 0);
end

function test_missing_required_columns_returns_empty_report_not_error(tc)
tbl = table({'A'; 'B'}, [1; 2], 'VariableNames', {'Metric', 'SomeOtherCol'});
report = compute_chi_squared_report(tbl, 7);
verifyEqual(tc, report.n_obs, 0);
end

function test_report_is_always_struct2table_compatible(tc)
% Direct regression test for the crash of 2026-08-29: struct2table requires
% every field of a scalar struct to represent exactly one row, so a field
% that is sometimes '' (0x0 char) and sometimes a non-empty string (1xN
% char) breaks it with "different numbers of rows". A completed 6-start,
% ~3.4-hour calibration run was lost to exactly this, because the failure
% happened in the final reporting step, downstream of a successful
% optimisation (export_chi_squared_summary in validation_report.m calls
% struct2table on this exact output). If struct2table throws below, the
% test framework fails this test with that error -- no wrapper needed.
names = {'A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'};
tbl = make_table(ones(9,1), true(9,1), names);

struct2table(compute_chi_squared_report(tbl, 7));   % dof=2, "insufficient_dof" path
struct2table(compute_chi_squared_report(tbl, 3));   % dof=6, "sufficient_dof" path
struct2table(compute_chi_squared_report(table(), 7)); % empty early-return path

verifyTrue(tc, true); % reaching here means none of the calls above threw
end
