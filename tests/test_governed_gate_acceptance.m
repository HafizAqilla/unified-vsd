function tests = test_governed_gate_acceptance()
% TEST_GOVERNED_GATE_ACCEPTANCE
% -----------------------------------------------------------------------
% ACCEPT must cover every metric the governed RMSE is computed over, not
% only the five selected primaries.
%
% Motivating case: the 2026-08-28 fair-prior Zhang baseline was classified
% ACCEPT with primary10_fail=0 while PAP_min sat at 23.83% and two of ten
% governed metrics were outside the band.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 1)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function report = report_fixture(governed_errors)
% REPORT_FIXTURE - report whose five primaries all pass, with the governed
% set's extra rows supplied by the caller.
primary = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
primary_err = [1.4; 3.3; 2.6; 0.2; 9.0];

metric = [primary(:); governed_errors.names(:)];
err = [primary_err; governed_errors.values(:)];
n = numel(metric);

report.primary_metrics = primary;
report.primary_gate = table(primary(:), nan(5,1), nan(5,1), primary_err, ...
    primary_err <= 5, primary_err <= 10, repmat({'accepted_10pct'},5,1), ...
    'VariableNames', {'Metric','Clinical','Model','AbsError_pct', ...
    'Pass_5pct','Pass_10pct','FitBand'});

report.full_metric_gate = struct();
report.full_metric_gate.table = table(metric, repmat({'mmHg'}, n, 1), ...
    repmat({'hard'}, n, 1), ones(n,1), ones(n,1), err, err, ...
    true(n,1), true(n,1), err <= 10, err <= 5, repmat({'none'}, n, 1), ...
    'VariableNames', {'Metric','Unit','Tier','Clinical','Model', ...
    'Error_pct','AbsError_pct','InObjective','InPrimaryRMSE', ...
    'WithinGate','WithinExcellent','Flag'});

report.rmse_baseline = 0.30;
report.rmse_cal = 0.10;
report.rmse_primary_baseline = 0.30;
report.rmse_primary_cal = 0.10;
end

function plausibility = clean_plausibility()
plausibility = struct('n_warning', 0, 'n_fail', 0, 'warning_fraction', 0);
end

function test_accept_when_whole_governed_set_passes(tc)
report = report_fixture(struct('names', {{'PAP_min','SAP_min'}}, ...
    'values', [7.0; 6.5]));
status = classify_calibration_run(report, clean_plausibility());

verifyEqual(tc, status.n_governed_fail, 0);
verifyTrue(tc, status.governed_gate_ok);
verifyEqual(tc, status.label, 'ACCEPT');
end

function test_reject_accept_when_a_governed_metric_fails(tc)
% The regression this test exists for.
report = report_fixture(struct('names', {{'PAP_min','SAP_min'}}, ...
    'values', [23.83; 13.98]));
status = classify_calibration_run(report, clean_plausibility());

verifyEqual(tc, status.n_primary_fail, 0, ...
    'All five selected primaries still pass in this fixture.');
verifyEqual(tc, status.n_governed_fail, 2);
verifyFalse(tc, status.governed_gate_ok);
verifyNotEqual(tc, status.label, 'ACCEPT', ...
    'A run failing governed metrics must not be labelled ACCEPT.');
end

function test_summary_names_the_failing_metrics(tc)
report = report_fixture(struct('names', {{'PAP_min','SAP_min'}}, ...
    'values', [23.83; 13.98]));
status = classify_calibration_run(report, clean_plausibility());

verifySubstring(tc, status.summary, 'governed_gate=5/7');
verifySubstring(tc, status.summary, 'PAP_min');
verifySubstring(tc, status.summary, 'SAP_min');
end

function test_reports_counts_over_the_governed_mask(tc)
report = report_fixture(struct('names', {{'PAP_min','SAP_min','PAP_max'}}, ...
    'values', [23.83; 13.98; 8.0]));
status = classify_calibration_run(report, clean_plausibility());

verifyEqual(tc, status.governed_gate_total, 8);
verifyEqual(tc, status.n_governed_fail, 2);
verifyEqual(tc, sort(status.governed_gate_failures), {'PAP_min','SAP_min'});
end

function test_report_without_gate_keeps_legacy_behaviour(tc)
% Backward compatibility: a bare report must neither pass nor block.
report = report_fixture(struct('names', {{}}, 'values', []));
report = rmfield(report, 'full_metric_gate');
status = classify_calibration_run(report, clean_plausibility());

verifyTrue(tc, status.governed_gate_ok);
verifyEqual(tc, status.n_governed_fail, 0);
verifyEqual(tc, status.governed_gate_total, 0);
verifyEqual(tc, status.label, 'ACCEPT');
end

function test_chi2_is_recorded_but_does_not_gate_accept(tc)
% PRD reyna_statistical_calibration_v1 Phase 2: chi2/N is recorded on status
% for visibility, but must NOT influence the ACCEPT/REJECT label in this
% phase. A run whose percentage gate and plausibility both pass must stay
% ACCEPT even when chi2/N reports 'underfit', and the chi2 fields must still
% be populated so the number is visible in status.summary.
report = report_fixture(struct('names', {{'PAP_min','SAP_min'}}, ...
    'values', [7.0; 6.5]));
report.chi_squared = struct('chi2_per_obs', 9.5, 'interpretation', 'underfit');

status = classify_calibration_run(report, clean_plausibility());

verifyEqual(tc, status.label, 'ACCEPT', ...
    'chi2/N must not gate ACCEPT in Phase 2; it is reporting-only.');
verifyEqual(tc, status.chi2_per_obs, 9.5, 'AbsTol', 1e-9);
verifyEqual(tc, status.chi2_interpretation, 'underfit');
verifySubstring(tc, status.summary, 'chi2_per_obs=9.50 (underfit)');
end

function test_missing_chi_squared_field_is_reported_as_unavailable(tc)
report = report_fixture(struct('names', {{}}, 'values', []));
status = classify_calibration_run(report, clean_plausibility());

verifyTrue(tc, isnan(status.chi2_per_obs));
verifyEqual(tc, status.chi2_interpretation, 'unavailable');
verifyFalse(tc, contains(status.summary, 'chi2_per_obs'), ...
    'Summary must not print a chi2 clause when chi2 is unavailable.');
end
