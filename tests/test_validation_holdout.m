function tests = test_validation_holdout()
% TEST_VALIDATION_HOLDOUT
% -----------------------------------------------------------------------
% Contract tests for print_validation_holdout, added by PRD
% reyna_statistical_calibration_v1 Phase 5.
%
% Deliberately does NOT test SVR as a holdout candidate: see
% tests/test_no_ungoverned_calibration_targets.m
% test_svr_is_excluded_from_calibration_and_primary_rmse for why the PRD's
% original suggestion to relabel SVR validation_holdout was not carried
% out (SVR = (SAP_mean-RAP_mean)/CO_Lmin is algebra over three already-
% fitted targets, not an independent measurement).
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md (Phase 5)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-29
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function tbl = synthetic_gate_table()
metric = {'RAP_mean'; 'PAP_mean'; 'HeldOutMetric'};
tier = {'hard'; 'hard'; 'validation_holdout'};
clinical = [5.0; 15.0; 42.0];
model = [5.1; 14.8; 40.0];
err_pct = 100*(model-clinical)./clinical;
abs_err = abs(err_pct);
within_gate = abs_err <= 10;
tbl = table(metric, tier, clinical, model, err_pct, abs_err, within_gate, ...
    'VariableNames', {'Metric','Tier','Clinical','Model','Error_pct', ...
    'AbsError_pct','WithinGate'});
end

function test_no_holdout_configured_returns_empty_report(tc)
tbl = synthetic_gate_table();
tbl.Tier = repmat({'hard'}, height(tbl), 1); % no holdout rows
out = print_validation_holdout(tbl);
verifyEqual(tc, out.n_total, 0);
verifyEqual(tc, height(out.table), 0);
end

function test_holdout_row_is_isolated_and_counted(tc)
tbl = synthetic_gate_table();
out = print_validation_holdout(tbl);
verifyEqual(tc, out.n_total, 1);
verifyEqual(tc, out.table.Metric{1}, 'HeldOutMetric');
end

function test_within_gate_count_reflects_the_band(tc)
% HeldOutMetric: clinical=42, model=40 -> |error| = 4.76%, within 10% band.
tbl = synthetic_gate_table();
out = print_validation_holdout(tbl);
verifyEqual(tc, out.n_within_gate, 1);
end

function test_multiple_holdout_rows_are_all_reported(tc)
tbl = synthetic_gate_table();
tbl.Tier{2} = 'validation_holdout'; % PAP_mean also held out in this synthetic case
out = print_validation_holdout(tbl);
verifyEqual(tc, out.n_total, 2);
verifyEqual(tc, sort(out.table.Metric), sort({'PAP_mean'; 'HeldOutMetric'}));
end

function test_missing_tier_column_returns_empty_without_error(tc)
tbl = table({'A'}, 1, 'VariableNames', {'Metric', 'SomeOtherCol'});
out = print_validation_holdout(tbl);
verifyEqual(tc, out.n_total, 0);
end

function test_missing_within_gate_column_does_not_error(tc)
tbl = synthetic_gate_table();
tbl = removevars(tbl, 'WithinGate');
out = print_validation_holdout(tbl);
verifyEqual(tc, out.n_total, 1);
verifyEqual(tc, out.n_within_gate, 0, ...
    'Without a WithinGate column the count must default to 0, not error.');
end

function test_current_reyna_recipe_has_no_holdout(tc)
% End-to-end confirmation of the Phase 5 finding: nothing in the real
% Reyna pre-surgery governed set is currently tiered validation_holdout.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

out = print_validation_holdout(table(tbl.Metric, tbl.Tier, ...
    'VariableNames', {'Metric','Tier'}));
verifyEqual(tc, out.n_total, 0);
end
