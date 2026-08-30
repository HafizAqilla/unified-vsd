function tests = test_full_metric_gate_export()
% TEST_FULL_METRIC_GATE_EXPORT
% -----------------------------------------------------------------------
% Contract tests for export_full_metric_gate: disclosed denominator,
% correct gate flags, worst-metric identification, and detection of
% metrics that are graded but never fitted.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function report = synthetic_report()
% SYNTHETIC_REPORT - fixed table with known errors and governance flags.
metric = {'A_ok'; 'B_edge'; 'C_fail'; 'D_nocomp'; 'E_graded_unfitted'};
unit   = {'mmHg'; 'mmHg'; 'L/min'; 'mL'; 'mmHg'};
tier   = {'hard'; 'soft'; 'hard'; 'unavailable'; 'validation_only'};
clin   = [100; 100; 100; NaN; 100];
cal    = [102; 110; 125; NaN; 94];
err    = [2; 10; 25; NaN; -6];
in_cal = [true; true; true; false; false];
in_pri = [true; true; true; false; true];
flag   = {'none'; 'none'; 'none'; 'none'; 'none'};

report.table_cal = table(metric, unit, tier, clin, cal, err, ...
    in_cal, in_pri, flag, ...
    'VariableNames', {'Metric','Unit','Tier','Clinical','Calibrated', ...
    'Error_pct','IncludedInCalibration','IncludedInPrimaryRMSE','Flag'});
report.table_baseline = report.table_cal;
end

function test_excludes_rows_without_clinical_comparator(tc)
gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', '');
verifyEqual(tc, gate.n_total, 4, ...
    'Rows with a NaN clinical comparator must be excluded.');
verifyFalse(tc, any(strcmp(gate.table.Metric, 'D_nocomp')), ...
    'A row without a patient comparator must not appear in the gate table.');
end

function test_gate_flags_are_inclusive_at_the_band(tc)
gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', '');

idx_edge = strcmp(gate.table.Metric, 'B_edge');
verifyTrue(tc, gate.table.WithinGate(idx_edge), ...
    'An error exactly at the 10%% band must count as within the gate.');

idx_fail = strcmp(gate.table.Metric, 'C_fail');
verifyFalse(tc, gate.table.WithinGate(idx_fail), ...
    'An error above the band must fall outside the gate.');

verifyEqual(tc, gate.n_within_gate, 3);
verifyEqual(tc, gate.n_within_excellent, 1, ...
    'Only the 2%% row is inside the 5%% excellent band.');
end

function test_primary_mask_counts_are_separate_from_totals(tc)
gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', '');
verifyEqual(tc, gate.n_primary_total, 4);
verifyEqual(tc, gate.n_primary_within, 3);
end

function test_detects_metric_graded_but_not_fitted(tc)
gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', '');
verifyEqual(tc, gate.n_graded_not_fitted, 1, ...
    'A metric inside the primary RMSE but outside the objective must be flagged.');
end

function test_worst_metric_and_ordering(tc)
gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', '');
verifyEqual(tc, gate.worst_metric, 'C_fail');
verifyEqual(tc, gate.worst_abs_pct, 25, 'AbsTol', 1e-12);
verifyEqual(tc, gate.table.Metric{1}, 'C_fail', ...
    'Table must be sorted worst-first.');
verifyEqual(tc, gate.outside_gate_metrics, {'C_fail'});
end

function test_custom_gate_band_is_honoured(tc)
gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', '', 5);
verifyEqual(tc, gate.gate_pct, 5);
verifyEqual(tc, gate.n_within_gate, 1, ...
    'At a 5%% band only the 2%% row passes.');
end

function test_csv_is_written_when_directory_supplied(tc)
out_dir = fullfile(tempdir, sprintf('vsd_gate_%s', ...
    datestr(now, 'yyyymmddHHMMSSFFF'))); %#ok<TNOW1,DATST>
cleanup = onCleanup(@() rmdir_if_present(out_dir));

gate = export_full_metric_gate(synthetic_report(), 'pre_surgery', out_dir);
verifyTrue(tc, isfile(gate.csv_path), 'Gate CSV must be written.');

written = readtable(gate.csv_path);
verifyEqual(tc, height(written), 4);
verifyTrue(tc, ismember('AbsError_pct', written.Properties.VariableNames));
end

function test_empty_report_is_handled(tc)
gate = export_full_metric_gate(struct(), 'pre_surgery', '');
verifyEqual(tc, gate.n_total, 0);
verifyEmpty(tc, gate.worst_metric);
end

function rmdir_if_present(path_str)
if exist(path_str, 'dir')
    rmdir(path_str, 's');
end
end
