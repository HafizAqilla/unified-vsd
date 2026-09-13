function result = run_reyna_scaling_experiment(experiment_ref, varargin)
% RUN_REYNA_SCALING_EXPERIMENT
% -----------------------------------------------------------------------
% Execute the locked Reyna scaling experiment contract through main_run.
%
% This runner is intentionally orchestration-only: it never changes the
% clinical profile, target governance, or acceptance thresholds.  Each arm
% is run with explicit environment metadata and every repeat is discovered
% from its own dated run folder.  Failed arms remain in the exported table.
%
% USAGE:
%   result = run_reyna_scaling_experiment( ...
%       'config/experiments/reyna_p1_scaling_v1.m');
%
% OPTIONS:
%   DryRun       - write the contract without launching MATLAB simulation
%   Arms         - cellstr subset of arm IDs; default all four
%   Repeats      - repeat count; default contract value (3)
%   DoGSA        - run PCE GSA; default true
%   GsaN         - PCE training samples; default locked contract value (128)
%   MaxFunctionEvaluations - per calibration stage budget
%   MaxIterations           - per calibration stage budget
%   FastCalibration         - use the existing fast calibration schedule
%   ScreeningMode           - skip optional polish stages during screening
%   DoPlots / DoOverlay      - plotting controls; defaults false
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

if nargin < 1 || isempty(experiment_ref)
    experiment_ref = fullfile(root, 'config', 'experiments', ...
        'reyna_p1_scaling_v1.m');
end
experiment = load_experiment_contract(experiment_ref, root);
opts = parse_options(experiment, varargin{:});

experiment_run_id = sprintf('%s_%s', experiment.id, datestr(now, 'yyyymmdd_HHMMSS'));
experiment_dir = fullfile(root, 'results', 'luna_experiments', experiment_run_id);
if ~exist(experiment_dir, 'dir')
    mkdir(experiment_dir);
end

contract_file = fullfile(experiment_dir, 'experiment_contract.txt');
write_contract(contract_file, experiment, opts);

if opts.DryRun
    result = dry_run_result(experiment, opts, experiment_dir, contract_file);
    write_result_exports(result, experiment_dir);
    fprintf('[run_reyna_scaling_experiment] Dry run only; no simulation launched.\n');
    return;
end

env_state = capture_environment();
cleanup_env = onCleanup(@() restore_environment(env_state)); %#ok<NASGU>

rows = cell(0, 24);
run_records = cell(0, 1);
for arm_idx = 1:numel(experiment.arms)
    arm = experiment.arms(arm_idx);
    if ~ismember(arm.id, opts.arms)
        continue;
    end

    for repeat_idx = 1:opts.repeats
        fprintf('\n============================================================\n');
        fprintf('[run_reyna_scaling_experiment] Arm %s | repeat %d/%d\n', ...
            arm.id, repeat_idx, opts.repeats);
        fprintf('============================================================\n');

        configure_environment(arm, experiment, opts);
        before_names = run_folder_names(root, experiment.scenario, experiment.patient_label);
        record = make_run_record(arm, repeat_idx, opts);

        try
            main_run(experiment.scenario, patient_reyna());
            after_names = run_folder_names(root, experiment.scenario, experiment.patient_label);
            new_names = setdiff(after_names, before_names, 'stable');
            if isempty(new_names)
                error('run_reyna_scaling_experiment:noRunFolder', ...
                    'main_run completed without creating a new dated run folder.');
            end
            record.run_dir = fullfile(root, 'results', 'runs', new_names{end});
            record = enrich_record_from_run(record, record.run_dir, experiment);
        catch ME
            record.status = 'failed';
            record.error_identifier = ME.identifier;
            record.error_message = ME.message;
            fprintf(2, '[run_reyna_scaling_experiment] FAILED %s repeat %d: %s\n', ...
                arm.id, repeat_idx, ME.message);
        end

        rows(end + 1, :) = record_to_row(record); %#ok<AGROW>
        run_records{end + 1, 1} = record; %#ok<AGROW>
    end
end

result = struct();
result.experiment = experiment;
result.options = opts;
result.experiment_dir = experiment_dir;
result.contract_file = contract_file;
result.records = run_records;
result.summary_table = records_to_table(rows);
result.decision = build_decision_summary(result.summary_table, experiment);

write_result_exports(result, experiment_dir);
fprintf('\n[run_reyna_scaling_experiment] Experiment exports:\n  %s\n', experiment_dir);
disp(result.summary_table);
fprintf('[run_reyna_scaling_experiment] Decision: %s\n', result.decision.label);
end

function experiment = load_experiment_contract(experiment_ref, root)
if isstruct(experiment_ref)
    experiment = experiment_ref;
    return;
end
experiment_ref = char(experiment_ref);
if contains(experiment_ref, filesep) || endsWith(experiment_ref, '.m')
    if ~isabsolute_path(experiment_ref)
        experiment_ref = fullfile(root, experiment_ref);
    end
    [folder, stem] = fileparts(experiment_ref);
    addpath(folder);
    factory = str2func(stem);
    experiment = factory();
else
    factory = str2func(experiment_ref);
    experiment = factory();
end
validate_contract(experiment);
end

function validate_contract(experiment)
required = {'id','version','scenario','primary_metrics','gsa_training_samples','arms'};
for idx = 1:numel(required)
    if ~isfield(experiment, required{idx})
        error('run_reyna_scaling_experiment:invalidContract', ...
            'Experiment contract is missing field %s.', required{idx});
    end
end
if ~strcmp(experiment.scenario, 'pre_surgery')
    error('run_reyna_scaling_experiment:invalidScenario', ...
        'This runner is locked to Reyna pre_surgery.');
end
if numel(experiment.primary_metrics) ~= 5 || ...
        ~isequal(experiment.primary_metrics, {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'})
    error('run_reyna_scaling_experiment:primaryMetricsChanged', ...
        'The Reyna primary metric contract must remain unchanged.');
end
arm_ids = {experiment.arms.id};
if numel(unique(arm_ids)) ~= numel(arm_ids)
    error('run_reyna_scaling_experiment:duplicateArm', ...
        'Experiment arm IDs must be unique.');
end
end

function opts = parse_options(experiment, varargin)
parser = inputParser();
addParameter(parser, 'DryRun', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'Arms', {experiment.arms.id}, @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(parser, 'Repeats', experiment.default_repeats, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'DoGSA', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'GsaN', experiment.gsa_training_samples, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'MaxFunctionEvaluations', 300, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'MaxIterations', 40, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(parser, 'FastCalibration', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'ScreeningMode', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'DoPlots', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'DoOverlay', false, @(x) islogical(x) || isnumeric(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.DryRun = logical(opts.DryRun);
opts.DoGSA = logical(opts.DoGSA);
opts.FastCalibration = logical(opts.FastCalibration);
opts.ScreeningMode = logical(opts.ScreeningMode);
opts.DoPlots = logical(opts.DoPlots);
opts.DoOverlay = logical(opts.DoOverlay);
opts.repeats = round(opts.Repeats);
opts.gsa_n = round(opts.GsaN);
opts.max_function_evaluations = round(opts.MaxFunctionEvaluations);
opts.max_iterations = round(opts.MaxIterations);
if ischar(opts.Arms) || isstring(opts.Arms)
    opts.arms = cellstr(opts.Arms);
else
    opts.arms = cellstr(string(opts.Arms));
end
valid_ids = {experiment.arms.id};
unknown = setdiff(opts.arms, valid_ids);
if ~isempty(unknown)
    error('run_reyna_scaling_experiment:unknownArm', ...
        'Unknown experiment arm(s): %s.', strjoin(unknown, ', '));
end
end

function configure_environment(arm, experiment, opts)
setenv('UNIFIED_VSD_RUN_MODE', 'exploratory');
setenv('UNIFIED_VSD_SCALING_MODE', arm.scaling_mode);
setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS', bool_text(arm.historical_seeds_disabled));
setenv('UNIFIED_VSD_FREEZE_CLINICAL_PROFILE', '1');
setenv('UNIFIED_VSD_DO_GSA', bool_text(opts.DoGSA));
setenv('UNIFIED_VSD_GSA_PCE_N', num2str(opts.gsa_n));
setenv('UNIFIED_VSD_FAST_CALIBRATION', bool_text(opts.FastCalibration));
setenv('UNIFIED_VSD_SCREENING_MODE', bool_text(opts.ScreeningMode));
setenv('UNIFIED_VSD_DO_PLOTS', bool_text(opts.DoPlots));
setenv('UNIFIED_VSD_DO_OVERLAY', bool_text(opts.DoOverlay));
setenv('UNIFIED_VSD_MAX_FUN_EVALS', num2str(opts.max_function_evaluations));
setenv('UNIFIED_VSD_MAX_ITERATIONS', num2str(opts.max_iterations));
setenv('UNIFIED_VSD_FMINCON_PARALLEL', '0');
setenv('UNIFIED_VSD_EXPERIMENT_ID', experiment.id);
end

function record = make_run_record(arm, repeat_idx, opts)
record = struct();
record.arm_id = arm.id;
record.arm_label = arm.label;
record.repeat = repeat_idx;
record.scaling_mode = arm.scaling_mode;
record.historical_seeds_disabled = arm.historical_seeds_disabled;
record.expected_seed_policy = arm.expected_seed_policy;
record.gsa_n = opts.gsa_n;
record.max_function_evaluations = opts.max_function_evaluations;
record.max_iterations = opts.max_iterations;
record.status = 'not_started';
record.run_dir = '';
record.manifest_path = '';
record.primary_rmse = NaN;
record.full_rmse = NaN;
record.hard_rmse = NaN;
record.soft_rmse = NaN;
record.baseline_primary_rmse = NaN;
record.metrics_within_gate = NaN;
record.metrics_total = NaN;
record.governed_within_gate = NaN;
record.governed_total = NaN;
record.worst_metric = '';
record.worst_abs_error_pct = NaN;
record.calibration_status = '';
record.calibration_accepted = false;
record.manifest_scaling_mode = '';
record.manifest_seed_policy = '';
record.clinical_profile_policy = '';
record.manifest_target_governance = '';
record.error_identifier = '';
record.error_message = '';
end

function record = enrich_record_from_run(record, run_dir, experiment)
record.status = 'completed';
record.manifest_path = fullfile(run_dir, 'run_manifest.txt');
if exist(record.manifest_path, 'file') ~= 2
    error('run_reyna_scaling_experiment:missingManifest', ...
        'Run folder is missing run_manifest.txt: %s', run_dir);
end
manifest = fileread(record.manifest_path);
record.calibration_status = manifest_value(manifest, 'CalibrationStatus');
record.calibration_accepted = strcmpi(record.calibration_status, 'ACCEPT');
record.manifest_scaling_mode = manifest_value(manifest, 'ScalingMode');
record.manifest_seed_policy = manifest_value(manifest, 'HistoricalSeedPolicy');
record.clinical_profile_policy = manifest_value(manifest, 'ClinicalProfilePolicy');
record.manifest_target_governance = manifest_value(manifest, 'TargetGovernance');
if ~strcmpi(record.manifest_scaling_mode, record.scaling_mode)
    error('run_reyna_scaling_experiment:scalingMismatch', ...
        'Manifest scaling mode %s does not match requested arm mode %s.', ...
        record.manifest_scaling_mode, record.scaling_mode);
end
if record.historical_seeds_disabled && ...
        ~strcmp(record.manifest_seed_policy, 'disabled_fair_prior')
    error('run_reyna_scaling_experiment:seedPolicyMismatch', ...
        'Fair-prior arm did not report disabled_fair_prior.');
end
if ~strcmp(record.clinical_profile_policy, 'frozen_input_profile')
    error('run_reyna_scaling_experiment:clinicalProfileNotFrozen', ...
        'Experiment run did not report frozen_input_profile.');
end
verify_clinical_profile_lock(run_dir);
verify_target_governance(run_dir, experiment.primary_metrics);

rmse_path = fullfile(run_dir, 'tables', ...
    'validation_rmse_summary_pre_surgery.csv');
if exist(rmse_path, 'file') ~= 2
    error('run_reyna_scaling_experiment:missingRmseSummary', ...
        'Run folder is missing validation RMSE summary: %s', rmse_path);
end
rmse_table = readtable(rmse_path);
record.primary_rmse = rmse_value(rmse_table, 'primary_governed', 'Calibrated');
record.full_rmse = rmse_value(rmse_table, 'full_transparent', 'Calibrated');
record.hard_rmse = rmse_value(rmse_table, 'hard_only', 'Calibrated');
record.soft_rmse = rmse_value(rmse_table, 'soft_only', 'Calibrated');
record.baseline_primary_rmse = rmse_value(rmse_table, 'primary_governed', 'Baseline');
if ~isfinite(record.primary_rmse)
    error('run_reyna_scaling_experiment:nonfinitePrimaryRmse', ...
        'Calibrated primary RMSE is not finite for %s.', run_dir);
end

record = enrich_record_from_metric_gate(record, run_dir);

    record.status = 'validated';
end

function record = enrich_record_from_metric_gate(record, run_dir)
% ENRICH_RECORD_FROM_METRIC_GATE - carry the n-of-N acceptance count into the
% experiment summary, so an arm can never be compared on RMSE alone.
gate_path = fullfile(run_dir, 'tables', 'full_metric_gate_pre_surgery.csv');
if exist(gate_path, 'file') ~= 2
    error('run_reyna_scaling_experiment:missingMetricGate', ...
        ['Run folder is missing the per-metric acceptance table: %s. ', ...
         'Every reported acceptance claim requires it.'], gate_path);
end

gate_tbl = readtable(gate_path);
record.metrics_total = height(gate_tbl);
record.metrics_within_gate = nnz(logical(gate_tbl.WithinGate));

governed = logical(gate_tbl.InPrimaryRMSE);
record.governed_total = nnz(governed);
record.governed_within_gate = nnz(governed & logical(gate_tbl.WithinGate));

[worst_pct, worst_ix] = max(gate_tbl.AbsError_pct);
if ~isempty(worst_ix) && isfinite(worst_pct)
    record.worst_metric = gate_tbl.Metric{worst_ix};
    record.worst_abs_error_pct = worst_pct;
end
end

function names = run_folder_names(root, scenario, patient_label)
pattern = ['*_', patient_label, '_', scenario];
entries = dir(fullfile(root, 'results', 'runs', pattern));
entries = entries([entries.isdir]);
names = {entries.name};
names = sort(names);
end

function value = rmse_value(tbl, rmse_type, column_name)
value = NaN;
if isempty(tbl) || ~ismember('RMSE_Type', tbl.Properties.VariableNames) || ...
        ~ismember(column_name, tbl.Properties.VariableNames)
    return;
end
row = strcmp(string(tbl.RMSE_Type), rmse_type);
if any(row)
    candidate = tbl.(column_name)(find(row, 1, 'first'));
    if isnumeric(candidate)
        value = candidate;
    end
end
end

function value = manifest_value(manifest, key)
value = '';
tokens = regexp(manifest, [key ':\s*([^\r\n]+)'], 'tokens', 'once');
if ~isempty(tokens)
    value = strtrim(tokens{1});
end
end

function verify_clinical_profile_lock(run_dir)
package_files = dir(fullfile(run_dir, 'mat', 'run_package_pre_surgery_*.mat'));
if isempty(package_files)
    error('run_reyna_scaling_experiment:missingRunPackage', ...
        'Run package is missing from %s.', run_dir);
end
[~, idx] = max([package_files.datenum]);
package = load(fullfile(package_files(idx).folder, package_files(idx).name), 'run_package');
if ~isfield(package, 'run_package') || ~isfield(package.run_package, 'clinical')
    error('run_reyna_scaling_experiment:missingClinicalSnapshot', ...
        'Run package does not contain a clinical snapshot.');
end
input_clinical = patient_reyna();
stored_clinical = package.run_package.clinical;
sections = {'common','pre_surgery'};
for section_idx = 1:numel(sections)
    section = sections{section_idx};
    if ~isfield(stored_clinical, section) || ~isfield(input_clinical, section)
        error('run_reyna_scaling_experiment:clinicalSectionMismatch', ...
            'Clinical section %s is missing from the run snapshot.', section);
    end
    names = fieldnames(input_clinical.(section));
    for field_idx = 1:numel(names)
        name = names{field_idx};
        if ~isfield(stored_clinical.(section), name) || ...
                ~isequaln(stored_clinical.(section).(name), input_clinical.(section).(name))
            error('run_reyna_scaling_experiment:clinicalValueChanged', ...
                'Clinical value changed in %s.%s.', section, name);
        end
    end
end
end

function verify_target_governance(run_dir, primary_metrics)
path_name = fullfile(run_dir, 'tables', 'validation_target_tiers_pre_surgery.csv');
if exist(path_name, 'file') ~= 2
    error('run_reyna_scaling_experiment:missingTargetTiers', ...
        'Target-tier export is missing from %s.', run_dir);
end
tbl = readtable(path_name);
for idx = 1:numel(primary_metrics)
    row = find(strcmp(string(tbl.Metric), primary_metrics{idx}), 1, 'first');
    if isempty(row) || ~strcmp(string(tbl.Tier(row)), 'hard') || ...
            ~logical(tbl.IncludedInPrimaryRMSE(row))
        error('run_reyna_scaling_experiment:targetGovernanceChanged', ...
            'Primary metric %s is not governed as a hard primary target.', primary_metrics{idx});
    end
end
pvr_row = find(strcmp(string(tbl.Metric), 'PVR'), 1, 'first');
if ~isempty(pvr_row) && logical(tbl.IncludedInPrimaryRMSE(pvr_row))
    error('run_reyna_scaling_experiment:derivedTargetPromoted', ...
        'Derived PVR was promoted into primary RMSE.');
end
end

function rows = record_to_row(record)
rows = {record.arm_id, record.arm_label, record.repeat, record.scaling_mode, ...
    record.historical_seeds_disabled, record.expected_seed_policy, record.gsa_n, ...
    record.max_function_evaluations, record.max_iterations, record.status, ...
    record.primary_rmse, record.full_rmse, record.hard_rmse, record.soft_rmse, ...
    record.baseline_primary_rmse, record.governed_within_gate, ...
    record.governed_total, record.metrics_within_gate, record.metrics_total, ...
    record.worst_metric, record.worst_abs_error_pct, record.calibration_status, ...
    record.calibration_accepted, record.run_dir};
end

function tbl = records_to_table(rows)
names = {'ArmId','ArmLabel','Repeat','ScalingMode','HistoricalSeedsDisabled', ...
    'ExpectedSeedPolicy','GsaN','MaxFunctionEvaluations','MaxIterations','Status', ...
    'PrimaryRMSE','FullRMSE','HardRMSE','SoftRMSE','BaselinePrimaryRMSE', ...
    'GovernedWithinGate','GovernedTotal','MetricsWithinGate','MetricsTotal', ...
    'WorstMetric','WorstAbsErrorPct','CalibrationStatus','CalibrationAccepted','RunDir'};
if isempty(rows)
    tbl = table(strings(0, 1), strings(0, 1), zeros(0, 1), strings(0, 1), ...
        false(0, 1), strings(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        strings(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), ...
        nan(0, 1), nan(0, 1), nan(0, 1), nan(0, 1), strings(0, 1), nan(0, 1), ...
        strings(0, 1), false(0, 1), strings(0, 1), 'VariableNames', names);
    return;
end
tbl = table( ...
    string(rows(:, 1)), string(rows(:, 2)), cell2mat(rows(:, 3)), ...
    string(rows(:, 4)), logical(cell2mat(rows(:, 5))), string(rows(:, 6)), ...
    cell2mat(rows(:, 7)), cell2mat(rows(:, 8)), cell2mat(rows(:, 9)), ...
    string(rows(:, 10)), cell2mat(rows(:, 11)), cell2mat(rows(:, 12)), ...
    cell2mat(rows(:, 13)), cell2mat(rows(:, 14)), cell2mat(rows(:, 15)), ...
    cell2mat(rows(:, 16)), cell2mat(rows(:, 17)), cell2mat(rows(:, 18)), ...
    cell2mat(rows(:, 19)), string(rows(:, 20)), cell2mat(rows(:, 21)), ...
    string(rows(:, 22)), logical(cell2mat(rows(:, 23))), string(rows(:, 24)), ...
    'VariableNames', names);
end

function decision = build_decision_summary(tbl, experiment)
decision = struct('label', 'insufficient_evidence', 'reason', '', ...
    'fair_winner', '', 'fair_selection_basis', '', ...
    'operational_winner', '', 'operational_selection_basis', '');
if isempty(tbl) || ~ismember('Status', tbl.Properties.VariableNames)
    decision.reason = 'No experiment records were produced.';
    return;
end
valid = strcmp(tbl.Status, 'validated') & isfinite(tbl.PrimaryRMSE);
valid_tbl = tbl(valid, :);
if isempty(valid_tbl)
    decision.reason = 'No validated primary RMSE records were produced.';
    return;
end
fair = valid_tbl(valid_tbl.HistoricalSeedsDisabled == true, :);
operational = valid_tbl(valid_tbl.HistoricalSeedsDisabled == false, :);
if ~isempty(fair)
    [decision.fair_winner, decision.fair_selection_basis] = select_family_winner(fair);
end
if ~isempty(operational)
    [decision.operational_winner, decision.operational_selection_basis] = ...
        select_family_winner(operational);
end
if ~isempty(decision.fair_winner) && ~isempty(decision.operational_winner)
    decision.label = 'report_fair_and_operational_results_separately';
    decision.reason = ['The runner reports fair-prior and operational winners ', ...
        'separately. Within each family, ACCEPT candidates are preferred; ', ...
        'validated RMSE is only a fallback when no ACCEPT candidate exists.'];
else
    decision.reason = 'One or both comparison families are incomplete.';
end
decision.validated_records = height(valid_tbl);
decision.contract_primary_metrics = experiment.primary_metrics;
end

function [winner, basis] = select_family_winner(tbl)
accepted = strcmpi(string(tbl.CalibrationStatus), 'ACCEPT');
if any(accepted)
    candidate_tbl = tbl(accepted, :);
    basis = 'accepted_calibration_status';
else
    candidate_tbl = tbl;
    basis = 'validated_rmse_fallback_no_accept';
end
means = group_mean_primary(candidate_tbl);
if isempty(means)
    winner = '';
else
    [~, idx] = min(means.MeanPrimaryRMSE);
    winner = char(means.ArmId(idx));
end
end

function summary = group_mean_primary(tbl)
ids = unique(tbl.ArmId, 'stable');
summary = table(strings(numel(ids), 1), nan(numel(ids), 1), ...
    'VariableNames', {'ArmId','MeanPrimaryRMSE'});
for idx = 1:numel(ids)
    summary.ArmId(idx) = ids(idx);
    values = tbl.PrimaryRMSE(strcmp(tbl.ArmId, ids(idx)));
    summary.MeanPrimaryRMSE(idx) = mean(values, 'omitnan');
end
end

function result = dry_run_result(experiment, opts, experiment_dir, contract_file)
result = struct();
result.experiment = experiment;
result.options = opts;
result.experiment_dir = experiment_dir;
result.contract_file = contract_file;
result.records = {};
result.summary_table = records_to_table({});
result.decision = struct('label', 'dry_run', 'reason', 'No simulation launched.');
end

function write_contract(path_name, experiment, opts)
fid = fopen(path_name, 'w');
if fid < 0
    error('run_reyna_scaling_experiment:contractOpenFailed', ...
        'Unable to write experiment contract: %s', path_name);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Luna Reyna Scaling Experiment Contract\n');
fprintf(fid, '======================================\n');
fprintf(fid, 'ExperimentId: %s\n', experiment.id);
fprintf(fid, 'Version: %s\n', experiment.version);
fprintf(fid, 'Scenario: %s\n', experiment.scenario);
fprintf(fid, 'PrimaryMetrics: %s\n', strjoin(experiment.primary_metrics, ','));
fprintf(fid, 'GsaTrainingSamples: %d\n', experiment.gsa_training_samples);
fprintf(fid, 'SelectedArms: %s\n', strjoin(opts.arms, ','));
fprintf(fid, 'Repeats: %d\n', opts.repeats);
fprintf(fid, 'DoGSA: %d\n', opts.DoGSA);
fprintf(fid, 'GsaN: %d\n', opts.gsa_n);
fprintf(fid, 'MaxFunctionEvaluations: %d\n', opts.max_function_evaluations);
fprintf(fid, 'MaxIterations: %d\n', opts.max_iterations);
fprintf(fid, 'ScreeningMode: %d\n', opts.ScreeningMode);
fprintf(fid, 'TargetGovernance: %s\n', experiment.target_governance);
fprintf(fid, 'ForbiddenChanges: %s\n', strjoin(experiment.forbidden_changes, '; '));
for idx = 1:numel(experiment.arms)
    arm = experiment.arms(idx);
    fprintf(fid, 'Arm_%s: mode=%s; historical_seeds_disabled=%d; policy=%s\n', ...
        arm.id, arm.scaling_mode, arm.historical_seeds_disabled, arm.expected_seed_policy);
end
end

function write_result_exports(result, experiment_dir)
summary_csv = fullfile(experiment_dir, 'reyna_scaling_experiment_summary.csv');
writetable(result.summary_table, summary_csv);
decision_file = fullfile(experiment_dir, 'decision_summary.txt');
fid = fopen(decision_file, 'w');
if fid < 0
    error('run_reyna_scaling_experiment:decisionOpenFailed', ...
        'Unable to write decision summary: %s', decision_file);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'DecisionLabel: %s\n', result.decision.label);
if isfield(result.decision, 'reason')
    fprintf(fid, 'Reason: %s\n', result.decision.reason);
end
if isfield(result.decision, 'fair_winner')
    fprintf(fid, 'FairWinner: %s\n', result.decision.fair_winner);
end
if isfield(result.decision, 'fair_selection_basis')
    fprintf(fid, 'FairSelectionBasis: %s\n', result.decision.fair_selection_basis);
end
if isfield(result.decision, 'operational_winner')
    fprintf(fid, 'OperationalWinner: %s\n', result.decision.operational_winner);
end
if isfield(result.decision, 'operational_selection_basis')
    fprintf(fid, 'OperationalSelectionBasis: %s\n', result.decision.operational_selection_basis);
end
if isfield(result.decision, 'validated_records')
    fprintf(fid, 'ValidatedRecords: %d\n', result.decision.validated_records);
end
fprintf(fid, 'PrimaryMetrics: %s\n', strjoin(result.experiment.primary_metrics, ','));
end

function state = capture_environment()
names = {'UNIFIED_VSD_RUN_MODE','UNIFIED_VSD_SCALING_MODE', ...
    'UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS','UNIFIED_VSD_FREEZE_CLINICAL_PROFILE', ...
    'UNIFIED_VSD_DO_GSA', ...
    'UNIFIED_VSD_UQLAB_PATH','UNIFIED_VSD_SOBIOS_PATH', ...
    'UNIFIED_VSD_GSA_PCE_N','UNIFIED_VSD_FAST_CALIBRATION', ...
    'UNIFIED_VSD_SCREENING_MODE', ...
    'UNIFIED_VSD_DO_PLOTS','UNIFIED_VSD_DO_OVERLAY', ...
    'UNIFIED_VSD_MAX_FUN_EVALS','UNIFIED_VSD_MAX_ITERATIONS', ...
    'UNIFIED_VSD_FMINCON_PARALLEL','UNIFIED_VSD_EXPERIMENT_ID'};
state.names = names;
state.values = cellfun(@getenv, names, 'UniformOutput', false);
end

function restore_environment(state)
for idx = 1:numel(state.names)
    setenv(state.names{idx}, state.values{idx});
end
end

function text = bool_text(value)
if value
    text = '1';
else
    text = '0';
end
end

function tf = isabsolute_path(path_name)
tf = ~isempty(regexp(path_name, '^[A-Za-z]:[\\/]', 'once')) || ...
    startsWith(path_name, filesep);
end

function project_path = build_clean_project_path(root)
root_paths = strsplit(genpath(root), pathsep);
root_paths = root_paths(~cellfun('isempty', root_paths));
is_shadow = contains(root_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
    contains(root_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
    contains(root_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, root_paths);
project_path = strjoin(root_paths(~is_shadow & is_existing), pathsep);
end
