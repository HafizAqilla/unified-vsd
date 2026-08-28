function result = compare_scaling_primary_gap(varargin)
% COMPARE_SCALING_PRIMARY_GAP
% -----------------------------------------------------------------------
% Diagnose why Zhang can be worse than Lundquist on governed Primary RMSE.
% This script separates prior defensibility from empirical calibrated fit.
% -----------------------------------------------------------------------

root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root_dir));
opts = parse_options(varargin{:});

result = struct();
result.PatientID = opts.PatientID;
result.Mode = opts.Mode;
result.Classification = 'mixed_or_inconclusive';
result.SameGovernance = false;
result.Notes = {};

fprintf('[compare_scaling_primary_gap] Patient: %s\n', opts.PatientID);
fprintf('[compare_scaling_primary_gap] Mode: %s\n', opts.Mode);

if strcmp(opts.Mode, 'dry_run')
    result.Notes{end + 1} = 'Dry run only; no run folders compared.';
    fprintf('[compare_scaling_primary_gap] Dry run only. Provide ZhangRunDir and LundquistRunDir for comparison.\n');
    return;
end

if isempty(opts.ZhangRunDir) || isempty(opts.LundquistRunDir)
    error('compare_scaling_primary_gap:missingRunDirs', ...
        'ZhangRunDir and LundquistRunDir are required unless Mode=dry_run.');
end

zhang = load_run_validation(opts.ZhangRunDir);
lundquist = load_run_validation(opts.LundquistRunDir);

result.Zhang = zhang;
result.Lundquist = lundquist;
result.PrimaryErrorDecomposition = build_error_decomposition(zhang, lundquist);
result.Classification = classify_gap(result.PrimaryErrorDecomposition);
result.SameGovernance = manifests_same_governance(zhang.Manifest, lundquist.Manifest);

out_dir = opts.OutputDir;
if isempty(out_dir)
    out_dir = fullfile(root_dir, 'results', 'scaling_primary_gap', opts.PatientID);
end
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
writetable(result.PrimaryErrorDecomposition, ...
    fullfile(out_dir, 'scaling_primary_error_decomposition.csv'));
write_gap_markdown(fullfile(out_dir, 'scaling_primary_gap_diagnostic.md'), result);
fprintf('[compare_scaling_primary_gap] Classification: %s\n', result.Classification);
fprintf('[compare_scaling_primary_gap] Output: %s\n', out_dir);
end

function opts = parse_options(varargin)
parser = inputParser();
addParameter(parser, 'PatientID', 'P1', @(x) ischar(x) || isstring(x));
addParameter(parser, 'Mode', 'compare', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ZhangRunDir', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'LundquistRunDir', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputDir', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.PatientID = char(opts.PatientID);
opts.Mode = lower(strtrim(char(opts.Mode)));
opts.ZhangRunDir = char(opts.ZhangRunDir);
opts.LundquistRunDir = char(opts.LundquistRunDir);
opts.OutputDir = char(opts.OutputDir);
end

function run = load_run_validation(run_dir)
run = struct();
run.RunDir = run_dir;
run.Manifest = parse_manifest(fullfile(run_dir, 'run_manifest.txt'));
candidate_files = {
    fullfile(run_dir, 'tables', 'validation_best_candidate_pre_surgery.csv')
    fullfile(run_dir, 'tables', 'validation_calibrated_pre_surgery.csv')
    fullfile(run_dir, 'tables', 'validation_baseline_pre_surgery.csv')
    };
validation_file = '';
for idx = 1:numel(candidate_files)
    if isfile(candidate_files{idx})
        validation_file = candidate_files{idx};
        break;
    end
end
if isempty(validation_file)
    error('compare_scaling_primary_gap:noValidationTable', ...
        'No validation table found under %s.', run_dir);
end
run.ValidationFile = validation_file;
run.Validation = readtable(validation_file, 'TextType', 'string');
end

function manifest = parse_manifest(path_value)
manifest = struct();
if ~isfile(path_value)
    return;
end
lines = regexp(fileread(path_value), '\r?\n', 'split');
for idx = 1:numel(lines)
    line = strtrim(lines{idx});
    sep = strfind(line, ':');
    if isempty(sep)
        continue;
    end
    key = matlab.lang.makeValidName(strtrim(line(1:sep(1)-1)));
    value = strtrim(line(sep(1)+1:end));
    manifest.(key) = value;
end
end

function tbl = build_error_decomposition(zhang, lundquist)
primary = ["RAP_mean","PAP_mean","SAP_mean","QpQs","CO_Lmin"];
z = normalize_validation_table(zhang.Validation);
l = normalize_validation_table(lundquist.Validation);

metric_col = {};
zhang_err = [];
lundquist_err = [];
delta_abs = [];
driver_col = {};
for idx = 1:numel(primary)
    metric = primary(idx);
    z_err = lookup_error(z, metric);
    l_err = lookup_error(l, metric);
    metric_col{end + 1, 1} = char(metric); %#ok<AGROW>
    zhang_err(end + 1, 1) = z_err; %#ok<AGROW>
    lundquist_err(end + 1, 1) = l_err; %#ok<AGROW>
    delta_abs(end + 1, 1) = abs(z_err) - abs(l_err); %#ok<AGROW>
    if isfinite(z_err) && isfinite(l_err) && abs(z_err) > abs(l_err)
        driver_col{end + 1, 1} = 'zhang_worse'; %#ok<AGROW>
    elseif isfinite(z_err) && isfinite(l_err)
        driver_col{end + 1, 1} = 'lundquist_worse_or_tied'; %#ok<AGROW>
    else
        driver_col{end + 1, 1} = 'missing'; %#ok<AGROW>
    end
end

tbl = table(metric_col, zhang_err, lundquist_err, delta_abs, driver_col, ...
    'VariableNames', {'Metric','ZhangError_pct','LundquistError_pct', ...
    'AbsErrorGap_ZhangMinusLundquist_pct','Driver'});
end

function tbl = normalize_validation_table(tbl)
if ismember('Calibrated', tbl.Properties.VariableNames)
    value_col = 'Calibrated';
elseif ismember('Baseline', tbl.Properties.VariableNames)
    value_col = 'Baseline';
else
    value_col = '';
end
if ~ismember('Error_pct', tbl.Properties.VariableNames) && ~isempty(value_col) && ...
        ismember('Clinical', tbl.Properties.VariableNames)
    tbl.Error_pct = 100 * (tbl.(value_col) - tbl.Clinical) ./ max(abs(tbl.Clinical), 1e-9);
end
end

function value = lookup_error(tbl, metric)
value = NaN;
if ~ismember('Metric', tbl.Properties.VariableNames) || ...
        ~ismember('Error_pct', tbl.Properties.VariableNames)
    return;
end
idx = find(strcmp(string(tbl.Metric), metric), 1);
if ~isempty(idx)
    value = tbl.Error_pct(idx);
end
end

function classification = classify_gap(tbl)
drivers = string(tbl.Metric(strcmp(tbl.Driver, 'zhang_worse')));
if any(ismember(drivers, ["CO_Lmin","QpQs"]))
    classification = 'target_conflict_gap';
elseif height(tbl) > 0 && any(strcmp(tbl.Driver, 'zhang_worse'))
    classification = 'scaling_prior_gap';
else
    classification = 'mixed_or_inconclusive';
end
end

function tf = manifests_same_governance(a, b)
keys = {'CalibrationCaseMode','CalibrationRecipeId','GSAEnabled'};
tf = true;
for idx = 1:numel(keys)
    key = keys{idx};
    if isfield(a, key) && isfield(b, key)
        tf = tf && strcmp(a.(key), b.(key));
    else
        tf = false;
    end
end
end

function write_gap_markdown(path_value, result)
fid = fopen(path_value, 'w');
if fid < 0
    error('compare_scaling_primary_gap:writeFailed', ...
        'Unable to write %s.', path_value);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '# Scaling Primary Gap Diagnostic\n\n');
fprintf(fid, '- PatientID: `%s`\n', result.PatientID);
fprintf(fid, '- Classification: `%s`\n', result.Classification);
fprintf(fid, '- Same governance: `%d`\n\n', result.SameGovernance);
fprintf(fid, ['This diagnostic separates Zhang as a publication prior from ' ...
    'Lundquist-BSA as an empirical comparator. A worse Zhang Primary RMSE ' ...
    'should be explained per metric instead of treated as proof that Zhang ' ...
    'is invalid.\n']);
end

