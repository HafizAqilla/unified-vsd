function summary_table = compare_reyna_pre_surgery_scaling_recipe()
% COMPARE_REYNA_PRE_SURGERY_SCALING_RECIPE
% -----------------------------------------------------------------------
% Runs Reyna pre-surgery with the explicit calibration recipe under each
% supported scaling mode and writes a compact comparison CSV.
%
% OUTPUTS:
%   summary_table - comparison table for Lundquist and Zhang runs        [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-24
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
cd(root);
addpath(build_clean_project_path(root));

recipe = reyna_pre_surgery();
scaling_modes = recipe.scaling_modes(:)';

previous_scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
previous_gsa = getenv('UNIFIED_VSD_DO_GSA');
previous_fast = getenv('UNIFIED_VSD_FAST_CALIBRATION');
cleanup_env = onCleanup(@() restore_env(previous_scaling_mode, previous_gsa, previous_fast)); %#ok<NASGU>

rows = repmat(init_row(), numel(scaling_modes), 1);
for idx = 1:numel(scaling_modes)
    scaling_mode = scaling_modes{idx};
    setenv('UNIFIED_VSD_SCALING_MODE', scaling_mode);
    setenv('UNIFIED_VSD_DO_GSA', '0');
    setenv('UNIFIED_VSD_FAST_CALIBRATION', '0');

    before_folder = latest_run_folder(root, 'reyna', 'pre_surgery');
    clinical = patient_reyna();
    main_run('pre_surgery', clinical);
    run_folder = latest_run_folder(root, 'reyna', 'pre_surgery');
    if strcmp(run_folder, before_folder)
        error('compare_reyna_scaling:noNewRun', ...
            'No new run folder was created for scaling mode %s.', scaling_mode);
    end

    rows(idx) = summarize_run(scaling_mode, run_folder);
end

summary_table = struct2table(rows);
out_dir = fullfile(root, 'results', 'scaling_recipe_comparison');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
writetable(summary_table, fullfile(out_dir, ...
    sprintf('reyna_pre_surgery_recipe_comparison_%s.csv', datestr(now, 'yyyymmdd_HHMMSS'))));
disp(summary_table);
end

function row = init_row()
row = struct( ...
    'scaling_mode', '', ...
    'run_folder', '', ...
    'recipe_id', '', ...
    'calibration_status', '', ...
    'rollback_applied', NaN, ...
    'primary_governed_rmse', NaN, ...
    'full_transparent_rmse', NaN, ...
    'co_error_pct', NaN, ...
    'primary_fail_count', NaN);
end

function row = summarize_run(scaling_mode, run_folder)
row = init_row();
row.scaling_mode = scaling_mode;
row.run_folder = run_folder;

manifest = parse_manifest(fullfile(run_folder, 'run_manifest.txt'));
row.recipe_id = lookup_manifest(manifest, 'CalibrationRecipeId');
row.calibration_status = lookup_manifest(manifest, 'CalibrationStatus');
row.rollback_applied = str2double(lookup_manifest(manifest, 'RollbackApplied'));

rmse_table = readtable(fullfile(run_folder, 'tables', ...
    'validation_rmse_summary_pre_surgery.csv'), 'TextType', 'string');
primary_idx = strcmp(rmse_table.RMSE_Type, 'primary_governed');
full_idx = strcmp(rmse_table.RMSE_Type, 'full_transparent');
row.primary_governed_rmse = rmse_table.Calibrated(primary_idx);
row.full_transparent_rmse = rmse_table.Calibrated(full_idx);

gate_table = readtable(fullfile(run_folder, 'tables', ...
    'validation_primary_gate_pre_surgery.csv'), 'TextType', 'string');
co_idx = strcmp(gate_table.Metric, 'CO_Lmin');
row.co_error_pct = gate_table.AbsError_pct(co_idx);
row.primary_fail_count = sum(~gate_table.Pass_10pct);
end

function manifest = parse_manifest(file_path)
manifest = containers.Map('KeyType', 'char', 'ValueType', 'char');
lines = splitlines(string(fileread(file_path)));
for idx = 1:numel(lines)
    line = char(lines(idx));
    sep = strfind(line, ':');
    if isempty(sep)
        continue;
    end
    key = strtrim(line(1:sep(1)-1));
    val = strtrim(line(sep(1)+1:end));
    manifest(key) = val;
end
end

function value = lookup_manifest(manifest, key)
if manifest.isKey(key)
    value = manifest(key);
else
    value = '';
end
end

function folder = latest_run_folder(root, patient_label, scenario)
runs_dir = fullfile(root, 'results', 'runs');
folder = '';
if ~exist(runs_dir, 'dir')
    return;
end
pattern = sprintf('*_%s_%s', patient_label, scenario);
listing = dir(fullfile(runs_dir, pattern));
listing = listing([listing.isdir]);
if isempty(listing)
    return;
end
[~, order] = sort([listing.datenum], 'descend');
folder = fullfile(listing(order(1)).folder, listing(order(1)).name);
end

function restore_env(scaling_mode, gsa, fast)
restore_one('UNIFIED_VSD_SCALING_MODE', scaling_mode);
restore_one('UNIFIED_VSD_DO_GSA', gsa);
restore_one('UNIFIED_VSD_FAST_CALIBRATION', fast);
end

function restore_one(name, value)
if isempty(value)
    setenv(name, '');
else
    setenv(name, value);
end
end
