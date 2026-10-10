%% test_supplement_export.m
% The raw-data export must account for every protocol entry (no entry
% silently ignored), compare fitted rows against the exact calibration
% targets, and write every supplement file.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));
set(groot, 'DefaultFigureVisible', 'off');

run_dir = fullfile(project_root, 'results', 'archive_2026-09-19', 'runs', ...
    '20260906_225713_reyna_pre_surgery');
if ~exist(run_dir, 'dir')
    fprintf('  SKIP: archived run not present locally.\n');
    return;
end

clinical = patient_reyna();
out = fullfile(tempdir, 'test_supplement_export');
if exist(out, 'dir'), rmdir(out, 's'); end
export_run_supplement(run_dir, clinical, 'OutDir', out);

expected = {'protocol_comparison_pre','protocol_comparison_post','inputs', ...
    'parameters_calibration','parameters_full_baseline','parameters_full_calibrated', ...
    'parameters_full_post','gsa_indices_initial','gsa_indices_final', ...
    'gsa_surrogate_quality','gsa_mask','multistart','identifiability_columns', ...
    'identifiability_pairs','identifiability_reduction','timeseries_pre_baseline', ...
    'timeseries_pre_calibrated','timeseries_post_prediction'};
for k = 1:numel(expected)
    assert(exist(fullfile(out, [expected{k} '.csv']), 'file') == 2, ...
        'Missing supplement file: %s', expected{k});
end
assert(exist(fullfile(out, 'run_metadata.json'), 'file') == 2);

% Every protocol entry appears exactly once, with a role.
for form = {'pre', 'post'}
    T = readtable(fullfile(out, ['protocol_comparison_' form{1} '.csv']), 'TextType', 'char');
    rows = T(~isnan(T.ProtocolRow), :);
    assert(height(rows) == numel(clinical.protocol.(form{1})), ...
        '%s: every protocol entry must appear in the comparison table.', form{1});
    assert(all(~cellfun(@isempty, T.Role)), '%s: every row needs a role.', form{1});
end

% Fitted rows are compared against the exact calibration target.
T = readtable(fullfile(out, 'protocol_comparison_pre.csv'), 'TextType', 'char');
map = T(strcmp(T.Variable, 'Mean arterial pressure (RFA, protocol formula)'), :);
assert(abs(map.ComparedValue - 71.3333333) < 1e-5, 'Pre MAP target must be 71.33.');
P = readtable(fullfile(out, 'protocol_comparison_post.csv'), 'TextType', 'char');
hr = readtable(fullfile(out, 'inputs.csv'), 'TextType', 'char');
assert(hr.ValueUsedByModel(strcmp(hr.Input, 'Heart rate, post-closure')) == 108);
assert(any(strcmp(P.ModelMetric, 'LVEF')), 'Post LVEF must be evaluated.');

fprintf('  PASS: all %d supplement files written; every protocol entry has a role.\n', numel(expected) + 1);
