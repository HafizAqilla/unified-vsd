function fig_paths = plot_calibration_error_comparison(report, scenario, figures_dir, gate_pct)
% PLOT_CALIBRATION_ERROR_COMPARISON
% -----------------------------------------------------------------------
% Grouped baseline-versus-calibrated absolute error chart with the patient
% acceptance band drawn as a reference line, plus a paired before/after
% slope chart showing the direction each metric moved.
%
% Replaces the dense per-metric error tables in the manuscript: the same
% information reads faster as a chart, and the acceptance band becomes a
% visible threshold rather than a number the reader has to compare against
% mentally.
%
% INPUTS:
%   report      - struct from validation_report                          [-]
%   scenario    - 'pre_surgery' | 'post_surgery'                         [-]
%   figures_dir - output directory for exported figures            [char]
%   gate_pct    - acceptance band, default 10                          [%]
%
% OUTPUTS:
%   fig_paths   - cellstr of written figure files
%
% ASSUMPTIONS:
%   - Only rows with a finite clinical comparator are plotted; a metric with
%     no patient measurement has no error to show.
%   - Metrics are ordered by calibrated error so the worst case reads first.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 4 || isempty(gate_pct) || ~isfinite(gate_pct)
    gate_pct = 10;
end
fig_paths = {};

if ~isstruct(report) || ~isfield(report, 'table_delta') || isempty(report.table_delta)
    return;
end

tbl = report.table_delta;
valid = ~isnan(tbl.BaseErr_pct) & ~isnan(tbl.CalErr_pct);
if ~any(valid)
    return;
end
tbl = tbl(valid, :);

base_abs = abs(tbl.BaseErr_pct);
cal_abs = abs(tbl.CalErr_pct);
[~, order] = sort(cal_abs, 'descend');
tbl = tbl(order, :);
base_abs = base_abs(order);
cal_abs = cal_abs(order);
labels = strrep(tbl.Metric, '_', '\_');

if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

fig_paths{end+1} = plot_grouped_bars(labels, base_abs, cal_abs, gate_pct, ...
    scenario, figures_dir);
fig_paths{end+1} = plot_paired_slopes(labels, base_abs, cal_abs, gate_pct, ...
    scenario, figures_dir);

fprintf('[plot_calibration_error_comparison] Figures written:\n');
for idx = 1:numel(fig_paths)
    fprintf('  %s\n', fig_paths{idx});
end
end

% =========================================================================
function out_path = plot_grouped_bars(labels, base_abs, cal_abs, gate_pct, scenario, figures_dir)
% PLOT_GROUPED_BARS - baseline vs calibrated absolute error per metric.
fig = figure('Visible', 'off', 'Position', [100, 100, 1100, 520]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

bars = bar(ax, [base_abs(:), cal_abs(:)], 'grouped');
bars(1).FaceColor = [0.62 0.68 0.72];
bars(2).FaceColor = [0.06 0.42 0.46];
bars(1).EdgeColor = 'none';
bars(2).EdgeColor = 'none';

yline(ax, gate_pct, '--', sprintf('%g%% acceptance band', gate_pct), ...
    'Color', [0.62 0.23 0.18], 'LineWidth', 1.4, ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 10);

set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
    'XTickLabelRotation', 40, 'TickDir', 'out', 'Box', 'off');
ylabel(ax, 'Absolute error (%)');
title(ax, sprintf('Calibration effect per metric — %s', ...
    strrep(scenario, '_', ' ')));
legend(ax, {'Baseline', 'Calibrated'}, 'Location', 'northeast', 'Box', 'off');
grid(ax, 'on');
ax.YGrid = 'on';
ax.XGrid = 'off';
hold(ax, 'off');

out_path = fullfile(figures_dir, ...
    sprintf('calibration_error_comparison_%s.png', scenario));
exportgraphics(fig, out_path, 'Resolution', 200);
close(fig);
end

% =========================================================================
function out_path = plot_paired_slopes(labels, base_abs, cal_abs, gate_pct, scenario, figures_dir)
% PLOT_PAIRED_SLOPES - before/after slope chart; direction of change per metric.
fig = figure('Visible', 'off', 'Position', [100, 100, 760, 620]);
ax = axes(fig); %#ok<LAXES>
hold(ax, 'on');

n = numel(labels);
improved = [0.10 0.45 0.32];
worsened = [0.62 0.23 0.18];

for idx = 1:n
    if cal_abs(idx) <= base_abs(idx)
        line_color = improved;
    else
        line_color = worsened;
    end
    plot(ax, [1, 2], [base_abs(idx), cal_abs(idx)], '-o', ...
        'Color', line_color, 'MarkerFaceColor', line_color, ...
        'MarkerSize', 5, 'LineWidth', 1.3);
    text(ax, 2.04, cal_abs(idx), labels{idx}, 'FontSize', 9, ...
        'VerticalAlignment', 'middle');
end

yline(ax, gate_pct, '--', sprintf('%g%% band', gate_pct), ...
    'Color', [0.62 0.23 0.18], 'LineWidth', 1.4, 'FontSize', 10);

xlim(ax, [0.85, 2.6]);
set(ax, 'XTick', [1, 2], 'XTickLabel', {'Baseline', 'Calibrated'}, ...
    'TickDir', 'out', 'Box', 'off');
ylabel(ax, 'Absolute error (%)');
title(ax, sprintf('Before/after calibration — %s', ...
    strrep(scenario, '_', ' ')));
grid(ax, 'on');
ax.YGrid = 'on';
ax.XGrid = 'off';
hold(ax, 'off');

out_path = fullfile(figures_dir, ...
    sprintf('calibration_error_slopes_%s.png', scenario));
exportgraphics(fig, out_path, 'Resolution', 200);
close(fig);
end
