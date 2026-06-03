function make_gsa_matrix_table(gsa_out, highlight_thresh, save_fig, varargin)
% MAKE_GSA_MATRIX_TABLE
% -----------------------------------------------------------------------
% Render a publication-style GSA heatmap:
%   - Rows    = uncertain parameters, sorted by strongest total effect
%   - Columns = output metrics
%   - Cells   = Sobol total-order index (S_T)
%
% Color carries the main interpretation, while numeric labels are reserved
% for the strongest cells so the figure stays readable in a thesis or paper.
%
% USAGE:
%   make_gsa_matrix_table(gsa_out)
%   make_gsa_matrix_table(gsa_out, 0.1)
%   make_gsa_matrix_table(gsa_out, 0.1, true, 'ResultsDir', out_dir)
%
% INPUTS:
%   gsa_out          - struct from gsa_run_pce/gsa_run_sobol
%   highlight_thresh - S_T threshold for cell labels [default 0.1]
%   save_fig         - true/false to save PNG/PDF/TIFF [default false]
% -----------------------------------------------------------------------

if nargin < 2 || isempty(highlight_thresh), highlight_thresh = 0.1; end
if nargin < 3 || isempty(save_fig), save_fig = false; end
opts = parse_matrix_table_options(varargin{:});

cfg = gsa_out.cfg;
metrics = cfg.all_metrics;
params = cfg.names;
nM = numel(metrics);
nP = numel(params);

ST_mat = nan(nP, nM);
for mi = 1:nM
    mf = metrics{mi};
    if isfield(gsa_out, mf) && isfield(gsa_out.(mf), 'ST')
        st_vec = gsa_out.(mf).ST;
        if numel(st_vec) == nP
            ST_mat(:, mi) = st_vec;
        end
    end
end

row_max = max(ST_mat, [], 2, 'omitnan');
row_max(isnan(row_max)) = -Inf;
[~, row_order] = sort(row_max, 'descend');

ST_plot = ST_mat(row_order, :);
params_plot = params(row_order);
row_max_plot = row_max(row_order);

param_labels = format_param_labels(params_plot);
metric_labels = format_metric_labels(metrics);

fig_w = max(9.5, min(18, 5.6 + 0.40 * nM));
fig_h = max(7.5, min(15, 3.8 + 0.36 * nP));
hfig = figure('Name', sprintf('GSA S_T Heatmap - %s', cfg.scenario), ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1, 1, fig_w, fig_h]);

bottom_margin = max(0.18, min(0.30, 0.12 + 0.007 * nM));
ax = axes(hfig, 'Position', [0.16, bottom_margin, 0.70, 0.84 - bottom_margin]);
imagesc(ax, ST_plot);
set(ax, 'YDir', 'reverse');
colormap(ax, scientific_colormap(256));

cmax = max(ST_plot(:), [], 'omitnan');
if isempty(cmax) || isnan(cmax) || cmax <= 0
    cmax = 1;
end
clim(ax, [0, min(1, max(cmax, highlight_thresh))]);

cb = colorbar(ax);
cb.Label.String = 'Sobol total-order index, S_T';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 9;
cb.TickDirection = 'out';
cb.FontSize = 8.5;

ax.XTick = 1:nM;
ax.XTickLabel = metric_labels;
ax.XTickLabelRotation = 45;
ax.YTick = 1:nP;
ax.YTickLabel = param_labels;
ax.TickLength = [0 0];
ax.FontName = 'Arial';
ax.FontSize = 8.5;
ax.LineWidth = 0.8;
ax.Box = 'on';
ax.TickLabelInterpreter = 'none';
if exist('axtoolbar', 'file') == 2
    axtoolbar(ax, {});
end

xlabel(ax, 'Model output metric', 'FontWeight', 'bold', 'FontSize', 9.5);
ylabel(ax, 'Uncertain parameter', 'FontWeight', 'bold', 'FontSize', 9.5);
title(ax, sprintf('Sobol total-effect sensitivity (%s)', ...
    strrep(cfg.scenario, '_', ' ')), ...
    'FontWeight', 'bold', 'FontSize', 10.5);

hold(ax, 'on');
draw_cell_grid(ax, nM, nP);
annotate_significant_cells(ax, ST_plot, highlight_thresh, row_max_plot);

if save_fig
    out_dir = opts.ResultsDir;
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    base = fullfile(out_dir, sprintf('gsa_st_heatmap_%s', cfg.scenario));
    exportgraphics(hfig, [base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(hfig, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    print(hfig, [base '.tiff'], '-dtiff', '-r600');
    fprintf('[make_gsa_matrix_table] Saved: %s.[png|pdf|tiff]\n', base);
end

end

function opts = parse_matrix_table_options(varargin)
default_dir = getenv('UNIFIED_VSD_GSA_DIR');
if isempty(default_dir)
    default_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'gsa');
end
parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'ResultsDir', default_dir, @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.ResultsDir = char(opts.ResultsDir);
end

function draw_cell_grid(ax, nM, nP)
for x = 0.5:1:(nM + 0.5)
    plot(ax, [x x], [0.5 nP + 0.5], '-', 'Color', [1 1 1], 'LineWidth', 0.5);
end
for y = 0.5:1:(nP + 0.5)
    plot(ax, [0.5 nM + 0.5], [y y], '-', 'Color', [1 1 1], 'LineWidth', 0.5);
end
end

function annotate_significant_cells(ax, ST_plot, highlight_thresh, row_max)
[nP, nM] = size(ST_plot);
for pi = 1:nP
    row_best = row_max(pi);
    for mi = 1:nM
        val = ST_plot(pi, mi);
        if isnan(val)
            continue;
        end
        is_row_best = abs(val - row_best) <= 1e-12 && row_best > 0;
        if val >= highlight_thresh || is_row_best
            if val >= 0.50
                txt_col = [1 1 1];
            else
                txt_col = [0.08 0.08 0.08];
            end
            text(ax, mi, pi, compact_value(val), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontName', 'Arial', ...
                'FontSize', 7.5, ...
                'FontWeight', 'bold', ...
                'Color', txt_col);
        end
    end
end
end

function s = compact_value(v)
if v >= 0.01
    s = sprintf('%.2f', v);
elseif v > 0
    s = sprintf('%.1e', v);
else
    s = '0';
end
end

function cmap = scientific_colormap(n)
if nargin < 1, n = 256; end
anchors = [
    0.95 0.97 1.00
    0.78 0.88 0.96
    0.48 0.71 0.86
    0.20 0.53 0.74
    0.08 0.34 0.56
    0.04 0.18 0.36
];
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'pchip');
cmap = max(0, min(1, cmap));
end

function labels = format_param_labels(names)
n = numel(names);
labels = cell(n, 1);
for i = 1:n
    nm = names{i};
    parts = strsplit(nm, '.');
    switch parts{1}
        case 'E'
            labels{i} = sprintf('E%s_%s', erase(parts{3}, 'E'), parts{2});
        case 'V0'
            labels{i} = sprintf('V0_%s', parts{2});
        otherwise
            labels{i} = strrep(nm, '.', '_');
    end
end
end

function labels = format_metric_labels(metrics)
n = numel(metrics);
labels = cell(n, 1);
for i = 1:n
    m = metrics{i};
    switch m
        case 'RAP_mean'; labels{i} = 'RAP mean';
        case 'LAP_mean'; labels{i} = 'LAP mean';
        case 'PAP_min'; labels{i} = 'PAP min';
        case 'PAP_max'; labels{i} = 'PAP max';
        case 'PAP_mean'; labels{i} = 'PAP mean';
        case 'PVP_mean'; labels{i} = 'PVP mean';
        case 'RVP_min'; labels{i} = 'RVP min';
        case 'RVP_max'; labels{i} = 'RVP max';
        case 'RVP_mean'; labels{i} = 'RVP mean';
        case 'LVP_min'; labels{i} = 'LVP min';
        case 'LVP_max'; labels{i} = 'LVP max';
        case 'LVP_mean'; labels{i} = 'LVP mean';
        case 'SAP_min'; labels{i} = 'SAP min';
        case 'SAP_max'; labels{i} = 'SAP max';
        case 'SAP_mean'; labels{i} = 'SAP mean';
        case 'QpQs'; labels{i} = 'Qp/Qs';
        otherwise; labels{i} = strrep(m, '_', ' ');
    end
end
end
