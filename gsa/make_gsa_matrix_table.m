function make_gsa_matrix_table(gsa_out, highlight_thresh, save_fig)
% MAKE_GSA_MATRIX_TABLE
% -----------------------------------------------------------------------
% Render a matrix-style GSA table:
%   - Rows    = uncertain parameters
%   - Columns = output metrics
%   - Cells   = Sobol total-order index (ST), displayed in scientific notation
%   - Yellow background on cells where ST >= highlight_thresh
%
% This replicates the table style from the reference paper figure.
%
% USAGE:
%   make_gsa_matrix_table(gsa_out)
%   make_gsa_matrix_table(gsa_out, 0.1)          % custom threshold
%   make_gsa_matrix_table(gsa_out, 0.1, true)    % also save PNG
%
% INPUTS:
%   gsa_out          - struct from gsa_run_pce (or loaded from .mat)
%   highlight_thresh - (optional) ST threshold for yellow highlight [default 0.1]
%   save_fig         - (optional) true/false to save figure as PNG [default false]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-12
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(highlight_thresh), highlight_thresh = 0.1; end
if nargin < 3 || isempty(save_fig),         save_fig         = false; end

%% =====================================================================
%  1. Collect metrics and parameters from gsa_out
%% =====================================================================
cfg = gsa_out.cfg;

metrics = cfg.all_metrics;   % ordered list of output metrics (columns)
params  = cfg.names;         % ordered list of parameter names (rows)
nM      = numel(metrics);
nP      = numel(params);

% Build ST matrix  [nP x nM]
ST_mat = nan(nP, nM);
for mi = 1:nM
    mf = metrics{mi};
    if isfield(gsa_out, mf) && isfield(gsa_out.(mf), 'ST')
        st_vec = gsa_out.(mf).ST;   % [nP x 1]  same order as cfg.names
        if numel(st_vec) == nP
            ST_mat(:, mi) = st_vec;
        end
    end
end

%% =====================================================================
%  2. Build display labels
%% =====================================================================
param_labels  = format_param_labels(params);
metric_labels = format_metric_labels(metrics);

%% =====================================================================
%  3. Draw the table as a MATLAB figure
%% =====================================================================
% Data-space sizing: keep these small so the axis coordinates are manageable
col_w   = 1.2;    % width per column in data units
row_h   = 0.55;   % height per row in data units
left_w  = 2.4;    % row-label column width
top_h   = 1.2;    % column-header row height

x_total = left_w + nM * col_w;
y_total = top_h  + nP * row_h;

% Open figure maximized so ALL columns and rows are visible immediately
hfig = figure('Name', sprintf('GSA Matrix Table — %s', cfg.scenario), ...
              'NumberTitle', 'off', ...
              'Color', 'w', ...
              'Units', 'normalized', ...
              'OuterPosition', [0 0 1 1]);   % <-- full-screen

% Axes covers full figure with small margins for the title
ax = axes(hfig, ...
    'Units', 'normalized', ...
    'Position', [0.01, 0.06, 0.98, 0.88], ...
    'XLim', [0, x_total], ...
    'YLim', [0, y_total], ...
    'YDir', 'reverse', ...
    'Visible', 'off');
hold(ax, 'on');

% Colors
col_header_bg  = [0.92 0.92 0.92];
row_label_bg   = [1.00 1.00 1.00];
cell_bg_normal = [1.00 1.00 1.00];
cell_bg_hi     = [1.00 1.00 0.50];   % yellow
border_col     = [0.55 0.55 0.55];
text_col       = [0.05 0.05 0.05];

% ---- Column headers (metrics) ----
for mi = 1:nM
    x0c = left_w + (mi-1)*col_w;
    y0c = 0;
    fill_rect(ax, x0c, y0c, col_w, top_h, col_header_bg, border_col);
    text(ax, x0c + col_w/2, top_h/2, metric_labels{mi}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 6.5, 'FontWeight', 'bold', 'Color', text_col, ...
        'Interpreter', 'tex');
end

% ---- Row labels (parameters) ----
for pi = 1:nP
    y0r = top_h + (pi-1)*row_h;
    fill_rect(ax, 0, y0r, left_w, row_h, row_label_bg, border_col);
    text(ax, left_w - 0.12, y0r + row_h/2, param_labels{pi}, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', 7.5, 'Color', text_col, 'Interpreter', 'tex');
end

% ---- Data cells ----
for pi = 1:nP
    for mi = 1:nM
        x0c = left_w + (mi-1)*col_w;
        y0c = top_h  + (pi-1)*row_h;
        val = ST_mat(pi, mi);

        if isnan(val)
            bg    = cell_bg_normal;
            label = '—';
        elseif val >= highlight_thresh
            bg    = cell_bg_hi;
            label = sci_str(val);
        else
            bg    = cell_bg_normal;
            label = sci_str(val);
        end

        fill_rect(ax, x0c, y0c, col_w, row_h, bg, border_col);
        text(ax, x0c + col_w/2, y0c + row_h/2, label, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 6.5, 'Color', text_col, 'Interpreter', 'tex');
    end
end

% Outer border
x_total = left_w + nM*col_w;
y_total = top_h  + nP*row_h;
rectangle(ax, 'Position', [0, 0, x_total, y_total], ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);

title(ax, sprintf('Sobol Total-Order Indices (S_T) — %s', ...
      strrep(cfg.scenario, '_', '\_')), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Visible', 'on');
ax.Visible = 'off';

% ---- Legend annotation ----
annotation(hfig, 'textbox', [0.01 0.01 0.40 0.04], ...
    'String', sprintf('Yellow = S_T \\geq %.2g   |   Values: Sobol Total-Order Index', highlight_thresh), ...
    'FitBoxToText', 'off', 'EdgeColor', 'none', ...
    'FontSize', 7, 'BackgroundColor', 'w');

%% =====================================================================
%  4. Optional save  (high-res PNG — captures all rows + columns)
%% =====================================================================
if save_fig
    out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'gsa');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    fname = fullfile(out_dir, sprintf('gsa_matrix_table_%s.png', cfg.scenario));
    % ContentType='vector' gives exact patch boundaries; Resolution for raster fallback
    exportgraphics(hfig, fname, 'Resolution', 250, 'BackgroundColor', 'white');
    fprintf('[make_gsa_matrix_table] Saved: %s\n', fname);
    fprintf('                        Open the PNG for a full-resolution view.\n');
end

end % make_gsa_matrix_table


%% =========================================================================
%  INTERNAL HELPERS
%% =========================================================================

function fill_rect(ax, x, y, w, h, fc, ec)
% Draw a filled rectangle (patch) in axis ax
patch(ax, [x x+w x+w x], [y y y+h y+h], fc, ...
    'EdgeColor', ec, 'LineWidth', 0.4);
end

function s = sci_str(v)
% Format a scalar as compact scientific notation for TeX: "a.b·10^{n}"
if v == 0
    s = '0';
    return;
end
exp_val = floor(log10(abs(v)));
mant    = v / 10^exp_val;
% Round mantissa to 1 decimal place
mant_r  = round(mant, 1);
if mant_r >= 10  % rounding edge case
    mant_r  = mant_r / 10;
    exp_val = exp_val + 1;
end
% Format: e.g. "4.2\cdot10^{-1}"
s = sprintf('%.1f{\\cdot}10^{%d}', mant_r, exp_val);
end

function labels = format_param_labels(names)
% Convert cfg.names like 'E.LV.EA' → TeX label 'E^{A}_{LV}'
n      = numel(names);
labels = cell(n, 1);
for i = 1:n
    nm = names{i};
    parts = strsplit(nm, '.');
    switch parts{1}
        case 'E'
            % E.LV.EA → E^{A}_{LV}
            chamber = parts{2};   % LV or RV or LA or RA
            subtype = parts{3};   % EA or EB
            sup = strrep(subtype, 'E', '');   % 'A' or 'B'
            labels{i} = sprintf('E^{%s}_{%s}', sup, chamber);
        case 'V0'
            % V0.LV → V^{0}_{LV}
            chamber = parts{2};
            labels{i} = sprintf('V^{0}_{%s}', chamber);
        case 'R'
            sub = strjoin(parts(2:end), '_{');
            close_braces = repmat('}', 1, numel(parts)-2);
            labels{i} = sprintf('R_{%s%s}', sub, close_braces);
        case 'C'
            sub = strjoin(parts(2:end), '_{');
            close_braces = repmat('}', 1, numel(parts)-2);
            labels{i} = sprintf('C_{%s%s}', sub, close_braces);
        otherwise
            labels{i} = strrep(nm, '.', '_{');
    end
end
end

function labels = format_metric_labels(metrics)
% Convert metric names to compact TeX column headers
n      = numel(metrics);
labels = cell(n, 1);
for i = 1:n
    m = metrics{i};
    switch m
        case 'RAP_mean';  labels{i} = 'P^{mean}_{RA}';
        case 'LAP_mean';  labels{i} = 'P^{mean}_{LA}';
        case 'PAP_min';   labels{i} = 'P^{min}_{PA}';
        case 'PAP_max';   labels{i} = 'P^{max}_{PA}';
        case 'PAP_mean';  labels{i} = 'P^{mean}_{PA}';
        case 'PVP_mean';  labels{i} = 'P^{PUL}_{VEN}';
        case 'RVP_min';   labels{i} = 'P^{min}_{RV}';
        case 'RVP_max';   labels{i} = 'P^{max}_{RV}';
        case 'RVP_mean';  labels{i} = 'P^{mean}_{RV}';
        case 'LVP_min';   labels{i} = 'P^{min}_{LV}';
        case 'LVP_max';   labels{i} = 'P^{max}_{LV}';
        case 'LVP_mean';  labels{i} = 'P^{mean}_{LV}';
        case 'SAP_min';   labels{i} = 'P^{min}_{AO}';
        case 'SAP_max';   labels{i} = 'P^{max}_{AO}';
        case 'SAP_mean';  labels{i} = 'P^{mean}_{AO}';
        case 'SVR';       labels{i} = 'R^{SYS}_{AR}';
        case 'PVR';       labels{i} = 'R^{PUL}_{AR}';
        case 'QpQs';      labels{i} = 'Q_{p}/Q_{s}';
        case 'LVEDV';     labels{i} = 'V^{ED}_{LV}';
        case 'LVESV';     labels{i} = 'V^{ES}_{LV}';
        case 'RVEDV';     labels{i} = 'V^{ED}_{RV}';
        case 'RVESV';     labels{i} = 'V^{ES}_{RV}';
        case 'LVEF';      labels{i} = 'EF_{LV}';
        case 'RVEF';      labels{i} = 'EF_{RV}';
        otherwise;        labels{i} = strrep(m, '_', '\_');
    end
end
end
