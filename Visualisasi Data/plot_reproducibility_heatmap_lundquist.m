% =========================================================================
% PLOT_REPRODUCIBILITY_HEATMAP_LUNDQUIST
% =========================================================================
% Generates a reproducibility heatmap of post-surgery calibrated parameters
% across 5 simulations (Sim_1_L to Sim_5_L) for patient Reyna (Lundquist).
%
% HOW TO RUN:
%   cd('C:\Users\asus\Documents\VSD Main\unified-vsd-main\Visualisasi Data')
%   plot_reproducibility_heatmap_lundquist
% =========================================================================

clear; clc; close all;

%% -----------------------------------------------------------------------
%  SECTION 0 — Paths
%% -----------------------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
data_root  = fullfile(script_dir, 'Hasil Fix', 'Reyna', 'Lundquist');
output_dir = fullfile(script_dir, 'output_figures');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% -----------------------------------------------------------------------
%  SECTION 1 — Define parameters and simulation labels
%% -----------------------------------------------------------------------
target_params = {
    'group.R\_sys\_scale'; ...
    'group.R\_pul\_scale'; ...
    'R.SVEN';              ...
    'C.SAR';               ...
    'E.LV.EA';             ...
    'E.LV.EB';             ...
    'E.RV.EA';             ...
    'E.RV.EB';             ...
    'E.LA.EA';             ...
    'E.RA.EA';             ...
    'V0.LV';               ...
    'V0.RV'                ...
};

% Raw names for CSV lookup (no TeX escaping)
param_keys = strrep(target_params, '\_', '_');

n_params    = numel(target_params);
sim_labels  = {'Sim 1','Sim 2','Sim 3','Sim 4','Sim 5'};
n_sims      = numel(sim_labels);

% Decimal places per parameter (to match reference image style)
fmt_dec = [4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 3];

%% -----------------------------------------------------------------------
%  SECTION 2 — Read calibrated values from each simulation folder
%% -----------------------------------------------------------------------
data_matrix = NaN(n_params, n_sims);

sim_folders = dir(fullfile(data_root, 'Sim_*_L'));
sim_folders = sim_folders([sim_folders.isdir]);
sim_names   = sort({sim_folders.name});

for s = 1:n_sims
    if s > numel(sim_names); break; end
    sim_dir   = fullfile(data_root, sim_names{s});
    post_runs = dir(fullfile(sim_dir, '*_post_surgery'));
    post_runs = post_runs([post_runs.isdir]);
    if isempty(post_runs); continue; end
    tbl_dir   = fullfile(sim_dir, post_runs(1).name, 'tables');

    pvals = containers.Map();

    % Primary: model_derived_parameter_findings
    f1 = fullfile(tbl_dir, 'model_derived_parameter_findings_post_surgery.csv');
    if isfile(f1)
        T = readtable(f1, 'TextType', 'string');
        for r = 1:height(T)
            k = char(T.Parameter(r));
            if ismember(k, param_keys)
                pvals(k) = double(T.FittedValue(r));
            end
        end
    end

    % Fallback: param_comparison_detailed
    f2 = fullfile(tbl_dir, 'param_comparison_detailed_post_surgery.csv');
    if isfile(f2)
        T2 = readtable(f2, 'TextType', 'string');
        for r = 1:height(T2)
            k = char(T2.Parameter(r));
            if ismember(k, param_keys) && ~isKey(pvals, k)
                v = str2double(char(T2.PostOp_Calibrated(r)));
                if ~isnan(v); pvals(k) = v; end
            end
        end
    end

    for p = 1:n_params
        if isKey(pvals, param_keys{p})
            data_matrix(p, s) = pvals(param_keys{p});
        end
    end
end

%% -----------------------------------------------------------------------
%  SECTION 3 — Compute SD column
%% -----------------------------------------------------------------------
sd_col = std(data_matrix, 0, 2, 'omitnan');   % 1 SD per row

%% -----------------------------------------------------------------------
%  SECTION 4 — Build combined matrix  [data | SD]
%% -----------------------------------------------------------------------
full_matrix = [data_matrix, sd_col];
n_cols      = n_sims + 1;                     % 5 sims + 1 SD

%% -----------------------------------------------------------------------
%  SECTION 5 — Draw heatmap using imagesc + text overlay
%% -----------------------------------------------------------------------
fig = figure('Color','w', 'Units','centimeters', 'Position',[2 2 26 14]);

% ---- colour map: white -> amber/gold (for Sim columns) -----------------
n_cmap = 256;
cmap_sim = [linspace(0.98,0.14,n_cmap)', ...
            linspace(0.94,0.52,n_cmap)', ...
            linspace(0.82,0.07,n_cmap)'];

% ---- axes with generous left margin for parameter labels ---------------
ax = axes('Parent', fig, ...
          'Units', 'normalized', ...
          'Position', [0.22 0.08 0.74 0.72]);

% Flip colour limits slightly to ensure visible gradient
clo = min(data_matrix(:), [], 'omitnan');
chi = max(data_matrix(:), [], 'omitnan');
if chi == clo; chi = clo + 1; end

% ----- draw each cell manually ----------------------------------------
hold(ax, 'on');

for col = 1:n_cols
    for row = 1:n_params
        val = full_matrix(row, col);

        % Cell background colour
        if col <= n_sims
            % Sim column: amber gradient
            frac = (val - clo) / (chi - clo);
            frac = max(0, min(1, frac));
            ci   = max(1, round(frac*(n_cmap-1)) + 1);
            fc   = cmap_sim(ci, :);
            tc   = [0.10 0.10 0.10];   % dark text
        else
            % SD column: light grey
            fc   = [0.88 0.88 0.88];
            tc   = [0.25 0.25 0.25];
        end

        % Rectangle  (x=col-1 .. col,  y=row-1 .. row)
        rectangle(ax, 'Position', [col-1, row-1, 1, 1], ...
                  'FaceColor', fc, ...
                  'EdgeColor', [0.65 0.65 0.65], ...
                  'LineWidth',  0.5);

        % Text label
        dec = fmt_dec(row);
        if col <= n_sims
            txt = sprintf(['%.' num2str(dec) 'f'], val);
        else
            txt = sprintf('%.6f', val);   % SD always 6 dp
        end
        text(ax, col-0.5, row-0.5, txt, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment',   'middle', ...
             'FontSize',  8.5, ...
             'FontWeight','bold', ...
             'Color', tc);
    end
end

% ---- axes limits & orientation -----------------------------------------
set(ax, 'XLim',         [0, n_cols], ...
        'YLim',         [0, n_params], ...
        'YDir',         'reverse', ...          % row 1 at top
        'TickLength',   [0 0], ...
        'XAxisLocation','top', ...
        'FontSize',      9, ...
        'Box',           'on', ...
        'XColor',       [0.3 0.3 0.3], ...
        'YColor',       [0.3 0.3 0.3]);

% ---- column header labels (X tick labels) ------------------------------
col_labels = [sim_labels, {'SD'}];
set(ax, 'XTick',       (0.5 : 1 : n_cols-0.5), ...
        'XTickLabel',  col_labels, ...
        'XTickLabelRotation', 0);

% Style SD header differently
% (done via drawing a text over the tick label position)
text(ax, n_cols-0.5, -0.3, 'SD', ...
     'HorizontalAlignment','center', ...
     'FontSize', 9, 'FontWeight','bold', 'Color',[0.35 0.35 0.35]);

% ---- row labels (Y tick labels = parameter names) ----------------------
set(ax, 'YTick',       (0.5 : 1 : n_params-0.5), ...
        'YTickLabel',  target_params, ...
        'TickLabelInterpreter', 'tex');

% ---- title & subtitle --------------------------------------------------
title_str = '\bfHasil Parameter Terkalibrasi ';


title(ax, title_str, 'Interpreter','tex', 'FontSize',11, 'Units','normalized', ...
      'Position',[0.5 1.07 0]);
annotation(fig, 'textbox', [0 0.845 1 0.04], ...
           'String', sub_str, ...
           'HorizontalAlignment','center', ...
           'VerticalAlignment','middle', ...
           'FontSize', 7.5, 'Color',[0.45 0.45 0.45], ...
           'EdgeColor','none', 'Interpreter','none');

hold(ax, 'off');
ax.Toolbar.Visible = 'off';   % suppress toolbar from exported image

%% -----------------------------------------------------------------------
%  SECTION 6 — Save
%% -----------------------------------------------------------------------
out_png = fullfile(output_dir, 'reproducibility_heatmap_lundquist.png');
out_pdf = fullfile(output_dir, 'reproducibility_heatmap_lundquist.pdf');

exportgraphics(fig, out_png, 'Resolution', 300);
exportgraphics(fig, out_pdf, 'ContentType', 'vector');

fprintf('\n[Done] Heatmap saved:\n  %s\n  %s\n', out_png, out_pdf);
