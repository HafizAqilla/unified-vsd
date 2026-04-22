function plot_post_surgery_validation(sim, params, clinical)
% PLOT_POST_SURGERY_VALIDATION
% -----------------------------------------------------------------------
% Overlays the model's post-surgery prediction against the actual
% post-operative clinical measurements for visual validation.
%
% For each measurable haemodynamic quantity, the plot shows:
%   - Solid coloured line  : model-predicted waveform over one cardiac cycle
%   - Dashed black line(s) : clinical measured value(s) (scalar references)
%
% INPUTS:
%   sim      - struct from integrate_system.m  (.t, .V)
%   params   - parameter struct (post-scaling, post warm-start)
%   clinical - unified clinical struct (from run_post_surgery.m)
%              clinical.post_surgery fields used for reference overlays
%
% PLOTS GENERATED:
%   Figure 1 — Systemic arterial pressure  (P_SAR) vs SAP_sys / SAP_dia / MAP
%   Figure 2 — Pulmonary artery pressure   (P_PAR) vs PAP_sys / PAP_dia / PAP_mean
%   Figure 3 — LV pressure-volume loop     vs LVEDV, LVESV, LVEF
%   Figure 4 — RV pressure-volume loop     vs RVEDV, RVESV, RVEF
%   Figure 5 — Clinical summary bar chart  (Predicted vs Clinical, all metrics)
%
% OUTPUT FILES (vector PDF, Guardrail §11.3):
%   Each run creates a unique timestamped subfolder to prevent overwriting:
%   results/figures/post_surgery_<YYYYMMDD_HHMMSS>/PostSurg_SAP_validation.pdf
%   results/figures/post_surgery_<YYYYMMDD_HHMMSS>/PostSurg_PAP_validation.pdf
%   results/figures/post_surgery_<YYYYMMDD_HHMMSS>/PostSurg_LV_PVloop.pdf
%   results/figures/post_surgery_<YYYYMMDD_HHMMSS>/PostSurg_RV_PVloop.pdf
%   results/figures/post_surgery_<YYYYMMDD_HHMMSS>/PostSurg_SummaryBar.pdf
%
% SIGN CONVENTIONS:
%   Positive Q_VSD = LV->RV (physiologically forward for unrepaired VSD).
%   Post-surgery: R.vsd = 1e6, so Q_VSD ~= 0 mL/s by construction.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis -- plotting conventions, Guardrail S11.
%   [2] compute_clinical_indices.m -- metric definitions and units.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-14
% VERSION:  1.1  (all helpers converted to subfunctions -- no nested fns)
% -----------------------------------------------------------------------

post = clinical.post_surgery;   % post-surgery clinical reference struct

%% Setup: extract last cycle
t    = sim.t(:);
XV   = sim.V;

T_HB      = 60 / params.HR;   % [s]  cardiac period
T1        = t(end);
T0        = T1 - T_HB;
time_mask = (t >= T0) & (t <= T1);
tc        = t(time_mask) - T0;   % [s]  normalised to cycle start at 0

%% Rebuild pressure signals for the trimmed cycle
P = pvsv_signals(tc, XV(time_mask, :), params);

sidx   = params.idx;              % state index struct (Guardrail §7.1)
V_LV_c = XV(time_mask, sidx.V_LV);   % [mL]
V_RV_c = XV(time_mask, sidx.V_RV);   % [mL]

%% Create a timestamped subfolder so each run is preserved separately
%  Format: results/figures/post_surgery_YYYYMMDD_HHMMSS/
%  fileparts x3 from src/utils/ -> src/ -> project root
root_dir  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
fig_dir   = fullfile(root_dir, 'results', 'figures', ...
                     sprintf('post_surgery_%s', timestamp));
mkdir(fig_dir);
fprintf('[plot_post_surgery_validation] Saving all figures to:\n  %s\n', fig_dir);

%% Colour / style constants
lw_model = 2.0;
lw_clin  = 1.4;
fs_label = 11;
fs_title = 12;

col_model = [0.18 0.45 0.72];   % steel blue  -- model prediction
col_clin  = [0.00 0.00 0.00];   % black       -- clinical reference lines
col_clin2 = [0.80 0.20 0.20];   % red         -- secondary clinical lines

%% =========================================================================
%  FIGURE 1 -- Systemic arterial pressure validation
%% =========================================================================
fig1 = make_fig('PostSurg_SAP_validation');

plot(tc, P.SAR, 'Color', col_model, 'LineWidth', lw_model, ...
     'DisplayName', 'Model: P_{SAR}(t)');
hold on;

if ~isnan(safe_get(post, 'SAP_sys_mmHg'))
    yline(post.SAP_sys_mmHg, '--', 'Color', col_clin, 'LineWidth', lw_clin, ...
          'DisplayName', sprintf('Clinical SAP_{sys} = %.1f mmHg', post.SAP_sys_mmHg), ...
          'LabelHorizontalAlignment', 'left', ...
          'Label', sprintf('Clin. SBP = %.1f mmHg', post.SAP_sys_mmHg));
end
if ~isnan(safe_get(post, 'SAP_dia_mmHg'))
    yline(post.SAP_dia_mmHg, '-.', 'Color', col_clin2, 'LineWidth', lw_clin, ...
          'DisplayName', sprintf('Clinical SAP_{dia} = %.1f mmHg', post.SAP_dia_mmHg), ...
          'LabelHorizontalAlignment', 'left', ...
          'Label', sprintf('Clin. DBP = %.1f mmHg', post.SAP_dia_mmHg));
end
if ~isnan(safe_get(post, 'MAP_mmHg'))
    yline(post.MAP_mmHg, ':', 'Color', col_clin, 'LineWidth', lw_clin + 0.4, ...
          'DisplayName', sprintf('Clinical MAP = %.1f mmHg', post.MAP_mmHg), ...
          'LabelHorizontalAlignment', 'right', ...
          'Label', sprintf('Clin. MAP = %.1f mmHg', post.MAP_mmHg));
end

xlabel('Time in cycle (s)',  'FontSize', fs_label, 'FontName', 'Arial');
ylabel('Pressure (mmHg)',    'FontSize', fs_label, 'FontName', 'Arial');
title('Systemic Arterial Pressure -- Model vs Clinical', ...
      'FontSize', fs_title, 'FontName', 'Arial');
legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Arial');
style_ax(gca); grid on;
export_pdf(fig1, fig_dir, 'PostSurg_SAP_validation');

%% =========================================================================
%  FIGURE 2 -- Pulmonary artery pressure validation
%% =========================================================================
fig2 = make_fig('PostSurg_PAP_validation');

plot(tc, P.PAR, 'Color', [0.80 0.20 0.20], 'LineWidth', lw_model, ...
     'DisplayName', 'Model: P_{PAR}(t)');
hold on;

if ~isnan(safe_get(post, 'PAP_sys_mmHg'))
    yline(post.PAP_sys_mmHg, '--', 'Color', col_clin, 'LineWidth', lw_clin, ...
          'DisplayName', sprintf('Clinical PAP_{sys} = %.1f mmHg', post.PAP_sys_mmHg), ...
          'LabelHorizontalAlignment', 'left', ...
          'Label', sprintf('Clin. PAsys = %.1f mmHg', post.PAP_sys_mmHg));
end
if ~isnan(safe_get(post, 'PAP_dia_mmHg'))
    yline(post.PAP_dia_mmHg, '-.', 'Color', col_clin2, 'LineWidth', lw_clin, ...
          'DisplayName', sprintf('Clinical PAP_{dia} = %.2f mmHg', post.PAP_dia_mmHg), ...
          'LabelHorizontalAlignment', 'left', ...
          'Label', sprintf('Clin. PAdia = %.1f mmHg', post.PAP_dia_mmHg));
end
if ~isnan(safe_get(post, 'PAP_mean_mmHg'))
    yline(post.PAP_mean_mmHg, ':', 'Color', col_clin, 'LineWidth', lw_clin + 0.4, ...
          'DisplayName', sprintf('Clinical PAP_{mean} = %.1f mmHg', post.PAP_mean_mmHg), ...
          'LabelHorizontalAlignment', 'right', ...
          'Label', sprintf('Clin. PAP_{mean} = %.1f mmHg', post.PAP_mean_mmHg));
end

xlabel('Time in cycle (s)',  'FontSize', fs_label, 'FontName', 'Arial');
ylabel('Pressure (mmHg)',    'FontSize', fs_label, 'FontName', 'Arial');
title('Pulmonary Artery Pressure -- Model vs Clinical', ...
      'FontSize', fs_title, 'FontName', 'Arial');
legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Arial');
style_ax(gca); grid on;
export_pdf(fig2, fig_dir, 'PostSurg_PAP_validation');

%% =========================================================================
%  FIGURE 3 -- LV Pressure-Volume loop: Model prediction vs Clinical rectangle
%
%  The continuous model PV loop is overlaid with a dashed clinical rectangle
%  constructed from the four known clinical corner points:
%    ED corner: (LVEDV,  P_lv_ed)   where P_lv_ed  ≈ LAP_mean (or 8 mmHg default)
%    ES corner: (LVESV,  P_lv_es)   where P_lv_es  ≈ SAP_sys  (LV peak ≈ systolic BP)
%  The rectangle connects ED and ES via iso-volume lines, approximating the
%  clinical loop geometry from scalar echo/cath measurements alone.
%  Source: Kass DA (1988) Circ 77:1422; Sagawa K (1981) Circ Res 48:451.
%% =========================================================================
fig3 = make_fig('PostSurg_LV_PVloop');

% --- Model predicted loop (continuous) ------------------------------------
plot(V_LV_c, P.LV, 'Color', col_model, 'LineWidth', lw_model, ...
     'DisplayName', 'Model: LV PV loop');
hold on;

% --- Clinical PV rectangle (dashed) --------------------------------------
% Pressure at end-diastole: prefer LAP_mean (filling pressure), else default
if ~isnan(safe_get(post, 'LAP_mean_mmHg'))
    P_lv_ed_clin = post.LAP_mean_mmHg;   % [mmHg]  filling pressure proxy
else
    P_lv_ed_clin = 8;                    % [mmHg]  typical post-op default
end
% Pressure at end-systole: LV peak ≈ systolic arterial pressure
if ~isnan(safe_get(post, 'SAP_sys_mmHg'))
    P_lv_es_clin = post.SAP_sys_mmHg;    % [mmHg]
else
    P_lv_es_clin = NaN;
end

clin_LVEDV = safe_get(post, 'LVEDV_mL');
clin_LVESV = safe_get(post, 'LVESV_mL');

if ~isnan(clin_LVEDV) && ~isnan(clin_LVESV) && ~isnan(P_lv_es_clin)
    % Rectangle corners (clockwise from ED):
    %   (LVEDV, P_ed) -> (LVEDV, P_es) -> (LVESV, P_es) -> (LVESV, P_ed) -> close
    rect_V = [clin_LVEDV, clin_LVEDV, clin_LVESV, clin_LVESV, clin_LVEDV];
    rect_P = [P_lv_ed_clin, P_lv_es_clin, P_lv_es_clin, P_lv_ed_clin, P_lv_ed_clin];
    plot(rect_V, rect_P, 'k--', 'LineWidth', lw_clin + 0.2, ...
         'DisplayName', 'Clinical PV rectangle');

    % ED and ES corner markers
    plot(clin_LVEDV, P_lv_ed_clin, 'ko', 'MarkerFaceColor', 'k', ...
         'MarkerSize', 7, 'DisplayName', ...
         sprintf('Clinical ED (%.0f mL, %.0f mmHg)', clin_LVEDV, P_lv_ed_clin));
    plot(clin_LVESV, P_lv_es_clin, 'rs', 'MarkerFaceColor', col_clin2, ...
         'MarkerSize', 7, 'DisplayName', ...
         sprintf('Clinical ES (%.0f mL, %.0f mmHg)', clin_LVESV, P_lv_es_clin));
end

% EF annotation box
model_LVEF = (max(V_LV_c) - min(V_LV_c)) / max(max(V_LV_c), 1e-6);
str_ef = sprintf('Model LVEF = %.1f%%', model_LVEF * 100);
if ~isnan(safe_get(post, 'LVEF'))
    str_ef = sprintf('%s\nClinical LVEF = %.1f%%', str_ef, post.LVEF * 100);
end
xl = xlim; yl = ylim;
text(xl(1) + 0.60*(xl(2)-xl(1)), yl(1) + 0.10*(yl(2)-yl(1)), str_ef, ...
     'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
     'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.6 0.6 0.6]);

xlabel('V_{LV} (mL)',   'FontSize', fs_label, 'FontName', 'Arial');
ylabel('P_{LV} (mmHg)', 'FontSize', fs_label, 'FontName', 'Arial');
title('LV PV Loop -- Model vs Clinical Rectangle', ...
      'FontSize', fs_title, 'FontName', 'Arial');
legend('Location', 'northwest', 'FontSize', 9, 'FontName', 'Arial');
style_ax(gca); grid on;
export_pdf(fig3, fig_dir, 'PostSurg_LV_PVloop');

%% =========================================================================
%  FIGURE 4 -- RV Pressure-Volume loop: Model prediction vs Clinical rectangle
%
%  Same approach as Figure 3 for the right ventricle:
%    P_rv_ed ≈ RAP_mean  (RV end-diastolic pressure ≈ right atrial pressure)
%    P_rv_es ≈ PAP_sys   (RV peak pressure ≈ PA systolic)
%% =========================================================================
fig4 = make_fig('PostSurg_RV_PVloop');

% --- Model predicted loop (continuous) ------------------------------------
plot(V_RV_c, P.RV, 'Color', [0.55 0.20 0.65], 'LineWidth', lw_model, ...
     'DisplayName', 'Model: RV PV loop');
hold on;

% --- Clinical PV rectangle (dashed) --------------------------------------
if ~isnan(safe_get(post, 'RAP_mean_mmHg'))
    P_rv_ed_clin = post.RAP_mean_mmHg;   % [mmHg]  RVEDP ≈ RAP
else
    P_rv_ed_clin = 5;                    % [mmHg]  typical default
end
if ~isnan(safe_get(post, 'PAP_sys_mmHg'))
    P_rv_es_clin = post.PAP_sys_mmHg;    % [mmHg]  RV peak ≈ PA systolic
else
    P_rv_es_clin = NaN;
end

clin_RVEDV = safe_get(post, 'RVEDV_mL');
clin_RVESV = safe_get(post, 'RVESV_mL');

if ~isnan(clin_RVEDV) && ~isnan(clin_RVESV) && ~isnan(P_rv_es_clin)
    rect_V = [clin_RVEDV, clin_RVEDV, clin_RVESV, clin_RVESV, clin_RVEDV];
    rect_P = [P_rv_ed_clin, P_rv_es_clin, P_rv_es_clin, P_rv_ed_clin, P_rv_ed_clin];
    plot(rect_V, rect_P, 'k--', 'LineWidth', lw_clin + 0.2, ...
         'DisplayName', 'Clinical PV rectangle');

    plot(clin_RVEDV, P_rv_ed_clin, 'ko', 'MarkerFaceColor', 'k', ...
         'MarkerSize', 7, 'DisplayName', ...
         sprintf('Clinical ED (%.0f mL, %.0f mmHg)', clin_RVEDV, P_rv_ed_clin));
    plot(clin_RVESV, P_rv_es_clin, 'rs', 'MarkerFaceColor', col_clin2, ...
         'MarkerSize', 7, 'DisplayName', ...
         sprintf('Clinical ES (%.0f mL, %.0f mmHg)', clin_RVESV, P_rv_es_clin));
end

% EF annotation box
model_RVEF = (max(V_RV_c) - min(V_RV_c)) / max(max(V_RV_c), 1e-6);
str_ef = sprintf('Model RVEF = %.1f%%', model_RVEF * 100);
if ~isnan(safe_get(post, 'RVEF'))
    str_ef = sprintf('%s\nClinical RVEF = %.1f%%', str_ef, post.RVEF * 100);
end
xl = xlim; yl = ylim;
text(xl(1) + 0.60*(xl(2)-xl(1)), yl(1) + 0.10*(yl(2)-yl(1)), str_ef, ...
     'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold', ...
     'BackgroundColor', [1 1 1 0.7], 'EdgeColor', [0.6 0.6 0.6]);

xlabel('V_{RV} (mL)',   'FontSize', fs_label, 'FontName', 'Arial');
ylabel('P_{RV} (mmHg)', 'FontSize', fs_label, 'FontName', 'Arial');
title('RV PV Loop -- Model vs Clinical Rectangle', ...
      'FontSize', fs_title, 'FontName', 'Arial');
legend('Location', 'northwest', 'FontSize', 9, 'FontName', 'Arial');
style_ax(gca); grid on;
export_pdf(fig4, fig_dir, 'PostSurg_RV_PVloop');

%% =========================================================================
%  FIGURE 5 -- Summary bar chart: Predicted vs Clinical (all scalar metrics)
%% =========================================================================
mean_t = @(y) trapz(tc, y) / max(tc(end) - tc(1), 1e-9);

% Compute model scalar metrics
model_SAP_max  = max(P.SAR);           % [mmHg]
model_SAP_dia  = min(P.SAR);           % [mmHg]
model_MAP      = mean_t(P.SAR);        % [mmHg]
model_PAP_sys  = max(P.PAR);           % [mmHg]
model_PAP_dia  = min(P.PAR);           % [mmHg]
model_PAP_mean = mean_t(P.PAR);        % [mmHg]
model_RAP      = mean_t(P.RA);         % [mmHg]
model_LVEF_pct = model_LVEF * 100;     % [%]
model_RVEF_pct = model_RVEF * 100;     % [%]
model_LVEDV    = max(V_LV_c);          % [mL]
model_LVESV    = min(V_LV_c);          % [mL]
model_RVEDV    = max(V_RV_c);          % [mL]
model_RVESV    = min(V_RV_c);          % [mL]

% Build paired lists — skip any metric where clinical value is NaN
metric_labels = {};
model_vals    = [];
clin_vals     = [];

[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'SAP_{sys}',  model_SAP_max,  safe_get(post, 'SAP_sys_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'SAP_{dia}',  model_SAP_dia,  safe_get(post, 'SAP_dia_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'MAP',        model_MAP,      safe_get(post, 'MAP_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'PAP_{sys}',  model_PAP_sys,  safe_get(post, 'PAP_sys_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'PAP_{dia}',  model_PAP_dia,  safe_get(post, 'PAP_dia_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'PAP_{mean}', model_PAP_mean, safe_get(post, 'PAP_mean_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'RAP_{mean}', model_RAP,      safe_get(post, 'RAP_mean_mmHg'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'LVEF (%)',   model_LVEF_pct, safe_get_pct(post, 'LVEF'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'RVEF (%)',   model_RVEF_pct, safe_get_pct(post, 'RVEF'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'LVEDV',      model_LVEDV,   safe_get(post, 'LVEDV_mL'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'LVESV',      model_LVESV,   safe_get(post, 'LVESV_mL'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'RVEDV',      model_RVEDV,   safe_get(post, 'RVEDV_mL'));
[metric_labels, model_vals, clin_vals] = add_row(metric_labels, model_vals, clin_vals, ...
    'RVESV',      model_RVESV,   safe_get(post, 'RVESV_mL'));

if ~isempty(metric_labels)
    n_metrics = numel(metric_labels);
    x = 1:n_metrics;

    fig5 = make_fig('PostSurg_SummaryBar');
    set(fig5, 'Position', [2 2 max(18, n_metrics * 1.8) 12]);

    bar(x - 0.2, model_vals, 0.35, 'FaceColor', [0.18 0.45 0.72], ...
        'DisplayName', 'Model Prediction');
    hold on;
    bar(x + 0.2, clin_vals,  0.35, 'FaceColor', [0.30 0.68 0.33], ...
        'DisplayName', 'Clinical Measurement');

    % Colour-coded % error labels
    for i = 1:n_metrics
        pct_err = 100 * (model_vals(i) - clin_vals(i)) / max(abs(clin_vals(i)), 1e-9);
        if     abs(pct_err) <= 5,  col_err = [0 0.55 0];    % green
        elseif abs(pct_err) <= 15, col_err = [0 0 0];        % black
        else,                      col_err = [0.80 0 0];     % red
        end
        top_y = max(model_vals(i), clin_vals(i));
        text(x(i), top_y * 1.04, sprintf('%.1f%%', pct_err), ...
             'HorizontalAlignment', 'center', 'FontSize', 8.5, ...
             'FontName', 'Arial', 'Color', col_err, 'FontWeight', 'bold');
    end

    xticks(x);
    xticklabels(metric_labels);
    ylabel('Value (mixed units)', 'FontSize', fs_label, 'FontName', 'Arial');
    title('Post-Surgery: Model Prediction vs Clinical Measurements', ...
          'FontSize', fs_title, 'FontName', 'Arial', 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Arial');
    style_ax(gca); grid on; grid minor;

    annotation(fig5, 'textbox', [0.13 0.01 0.75 0.05], ...
        'String', '% error: green <= 5% | black <= 15% | red > 15%', ...
        'FontSize', 9, 'FontName', 'Arial', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center');

    export_pdf(fig5, fig_dir, 'PostSurg_SummaryBar');
else
    fprintf('[plot_post_surgery_validation] No clinical targets provided -- summary bar skipped.\n');
end

fprintf('[plot_post_surgery_validation] All plots complete.\n');
fprintf('  Folder: %s\n', fig_dir);

end  % plot_post_surgery_validation

% =========================================================================
%  SUBFUNCTIONS  (must be after the main function end — no nesting allowed)
% =========================================================================

function style_ax(ax)
% STYLE_AX — apply publication axes style (Guardrail §11.1)
set(ax, 'FontSize', 10, 'FontName', 'Arial', 'Box', 'on', ...
        'TickDir', 'out', 'LineWidth', 0.8);
end

function h = make_fig(name)
% MAKE_FIG — create a standard publication-size figure window [18x12 cm]
h = figure('Name', name, 'Color', 'w', 'NumberTitle', 'off', ...
           'Units', 'centimeters', 'Position', [2 2 18 12]);
end

function export_pdf(fig_handle, fig_dir, fname)
% EXPORT_PDF — save figure as vector PDF (Guardrail §11.3)
out_path = fullfile(fig_dir, [fname '.pdf']);
exportgraphics(fig_handle, out_path, 'ContentType', 'vector', 'Resolution', 300);
fprintf('[plot_post_surgery_validation] Saved: %s\n', out_path);
end

function [ml, mv, cv] = add_row(ml, mv, cv, label, model_v, clin_v)
% ADD_ROW — append a metric pair only if the clinical value is not NaN
if ~isnan(clin_v)
    ml{end+1} = label;    %#ok<AGROW>
    mv(end+1) = model_v;  %#ok<AGROW>
    cv(end+1) = clin_v;   %#ok<AGROW>
end
end

function v = safe_get(s, field)
% SAFE_GET — return NaN if field is missing or already NaN
if isfield(s, field) && ~isnan(s.(field))
    v = s.(field);
else
    v = NaN;
end
end

function v = safe_get_pct(s, field)
% SAFE_GET_PCT — read fraction field and convert to % for display (Guardrail §6.4)
if isfield(s, field) && ~isnan(s.(field))
    v = s.(field) * 100;   % fraction -> %
else
    v = NaN;
end
end

function P = pvsv_signals(tc, XV, params)
% PVSV_SIGNALS — rebuild pressure struct P from a trimmed state window
%   P fields [mmHg]: RA, RV, LA, LV, SAR, SVEN, PAR, PVEN

n    = numel(tc);
sidx = params.idx;   % state index struct (Guardrail §7.1)

P.RA   = zeros(n,1);  P.RV   = zeros(n,1);
P.LA   = zeros(n,1);  P.LV   = zeros(n,1);
P.SAR  = zeros(n,1);  P.SVEN = zeros(n,1);
P.PAR  = zeros(n,1);  P.PVEN = zeros(n,1);

for i = 1:n
    xi     = XV(i,:)';
    V_RA   = xi(sidx.V_RA);    % [mL]
    V_RV   = xi(sidx.V_RV);    % [mL]
    V_LA   = xi(sidx.V_LA);    % [mL]
    V_LV   = xi(sidx.V_LV);    % [mL]
    V_SAR  = xi(sidx.V_SAR);   % [mL]
    V_SVEN = xi(sidx.V_SVEN);  % [mL]
    V_PAR  = xi(sidx.V_PAR);   % [mL]
    V_PVEN = xi(sidx.V_PVEN);  % [mL]

    [E_LV_i, E_RV_i, E_LA_i, E_RA_i] = elastance_model(tc(i), params);

    P.RA(i)   = max(E_RA_i * (V_RA  - params.V0.RA),  -5);           % [mmHg]
    P.RV(i)   = E_RV_i     * (V_RV  - params.V0.RV);                 % [mmHg]
    P.LA(i)   = max(E_LA_i * (V_LA  - params.V0.LA),  -5);           % [mmHg]
    P.LV(i)   = E_LV_i     * (V_LV  - params.V0.LV);                 % [mmHg]
    P.SAR(i)  = (V_SAR  - params.V0.SAR)  / params.C.SAR;            % [mmHg]
    P.SVEN(i) = max((V_SVEN - params.V0.SVEN) / params.C.SVEN, -5);  % [mmHg]
    P.PAR(i)  = (V_PAR  - params.V0.PAR)  / params.C.PAR;            % [mmHg]
    P.PVEN(i) = max((V_PVEN - params.V0.PVEN) / params.C.PVEN, -5);  % [mmHg]

    V_SVEN; V_PVEN; %#ok<VUNUS>
end
end  % pvsv_signals
