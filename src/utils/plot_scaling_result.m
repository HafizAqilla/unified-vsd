function plot_scaling_result(params_ref, params_scaled)
% PLOT_SCALING_RESULT
% -----------------------------------------------------------------------
% Visualises the effect of apply_scaling on all scaled parameter groups.
%
% Produces a 6-panel figure showing the reference (adult) vs scaled
% (patient) values for:
%   Panel 1 — Vascular Resistances  R  [mmHg·s/mL]
%   Panel 2 — Vascular Compliances  C  [mL/mmHg]
%   Panel 3 — Chamber Elastances    E  (EA and EB)  [mmHg/mL]
%   Panel 4 — Unstressed Volumes    V0  [mL]
%   Panel 5 — Heart Rate & Blood Volume summary  (scalar quantities)
%   Panel 6 — Initial Condition vector  [mL / mL/s / mmHg]
%
% INPUTS:
%   params_ref    - adult reference struct from default_parameters()
%   params_scaled - patient-specific struct from apply_scaling()
%
% USAGE EXAMPLE:
%   params_ref    = default_parameters();
%   patient.age_years = 2;
%   patient.weight_kg = 12;
%   patient.height_cm = 87;
%   patient.sex       = 'M';
%   params_scaled = apply_scaling(params_ref, patient);
%   plot_scaling_result(params_ref, params_scaled);
%
% OUTPUTS:
%   Figure with 6 panels (displayed; not saved automatically).
%
% SIGN CONVENTIONS / NOTES:
%   - All values are in the project internal unit system (mmHg, mL, s).
%   - EF and HR are scalars; they appear in Panel 5.
%   - The scale factor s and patient info are printed in the figure title.
%
% REFERENCES:
%   [1] apply_scaling.m  — produces params_scaled used here.
%   [2] VIBECODING_GUARDRAILS.md, §11 — Plotting Standards.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-31
% VERSION:  1.0
% -----------------------------------------------------------------------

%% --- shared style -------------------------------------------------------
clr_ref    = [0.18 0.44 0.70];   % blue  — reference adult
clr_scaled = [0.84 0.22 0.22];   % red   — scaled patient
bar_w      = 0.35;               % bar width fraction
font_sz    = 9;

sc  = params_scaled.scaling;
s   = sc.s;
BSA = sc.BSA_patient;
BV  = sc.BV_patient;
age = sc.patient.age_years;
wt  = sc.patient.weight_kg;
ht  = sc.patient.height_cm;

fig_title = sprintf( ...
    'apply\\_scaling result  |  Age %.1f yr, %.1f kg, %.0f cm  |  BSA = %.3f m²  |  s = %.3f  |  BV = %.0f mL', ...
    age, wt, ht, BSA, s);

fig = figure('Name', 'Scaling Result', 'NumberTitle', 'off', ...
             'Color', 'w', 'Position', [60 60 1400 900]);
sgtitle(fig_title, 'FontSize', 11, 'FontWeight', 'bold', 'Interpreter', 'tex');

%% =========================================================================
%  Panel 1 — Vascular Resistances
%% =========================================================================
ax1 = subplot(3, 2, 1);

R_fields = {'SAR','SC','SVEN','PAR','PCOX','PCNO','PVEN'};
R_ref_vals    = arrayfun(@(f) params_ref.R.(f{1}),    R_fields);
R_scaled_vals = arrayfun(@(f) params_scaled.R.(f{1}), R_fields);

helper_grouped_bar(ax1, R_fields, R_ref_vals, R_scaled_vals, ...
    clr_ref, clr_scaled, bar_w);
ylabel(ax1, 'R  [mmHg·s/mL]', 'FontSize', font_sz);
title(ax1,  'Vascular Resistances', 'FontSize', font_sz, 'FontWeight', 'bold');
set(ax1, 'YScale', 'log');   % log scale — spans 4 orders of magnitude

%% =========================================================================
%  Panel 2 — Vascular Compliances
%% =========================================================================
ax2 = subplot(3, 2, 2);

C_fields = {'SAR','SC','SVEN','PAR','PCOX','PCNO','PVEN'};
C_ref_vals    = arrayfun(@(f) params_ref.C.(f{1}),    C_fields);
C_scaled_vals = arrayfun(@(f) params_scaled.C.(f{1}), C_fields);

helper_grouped_bar(ax2, C_fields, C_ref_vals, C_scaled_vals, ...
    clr_ref, clr_scaled, bar_w);
ylabel(ax2, 'C  [mL/mmHg]', 'FontSize', font_sz);
title(ax2,  'Vascular Compliances', 'FontSize', font_sz, 'FontWeight', 'bold');
set(ax2, 'YScale', 'log');

%% =========================================================================
%  Panel 3 — Chamber Elastances  (EA — active; EB — passive)
%% =========================================================================
ax3 = subplot(3, 2, 3);

chambers = {'LV','RV','LA','RA'};
nC = numel(chambers);
E_ref_EA    = zeros(1, nC);  E_ref_EB    = zeros(1, nC);
E_scaled_EA = zeros(1, nC);  E_scaled_EB = zeros(1, nC);
for k = 1:nC
    ch = chambers{k};
    E_ref_EA(k)    = params_ref.E.(ch).EA;
    E_ref_EB(k)    = params_ref.E.(ch).EB;
    E_scaled_EA(k) = params_scaled.E.(ch).EA;
    E_scaled_EB(k) = params_scaled.E.(ch).EB;
end

x = 1:nC;
bw = 0.2;
b1 = bar(ax3, x - bw*1.5, E_ref_EA,    bw, 'FaceColor', clr_ref,    'EdgeColor', 'none');
hold(ax3, 'on');
     bar(ax3, x - bw*0.5, E_scaled_EA, bw, 'FaceColor', clr_scaled,  'EdgeColor', 'none');
     bar(ax3, x + bw*0.5, E_ref_EB,    bw, 'FaceColor', clr_ref,    'EdgeColor', 'none', 'FaceAlpha', 0.45);
     bar(ax3, x + bw*1.5, E_scaled_EB, bw, 'FaceColor', clr_scaled,  'EdgeColor', 'none', 'FaceAlpha', 0.45);
hold(ax3, 'off');

xticks(ax3, x);  xticklabels(ax3, chambers);
ylabel(ax3, 'E  [mmHg/mL]', 'FontSize', font_sz);
title(ax3,  'Chamber Elastances  (solid = EA, faded = EB)', ...
    'FontSize', font_sz, 'FontWeight', 'bold');
legend(ax3, {'Ref EA','Scaled EA','Ref EB','Scaled EB'}, ...
    'Location', 'northwest', 'FontSize', 7);
grid(ax3, 'on');  box(ax3, 'off');

%% =========================================================================
%  Panel 4 — Unstressed Volumes V0
%% =========================================================================
ax4 = subplot(3, 2, 4);

V0_fields = {'RA','RV','LA','LV','SAR','SC','SVEN','PAR','PVEN'};
V0_ref_vals    = arrayfun(@(f) params_ref.V0.(f{1}),    V0_fields);
V0_scaled_vals = arrayfun(@(f) params_scaled.V0.(f{1}), V0_fields);

helper_grouped_bar(ax4, V0_fields, V0_ref_vals, V0_scaled_vals, ...
    clr_ref, clr_scaled, bar_w);
ylabel(ax4, 'V_0  [mL]', 'FontSize', font_sz);
title(ax4,  'Unstressed Volumes  (BV-ratio scaling)', ...
    'FontSize', font_sz, 'FontWeight', 'bold');

%% =========================================================================
%  Panel 5 — Scalar Quantities  (HR, BV, BV_scale, s)
%% =========================================================================
ax5 = subplot(3, 2, 5);

scalar_labels = {'HR [bpm]', 'BV [mL]', 'BV_{scale} [×]', 's [×]'};
ref_scalars    = [params_ref.HR,  4900,               1,              1];
scaled_scalars = [params_scaled.HR, BV, sc.BV_scale, s];

helper_grouped_bar(ax5, scalar_labels, ref_scalars, scaled_scalars, ...
    clr_ref, clr_scaled, bar_w);
title(ax5, 'Scalar Quantities', 'FontSize', font_sz, 'FontWeight', 'bold');
ylabel(ax5, 'Value', 'FontSize', font_sz);

% Annotate the ratio above each pair
x_pos = (1:numel(scalar_labels));
ratios = scaled_scalars ./ ref_scalars;
for k = 1:numel(ratios)
    txt = sprintf('×%.3f', ratios(k));
    text(ax5, x_pos(k), max(ref_scalars(k), scaled_scalars(k)) * 1.05, txt, ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.3 0.3 0.3]);
end

%% =========================================================================
%  Panel 6 — Initial Condition Vector
%% =========================================================================
ax6 = subplot(3, 2, 6);

ic_ref_row    = params_ref.ic.V(:);      % column vector [14×1]
ic_scaled_row = params_scaled.ic.V(:);

n_states = numel(ic_ref_row);
x6 = 1:n_states;

% State labels — matches default_parameters.m idx layout
state_labels = {'V_{RA}','V_{RV}','V_{LA}','V_{LV}', ...
                'V_{SAR}','Q_{SAR}','V_{SC}','V_{SVEN}','Q_{SVEN}', ...
                'V_{PAR}','Q_{PAR}','P_{PC}','V_{PVEN}','Q_{PVEN}'};

bb1 = bar(ax6, x6 - bar_w/2, ic_ref_row,    bar_w, 'FaceColor', clr_ref,   'EdgeColor','none');
hold(ax6, 'on');
      bar(ax6, x6 + bar_w/2, ic_scaled_row, bar_w, 'FaceColor', clr_scaled,'EdgeColor','none');
hold(ax6, 'off');

xticks(ax6, x6);
xticklabels(ax6, state_labels);
ax6.XTickLabelRotation = 40;
ylabel(ax6, 'mL | mL/s | mmHg', 'FontSize', font_sz);
title(ax6, 'Initial Condition Vector  X_0', ...
    'FontSize', font_sz, 'FontWeight', 'bold');
legend(ax6, {'Adult ref','Scaled patient'}, 'Location', 'northeast', 'FontSize', 7);
grid(ax6, 'on');  box(ax6, 'off');

%% --- global legend (Panel 1 / 2) ----------------------------------------
legend(ax1, {'Adult ref','Scaled patient'}, 'Location', 'best', 'FontSize', 7);
legend(ax2, {'Adult ref','Scaled patient'}, 'Location', 'best', 'FontSize', 7);
legend(ax4, {'Adult ref','Scaled patient'}, 'Location', 'best', 'FontSize', 7);
legend(ax5, {'Adult ref','Scaled patient'}, 'Location', 'best', 'FontSize', 7);

fprintf('[plot_scaling_result] Figure rendered: %d panels, %d IC states.\n', ...
    6, n_states);

end  % plot_scaling_result


% =========================================================================
%  LOCAL HELPER — helper_grouped_bar
%  Draws a grouped bar chart on ax for two vectors (ref, scaled).
%  labels : cell array of strings (x-axis tick labels)
%  v_ref, v_sc : numeric row vectors of matching length
% =========================================================================
function helper_grouped_bar(ax, labels, v_ref, v_scaled, clr_ref, clr_sc, bw)
n = numel(labels);
x = 1:n;
bar(ax, x - bw/2, v_ref,    bw, 'FaceColor', clr_ref, 'EdgeColor', 'none');
hold(ax, 'on');
bar(ax, x + bw/2, v_scaled, bw, 'FaceColor', clr_sc,  'EdgeColor', 'none');
hold(ax, 'off');
xticks(ax, x);
xticklabels(ax, labels);
ax.XTickLabelRotation = 35;
grid(ax, 'on');
box(ax, 'off');
end
