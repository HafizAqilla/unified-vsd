function plot_overlay_comparison(sim_a, params_a, label_a, sim_b, params_b, label_b, scenario, varargin)
% PLOT_OVERLAY_COMPARISON
% -----------------------------------------------------------------------
% Plot direct-overlap comparison (two curves in one canvas) for key
% pressure traces over the last cardiac cycle.
%
% INPUTS:
%   sim_a, params_a, label_a - first simulation/parameter set/label
%   sim_b, params_b, label_b - second simulation/parameter set/label
%   scenario                 - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   <ResultsDir>/Overlay_Chambers_<scenario>_<label_a>_vs_<label_b>.pdf
%   <ResultsDir>/Overlay_Vascular_<scenario>_<label_a>_vs_<label_b>.pdf
%
% AUTHOR: Unified VSD Model
% DATE:   2026-04-08
% -----------------------------------------------------------------------

opts = parse_overlay_options(varargin{:});
fig_dir = opts.ResultsDir;
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

[ta, Pa] = last_cycle_pressures(sim_a, params_a);
[tb, Pb] = last_cycle_pressures(sim_b, params_b);

fa = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 12]);
subplot(2,2,1);
plot(ta, Pa.RA, 'b-', 'LineWidth', 1.8, 'DisplayName', ['RA ' label_a]); hold on;
plot(tb, Pb.RA, 'b--', 'LineWidth', 1.8, 'DisplayName', ['RA ' label_b]);
xlabel('Time (s)'); ylabel('P_{RA} (mmHg)'); title('Right Atrium'); grid on; legend('Location', 'best');

subplot(2,2,2);
plot(ta, Pa.RV, 'r-', 'LineWidth', 1.8, 'DisplayName', ['RV ' label_a]); hold on;
plot(tb, Pb.RV, 'r--', 'LineWidth', 1.8, 'DisplayName', ['RV ' label_b]);
xlabel('Time (s)'); ylabel('P_{RV} (mmHg)'); title('Right Ventricle'); grid on; legend('Location', 'best');

subplot(2,2,3);
plot(ta, Pa.LA, 'c-', 'LineWidth', 1.8, 'DisplayName', ['LA ' label_a]); hold on;
plot(tb, Pb.LA, 'c--', 'LineWidth', 1.8, 'DisplayName', ['LA ' label_b]);
xlabel('Time (s)'); ylabel('P_{LA} (mmHg)'); title('Left Atrium'); grid on; legend('Location', 'best');

subplot(2,2,4);
plot(ta, Pa.LV, 'm-', 'LineWidth', 1.8, 'DisplayName', ['LV ' label_a]); hold on;
plot(tb, Pb.LV, 'm--', 'LineWidth', 1.8, 'DisplayName', ['LV ' label_b]);
xlabel('Time (s)'); ylabel('P_{LV} (mmHg)'); title('Left Ventricle'); grid on; legend('Location', 'best');

sgtitle(sprintf('Overlay Chamber Pressures (%s): %s vs %s', scenario, label_a, label_b), ...
    'FontSize', 12, 'FontWeight', 'bold');
out1 = fullfile(fig_dir, sprintf('Overlay_Chambers_%s_%s_vs_%s.pdf', scenario, safe_label(label_a), safe_label(label_b)));
exportgraphics(fa, out1, 'ContentType', 'vector', 'Resolution', 300);

fv = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 16 12]);
plot(ta, Pa.SAR, 'k-', 'LineWidth', 1.8, 'DisplayName', ['SAR ' label_a]); hold on;
plot(tb, Pb.SAR, 'k--', 'LineWidth', 1.8, 'DisplayName', ['SAR ' label_b]);
plot(ta, Pa.PAR, 'r-', 'LineWidth', 1.8, 'DisplayName', ['PAR ' label_a]);
plot(tb, Pb.PAR, 'r--', 'LineWidth', 1.8, 'DisplayName', ['PAR ' label_b]);
plot(ta, Pa.SVEN, 'b-', 'LineWidth', 1.8, 'DisplayName', ['SVEN ' label_a]);
plot(tb, Pb.SVEN, 'b--', 'LineWidth', 1.8, 'DisplayName', ['SVEN ' label_b]);
plot(ta, Pa.PVEN, 'g-', 'LineWidth', 1.8, 'DisplayName', ['PVEN ' label_a]);
plot(tb, Pb.PVEN, 'g--', 'LineWidth', 1.8, 'DisplayName', ['PVEN ' label_b]);
xlabel('Time (s)'); ylabel('Pressure (mmHg)');
title(sprintf('Overlay Vascular Pressures (%s): %s vs %s', scenario, label_a, label_b));
grid on; legend('Location', 'best');
out2 = fullfile(fig_dir, sprintf('Overlay_Vascular_%s_%s_vs_%s.pdf', scenario, safe_label(label_a), safe_label(label_b)));
exportgraphics(fv, out2, 'ContentType', 'vector', 'Resolution', 300);

fprintf('[plot_overlay_comparison] Saved: %s\n', out1);
fprintf('[plot_overlay_comparison] Saved: %s\n', out2);

end

function opts = parse_overlay_options(varargin)
% PARSE_OVERLAY_OPTIONS â€” parse optional output directory.
root_dir = fileparts(fileparts(mfilename('fullpath')));
parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'ResultsDir', fullfile(root_dir, 'results', 'figures'), ...
    @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.ResultsDir = char(opts.ResultsDir);
end

function [tc, P] = last_cycle_pressures(sim, params)
t = sim.t(:);
X = sim.V;
idx = params.idx;
T = 60 / params.HR;
mask = t >= (t(end) - T);
tc = t(mask);

V_RA = X(mask, idx.V_RA);
V_RV = X(mask, idx.V_RV);
V_LA = X(mask, idx.V_LA);
V_LV = X(mask, idx.V_LV);
V_SAR = X(mask, idx.V_SAR);
V_SVEN = X(mask, idx.V_SVEN);
V_PAR = X(mask, idx.V_PAR);
V_PVEN = X(mask, idx.V_PVEN);

[E_LV, E_RV, E_LA, E_RA] = elastance_model(tc, params);

P.RA = max(E_RA .* (V_RA - params.V0.RA), -5);
P.RV = E_RV .* (V_RV - params.V0.RV);
P.LA = max(E_LA .* (V_LA - params.V0.LA), -5);
P.LV = E_LV .* (V_LV - params.V0.LV);
P.SAR = (V_SAR - params.V0.SAR) ./ params.C.SAR;
P.SVEN = max((V_SVEN - params.V0.SVEN) ./ params.C.SVEN, -5);
P.PAR = (V_PAR - params.V0.PAR) ./ params.C.PAR;
P.PVEN = max((V_PVEN - params.V0.PVEN) ./ params.C.PVEN, -5);
end

function s = safe_label(lbl)
s = lower(lbl);
s = regexprep(s, '\\s+', '_');
s = regexprep(s, '[^a-z0-9_]', '');
end
