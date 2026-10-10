function out = export_post_closure_calibration(run_dir)
% EXPORT_POST_CLOSURE_CALIBRATION
% -----------------------------------------------------------------------
% Writes the numbers behind the paper's post-closure calibration tables
% from one run of run_post_closure_calibration:
%   post_calibration_targets.csv  every post-closure measurement against
%                                 the calibrated model (Table post_calib)
%   post_calibration_extra.csv    model-only quantities (LV/RV systolic
%                                 pressure, Qp, LA pressure)
%   post_calibration_params.csv   fitted parameters with seeds and bounds
%
% The calibrated parameters are re-simulated with the defect closed at the
% post-closure heart rate (evaluate_post_closure_prediction), so the values
% are computed exactly as the prediction in Table 8 was.
%
% INPUTS:
%   run_dir - main_run output folder of the post-closure calibration     [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-10-07
% VERSION:  1.0
% -----------------------------------------------------------------------

clinical = patient_reyna();
cand = fullfile(run_dir, 'mat', 'params_accepted_candidate_post_surgery.mat');
out_dir = fullfile(run_dir, 'paper_export');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% The post_surgery scenario closes the defect with a linear resistance of
% 1e6 mmHg s/mL (leak about 1e-4 mL/s at 70 mmHg, i.e. < 1e-5 of Qs). The
% prediction closed it exactly (orifice law, A_VSD = 0); re-simulate the
% calibrated parameters the same way so both tables use one closure.
S = load(cand);
fn = fieldnames(S);
params = [];
for i = 1:numel(fn)
    v = S.(fn{i});
    if isstruct(v) && isfield(v, 'R') && isfield(v, 'vsd')
        params = v; break;
    elseif isstruct(v) && isfield(v, 'params') && isstruct(v.params)
        params = v.params; break;
    end
end
assert(~isempty(params), 'export_post_closure_calibration:noParams', ...
    'No parameter struct in %s', cand);
params.vsd.mode = 'orifice_bidirectional';
params.vsd.area_mm2 = 0;
closed_file = fullfile(out_dir, 'params_post_calibrated_exact_closure.mat');
save(closed_file, 'params');
ev = evaluate_post_closure_prediction(closed_file, clinical, true);

T = ev.table;
targets = get_calibration_targets('post_surgery', clinical);
sigma = nan(height(T), 1);
for i = 1:height(T)
    t = targets(strcmp({targets.Metric}, T.Metric{i}));
    sigma(i) = (T.Predicted(i) - T.Measured(i)) / T.ZScore(i);
    if isfield(t, 'UncertaintyAbs') && isfinite(t.UncertaintyAbs) && t.UncertaintyAbs > 0
        sigma(i) = t.UncertaintyAbs;
    end
end
T.Sigma = sigma;
T.Properties.VariableNames{'Predicted'} = 'Model';
writetable(T, fullfile(out_dir, 'post_calibration_targets.csv'));

m = ev.metrics;
extra = table({'LVP_max';'RVP_max';'Qp_Lmin';'Qs_Lmin';'LAP_mean';'SVR'}, ...
    [m.LVP_max; m.RVP_max; m.Qp_Lmin; m.CO_Lmin; m.LAP_mean; m.SVR], ...
    'VariableNames', {'Metric', 'Model'});
writetable(extra, fullfile(out_dir, 'post_calibration_extra.csv'));
disp(extra);

P = readtable(fullfile(run_dir, 'tables', ...
    'parameter_plausibility_accepted_candidate_post_surgery.csv'));
writetable(P, fullfile(out_dir, 'post_calibration_params.csv'));
disp(P(:, {'Parameter','SeededValue','FittedValue','RatioToBaseline', ...
    'NearLowerBound','NearUpperBound'}));

out = struct('targets', T, 'extra', extra, 'params', P, 'metrics', m);
end
