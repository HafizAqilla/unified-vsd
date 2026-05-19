function test_post_surgery_warm_start()
% TEST_POST_SURGERY_WARM_START
% -----------------------------------------------------------------------
% Verifies that post-op warm-start copies physiology but never copies R.vsd.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-18
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));

params_start = default_parameters();
params_pre = params_start;

params_pre.E.LV.EA = params_start.E.LV.EA * 1.25;  % [mmHg/mL]
params_pre.E.RV.EB = params_start.E.RV.EB * 0.90;  % [mmHg/mL]
params_pre.R.SAR   = params_start.R.SAR * 1.40;    % [mmHg*s/mL]
params_pre.C.PAR   = params_start.C.PAR * 0.75;    % [mL/mmHg]
params_pre.R.vsd   = 0.01;                         % [mmHg*s/mL]
params_start.R.vsd = 1e6;                          % [mmHg*s/mL]

clinical = patient_template();
clinical.pre_surgery.CalibParams = params_pre;

[params_post, warm_start] = apply_post_surgery_warm_start( ...
    params_start, clinical, 'post_surgery');

assert(warm_start.applied, 'Expected warm-start to be applied.');
assert(abs(params_post.E.LV.EA - params_pre.E.LV.EA) < 1e-12, ...
    'LV active elastance was not copied.');
assert(abs(params_post.E.RV.EB - params_pre.E.RV.EB) < 1e-12, ...
    'RV passive elastance was not copied.');
assert(abs(params_post.R.SAR - params_pre.R.SAR) < 1e-12, ...
    'Systemic arterial resistance was not copied.');
assert(abs(params_post.C.PAR - params_pre.C.PAR) < 1e-12, ...
    'Pulmonary arterial compliance was not copied.');
assert(params_post.R.vsd == params_start.R.vsd, ...
    'R.vsd must remain controlled by the post-surgery scenario, not pre-op calibration.');

[params_pre_scenario, warm_start_pre] = apply_post_surgery_warm_start( ...
    params_start, clinical, 'pre_surgery');
assert(~warm_start_pre.applied, 'Warm-start should only apply to post_surgery.');
assert(params_pre_scenario.R.SAR == params_start.R.SAR, ...
    'pre_surgery scenario should not be modified by post-op warm-start.');
end
