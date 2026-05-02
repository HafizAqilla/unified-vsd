function run_quick_regression_suite()
% RUN_QUICK_REGRESSION_SUITE
% -----------------------------------------------------------------------
% Run a lightweight no-plot, no-GSA regression sweep across benchmark
% patient profiles.
%
% This helper is intended for overnight or checkpoint verification after
% structural calibration changes.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-02
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));

setenv('UNIFIED_VSD_DO_PLOTS', '0');
setenv('UNIFIED_VSD_DO_OVERLAY', '0');
setenv('UNIFIED_VSD_DO_GSA', '0');

cases = {
    'pre_surgery', patient_reyna()
    'pre_surgery', patient_profile_Razka()
    'pre_surgery', patient_profile_A()
    };

for i = 1:size(cases, 1)
    scenario = cases{i, 1};
    clinical = cases{i, 2};
    fprintf('\n[run_quick_regression_suite] Case %d/%d: %s (%s)\n', ...
        i, size(cases, 1), clinical.common.patient_name, scenario);
    main_run(scenario, clinical);
end
end
