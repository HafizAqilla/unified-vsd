function recipe = reyna_post_surgery()
% REYNA_POST_SURGERY
% -----------------------------------------------------------------------
% Calibration recipe for Reyna's post-closure state (paper Section 2.10,
% "Post-closure calibration").
%
% This calibration is run AFTER the zero-shot post-closure prediction and
% does not change it. It applies the pre-closure procedure
% (reyna_pre_surgery.m) to the post-closure record, so that the parameters
% of the two states can be compared:
%   - free parameters: the GSA-selected pre-closure set without vsd.Cd,
%     which has no effect once the defect is closed;
%   - targets: the seven post-closure catheter pressures and Qs, with the
%     same hard/soft split as before closure;
%   - Qp/Qs (0.90, below the value of 1 that a closed defect imposes) and
%     the chamber volumes are compared but not fitted, as before closure;
%   - same bound scales, metric weights and staged multi-start optimizer.
%
% OUTPUTS:
%   recipe - calibration recipe struct                                  [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-10-07
% VERSION:  1.0
% -----------------------------------------------------------------------

recipe = struct();
recipe.id = 'reyna_post_surgery';
recipe.version = '2026-10-07';
recipe.patient_label = 'reyna';
recipe.scenario = 'post_surgery';
recipe.profile_mode = 'reyna_recipe';
recipe.description = ['Reyna post-closure calibration; the pre-closure ', ...
    'procedure applied to the post-closure record after the prediction.'];

% Same measurement-day anthropometry as reyna_pre_surgery.m.
recipe.demographics = struct( ...
    'weight_kg', 13.4, ...
    'height_cm', 95.0, ...
    'BSA', 0.588);

% Same target split as before closure, without the shunt targets.
recipe.primary_metrics = {'RAP_mean','PAP_mean','SAP_mean','CO_Lmin'};
recipe.soft_metrics = {'SAP_max','SAP_min','PAP_max','PAP_min'};
recipe.consistency_only = {'QpQs','LVEDV','LVESV','RVEDV','RVESV','LVEF','RVEF'};
recipe.derived_validation = {'SVR','PVR'};
recipe.validation_holdout = {};
recipe.primary_rmse_holdout = {};

% GSA-selected pre-closure set (results/runs/20260920_005419_reyna_pre_surgery)
% without vsd.Cd.
recipe.active_parameters = {'group.R_sys_scale','R.SVEN', ...
    'group.R_pul_scale','C.SAR','C.PAR','E.RV.EB'};
recipe.stage_c_parameters = {'group.R_pul_scale','C.PAR', ...
    'group.R_sys_scale','R.SVEN','C.SAR','E.RV.EB'};
recipe.systemic_polish_parameters = {'group.R_sys_scale','R.SVEN', ...
    'C.SAR','E.RV.EB'};
recipe.systemic_polish_metrics = {'SAP_mean','CO_Lmin','RAP_mean','PAP_mean'};

% No historical seeds: the multi-start begins from the scaled and
% clinically seeded vector, as before closure. The scaling_modes entry
% keeps these empty blocks inert for the Zhang scaling used in the paper.
recipe.initial_parameter_values = struct('names', {{}}, 'values', [], ...
    'scaling_modes', {{'lundquist_bsa'}});
recipe.accept_initial_seed_if_pass = false;
recipe.accept_initial_seed_rmse_max = 0;
recipe.accept_initial_seed_scaling_modes = {'lundquist_bsa'};
recipe.fixed_parameter_values = struct('names', {{}}, 'values', [], ...
    'scaling_modes', {{'lundquist_bsa'}});

% Same bound scales as reyna_pre_surgery.m.
recipe.bound_scale = struct();
recipe.bound_scale.names = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB', ...
    'E.LA.EA','E.RA.EA','V0.LV','V0.RV'};
recipe.bound_scale.lower = [0.25, 0.25, 0.45, 0.75, 0.70, 0.60, ...
    0.50, 0.55, 0.50, 0.50, 0.50, 0.20, 0.35];
recipe.bound_scale.upper = [2.80, 2.80, 2.80, 1.35, 1.45, 2.20, ...
    2.80, 2.70, 2.80, 1.80, 1.80, 1.90, 1.90];

% Same weights as reyna_pre_surgery.m, without the shunt targets.
recipe.metric_weight_overrides = struct( ...
    'PAP_mean', 0.90, ...
    'PAP_max', 0.50, ...
    'PAP_min', 0.45, ...
    'SAP_mean', 1.35, ...
    'RAP_mean', 1.20, ...
    'CO_Lmin', 1.80);

recipe.acceptance = struct( ...
    'primary_rmse_max', 0.09, ...
    'primary_gate_pct', 10, ...
    'excellent_gate_pct', 5, ...
    'secondary_gate_pct', 15);

recipe.scaling_modes = {'lundquist_bsa','zhang'};
recipe.preferred_scaling_mode = 'lundquist_bsa';
end
