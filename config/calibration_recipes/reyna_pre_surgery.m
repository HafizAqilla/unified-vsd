function recipe = reyna_pre_surgery()
% REYNA_PRE_SURGERY
% -----------------------------------------------------------------------
% Explicit calibration recipe for Reyna pre-surgery pressure-flow fitting.
%
% The recipe is the source of truth for patient-scenario calibration
% decisions that should not be inferred implicitly from clinical-field
% availability. It preserves direct hemodynamic targets while retaining
% inconsistent chamber evidence only for consistency checks and seeding.
%
% OUTPUTS:
%   recipe - calibration recipe struct                                  [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-24
% VERSION:  1.0
% -----------------------------------------------------------------------

recipe = struct();
recipe.id = 'reyna_pre_surgery';
recipe.version = '2026-05-24';
recipe.patient_label = 'reyna';
recipe.scenario = 'pre_surgery';
recipe.profile_mode = 'reyna_recipe';
recipe.description = ['Reyna pre-surgery explicit pressure-flow recipe; ', ...
    'direct hemodynamics drive calibration while inconsistent chamber rows ', ...
    'remain transparent consistency evidence.'];

recipe.demographics = struct( ...
    'weight_kg', 14.0, ...
    'height_cm', 98.0, ...
    'BSA', 0.6173419726);

pre = struct();
pre.VSD_diameter_mm = 3.025;       % [mm] accepted effective RV-side diameter
pre.LAP_mean_mmHg = 8;             % [mmHg] validation-only filling estimate
pre.LVEDP_mmHg = 8;                % [mmHg] validation-only filling estimate
pre.LVEDV_mL = 41.0;               % [mL] consistency-only Teichholz evidence
pre.LVESV_mL = 19.3;               % [mL] consistency-only Teichholz evidence
pre.RVEDV_mL = 30.5;               % [mL] consistency-only/audit evidence
pre.RVESV_mL = 12.0;               % [mL] soft guard evidence
pre.LVEF = 0.528;                  % [-] consistency-only derived EF
pre.EF = 0.528;                    % [-] legacy alias for target mapping
pre.override_IC = true;            % [-] seed from audited chamber evidence
recipe.pre_surgery_overrides = pre;

recipe.primary_metrics = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
recipe.soft_metrics = {'Q_shunt_Lmin','SAP_max','SAP_min','RVESV'};
recipe.consistency_only = {'LVEDV','LVESV','RVEDV','LVEF'};
recipe.derived_validation = {'SVR','PVR'};
recipe.validation_holdout = {};
recipe.primary_rmse_holdout = {'Q_shunt_Lmin'};

recipe.active_parameters = {'group.R_sys_scale','R.SVEN', ...
    'group.R_pul_scale','C.SAR','C.PAR','E.LV.EA','E.LV.EB', ...
    'E.RV.EA','E.RV.EB','E.LA.EA','E.RA.EA','V0.LV','V0.RV','vsd.Cd'};
recipe.stage_c_parameters = {'vsd.Cd','group.R_pul_scale','C.PAR', ...
    'group.R_sys_scale','R.SVEN','C.SAR','E.LV.EA','E.LV.EB', ...
    'E.RV.EA','E.RV.EB','E.RA.EA','V0.LV','V0.RV'};
recipe.systemic_polish_parameters = {'group.R_sys_scale','R.SVEN', ...
    'C.SAR','E.LV.EA','E.LV.EB','V0.LV','E.RV.EA','E.RV.EB','V0.RV'};
recipe.systemic_polish_metrics = {'SAP_mean','CO_Lmin','RAP_mean', ...
    'QpQs','PAP_mean'};

recipe.initial_parameter_values = struct();
recipe.initial_parameter_values.names = recipe.active_parameters;
recipe.initial_parameter_values.values = [0.457982333986, 0.172102628970, ...
    0.450000000000, 0.505791080716, 1.695014086620, 6.986140864970, ...
    0.160608281828, 2.079243509970, 0.108365846196, 0.583711160536, ...
    1.228992521310, 1.427896255390, 3.036185781620, 0.520742612671];
recipe.initial_parameter_values.source = ...
    'results/runs/20260521_214941_reyna_pre_surgery accepted_candidate';
recipe.accept_initial_seed_if_pass = true;
recipe.accept_initial_seed_rmse_max = 0.095;

recipe.fixed_parameter_values = struct();
recipe.fixed_parameter_values.names = {'V0.SVEN'};
recipe.fixed_parameter_values.values = 421.1753037549494;
recipe.fixed_parameter_values.source = ...
    'results/runs/20260521_214941_reyna_pre_surgery accepted_candidate';

recipe.initial_conditions = struct();
recipe.initial_conditions.V = [11.9210; 46.3451; 15.0976; 41.0000; ...
    80.1331; 57.0500; 53.2276; 571.5634; 57.0500; 121.4592; ...
    57.0500; 11.5000; 39.2530; 57.0500];
recipe.initial_conditions.source = ...
    'results/runs/20260521_214941_reyna_pre_surgery accepted_candidate';

recipe.bound_scale = struct();
recipe.bound_scale.names = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB', ...
    'E.LA.EA','E.RA.EA','V0.LV','V0.RV','R.vsd','vsd.Cd'};
recipe.bound_scale.lower = [0.25, 0.25, 0.45, 0.75, 0.70, 0.60, ...
    0.50, 0.55, 0.50, 0.50, 0.50, 0.20, 0.35, 0.25, 0.60];
recipe.bound_scale.upper = [2.80, 2.80, 2.80, 1.35, 1.45, 2.20, ...
    2.80, 2.70, 2.80, 1.80, 1.80, 1.90, 1.90, 4.00, 2.00];

recipe.metric_weight_overrides = struct( ...
    'PAP_mean', 0.90, ...
    'QpQs', 1.10, ...
    'SAP_mean', 1.35, ...
    'RAP_mean', 1.20, ...
    'CO_Lmin', 1.80, ...
    'Q_shunt_Lmin', 0.80);

recipe.acceptance = struct( ...
    'primary_rmse_max', 0.09, ...
    'primary_gate_pct', 10, ...
    'excellent_gate_pct', 5, ...
    'secondary_gate_pct', 15);

recipe.scaling_modes = {'lundquist_bsa','zhang'};
recipe.preferred_scaling_mode = 'lundquist_bsa';
end
