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

% Measurement-day anthropometry from the study "reyna" procedure log,
% 06/04/2026 07.53 — the same session that produced every pre- and
% post-closure pressure this recipe fits. Full provenance is kept out of
% this tracked file; see config/private/patient_provenance.local.m.
%
% These values are merged OVER clinical.common by
% apply_calibration_recipe_to_clinical.m:24, so they must stay in step with
% config/patient_reyna.m. The previous entry (14.0 kg, 98.0 cm, Mosteller BSA
% 0.6173) came from a 2026-05-11 revision recorded five weeks after this
% catheterisation, and silently re-scaled the model to a larger child than
% the one who was measured.
recipe.demographics = struct( ...
    'weight_kg', 13.4, ...
    'height_cm', 95.0, ...
    'BSA', 0.588);

pre = struct();
pre.VSD_diameter_mm = 3.025;       % [mm] accepted effective RV-side diameter
pre.LAP_mean_mmHg = NaN;           % [mmHg] not directly measured; exclude from RMSE
pre.LVEDP_mmHg = NaN;              % [mmHg] not directly measured; exclude from RMSE

% ---- Chamber volumes: removed as pre-surgery targets -------------------
% The LV/RV volume and EF block is H+1 POST-operative echo. It measures a
% different physiological state: with the VSD open the LV is volume-loaded
% by the left-to-right shunt, so pre-operative LVEDV is expected to exceed
% the post-closure value, not equal it. Fitting a pre-operative model to
% post-closure volumes would make the model wrong, not accurate.
%
% patient_reyna() already records these as NaN for exactly this reason. The
% recipe previously re-injected them as consistency-only evidence, which
% still let them drive the clinical consistency audit and (through
% override_IC) the initial chamber state. Both channels are now closed.
%
% The values are retained below as documented excluded evidence, and are
% recommended for relocation to clinical.post_surgery, where their timing is
% valid and where no clinical target currently exists at all.
pre.LVEDV_mL = NaN;                % [mL] post-operative evidence; not a pre-surgery target
pre.LVESV_mL = NaN;                % [mL] post-operative evidence; not a pre-surgery target
pre.RVEDV_mL = NaN;                % [mL] post-operative evidence; not a pre-surgery target
pre.RVESV_mL = NaN;                % [mL] post-operative evidence; not a pre-surgery target
pre.LVEF = NaN;                    % [-] derived from the same post-operative block
pre.EF = NaN;                      % [-] legacy alias for target mapping
pre.override_IC = false;           % [-] pre-surgery fit stays purely haemodynamic
recipe.pre_surgery_overrides = pre;

% Excluded evidence, retained for provenance and reporting only. These values
% are never mapped into clinical.pre_surgery targets; they document what was
% measured, when, and why it is not a pre-operative comparator.
recipe.excluded_evidence = struct( ...
    'timing', 'post_operative_H1', ...
    'reason', ['Chamber volumes and EF were recorded at H+1 after VSD ', ...
        'closure. Pre-operative LV loading differs (left-to-right shunt), ', ...
        'so these are not valid pre-surgery comparators.'], ...
    'recommended_scenario', 'post_surgery', ...
    'LVEDV_mL', 41.0, ...
    'LVESV_mL', 19.3, ...
    'RVEDV_mL', 30.5, ...
    'RVESV_mL', 12.0, ...
    'LVEF', 0.528);

% Provenance for recipe-supplied clinical overrides. Enforced by
% assert_evidence_timing_governance: nothing marked post_operative_H1 may
% carry a fitted tier in a pre_surgery run. Retained so that re-adding any of
% the excluded rows fails loudly rather than silently.
recipe.evidence_timing = struct( ...
    'LVEDV_mL', 'post_operative_H1', ...
    'LVESV_mL', 'post_operative_H1', ...
    'RVEDV_mL', 'post_operative_H1', ...
    'RVESV_mL', 'post_operative_H1', ...
    'LVEF', 'post_operative_H1', ...
    'EF', 'post_operative_H1');
recipe.allow_cross_timing_evidence = false;

recipe.primary_metrics = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
% PAP_max/PAP_min are directly measured catheter pressures (rows 16-17,
% three repeated measures each). They were previously graded inside the
% governed primary RMSE without ever entering the objective; fitting them
% closes that gap and constrains pulmonary compliance/resistance directly.
recipe.soft_metrics = {'Q_shunt_Lmin','SAP_max','SAP_min','PAP_max','PAP_min'};
% The chamber block is no longer supplied to pre_surgery at all (see
% pre_surgery_overrides above), so there is nothing left to demote. These
% names are retained so that re-adding any of them lands in a report-only
% tier rather than silently becoming a fitted target.
recipe.consistency_only = {'LVEDV','LVESV','RVEDV','RVESV','LVEF'};
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
recipe.initial_parameter_values.scaling_modes = {'lundquist_bsa'};
recipe.initial_parameter_values.representation = ...
    'calibration_vector_relative_to_lundquist_bsa_reference';
recipe.accept_initial_seed_if_pass = true;
recipe.accept_initial_seed_rmse_max = 0.095;
recipe.accept_initial_seed_scaling_modes = {'lundquist_bsa'};

recipe.fixed_parameter_values = struct();
recipe.fixed_parameter_values.names = {'V0.SVEN'};
recipe.fixed_parameter_values.values = 421.1753037549494;
recipe.fixed_parameter_values.source = ...
    'results/runs/20260521_214941_reyna_pre_surgery accepted_candidate';
recipe.fixed_parameter_values.scaling_modes = {'lundquist_bsa'};
recipe.fixed_parameter_values.representation = 'absolute_physical_parameter';

recipe.initial_conditions = struct();
recipe.initial_conditions.V = [11.9210; 46.3451; 15.0976; 41.0000; ...
    80.1331; 57.0500; 53.2276; 571.5634; 57.0500; 121.4592; ...
    57.0500; 11.5000; 39.2530; 57.0500];
recipe.initial_conditions.source = ...
    'results/runs/20260521_214941_reyna_pre_surgery accepted_candidate';
recipe.initial_conditions.scaling_modes = {'lundquist_bsa'};
recipe.initial_conditions.representation = 'absolute_state_vector';

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
    'PAP_max', 0.50, ...
    'PAP_min', 0.45, ...
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
