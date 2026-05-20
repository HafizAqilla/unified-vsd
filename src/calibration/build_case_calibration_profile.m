function profile = build_case_calibration_profile(clinical, scenario)
% BUILD_CASE_CALIBRATION_PROFILE
% -----------------------------------------------------------------------
% Defines evidence-aware calibration governance for a patient scenario.
%
% The profile separates rich real cases, sparse catheterisation cases, and
% synthetic benchmarks before optimisation. This prevents underidentified
% cases from using the same free-parameter set as full-data cases.
%
% INPUTS:
%   clinical  - unified clinical struct from config/ patient profiles    [-]
%   scenario  - scenario string: 'pre_surgery' | 'post_surgery'          [-]
%
% OUTPUTS:
%   profile   - calibration governance profile struct                    [-]
%
% ASSUMPTIONS:
%   - Direct catheter pressures and measured Qp/Qs are more informative
%     than derived resistance values when cardiac output is unavailable.
%   - Synthetic profiles are benchmarks, not patient-specific evidence.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%   [2] AGENTS.md, Sections 9-10: clinical reliability and validation.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-05
% VERSION:  1.1
% -----------------------------------------------------------------------

if nargin < 2 || isempty(scenario)
    scenario = 'pre_surgery';
end
if nargin < 1 || isempty(clinical)
    clinical = patient_template();
end

case_id = resolve_case_id(clinical);          % [char]
src = clinical.(scenario);                    % [-]
has_CO = has_finite_field(src, 'CO_Lmin');    % [-]
has_QpQs = has_finite_field(src, 'QpQs');     % [-]
has_pressures = has_finite_field(src, 'PAP_mean_mmHg') && ...
    has_finite_field(src, 'SAP_mean_mmHg');   % [-]
has_LV_volume = has_finite_field(src, 'LVEDV_mL') || ...
    has_finite_field(src, 'LVESV_mL');        % [-]
has_RV_volume = has_finite_field(src, 'RVEDV_mL') || ...
    has_finite_field(src, 'RVESV_mL');        % [-]
is_synthetic = contains(lower(case_id), 'patient_profile') || ...
    contains(lower(case_id), 'synthetic');    % [-]

profile = base_profile(case_id, scenario);
if is_synthetic
    profile.mode = 'synthetic_benchmark';
    profile.description = 'Synthetic benchmark; useful for stress testing, not patient-specific inference.';
elseif has_pressures && has_QpQs && ...
        ~(has_CO && has_QpQs && has_pressures && has_LV_volume && has_RV_volume)
    profile.mode = 'sparse_cath';
    profile.description = 'Sparse catheterisation case; use representative vascular and shunt knobs only.';
else
    profile.mode = 'adaptive_patient';
    profile.description = 'Real patient adaptive calibration; targets and free parameters follow available evidence.';
end

profile.evidence = struct( ...
    'has_CO_Lmin', has_CO, ...
    'has_QpQs', has_QpQs, ...
    'has_pressure_targets', has_pressures, ...
    'has_LV_volume_target', has_LV_volume, ...
    'has_RV_volume_target', has_RV_volume, ...
    'has_full_data_block', has_CO && has_QpQs && has_pressures && has_LV_volume && has_RV_volume);

switch profile.mode
    case 'adaptive_patient'
        profile = apply_adaptive_patient_governance(profile, src, scenario);
    case 'sparse_cath'
        profile.maxPrimaryMetrics = 5;        % [-] five direct hemodynamic anchors when CO is available
        profile.regLambda = 0.05;             % [-] soft prior against unsupported drift
        profile.allowedMetricFields = {'RAP_mean','PAP_min','PAP_max','PAP_mean', ...
            'SAP_min','SAP_max','SAP_mean','QpQs','Q_shunt_Lmin','CO_Lmin'};
        profile.preferredPrimaryMetrics = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
        % Batch B identifiability audit showed near-collinearity inside
        % serial systemic and pulmonary beds. Sparse cases therefore use
        % vascular/shunt knobs by default; chamber mechanics are enabled
        % only when matching pre-surgery chamber evidence is available.
        profile.allowedFreeParameters = sparse_cath_free_parameters(src, scenario);
        profile.stageCPreferredNames = sparse_cath_stage_c_names( ...
            profile.allowedFreeParameters);
        profile.boundScale = make_sparse_bound_scale();
        profile.metricWeightOverrides = struct('QpQs', 1.3, 'Q_shunt_Lmin', 1.5, ...
            'PAP_mean', 1.2, 'SAP_mean', 1.2, 'CO_Lmin', 1.3);
    otherwise
        profile.stageCPreferredNames = {'R.vsd','vsd.Cd','group.R_pul_scale', ...
            'group.R_sys_scale','C.SAR','C.PAR','E.LV.EA'};
        profile.boundScale = make_full_data_bound_scale();
        profile.systemicLoadWeights = make_full_data_systemic_load_weights();
end

profile = apply_target_tier_governance(profile, clinical, scenario);
profile = apply_age_validity_prior_adjustment(profile, clinical);
end

function profile = base_profile(case_id, scenario)
profile = struct();
profile.case_id = case_id;                    % [char]
profile.scenario = scenario;                  % [char]
profile.mode = 'full_data';                   % [char]
profile.description = '';                     % [char]
profile.allowedMetricFields = {};             % [cellstr]
profile.allowedFreeParameters = {};           % [cellstr]
profile.stageCPreferredNames = {};            % [cellstr]
profile.coupledParameterGroups = make_common_coupled_groups();
profile.metricWeightOverrides = struct();     % [-]
profile.boundScale = struct('names', {{}}, 'lower', [], 'upper', []);
profile.maxPrimaryMetrics = 5;                % [-]
profile.preferredPrimaryMetrics = {};         % [cellstr]
profile.validationHoldoutMetrics = {};        % [cellstr]
profile.regLambda = 0;                        % [-]
profile.useVascularRcCoupling = true;         % [-]
profile.preferredScalingMode = 'lundquist_bsa';
profile.ageValidityRegime = 'unspecified';    % [char]
profile.preschoolPriorAdjustmentApplied = false; % [-]
profile.paramPlausibilityLambdaScale = 1.0;    % [-]
profile.boundaryPlausibilityLambdaScale = 1.0; % [-]
profile.volumeFlowGuard = struct();           % [-]
profile.clinicalConsistencyAudit = struct();   % [-]
profile.targetTiers = struct();                % [-]
profile.useRVEDVForClinicalSeeding = true;     % [-]
end

function bound_scale = make_sparse_bound_scale()
bound_scale = struct();
bound_scale.names = {'group.R_sys_scale','R.SVEN','C.SAR','group.R_pul_scale', ...
    'C.PAR','E.LV.EA','E.LV.EB','E.RV.EA','E.RV.EB','V0.LV','V0.RV', ...
    'R.vsd','vsd.Cd'};
% Sparse catheter cases need enough local freedom to preserve a validated
% pressure-flow basin while adding derived shunt-flow consistency.
bound_scale.lower = [0.50, 0.40, 0.45, 0.50, 0.30, 0.42, 0.50, 0.50, ...
    0.50, 0.50, 0.35, 0.25, 0.60];
bound_scale.upper = [2.50, 2.50, 1.75, 2.50, 1.75, 2.80, 2.80, 2.80, ...
    2.80, 1.90, 1.90, 4.00, 2.00];
end

function names = sparse_cath_free_parameters(src, scenario)
% SPARSE_CATH_FREE_PARAMETERS - evidence-gated sparse fit parameters.
names = {'group.R_sys_scale','R.SVEN','C.SAR','group.R_pul_scale','C.PAR'};

if strcmp(scenario, 'pre_surgery') && has_finite_field(src, 'QpQs')
    if isfield(src, 'VSD_mode') && strcmpi(char(src.VSD_mode), 'orifice_bidirectional')
        names = [names, {'vsd.Cd'}];
    else
        names = [names, {'R.vsd'}];
    end
end

if has_finite_field(src, 'LVEDV_mL') || has_finite_field(src, 'LVESV_mL') || ...
        has_finite_field(src, 'EF')
    names = [names, {'E.LV.EA','E.LV.EB','V0.LV'}];
end
if has_finite_field(src, 'RVEDV_mL') || has_finite_field(src, 'RVESV_mL') || ...
        has_finite_field(src, 'RVEF')
    names = [names, {'E.RV.EA','E.RV.EB','V0.RV'}];
end

names = unique(names, 'stable');
end

function names = sparse_cath_stage_c_names(allowed_names)
% SPARSE_CATH_STAGE_C_NAMES - ordered stage-C subset from active evidence.
preferred = {'R.vsd','vsd.Cd','group.R_pul_scale','C.PAR', ...
    'group.R_sys_scale','R.SVEN','C.SAR','E.LV.EA','E.LV.EB', ...
    'E.RV.EA','E.RV.EB','V0.LV','V0.RV'};
names = preferred(ismember(preferred, allowed_names));
end

function bound_scale = make_full_data_bound_scale()
bound_scale = struct();
bound_scale.names = {'group.R_sys_scale','R.SVEN','group.R_pul_scale', ...
    'C.SAR','C.PAR','E.LV.EA','E.RV.EA','V0.LV','V0.RV','R.vsd','vsd.Cd'};
bound_scale.lower = [0.25, 0.25, 0.45, 0.75, 0.70, 0.60, 0.55, 0.75, 0.70, 0.25, 0.80];
bound_scale.upper = [2.80, 2.80, 2.80, 1.35, 1.45, 2.20, 2.40, 1.40, 1.35, 4.00, 1.20];
end

function weights = make_full_data_systemic_load_weights()
weights = struct();
weights.qs = 0.95;
weights.sap = 0.85;
weights.rap = 0.45;
weights.svr = 0.65;
weights.high_svr = 0.40;
weights.low_qs_high_svr = 0.25;
weights.co_floor = 0.90;
weights.svr_ceiling = 0.70;
weights.low_flow_svr_wall = 0.50;
weights.flow_balance = 0.20;
end

function profile = apply_adaptive_patient_governance(profile, src, scenario)
% APPLY_ADAPTIVE_PATIENT_GOVERNANCE - choose fit surface from evidence.
profile.maxPrimaryMetrics = 5;
profile.regLambda = 0.02;
profile.allowedMetricFields = adaptive_metric_fields(src, scenario);
profile.allowedFreeParameters = adaptive_free_parameters(src, scenario);
profile.stageCPreferredNames = adaptive_stage_c_names(src, scenario);
profile.boundScale = make_full_data_bound_scale();
profile.systemicLoadWeights = make_adaptive_systemic_load_weights(src);
profile.systemicPolishEnabled = strcmp(scenario, 'pre_surgery') && has_systemic_evidence(src, scenario);
profile.systemicPolishNames = adaptive_systemic_polish_names(src, scenario);
profile.systemicPolishMetrics = adaptive_systemic_polish_metrics(src, scenario);
profile.systemicPolishWeights = make_systemic_polish_weights(src);
profile.plausibilityPolishEnabled = true;
profile.plausibilityPolishParamLambda = 2.25;
profile.plausibilityPolishBoundaryLambda = 120.0;
profile.plausibilityPolishRmseTolerance = 0.005;
profile.plausibilityPolishSystemicPairTolerance = 0.05;
profile.metricWeightOverrides = make_adaptive_metric_weights(src);
profile.validationHoldoutMetrics = adaptive_holdout_metrics(src, scenario);
profile.targetGovernance = 'direct_measurements_fit_derived_values_check';
profile = apply_reyna_systemic_flow_profile(profile, src, scenario);
end

function profile = apply_reyna_systemic_flow_profile(profile, src, scenario)
% APPLY_REYNA_SYSTEMIC_FLOW_PROFILE - prioritize Reyna systemic flow/volume fit.
if ~strcmpi(profile.case_id, 'reyna') || ~strcmpi(char(scenario), 'pre_surgery')
    return;
end

has_systemic_flow_targets = has_finite_field(src, 'CO_Lmin') && ...
    has_finite_field(src, 'SAP_mean_mmHg') && has_finite_field(src, 'RAP_mean_mmHg');
has_volume_function_targets = has_finite_field(src, 'RVEDV_mL') && ...
    (has_finite_field(src, 'LVESV_mL') || has_finite_field(src, 'EF'));
if ~(has_systemic_flow_targets && has_volume_function_targets)
    return;
end

profile.preferredPrimaryMetrics = {'PAP_mean', 'QpQs', 'SAP_mean', 'CO_Lmin', 'RAP_mean'};
profile.validationHoldoutMetrics = setdiff(profile.validationHoldoutMetrics, {'CO_Lmin'}, 'stable');

reyna_weights = struct();
reyna_weights.PAP_mean = 0.90;
reyna_weights.QpQs = 1.10;
reyna_weights.SAP_mean = 1.35;
reyna_weights.RAP_mean = 1.20;
reyna_weights.CO_Lmin = 1.80;
reyna_weights.SVR = 1.10;
reyna_weights.LVESV = 1.35;
reyna_weights.LVEF = 1.45;
profile.metricWeightOverrides = merge_metric_weight_overrides( ...
    profile.metricWeightOverrides, reyna_weights);

profile.systemicLoadWeights = reyna_systemic_load_weights();
profile.systemicPolishWeights = reyna_systemic_polish_weights();
profile.systemicPolishMetrics = available_reyna_systemic_polish_metrics(src);
profile.systemicPolishNames = available_reyna_systemic_polish_names(src, scenario);
profile.volumeFlowGuard = reyna_volume_flow_guard();
profile.targetGovernance = sprintf('%s; Reyna systemic-flow polish: protocol anthropometry, PAP_mean-centered pulmonary target, stronger CO/SAP/RAP/LVESV constraints with EF retained as a derived consistency check.', ...
    profile.targetGovernance);
end

function merged = merge_metric_weight_overrides(existing, override)
% MERGE_METRIC_WEIGHT_OVERRIDES - override named metric weights by field.
merged = existing;
names = fieldnames(override);
for idx = 1:numel(names)
    merged.(names{idx}) = override.(names{idx});
end
end

function weights = reyna_systemic_load_weights()
% REYNA_SYSTEMIC_LOAD_WEIGHTS - coupled penalty emphasizing flow recovery.
weights = struct();
weights.qs = 1.45;
weights.sap = 1.35;
weights.rap = 0.65;
weights.svr = 0.95;
weights.high_svr = 1.00;
weights.low_qs_high_svr = 2.10;
weights.co_floor = 3.00;
weights.svr_ceiling = 1.80;
weights.low_flow_svr_wall = 3.00;
weights.flow_balance = 0.40;
end

function weights = reyna_systemic_polish_weights()
% REYNA_SYSTEMIC_POLISH_WEIGHTS - final pass weights for systemic consistency.
weights = struct();
weights.qs = 1.75;
weights.sap = 1.50;
weights.rap = 0.75;
weights.svr = 1.00;
weights.high_svr = 1.10;
weights.low_qs_high_svr = 2.40;
weights.co_floor = 3.50;
weights.svr_ceiling = 2.00;
weights.low_flow_svr_wall = 3.50;
weights.flow_balance = 0.45;
end

function guard = reyna_volume_flow_guard()
% REYNA_VOLUME_FLOW_GUARD - profile-specific penalties for residual errors.
guard = struct();
guard.enabled = true;
guard.RVEDV_upper_ratio = 1.10;
guard.RVEDV_upper_scale = 120;
guard.LVESV_lower_ratio = 0.90;
guard.LVESV_lower_scale = 110;
guard.LVEF_upper_ratio = 1.10;
guard.LVEF_upper_scale = 130;
guard.CO_low_rel_thresh = 0.12;
guard.CO_low_scale = 120;
end

function metrics = available_reyna_systemic_polish_metrics(src)
% AVAILABLE_REYNA_SYSTEMIC_POLISH_METRICS - keep only measured targets.
preferred = {'SAP_mean', 'CO_Lmin', 'RAP_mean', 'QpQs', 'PAP_mean', ...
    'RVEDV', 'LVESV', 'LVEF', 'LVEDV', 'RVESV'};
available = adaptive_metric_fields(src, 'pre_surgery');
metrics = preferred(ismember(preferred, available));
end

function names = available_reyna_systemic_polish_names(src, scenario)
% AVAILABLE_REYNA_SYSTEMIC_POLISH_NAMES - allow systemic, LV, and RV knobs.
preferred = {'group.R_sys_scale', 'R.SVEN', 'C.SAR', ...
    'E.LV.EA', 'E.LV.EB', 'V0.LV', ...
    'E.RV.EA', 'E.RV.EB', 'V0.RV'};
allowed = adaptive_free_parameters(src, scenario);
names = preferred(ismember(preferred, allowed));
end

function metrics = adaptive_metric_fields(src, scenario)
metrics = {};

metrics = add_if_available(metrics, src, 'RAP_mean', 'RAP_mean_mmHg');
metrics = add_if_available(metrics, src, 'LAP_mean', 'LAP_mean_mmHg');
metrics = add_if_available(metrics, src, 'PAP_min', 'PAP_dia_mmHg');
metrics = add_if_available(metrics, src, 'PAP_max', 'PAP_sys_mmHg');
metrics = add_if_available(metrics, src, 'PAP_mean', 'PAP_mean_mmHg');
metrics = add_if_available(metrics, src, 'SAP_min', 'SAP_dia_mmHg');
metrics = add_if_available(metrics, src, 'SAP_max', 'SAP_sys_mmHg');
if strcmp(scenario, 'post_surgery')
    metrics = add_if_available(metrics, src, 'SAP_mean', 'MAP_mmHg');
else
    metrics = add_if_available(metrics, src, 'SAP_mean', 'SAP_mean_mmHg');
end
metrics = add_if_available(metrics, src, 'QpQs', 'QpQs');
metrics = add_if_available(metrics, src, 'CO_Lmin', 'CO_Lmin');
metrics = add_if_available(metrics, src, 'LVEDV', 'LVEDV_mL');
metrics = add_if_available(metrics, src, 'LVESV', 'LVESV_mL');
metrics = add_if_available(metrics, src, 'RVEDV', 'RVEDV_mL');
metrics = add_if_available(metrics, src, 'RVESV', 'RVESV_mL');
metrics = add_if_available(metrics, src, 'LVEF', 'EF');
metrics = add_if_available(metrics, src, 'RVEF', 'RVEF');

% Derived resistances are held out when their source pressure/flow bundle is
% present, so calibration does not count the same physiology twice.
if has_finite_field(src, 'PVR_WU') && ~has_pulmonary_bundle(src)
    metrics{end + 1} = 'PVR';
end
if has_finite_field(src, 'SVR_WU') && ~has_systemic_bundle(src, scenario)
    metrics{end + 1} = 'SVR';
end
end

function names = adaptive_free_parameters(src, scenario)
names = {};
if has_systemic_evidence(src, scenario)
    names = [names, {'group.R_sys_scale','R.SVEN'}];
end
if has_pulmonary_evidence(src)
    names = [names, {'group.R_pul_scale'}];
end
if has_systemic_pulse_evidence(src, scenario)
    names = [names, {'C.SAR'}];
end
if has_pulmonary_pulse_evidence(src)
    names = [names, {'C.PAR'}];
end
if strcmp(scenario, 'pre_surgery') && has_finite_field(src, 'QpQs')
    if isfield(src, 'VSD_mode') && strcmpi(char(src.VSD_mode), 'orifice_bidirectional')
        names = [names, {'vsd.Cd'}];
    else
        names = [names, {'R.vsd'}];
    end
end
if has_finite_field(src, 'LVEDV_mL') || has_finite_field(src, 'LVESV_mL') || ...
        has_finite_field(src, 'EF')
    names = [names, {'E.LV.EA','E.LV.EB','V0.LV'}];
end
if has_finite_field(src, 'RVEDV_mL') || has_finite_field(src, 'RVESV_mL') || ...
        has_finite_field(src, 'RVEF')
    names = [names, {'E.RV.EA','E.RV.EB','V0.RV'}];
end
if has_finite_field(src, 'LAP_mean_mmHg')
    names = [names, {'E.LA.EA'}];
end
if has_finite_field(src, 'RAP_mean_mmHg')
    names = [names, {'E.RA.EA'}];
end
names = unique(names, 'stable');
end

function names = adaptive_stage_c_names(src, scenario)
names = {};
if strcmp(scenario, 'pre_surgery') && has_finite_field(src, 'QpQs')
    if isfield(src, 'VSD_mode') && strcmpi(char(src.VSD_mode), 'orifice_bidirectional')
        names = [names, {'vsd.Cd'}];
    else
        names = [names, {'R.vsd'}];
    end
end
if has_pulmonary_evidence(src)
    names = [names, {'group.R_pul_scale'}];
end
if has_systemic_evidence(src, scenario)
    names = [names, {'group.R_sys_scale','R.SVEN'}];
end
if has_systemic_pulse_evidence(src, scenario)
    names = [names, {'C.SAR'}];
end
if has_pulmonary_pulse_evidence(src)
    names = [names, {'C.PAR'}];
end
if has_finite_field(src, 'LVEDV_mL') || has_finite_field(src, 'EF')
    names = [names, {'E.LV.EA','V0.LV'}];
end
names = unique(names, 'stable');
end

function names = adaptive_systemic_polish_names(src, scenario)
% ADAPTIVE_SYSTEMIC_POLISH_NAMES - final pass knobs for pressure-flow fit.
if ~strcmp(scenario, 'pre_surgery') || ~has_systemic_evidence(src, scenario)
    names = {};
    return;
end

preferred = {'group.R_sys_scale','R.SVEN','C.SAR','E.LV.EA','E.LV.EB','V0.LV'};
allowed = adaptive_free_parameters(src, scenario);
names = preferred(ismember(preferred, allowed));
end

function metrics = adaptive_systemic_polish_metrics(src, scenario)
% ADAPTIVE_SYSTEMIC_POLISH_METRICS - coupled systemic targets plus guards.
metrics = {};

if strcmp(scenario, 'post_surgery')
    metrics = add_if_available(metrics, src, 'SAP_mean', 'MAP_mmHg');
else
    metrics = add_if_available(metrics, src, 'SAP_mean', 'SAP_mean_mmHg');
end
metrics = add_if_available(metrics, src, 'CO_Lmin', 'CO_Lmin');
if has_finite_field(src, 'SVR_WU')
    metrics{end + 1} = 'SVR';
end
metrics = add_if_available(metrics, src, 'SAP_max', 'SAP_sys_mmHg');
metrics = add_if_available(metrics, src, 'QpQs', 'QpQs');
metrics = add_if_available(metrics, src, 'PAP_mean', 'PAP_mean_mmHg');
metrics = add_if_available(metrics, src, 'LVEDV', 'LVEDV_mL');
metrics = add_if_available(metrics, src, 'LVESV', 'LVESV_mL');
metrics = add_if_available(metrics, src, 'RVEDV', 'RVEDV_mL');
metrics = add_if_available(metrics, src, 'RVESV', 'RVESV_mL');
metrics = add_if_available(metrics, src, 'LVEF', 'EF');
metrics = unique(metrics, 'stable');
end

function weights = make_adaptive_systemic_load_weights(src)
weights = make_full_data_systemic_load_weights();
if has_finite_field(src, 'CO_uncertainty_Lmin')
    weights.qs = 0.75;
    weights.svr = 0.50;
    weights.high_svr = 0.30;
    weights.low_qs_high_svr = 0.35;
    weights.co_floor = 1.20;
    weights.svr_ceiling = 0.90;
    weights.low_flow_svr_wall = 1.00;
end
end

function weights = make_systemic_polish_weights(src)
weights = make_full_data_systemic_load_weights();
weights.qs = 1.35;
weights.sap = 1.05;
weights.rap = 0.35;
weights.svr = 1.10;
weights.high_svr = 0.90;
weights.low_qs_high_svr = 1.30;
weights.co_floor = 2.50;
weights.svr_ceiling = 1.60;
weights.low_flow_svr_wall = 2.50;
weights.flow_balance = 0.35;

% CO is still guarded when documented as uncertain, but the polish pass
% needs enough leverage to test the systemic pressure-flow hypothesis.
if has_finite_field(src, 'CO_uncertainty_Lmin')
    weights.qs = 1.10;
    weights.svr = 0.90;
    weights.co_floor = 2.20;
    weights.svr_ceiling = 1.40;
    weights.low_flow_svr_wall = 2.20;
end
end

function weights = make_adaptive_metric_weights(src)
weights = struct();
if has_finite_field(src, 'QpQs'), weights.QpQs = 1.25; end
if has_finite_field(src, 'PAP_mean_mmHg'), weights.PAP_mean = 1.20; end
if has_finite_field(src, 'SAP_mean_mmHg') || has_finite_field(src, 'MAP_mmHg')
    weights.SAP_mean = 1.15;
end
if has_finite_field(src, 'CO_uncertainty_Lmin')
    weights.CO_Lmin = 0.85;
end
if has_finite_field(src, 'LVEDV_mL') || has_finite_field(src, 'LVESV_mL')
    weights.LVEDV = 0.85;
    weights.LVESV = 0.75;
end
if has_finite_field(src, 'RVEDV_mL') || has_finite_field(src, 'RVESV_mL')
    weights.RVEDV = 0.85;
    weights.RVESV = 0.75;
end
end

function holdout = adaptive_holdout_metrics(src, scenario)
holdout = {};
if has_finite_field(src, 'SVR_WU') && has_systemic_bundle(src, scenario)
    holdout{end + 1} = 'SVR';
end
if has_finite_field(src, 'PVR_WU') && has_pulmonary_bundle(src)
    holdout{end + 1} = 'PVR';
end
if has_finite_field(src, 'CO_uncertainty_Lmin')
    holdout{end + 1} = 'CO_Lmin';
end
end

function profile = apply_target_tier_governance(profile, clinical, scenario)
% APPLY_TARGET_TIER_GOVERNANCE - make target inclusion auditable.
audit = audit_clinical_consistency(clinical, scenario);
target_tiers = build_target_tiers(clinical, scenario, audit);

profile.clinicalConsistencyAudit = audit;
profile.targetTiers = target_tiers;
profile.allowedMetricFields = merge_allowed_metrics_with_tiers( ...
    profile.allowedMetricFields, target_tiers);

if isfield(profile, 'systemicPolishMetrics') && ~isempty(profile.systemicPolishMetrics)
    profile.systemicPolishMetrics = intersect( ...
        profile.systemicPolishMetrics(:)', ...
        target_tiers.included_in_calibration(:)', 'stable');
end
if isfield(profile, 'preferredPrimaryMetrics') && ~isempty(profile.preferredPrimaryMetrics)
    profile.preferredPrimaryMetrics = setdiff( ...
        profile.preferredPrimaryMetrics(:)', ...
        target_tiers.consistency_only(:)', 'stable');
end

profile.validationHoldoutMetrics = unique([ ...
    profile.validationHoldoutMetrics(:)', ...
    target_tiers.consistency_only(:)'], 'stable');

if isfield(profile.metricWeightOverrides, 'RVEDV') && ...
        ismember('RVEDV', target_tiers.consistency_only)
    profile.metricWeightOverrides = rmfield(profile.metricWeightOverrides, 'RVEDV');
end

if ismember('RVEDV', target_tiers.consistency_only)
    profile.useRVEDVForClinicalSeeding = false;
    if isfield(profile, 'volumeFlowGuard') && isstruct(profile.volumeFlowGuard)
        profile.volumeFlowGuard.RVEDV_upper_ratio = Inf;
        profile.volumeFlowGuard.RVEDV_upper_scale = 0;
    end
end

profile = append_governance_note(profile, sprintf( ...
    'target tiers=%s; consistency-only=%s; audit=%s', ...
    target_tiers.policy, strjoin_or_none(target_tiers.consistency_only), ...
    audit.severity));
end

function metrics = merge_allowed_metrics_with_tiers(existing_metrics, target_tiers)
included = target_tiers.included_in_calibration(:)';
if isempty(existing_metrics)
    metrics = included;
else
    metrics = unique([intersect(existing_metrics(:)', included, 'stable'), ...
        included], 'stable');
end
metrics = setdiff(metrics, target_tiers.consistency_only(:)', 'stable');
end

function profile = append_governance_note(profile, note)
if ~isfield(profile, 'targetGovernance') || isempty(profile.targetGovernance)
    profile.targetGovernance = note;
else
    profile.targetGovernance = sprintf('%s; %s', profile.targetGovernance, note);
end
end

function text = strjoin_or_none(values)
if isempty(values)
    text = 'none';
else
    text = strjoin(values(:)', ',');
end
end

function metrics = add_if_available(metrics, src, metric_name, field_name)
if has_finite_field(src, field_name)
    metrics{end + 1} = metric_name;
end
end

function tf = has_systemic_bundle(src, scenario)
if strcmp(scenario, 'post_surgery')
    has_map = has_finite_field(src, 'MAP_mmHg');
else
    has_map = has_finite_field(src, 'SAP_mean_mmHg');
end
tf = has_map && has_finite_field(src, 'RAP_mean_mmHg') && ...
    has_finite_field(src, 'CO_Lmin');
end

function tf = has_pulmonary_bundle(src)
tf = has_finite_field(src, 'PAP_mean_mmHg') && ...
    has_finite_field(src, 'LAP_mean_mmHg') && ...
    has_finite_field(src, 'CO_Lmin') && has_finite_field(src, 'QpQs');
end

function tf = has_systemic_evidence(src, scenario)
tf = has_systemic_bundle(src, scenario) || has_finite_field(src, 'SVR_WU') || ...
    has_finite_field(src, 'SAP_mean_mmHg') || has_finite_field(src, 'MAP_mmHg');
end

function tf = has_pulmonary_evidence(src)
tf = has_pulmonary_bundle(src) || has_finite_field(src, 'PVR_WU') || ...
    has_finite_field(src, 'PAP_mean_mmHg') || has_finite_field(src, 'QpQs');
end

function tf = has_systemic_pulse_evidence(src, scenario)
if strcmp(scenario, 'post_surgery')
    has_mean = has_finite_field(src, 'MAP_mmHg');
else
    has_mean = has_finite_field(src, 'SAP_mean_mmHg');
end
tf = has_mean && has_finite_field(src, 'SAP_sys_mmHg') && ...
    has_finite_field(src, 'SAP_dia_mmHg');
end

function tf = has_pulmonary_pulse_evidence(src)
tf = has_finite_field(src, 'PAP_mean_mmHg') && ...
    has_finite_field(src, 'PAP_sys_mmHg') && ...
    has_finite_field(src, 'PAP_dia_mmHg');
end

function profile = apply_age_validity_prior_adjustment(profile, clinical)
% APPLY_AGE_VALIDITY_PRIOR_ADJUSTMENT - record preschool extrapolation.
if ~isfield(clinical, 'common') || ~isfield(clinical.common, 'age_years') || ...
        ~isfinite(clinical.common.age_years)
    return;
end

age_years = clinical.common.age_years;        % [years]
if age_years < 2 || age_years > 6
    profile.ageValidityRegime = 'default_age_range';
    return;
end

profile.ageValidityRegime = 'preschool_extrapolated';
profile.preschoolPriorAdjustmentApplied = true;
profile.paramPlausibilityLambdaScale = 0.75;
profile.boundaryPlausibilityLambdaScale = 0.85;

if strcmp(profile.mode, 'full_data')
    profile.boundScale = widen_selected_bound_scales(profile.boundScale, ...
        {'group.R_sys_scale','R.SVEN','C.SAR','C.PAR'}, 0.85, 1.15);
end
end

function bound_scale = widen_selected_bound_scales(bound_scale, names, lower_factor, upper_factor)
% WIDEN_SELECTED_BOUND_SCALES - small documented relaxation for extrapolated ages.
for idx = 1:numel(names)
    hit = strcmp(bound_scale.names, names{idx});
    bound_scale.lower(hit) = bound_scale.lower(hit) * lower_factor;
    bound_scale.upper(hit) = bound_scale.upper(hit) * upper_factor;
end
end

function group_defs = make_common_coupled_groups()
group_defs = repmat(struct( ...
    'name', '', ...
    'memberNames', {{}}, ...
    'groupType', '', ...
    'found', false), 2, 1);

group_defs(1).name = 'group.R_sys_scale';
group_defs(1).memberNames = {'R.SAR','R.SC'};
group_defs(1).groupType = 'R';

group_defs(2).name = 'group.R_pul_scale';
group_defs(2).memberNames = {'R.PAR','R.PCOX','R.PVEN'};
group_defs(2).groupType = 'R';
end

function case_id = resolve_case_id(clinical)
case_id = 'patient_unknown';
if isfield(clinical, 'common')
    if isfield(clinical.common, 'patient_name') && ~isempty(clinical.common.patient_name)
        case_id = char(clinical.common.patient_name);
    elseif isfield(clinical.common, 'patient_id') && ~isempty(clinical.common.patient_id)
        case_id = char(clinical.common.patient_id);
    end
end
end

function tf = has_finite_field(src, field_name)
tf = isfield(src, field_name) && isnumeric(src.(field_name)) && ...
    isscalar(src.(field_name)) && isfinite(src.(field_name));
end
