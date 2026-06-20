function params = apply_scaling(params_ref, patient)
% APPLY_SCALING
% -----------------------------------------------------------------------
% Backward-compatible wrapper that applies physiological scaling,
% pediatric maturation, and initial-condition construction in sequence.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  3.0
% -----------------------------------------------------------------------

scaling_policy = resolve_scaling_policy_from_patient(patient);
scaling_mode = scaling_policy.ScalingMode;
params = apply_physiological_scaling(params_ref, patient, scaling_mode);

if isfield(patient, 'age_days') && ~isnan(patient.age_days)
    age_days = patient.age_days;
else
    age_days = 365.25 * patient.age_years;
end

maturation_mode = 'normal';
if isfield(patient, 'maturation_mode') && ~isempty(patient.maturation_mode)
    maturation_mode = patient.maturation_mode;
end
params = apply_maturation(params, age_days, maturation_mode);

% V0 separation: chamber V0 supports chamber mechanics; vascular V0 supports preload.
% Reconcile vascular V0 (especially systemic venous V0) with blood volume before IC build.
params = reconcile_vascular_v0(params, patient, []);
params.ic.V = build_initial_conditions(params, patient);
params.scaling.age_days = age_days;
params.scaling.maturation_mode = maturation_mode;
params.scaling.requested_mode = scaling_mode;
params.scaling.policy = scaling_policy;
params.scaling.role = scaling_policy.ScalingRole;
params.scaling.citation = scaling_policy.ScalingCitation;
params.scaling.implementation_variant = scaling_policy.ImplementationVariant;
params.scaling.deviation_from_citation = scaling_policy.DeviationFromCitation;

fprintf('[apply_scaling] weight=%.1f kg | age=%.0f days | mode=%s+%s\n', ...
    patient.weight_kg, age_days, scaling_mode, maturation_mode);

end

% RECONCILE_VASCULAR_V0 — reconcile vascular V0 with BV target [mL]
function params = reconcile_vascular_v0(params, patient, clinical)
if nargin < 3
    clinical = [];
end

BV_patient = patient.weight_kg * blood_volume_per_kg(patient.age_years);

[src, scenario] = resolve_clinical_source(clinical, 'pre_surgery');

P_nom_SAR = first_valid(src, {'SAP_mean_mmHg', 'MAP_mmHg'}, 65);
P_nom_SC = max(0.5 * P_nom_SAR, 15);
P_nom_SVEN = first_valid(src, {'RAP_mean_mmHg'}, 2);
P_nom_PAR = first_valid(src, {'PAP_mean_mmHg'}, 15);
P_nom_PVEN = first_valid(src, {'LAP_mean_mmHg', 'LVEDP_mmHg'}, 6);
P_nom_LV_ED = first_valid(src, {'LVEDP_mmHg', 'LAP_mean_mmHg'}, 8);
P_nom_RV_ED = first_valid(src, {'RVEDP_mmHg', 'RAP_mean_mmHg'}, 4);
P_nom_LA_ED = first_valid(src, {'LAP_mean_mmHg'}, 6);
P_nom_RA_ED = first_valid(src, {'RAP_mean_mmHg'}, 4);

V_RA = params.V0.RA + P_nom_RA_ED / max(params.E.RA.EB, 1e-6);
V_RV = params.V0.RV + P_nom_RV_ED / max(params.E.RV.EB, 1e-6);
V_LA = params.V0.LA + P_nom_LA_ED / max(params.E.LA.EB, 1e-6);
V_LV = params.V0.LV + P_nom_LV_ED / max(params.E.LV.EB, 1e-6);
V_SAR = params.V0.SAR + P_nom_SAR * params.C.SAR;
V_SC = params.V0.SC + P_nom_SC * params.C.SC;
V_PAR = params.V0.PAR + P_nom_PAR * params.C.PAR;
V_PVEN = params.V0.PVEN + P_nom_PVEN * params.C.PVEN;

V_other = V_RA + V_RV + V_LA + V_LV + V_SAR + V_SC + V_PAR + V_PVEN;
V0_sven_required = BV_patient - V_other - P_nom_SVEN * params.C.SVEN;

V0_sven_min = max(0.01 * BV_patient, 1.0);
V0_sven_max = 0.85 * BV_patient;
V0_sven_applied = min(max(V0_sven_required, V0_sven_min), V0_sven_max);
params.V0.SVEN = V0_sven_applied;

params.scaling.BV_patient = BV_patient;
params.scaling.V0_SVEN_required = V0_sven_required;
params.scaling.V0_SVEN_applied = V0_sven_applied;
params.scaling.V0_SVEN_bounds = [V0_sven_min V0_sven_max];
params.scaling.V0_SVEN_scenario = scenario;
params.scaling.V0_SVEN_reconciled = true;
end

function BV_per_kg = blood_volume_per_kg(age_years)
if age_years < 1
    BV_per_kg = 82;
else
    BV_per_kg = 70;
end
end

function [src, scenario] = resolve_clinical_source(clinical, scenario)
src = struct();
if isempty(clinical) || ~isstruct(clinical)
    return;
end
if isfield(clinical, scenario)
    src = clinical.(scenario);
elseif strcmp(scenario, 'pre_surgery') && isfield(clinical, 'pre_surgery')
    src = clinical.pre_surgery;
elseif strcmp(scenario, 'post_surgery') && isfield(clinical, 'post_surgery')
    src = clinical.post_surgery;
end
end

function value = first_valid(src, field_names, fallback)
value = fallback;
for k = 1:numel(field_names)
    fn = field_names{k};
    if isfield(src, fn) && ~isnan(src.(fn))
        value = src.(fn);
        return;
    end
end
end

% RESOLVE_SCALING_POLICY_FROM_PATIENT - choose patient/env scaling policy.
function scaling_policy = resolve_scaling_policy_from_patient(patient)
scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
if isfield(patient, 'scaling_mode') && ~isempty(patient.scaling_mode)
    scaling_mode = patient.scaling_mode;
end

run_mode = getenv('UNIFIED_VSD_RUN_MODE');
if isfield(patient, 'run_mode') && ~isempty(patient.run_mode)
    run_mode = patient.run_mode;
end

scaling_role = '';
if isfield(patient, 'scaling_role') && ~isempty(patient.scaling_role)
    scaling_role = patient.scaling_role;
end

override_rationale = '';
if isfield(patient, 'scaling_override_rationale') && ~isempty(patient.scaling_override_rationale)
    override_rationale = patient.scaling_override_rationale;
end

scaling_policy = resolve_scaling_policy(scaling_mode, run_mode, ...
    'ScalingRole', scaling_role, ...
    'OverrideRationale', override_rationale);
end
