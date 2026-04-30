function params = params_from_clinical(params, clinical, scenario)
% PARAMS_FROM_CLINICAL
% -----------------------------------------------------------------------
% Maps scenario-specific clinical data into model parameters and rebuilds
% the initial-condition vector. Vascular V0 (especially V0.SVEN) is
% reconciled for blood-volume/preload consistency without mutating chamber V0.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  2.0
% -----------------------------------------------------------------------

R_VSD_CLOSED = 1e6;
common = clinical.common;
k = params.conv.WU_to_R;
Lmin_to_mLs = params.conv.Lmin_to_mLs;

switch scenario
    case 'pre_surgery'
        src = clinical.pre_surgery;
    case 'post_surgery'
        src = clinical.post_surgery;
    otherwise
        error('params_from_clinical:unknownScenario', ...
            'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

if isfield(common, 'HR') && ~isnan(common.HR)
    params.HR = common.HR;
    params = recompute_timing(params);
end

if isfield(src, 'SVR_WU') && ~isnan(src.SVR_WU)
    SVR_target_WU = select_resistance_target(src, 'SVR');
    SVR_mech = SVR_target_WU * k;
    SVR_ref = params.R.SAR + params.R.SC + params.R.SVEN;
    ratio = SVR_mech / max(SVR_ref, 1e-6);
    params.R.SAR = params.R.SAR * ratio;
    params.R.SC = params.R.SC * ratio;
    params.R.SVEN = params.R.SVEN * ratio;
    params.clinical_override.SVR_seed_WU = SVR_target_WU;
end

if isfield(src, 'PVR_WU') && ~isnan(src.PVR_WU)
    PVR_target_WU = select_resistance_target(src, 'PVR');
    PVR_mech = PVR_target_WU * k;
    R_cap_ref = 1 / max(1 / params.R.PCOX + 1 / params.R.PCNO, 1e-9);
    PVR_ref = params.R.PAR + R_cap_ref + params.R.PVEN;
    ratio = PVR_mech / max(PVR_ref, 1e-6);
    params.R.PAR = params.R.PAR * ratio;
    params.R.PCOX = params.R.PCOX * ratio;
    params.R.PCNO = params.R.PCNO * ratio;
    params.R.PVEN = params.R.PVEN * ratio;
    params.clinical_override.PVR_seed_WU = PVR_target_WU;
end

params = configure_vsd(params, src, scenario, R_VSD_CLOSED, Lmin_to_mLs);

if isfield(src, 'override_IC') && isequal(src.override_IC, true)
    params = apply_chamber_tuning_from_clinical(params, src);
end

patient = params.scaling.patient;
params = reconcile_vascular_v0(params, patient, clinical, scenario);
params.ic.V = build_initial_conditions(params, patient, clinical, scenario);

end

function params = configure_vsd(params, src, scenario, R_VSD_CLOSED, Lmin_to_mLs)
params.vsd.mode = 'linear_bidirectional';
if isfield(src, 'VSD_mode') && ~isempty(src.VSD_mode)
    params.vsd.mode = lower(char(src.VSD_mode));
elseif isfield(src, 'vsd_mode') && ~isempty(src.vsd_mode)
    params.vsd.mode = lower(char(src.vsd_mode));
end

if isfield(src, 'VSD_diameter_mm') && ~isnan(src.VSD_diameter_mm)
    params.vsd.diameter_mm = src.VSD_diameter_mm;
    params.vsd.area_mm2 = pi * (src.VSD_diameter_mm / 2)^2;
else
    params.vsd.diameter_mm = 0;
    params.vsd.area_mm2 = 0;
end

switch scenario
    case 'pre_surgery'
        R_vsd = params.R.vsd;
        has_gradient = isfield(src, 'VSD_gradient_mmHg') && ~isnan(src.VSD_gradient_mmHg);
        has_flow = isfield(src, 'Q_shunt_Lmin') && ~isnan(src.Q_shunt_Lmin);

        if has_gradient && has_flow
            delta_P = 0.5 * src.VSD_gradient_mmHg;
            Q_shunt = src.Q_shunt_Lmin * Lmin_to_mLs;
            if Q_shunt > 1e-4
                R_vsd = delta_P / Q_shunt;
            end
        elseif params.vsd.area_mm2 > 0
            R_vsd = estimate_linear_vsd_resistance(params, src);
        else
            R_vsd = max(params.R.vsd, 0.01);
        end

        params.R.vsd = R_vsd;

    case 'post_surgery'
        params.R.vsd = R_VSD_CLOSED;
        params.vsd.area_mm2 = 0;
        params.vsd.diameter_mm = 0;
        params.vsd.mode = 'linear_bidirectional';
end
end

% RECONCILE_VASCULAR_V0 — reconcile vascular V0 with BV target [mL]
function params = reconcile_vascular_v0(params, patient, clinical, scenario)
if nargin < 3
    clinical = [];
end
if nargin < 4 || isempty(scenario)
    scenario = 'pre_surgery';
end

[src, scenario] = resolve_clinical_source(clinical, scenario);
BV_patient = resolve_total_blood_volume(patient, src);

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

function target_WU = select_resistance_target(src, resistance_name)
doc_value = NaN;
derived_value = NaN;

switch upper(resistance_name)
    case 'SVR'
        if isfield(src, 'SVR_WU')
            doc_value = src.SVR_WU;
        end
        if isfield(src, 'SAP_mean_mmHg') && isfield(src, 'RAP_mean_mmHg') && ...
                isfield(src, 'CO_Lmin') && all(~isnan([src.SAP_mean_mmHg src.RAP_mean_mmHg src.CO_Lmin]))
            derived_value = (src.SAP_mean_mmHg - src.RAP_mean_mmHg) / max(src.CO_Lmin, 1e-6);
        end
    case 'PVR'
        if isfield(src, 'PVR_WU')
            doc_value = src.PVR_WU;
        end
        if isfield(src, 'PAP_mean_mmHg') && isfield(src, 'LAP_mean_mmHg') && ...
                isfield(src, 'CO_Lmin') && isfield(src, 'QpQs') && ...
                all(~isnan([src.PAP_mean_mmHg src.LAP_mean_mmHg src.CO_Lmin src.QpQs]))
            Qpul_Lmin = src.CO_Lmin * src.QpQs;
            derived_value = (src.PAP_mean_mmHg - src.LAP_mean_mmHg) / max(Qpul_Lmin, 1e-6);
        end
    otherwise
        error('select_resistance_target:unknownResistance', ...
            'Unknown resistance selector: %s', resistance_name);
end

target_WU = first_non_nan(doc_value, derived_value, 1.0);

    if strcmpi(resistance_name, 'SVR') && ~isnan(doc_value) && ~isnan(derived_value)
        relative_gap = abs(doc_value - derived_value) / max(abs(derived_value), 1e-6);
        if relative_gap > 0.20
            target_WU = derived_value;
        end
    end
end

function BV_patient = resolve_total_blood_volume(patient, src)
if isfield(src, 'BV_total_mL') && ~isnan(src.BV_total_mL)
    BV_patient = src.BV_total_mL;
else
    BV_patient = patient.weight_kg * blood_volume_per_kg(patient.age_years);
end
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

function R_vsd = estimate_linear_vsd_resistance(params, src)
if isfield(src, 'DeltaP_VSD_peak_mmHg') && ~isnan(src.DeltaP_VSD_peak_mmHg)
    DP_ref = 0.5 * src.DeltaP_VSD_peak_mmHg;
elseif isfield(src, 'VSD_gradient_mmHg') && ~isnan(src.VSD_gradient_mmHg)
    DP_ref = 0.5 * src.VSD_gradient_mmHg;
else
    DP_ref = params.vsd.reference_gradient_mmHg;
end

A_m2 = params.vsd.area_mm2 * 1e-6;
DP_Pa = DP_ref * params.conv.mmHg_to_Pa;
Q_m3s = params.vsd.Cd * A_m2 * sqrt(2 * DP_Pa / params.vsd.rho_blood);
Q_mLs = Q_m3s * params.conv.m3_to_mL;
R_vsd = DP_ref / max(Q_mLs, 1e-6);
end

function params = apply_chamber_tuning_from_clinical(params, src)
if ~(isfield(src, 'LVEDV_mL') && isfield(src, 'LVESV_mL') && ...
        isfield(src, 'RVEDV_mL') && isfield(src, 'RVESV_mL'))
    return;
end

if any(isnan([src.LVEDV_mL, src.LVESV_mL, src.RVEDV_mL, src.RVESV_mL]))
    return;
end

P_lv_ed = first_valid(src, {'LVEDP_mmHg', 'LAP_mean_mmHg'}, 8);
P_lv_es = first_valid(src, {'SAP_sys_mmHg', 'SAP_mean_mmHg', 'MAP_mmHg'}, 80);
P_rv_ed = first_valid(src, {'RVEDP_mmHg', 'RAP_mean_mmHg'}, 5);
P_rv_es = first_valid(src, {'PAP_sys_mmHg', 'PAP_mean_mmHg'}, 25);

f_min_LV = 0.03;
f_min_RV = 0.04;
LAP_fill = first_valid(src, {'LAP_mean_mmHg', 'LVEDP_mmHg'}, P_lv_ed);
RAP_fill = first_valid(src, {'RAP_mean_mmHg', 'RVEDP_mmHg'}, P_rv_ed);

A_LV = LAP_fill / max(f_min_LV * P_lv_es, 1e-6);
A_RV = RAP_fill / max(f_min_RV * P_rv_es, 1e-6);

denom_LV = 1 - A_LV;
V0_cand_LV = (src.LVEDV_mL - A_LV * src.LVESV_mL) / denom_LV;
if isfinite(V0_cand_LV) && V0_cand_LV > 0 && V0_cand_LV < src.LVESV_mL
    params.V0.LV = max(V0_cand_LV, 0.5);
end

denom_RV = 1 - A_RV;
V0_cand_RV = (src.RVEDV_mL - A_RV * src.RVESV_mL) / denom_RV;
if isfinite(V0_cand_RV) && V0_cand_RV > 0 && V0_cand_RV < src.RVESV_mL
    params.V0.RV = max(V0_cand_RV, 0.5);
end

Emax_LV = P_lv_es / max(src.LVESV_mL - params.V0.LV, 0.1);
Emax_RV = P_rv_es / max(src.RVESV_mL - params.V0.RV, 0.1);

params.E.LV.EB = max(f_min_LV * Emax_LV, 0.01);
params.E.LV.EA = max(Emax_LV - params.E.LV.EB, 0.1);
params.E.RV.EB = max(f_min_RV * Emax_RV, 0.01);
params.E.RV.EA = max(Emax_RV - params.E.RV.EB, 0.1);
params.clinical_override.chamber_tuning = 'echo_volume_elastance';

SV_LV = max(src.LVEDV_mL - src.LVESV_mL, 1);
SV_RV = max(src.RVEDV_mL - src.RVESV_mL, 1);
PP_sys = max(first_valid(src, {'SAP_sys_mmHg'}, 90) - first_valid(src, {'SAP_dia_mmHg'}, 55), 10);
PP_pul = max(first_valid(src, {'PAP_sys_mmHg'}, 30) - first_valid(src, {'PAP_dia_mmHg'}, 12), 6);

params.C.SAR = max(SV_LV / PP_sys, 0.05);
params.C.PAR = max(SV_RV / PP_pul, 0.05);
end

function params = recompute_timing(params)
T_HB = 60 / params.HR;
params.Tc_LV   = params.Tc_LV_frac   * T_HB;
params.Tr_LV   = params.Tr_LV_frac   * T_HB;
params.Tc_RV   = params.Tc_RV_frac   * T_HB;
params.Tr_RV   = params.Tr_RV_frac   * T_HB;
params.t_ac_LA = params.t_ac_LA_frac * T_HB;
params.Tc_LA   = params.Tc_LA_frac   * T_HB;
params.t_ar_LA = params.t_ac_LA + params.Tc_LA;
params.Tr_LA   = params.Tr_LA_frac   * T_HB;
params.t_ac_RA = params.t_ac_RA_frac * T_HB;
params.Tc_RA   = params.Tc_RA_frac   * T_HB;
params.t_ar_RA = params.t_ac_RA + params.Tc_RA;
params.Tr_RA   = params.Tr_RA_frac   * T_HB;
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

function value = first_non_nan(varargin)
value = NaN;
for k = 1:nargin
    candidate = varargin{k};
    if ~isnan(candidate)
        value = candidate;
        return;
    end
end
end
