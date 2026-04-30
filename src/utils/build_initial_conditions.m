function ic = build_initial_conditions(params, patient, clinical, scenario)
% BUILD_INITIAL_CONDITIONS
% -----------------------------------------------------------------------
% Builds a pressure-aware, blood-volume-consistent initial-condition
% vector without mutating any physiological parameter such as params.V0.*.
%
% INPUTS:
%   params    - parameter struct after scaling and clinical mapping   [-]
%   patient   - patient demographics struct                           [-]
%   clinical  - optional clinical struct                              [-]
%   scenario  - optional scenario string                              [-]
%
% OUTPUTS:
%   ic        - 14x1 initial condition vector                         [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3
    clinical = [];
end
if nargin < 4 || isempty(scenario)
    scenario = 'pre_surgery';
end

sidx = params.idx;
src = struct();
if ~isempty(clinical) && isstruct(clinical)
    if isfield(clinical, scenario)
        src = clinical.(scenario);
    elseif strcmp(scenario, 'pre_surgery') && isfield(clinical, 'pre_surgery')
        src = clinical.pre_surgery;
    elseif strcmp(scenario, 'post_surgery') && isfield(clinical, 'post_surgery')
        src = clinical.post_surgery;
    end
end

if patient.age_years < 1
    BV_per_kg = 82;
else
    BV_per_kg = 70;
end
BV_patient = patient.weight_kg * BV_per_kg;

P_nom_SAR = first_valid(src, {'SAP_mean_mmHg', 'MAP_mmHg'}, 65);
P_nom_SC = max(0.5 * P_nom_SAR, 15);
P_nom_SVEN = first_valid(src, {'RAP_mean_mmHg'}, 2);
P_nom_PAR = first_valid(src, {'PAP_mean_mmHg'}, 15);
P_nom_PVEN = first_valid(src, {'LAP_mean_mmHg', 'LVEDP_mmHg'}, 6);
P_nom_PC = 0.5 * (P_nom_PAR + P_nom_PVEN);
P_nom_LV_ED = first_valid(src, {'LVEDP_mmHg', 'LAP_mean_mmHg'}, 8);
P_nom_RV_ED = first_valid(src, {'RVEDP_mmHg', 'RAP_mean_mmHg'}, 4);
P_nom_LA_ED = first_valid(src, {'LAP_mean_mmHg'}, 6);
P_nom_RA_ED = first_valid(src, {'RAP_mean_mmHg'}, 4);

ic = zeros(14, 1);

ic(sidx.V_SAR) = params.V0.SAR + P_nom_SAR * params.C.SAR;
ic(sidx.V_SC) = params.V0.SC + P_nom_SC * params.C.SC;
ic(sidx.V_PAR) = params.V0.PAR + P_nom_PAR * params.C.PAR;
ic(sidx.V_PVEN) = params.V0.PVEN + P_nom_PVEN * params.C.PVEN;
ic(sidx.P_PC) = P_nom_PC;

if isfield(src, 'LVEDV_mL') && ~isnan(src.LVEDV_mL)
    ic(sidx.V_LV) = src.LVEDV_mL;
else
    ic(sidx.V_LV) = params.V0.LV + P_nom_LV_ED / max(params.E.LV.EB, 1e-6);
end
if isfield(src, 'RVEDV_mL') && ~isnan(src.RVEDV_mL)
    ic(sidx.V_RV) = src.RVEDV_mL;
else
    ic(sidx.V_RV) = params.V0.RV + P_nom_RV_ED / max(params.E.RV.EB, 1e-6);
end

ic(sidx.V_LA) = params.V0.LA + P_nom_LA_ED / max(params.E.LA.EB, 1e-6);
ic(sidx.V_RA) = params.V0.RA + P_nom_RA_ED / max(params.E.RA.EB, 1e-6);

Q_init = initial_flow_seed(params, src, ic(sidx.V_LV));
ic(sidx.Q_SAR) = Q_init;
ic(sidx.Q_SVEN) = Q_init;
ic(sidx.Q_PAR) = Q_init;
ic(sidx.Q_PVEN) = Q_init;

ic(sidx.V_SVEN) = params.V0.SVEN + P_nom_SVEN * params.C.SVEN;

vol_idx = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV ...
           sidx.V_SAR sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];
BV_current = sum(ic(vol_idx));
BV_residual = BV_patient - BV_current;

% Distribute any residual across vascular volumes by compliance to avoid
% dumping all mismatch into stressed venous volume alone.
vascular_idx = [sidx.V_SAR sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];
vascular_C = [params.C.SAR; params.C.SC; params.C.SVEN; params.C.PAR; params.C.PVEN];
if abs(BV_residual) > 1e-6
    if abs(BV_residual) > 0.05 * BV_patient
        ic(sidx.V_SVEN) = ic(sidx.V_SVEN) + BV_residual;
    else
        weights = vascular_C / sum(vascular_C);
        ic(vascular_idx) = ic(vascular_idx) + BV_residual * weights;
    end
end

end

function Q_init = initial_flow_seed(params, src, V_LV_ed)
if isfield(src, 'CO_Lmin') && ~isnan(src.CO_Lmin)
    Q_init = src.CO_Lmin * params.conv.Lmin_to_mLs;
    return;
end

SV_seed = max(0.25 * V_LV_ed, 5.0);
Q_init = SV_seed * params.HR / 60;
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
