function params = params_from_clinical(params, clinical, scenario)
% PARAMS_FROM_CLINICAL
% -----------------------------------------------------------------------
% Maps clinical measurements into model parameters.
% Handles both pre-surgery and post-surgery scenarios via the 'scenario'
% flag, which controls how R_VSD is assigned.
%
% INPUTS:
%   params    - parameter struct (from apply_scaling.m)
%   clinical  - unified clinical struct (from patient_template.m)
%               with sub-fields: .common, .pre_surgery, .post_surgery
%   scenario  - string: 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   params    - updated parameter struct with:
%               - HR from clinical.common
%               - SVR/PVR distributed across compartments
%               - R.vsd assigned per scenario
%               - absolute timing fields updated to match HR
%
% KEY SCENARIO BEHAVIOUR:
%   'pre_surgery':
%     R.vsd is COMPUTED from VSD clinical data (gradient + shunt flow, or
%     orifice equation from diameter).  Finite R.vsd enables the shunt.
%   'post_surgery':
%     R.vsd is set to R_VSD_CLOSED (1e6 mmHg·s/mL).  The ODE Q_VSD ≈ 0,
%     which is physically equivalent to surgical VSD closure.
%
% SIGN CONVENTIONS:
%   Q_VSD > 0 = LV→RV (left-to-right shunt, typical in unrepaired VSD).
%
% REFERENCES:
%   [1] Ohm's law for vascular resistances: R = ΔP / Q
%   [2] Gorlin orifice equation for R_VSD from defect diameter (see below).
%   [3] Valenti (2023). Thesis — SVR/PVR distribution across RLC segments.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

R_VSD_CLOSED = 1e6;   % [mmHg·s/mL]  effectively infinite resistance post-surgery

common = clinical.common;
k      = params.conv.WU_to_R;   % 0.06  mmHg·s/mL per Wood unit

%% =====================================================================
%  1. HEART RATE
%% =====================================================================
if isfield(common, 'HR') && ~isnan(common.HR)
    params.HR = common.HR;
    % Re-compute absolute timing with new HR
    T_HB = 60 / params.HR;
    params = recompute_timing(params, T_HB);
end

%% =====================================================================
%  2. SYSTEMIC VASCULAR RESISTANCE  →  distribute across SAR / SC / SVEN
%% =====================================================================
switch scenario
    case 'pre_surgery'
        svr_field = 'SVR_WU';
        pvr_field = 'PVR_WU';
        src = clinical.pre_surgery;
    case 'post_surgery'
        svr_field = 'SVR_WU';
        pvr_field = 'PVR_WU';
        src = clinical.post_surgery;
    otherwise
        error('params_from_clinical:unknownScenario', ...
              'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

if isfield(src, svr_field) && ~isnan(src.(svr_field))
    SVR_mech = src.(svr_field) * k;                    % [mmHg·s/mL]
    SVR_ref  = params.R.SAR + params.R.SC + params.R.SVEN;
    ratio    = SVR_mech / SVR_ref;
    params.R.SAR  = params.R.SAR  * ratio;
    params.R.SC   = params.R.SC   * ratio;
    params.R.SVEN = params.R.SVEN * ratio;
end

%% =====================================================================
%  3. PULMONARY VASCULAR RESISTANCE  →  distribute across PAR / PCOX||PCNO / PVEN
%% =====================================================================
if isfield(src, pvr_field) && ~isnan(src.(pvr_field))
    PVR_mech  = src.(pvr_field) * k;
    R_cap_ref = 1 / (1/params.R.PCOX + 1/params.R.PCNO);   % parallel combination
    PVR_ref   = params.R.PAR + R_cap_ref + params.R.PVEN;
    ratio     = PVR_mech / PVR_ref;
    params.R.PAR  = params.R.PAR  * ratio;
    params.R.PCOX = params.R.PCOX * ratio;
    params.R.PCNO = params.R.PCNO * ratio;
    params.R.PVEN = params.R.PVEN * ratio;
end

%% =====================================================================
%  4. VSD SHUNT RESISTANCE  —  scenario-specific assignment
%% =====================================================================
switch scenario

    case 'pre_surgery'
        pre = clinical.pre_surgery;
        R_vsd = R_VSD_CLOSED;   % default unless clinical data overrides

        % Option A: direct Ohm's law from gradient and shunt flow
        has_gradient = isfield(pre, 'VSD_gradient_mmHg') && ~isnan(pre.VSD_gradient_mmHg);
        has_flow     = isfield(pre, 'Q_shunt_Lmin')      && ~isnan(pre.Q_shunt_Lmin);

        if has_gradient && has_flow
            delta_P    = pre.VSD_gradient_mmHg;            % [mmHg]
            Q_shunt    = pre.Q_shunt_Lmin * 1000/60;       % [L/min] → [mL/s]
            if Q_shunt > 1e-4
                R_vsd = delta_P / Q_shunt;
                fprintf('[params_from_clinical] R_VSD from Ohm s law: %.4f mmHg·s/mL\n', R_vsd);
            end

        elseif isfield(pre, 'VSD_diameter_mm') && ~isnan(pre.VSD_diameter_mm)
            % Option B: Gorlin orifice equation
            %   Q = Cc * A * sqrt(2*ΔP/rho)  →  R = ΔP/Q = sqrt(ΔP/2*rho) / (Cc*A)
            %   Simplified to lumped resistance from defect area.
            %   Cc = 0.7 (contraction coefficient), rho = 1060 kg/m³ (blood)
            %   Pressure reference ~ mean LV-RV gradient ≈ 0.5 * VSD_gradient
            D_mm   = pre.VSD_diameter_mm;
            A_m2   = pi * (D_mm / 2000)^2;   % defect area [m²]
            Cc     = 0.7;
            rho    = 1060;                    % kg/m³
            DP_ref = 20;                      % mmHg  assumed mean gradient if not supplied
            DP_Pa  = DP_ref * 133.322;        % mmHg → Pa
            Q_m3s  = Cc * A_m2 * sqrt(2 * DP_Pa / rho);   % m³/s
            Q_mLs  = Q_m3s * 1e6;                          % m³/s → mL/s
            DP_ref_mLs = DP_ref / max(Q_mLs, 1e-6);
            R_vsd  = DP_ref_mLs;
            fprintf('[params_from_clinical] R_VSD from Gorlin orifice (D=%.1f mm): %.4f mmHg·s/mL\n', ...
                    D_mm, R_vsd);
        else
            warning('params_from_clinical:noVSDdata', ...
                ['No VSD gradient/flow/diameter provided for pre_surgery. ' ...
                 'R_VSD set to default (%.0f) — calibrate R.vsd.'], R_VSD_CLOSED);
        end

        params.R.vsd = R_vsd;

    case 'post_surgery'
        params.R.vsd = R_VSD_CLOSED;
        fprintf('[params_from_clinical] Post-surgery: R_VSD = %.0e mmHg·s/mL (closed)\n', ...
                R_VSD_CLOSED);
end

end  % params_from_clinical

% =========================================================================
%  LOCAL HELPER
% =========================================================================

function params = recompute_timing(params, T_HB)
% RECOMPUTE_TIMING — convert fractional phase parameters to absolute seconds
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
