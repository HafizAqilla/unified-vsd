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
%   [2] Gorlin orifice equation:
%         Gorlin R & Gorlin SG (1951). Am Heart J 41(1):1–29.
%         Baumgartner H et al. (2009). JASE 22(1):1–23.  (Cd = 0.7)
%         Levick JR (2010). Cardiovascular Physiology 5th ed. p. 13.  (rho_blood)
%   [3] Valenti (2023). Thesis — SVR/PVR distribution across RLC segments.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

R_VSD_CLOSED = 1e6;   % [mmHg·s/mL]  effectively infinite resistance post-surgery

common = clinical.common;
k      = params.conv.WU_to_R;         % [mmHg·s/mL per WU]  Wood unit conversion
Lmin_to_mLs = params.conv.Lmin_to_mLs;  % [mL/s per L/min]

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
            % Convert peak echo Doppler gradient to mean gradient via × 0.5.
            % Echo reports peak instantaneous gradient; haemodynamically relevant
            % mean gradient ≈ 0.5 × peak for sub-pulmonic VSD jets.
            % Source: Baumgartner H et al. (2009) JASE 22(1):1–23, §3.4;
            %         Hornberger LK et al. (1999) Am Heart J 138:179–86.
            mean_gradient_factor = 0.5;   % [-]  peak → mean gradient conversion
            delta_P    = pre.VSD_gradient_mmHg * mean_gradient_factor;  % [mmHg] mean
            Q_shunt    = pre.Q_shunt_Lmin * Lmin_to_mLs;               % [L/min] → [mL/s]
            if Q_shunt > 1e-4
                R_vsd = delta_P / Q_shunt;
                fprintf('[params_from_clinical] R_VSD from Ohm''s law (mean grad): %.4f mmHg·s/mL\n', R_vsd);
            end

        elseif isfield(pre, 'VSD_diameter_mm') && ~isnan(pre.VSD_diameter_mm)
            % Option B: Gorlin orifice equation
            %   Q = Cc * A * sqrt(2*ΔP/rho)  →  R = ΔP/Q
            %   Linearised to lumped resistance from defect area at a reference gradient.
            %   Cc = 0.7 (contraction / discharge coefficient)
            %     Source: Gorlin R & Gorlin SG (1951). Hydraulic formula for calculation
            %             of the area of the stenotic mitral valve. Am Heart J 41(1):1–29.
            %             Value Cd ≈ 0.7 is standard for small intracardiac orifices;
            %             see also Baumgartner H et al. (2009). Echocardiographic assessment
            %             of valve stenosis: EAE/ASE recommendations. JASE 22(1):1–23.
            %   rho_blood = 1060 kg/m³
            %     Source: Levick JR (2010). An Introduction to Cardiovascular Physiology,
            %             5th ed. Hodder Arnold, London. p. 13.
            D_mm      = pre.VSD_diameter_mm;                       % [mm]
            A_m2      = pi * (D_mm * params.conv.mm_to_m / 2)^2;  % [m²]  defect orifice area
            Cc        = 0.7;                                       % [-]   discharge coeff. (Gorlin 1951; Baumgartner 2009)
            rho_blood = 1060;                                      % [kg/m³] blood density (Levick 2010, p. 13)
            if isfield(pre, 'DeltaP_VSD_peak_mmHg') && ~isnan(pre.DeltaP_VSD_peak_mmHg) && pre.DeltaP_VSD_peak_mmHg ~= 0
                peak_to_mean_factor = 0.5;                         % approximation: mean DeltaP ~= 0.5 x peak for VSD jet
                % Peak-to-mean factor: Baumgartner et al. (1999), simplified VSD jet assumption
                DP_ref = pre.DeltaP_VSD_peak_mmHg * peak_to_mean_factor;  % [mmHg]
            else
                DP_ref = 20;                                       % [mmHg] default assumed mean LV-RV gradient
            end
            DP_Pa     = DP_ref * params.conv.mmHg_to_Pa;          % [Pa]
            Q_m3s     = Cc * A_m2 * sqrt(2 * DP_Pa / rho_blood);  % [m³/s]
            Q_mLs     = Q_m3s * params.conv.m3_to_mL;             % [mL/s]
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

%% =====================================================================
%  5. IC + ELASTANCE OVERRIDE  —  inject echo volumes & re-calibrate E
%     Bypasses BSA/BV allometric scaling for small paediatric patients.
%     Activated by src.override_IC == true (set in patient profile).
%
%  RATIONALE:
%  For infants with BSA 0.22–0.25 m² the allometric scale factor s ≈ 0.13
%  drives V0 to < 0.3 mL and E to > 20 mmHg/mL.  The resulting steady-
%  state LVEDV (< 5 mL) is non-physiological.  When echo volumes are known,
%  V0 and E are re-anchored to echo LVESV / RVEDV so the pressure-volume
%  relationship matches the clinical measurements directly.
%
%  PHYSICS (Valenti Eq 2.2 — elastance model):
%   P_ch = E_ch(t) * (V_ch - V0_ch)
%   At end-diastole (e=0): P_ed = EB * (VEDV - V0)  → EB = P_ed / (VEDV-V0)
%   At end-systole  (e=1): P_es = (EA+EB) * (VESV-V0) → EA = Emax - EB
%          where Emax = P_es / (VESV - V0)
%% =====================================================================
has_override = isfield(src, 'override_IC') && isequal(src.override_IC, true);
has_echo_vols = isfield(src,'LVEDV_mL') && ~isnan(src.LVEDV_mL) && ...
                isfield(src,'LVESV_mL') && ~isnan(src.LVESV_mL) && ...
                isfield(src,'RVEDV_mL') && ~isnan(src.RVEDV_mL) && ...
                isfield(src,'RVESV_mL') && ~isnan(src.RVESV_mL);

if has_override && has_echo_vols

    LVEDV = src.LVEDV_mL;   % [mL]  LV end-diastolic volume (echo)
    LVESV = src.LVESV_mL;   % [mL]  LV end-systolic volume  (echo)
    RVEDV = src.RVEDV_mL;   % [mL]  RV end-diastolic volume (echo)
    RVESV = src.RVESV_mL;   % [mL]  RV end-systolic volume  (echo)

    % Blood volume: prefer explicit field, fall back to weight-based estimate
    if isfield(src,'BV_total_mL') && ~isnan(src.BV_total_mL)
        BV = src.BV_total_mL;             % [mL]  from patient profile
    elseif isfield(common,'weight_kg') && ~isnan(common.weight_kg)
        BV = common.weight_kg * 85;       % [mL]  85 mL/kg infant <1yr (Lundquist 2025)
    else
        BV = LVEDV + RVEDV + 200;         % [mL]  rough fallback
        warning('params_from_clinical:noBV','BV estimated from chambers only — set BV_total_mL in profile');
    end

    % ---- Pressures for elastance calibration and V0 derivation --------
    %  Defined first so both the self-consistent V0 formula (below) and
    %  the elastance block can reference these variables.
    if isfield(src,'LVEDP_mmHg') && ~isnan(src.LVEDP_mmHg)
        P_lv_ed = src.LVEDP_mmHg;            % [mmHg]  from profile
    else
        P_lv_ed = 8;                          % [mmHg]  default: typical infant VSD
    end
    if isfield(src,'SAP_sys_mmHg') && ~isnan(src.SAP_sys_mmHg)
        P_lv_es = src.SAP_sys_mmHg;          % [mmHg]  LV peak ≈ systemic systolic
    elseif isfield(src,'SAP_mean_mmHg') && ~isnan(src.SAP_mean_mmHg)
        P_lv_es = src.SAP_mean_mmHg + 25;    % [mmHg]  MAP + 25 (typical infant PP/2)
    else
        P_lv_es = 80;                         % [mmHg]  default infant
    end
    if isfield(src,'RAP_mean_mmHg') && ~isnan(src.RAP_mean_mmHg)
        P_rv_ed = src.RAP_mean_mmHg;         % [mmHg]  RVEDP ≈ RAP
    else
        P_rv_ed = 5;                          % [mmHg]  default infant
    end
    if isfield(src,'PAP_sys_mmHg') && ~isnan(src.PAP_sys_mmHg)
        P_rv_es = src.PAP_sys_mmHg;          % [mmHg]  RV peak ≈ PA systolic
    elseif isfield(src,'PAP_mean_mmHg') && ~isnan(src.PAP_mean_mmHg)
        P_rv_es = src.PAP_mean_mmHg + 15;    % [mmHg]  PAP_mean + 15 approximation
    else
        P_rv_es = 25;                         % [mmHg]  default infant
    end

    % ---- Unstressed volumes: self-consistent with 0.03×Emax relationship
    %  DERIVATION:
    %   E_min = f_min × E_max = f_min × P_es / (VESV - V0)     ... (i)
    %   Diastolic equilibrium: VEDV = V0 + LAP / E_min           ... (ii)
    %   Substitute (i) into (ii) and solve for V0:
    %     A = LAP / (f_min × P_es)
    %     V0 = (VEDV - A × VESV) / (1 - A)
    %   This ensures LVEDV_eq at the clinical LAP equals the measured LVEDV
    %   without prescribing a fixed gamma fraction.
    %   Source: Valenti (2023) elastance model; Suga & Sagawa (1974).
    %
    %  f_min_LV = 0.03 (Kass 1988 Circ 77:1477; Lundquist 2025 Table 2)
    %  f_min_RV = 0.04 (slightly higher passive stiffness, Valenti 2023)
    f_min_LV = 0.03;    % [-]  LV ratio E_min/E_max
    f_min_RV = 0.04;    % [-]  RV ratio E_min/E_max

    % LAP for LV diastolic equilibrium
    if isfield(src,'LAP_mean_mmHg') && ~isnan(src.LAP_mean_mmHg)
        LAP_fill = src.LAP_mean_mmHg;          % [mmHg]  preferred
    else
        LAP_fill = P_lv_ed;                    % [mmHg]  fallback: LVEDP
    end
    % RAP for RV diastolic equilibrium
    RVfill = P_rv_ed;                          % [mmHg]  RVEDP ≈ RAP

    A_LV = LAP_fill / (f_min_LV * P_lv_es);   % dimensionless
    A_RV = RVfill   / (f_min_RV * P_rv_es);   % dimensionless

    % Guard: evaluate the candidate V0 directly.
    %   NOTE: A > 1 does NOT indicate failure — when A > 1 both numerator and
    %   denominator are negative, so the ratio stays positive and physiological.
    %   Only reject if the candidate V0 is non-finite, negative, or ≥ VESV
    %   (which would make dV_es ≤ 0 and Emax undefined).
    V0_cand_LV = (LVEDV - A_LV * LVESV) / (1 - A_LV);   % [mL]
    if ~isfinite(V0_cand_LV) || V0_cand_LV < 0 || V0_cand_LV >= LVESV
        params.V0.LV = max(0.10 * LVESV, 0.5);   % [mL]  fallback
    else
        params.V0.LV = max(V0_cand_LV, 0.5);     % [mL]
    end

    V0_cand_RV = (RVEDV - A_RV * RVESV) / (1 - A_RV);   % [mL]
    if ~isfinite(V0_cand_RV) || V0_cand_RV < 0 || V0_cand_RV >= RVESV
        params.V0.RV = max(0.10 * RVESV, 0.5);   % [mL]  fallback
    else
        params.V0.RV = max(V0_cand_RV, 0.5);     % [mL]
    end

    % ---- LV elastance re-calibration ----------------------------------
    %  E_max from end-systolic PV relation: Emax = P_ES / (VESV - V0)
    %  E_min (EB, diastolic) = 0.03 × Emax  (physiological ratio)
    %    Source: Suga & Sagawa (1974); Kass DA (1988) Circ 77:1477;
    %    validated in both adult and paediatric models (Lundquist 2025).
    %  NOT from diastolic pressure/volume slope — that approach over-
    %  estimates EB for small chambers where P_ED/(EDV-V0) ≈ 0.35 vs
    %  the true passive myocardial wall stiffness ≈ 0.03–0.05 × Emax.
    dV_LV_es        = max(LVESV - params.V0.LV, 0.1);  % [mL]
    Emax_LV         = P_lv_es / dV_LV_es;              % [mmHg/mL]  peak systolic
    params.E.LV.EB  = 0.03 * Emax_LV;                  % [mmHg/mL]  diastolic (3% of Emax)
    params.E.LV.EA  = max(Emax_LV - params.E.LV.EB, 0.1);  % [mmHg/mL]  active component

    % ---- RV elastance re-calibration ----------------------------------
    dV_RV_es        = max(RVESV - params.V0.RV, 0.1);  % [mL]
    Emax_RV         = P_rv_es / dV_RV_es;              % [mmHg/mL]
    params.E.RV.EB  = 0.04 * Emax_RV;                  % [mmHg/mL]  diastolic (4% of Emax)
    params.E.RV.EA  = max(Emax_RV - params.E.RV.EB, 0.1);  % [mmHg/mL]

    % ---- Clinical pressures for IC (with safe fallbacks) ---------------
    %  Use known catheter/echo pressures to seed each compartment at the
    %  correct operating pressure.  This is critical for small infants where
    %  the allometric-scaled nominal pressures (Section C of apply_scaling)
    %  differ substantially from the pathological steady-state.
    %
    %  P_PVEN_target drives LAP → LV filling pressure → LVEDV.
    %  P_PAR_target  drives PAP → pulmonary flow → PVEN pressure.
    %  Without matching these to clinical values the model converges to a
    %  low-volume limit cycle even when the IC echo volumes are correct.
    %
    %  PHYSICS: V_ic = V0_scaled + P_target × C_scaled   (vascular)
    %           V_ic = V0_ch     + P_target / E_ch.EB     (chambers)

    % Arterial pressure → systemic
    if isfield(src,'SAP_mean_mmHg') && ~isnan(src.SAP_mean_mmHg)
        P_SAR_ic = src.SAP_mean_mmHg;         % [mmHg]  MAP
    else
        P_SAR_ic = 65;                          % [mmHg]  healthy adult nominal
    end
    % Systemic capillary pressure ≈ midpoint between MAP and CVP
    P_SC_ic = (P_SAR_ic + P_rv_ed) / 2;       % [mmHg]

    % Pulmonary arterial pressure
    if isfield(src,'PAP_mean_mmHg') && ~isnan(src.PAP_mean_mmHg)
        P_PAR_ic = src.PAP_mean_mmHg;          % [mmHg]  mean PAP (clinical)
    else
        P_PAR_ic = 15;                          % [mmHg]  nominal healthy
    end

    % Pulmonary venous / LAP: prefer explicit LAP, fall back to LVEDP.
    % LAP drives V_PVEN_ic and hence LV filling via the mitral pathway.
    % Using LAP_mean (≈ PAWP) is more accurate than LVEDP alone when
    % there is mitral inflow augmentation (elevated LA pressure in VSD).
    if isfield(src,'LAP_mean_mmHg') && ~isnan(src.LAP_mean_mmHg)
        P_PVEN_ic = src.LAP_mean_mmHg;         % [mmHg]  LAP from profile
    else
        P_PVEN_ic = P_lv_ed;                   % [mmHg]  fallback: PAWP ≈ LVEDP
    end

    % Pulmonary capillary pressure state: between PAP_mean and P_PVEN
    P_PC_ic   = (P_PAR_ic + P_PVEN_ic) / 2;   % [mmHg]

    % ---- Windkessel arterial compliance from SV / pulse pressure --------
    %   C_SAR = SV_LV / PP_sys   (Frank-Starling Windkessel)
    %   C_PAR = SV_RV / PP_pul
    %   More accurate than allometric scaling for pathological states.
    %   Source: Segers P et al. (2008) Am J Physiol 295:H154–H163.
    SV_LV_wk = max(LVEDV - LVESV, 1);  % [mL]
    SV_RV_wk = max(RVEDV - RVESV, 1);  % [mL]
    PP_sys_wk = 35;  % [mmHg]  typical infant systemic pulse pressure
    PP_pul_wk = 20;  % [mmHg]  typical infant pulmonary pulse pressure
    params.C.SAR  = SV_LV_wk / PP_sys_wk;   % [mL/mmHg]  Windkessel SAR compliance
    params.C.PAR  = SV_RV_wk / PP_pul_wk;   % [mL/mmHg]  Windkessel PAR compliance

    % ---- Vascular IC volumes from P_target × C (pressure-based) --------
    %   V_ic = V0_scaled + P_target × C_scaled
    %   (V0_scaled already set by apply_scaling Section C using BV_scale)
    V_SAR_ic  = params.V0.SAR  + P_SAR_ic  * params.C.SAR;   % [mL]
    V_SC_ic   = params.V0.SC   + P_SC_ic   * params.C.SC;    % [mL]
    V_PAR_ic  = params.V0.PAR  + P_PAR_ic  * params.C.PAR;   % [mL]
    V_PVEN_ic = params.V0.PVEN + P_PVEN_ic * params.C.PVEN;  % [mL]

    % ---- Atrial IC volumes from P / E.EB (elastance-based) -------------
    %   V_ic = V0_scaled + P_atrium / E.EB_scaled
    V_LA_ic = params.V0.LA + P_lv_ed / params.E.LA.EB;   % [mL]  LAP ≈ LVEDP
    V_RA_ic = params.V0.RA + P_rv_ed / params.E.RA.EB;   % [mL]  RAP

    % ---- SVEN: remainder for exact blood conservation -------------------
    V_SVEN_ic = BV - LVEDV - RVEDV - V_LA_ic - V_RA_ic ...
                   - V_SAR_ic - V_SC_ic - V_PAR_ic - V_PVEN_ic;  % [mL]
    V_SVEN_ic = max(V_SVEN_ic, 0.05*BV);  % [mL]  floor: ≥ 5% BV safety guard

    % Adjust V0.SVEN for P_SVEN ≈ CVP ≈ 2 mmHg at this volume
    %   P_SVEN = (V_SVEN - V0.SVEN) / C_SVEN  →  V0.SVEN = V_SVEN - 2×C_SVEN
    params.V0.SVEN = max(V_SVEN_ic - 2 * params.C.SVEN, 0.0);  % [mL]

    % ---- Flow ICs (mean CO in mL/s) ------------------------------------
    SV_ic = max(LVEDV - LVESV, 1);     % [mL]
    Q_ic  = SV_ic * params.HR / 60;    % [mL/s]

    % ---- Assemble override IC vector ------------------------------------
    idx_s = params.idx;
    ic    = params.ic.V(:);            % start from BSA-scaled baseline
    ic(idx_s.V_LV)   = LVEDV;         % [mL]  LV end-diastolic (echo)
    ic(idx_s.V_RV)   = RVEDV;         % [mL]  RV end-diastolic (echo)
    ic(idx_s.V_LA)   = V_LA_ic;       % [mL]
    ic(idx_s.V_RA)   = V_RA_ic;       % [mL]
    ic(idx_s.V_SAR)  = V_SAR_ic;      % [mL]
    ic(idx_s.V_SC)   = V_SC_ic;       % [mL]
    ic(idx_s.V_SVEN) = V_SVEN_ic;     % [mL]
    ic(idx_s.V_PAR)  = V_PAR_ic;      % [mL]
    ic(idx_s.V_PVEN) = V_PVEN_ic;     % [mL]
    ic(idx_s.Q_SAR)  = Q_ic;          % [mL/s]
    ic(idx_s.Q_SVEN) = Q_ic;          % [mL/s]
    ic(idx_s.Q_PAR)  = Q_ic;          % [mL/s]
    ic(idx_s.Q_PVEN) = Q_ic;          % [mL/s]
    ic(idx_s.P_PC)   = P_PC_ic;       % [mmHg]
    params.ic.V = ic(:)';

    fprintf('[params_from_clinical] IC override (pressure-based): BV=%.0f mL\n', BV);
    fprintf('  Pressures used: P_SAR=%.0f  P_PAR=%.0f  P_PVEN=%.0f  P_SVEN=2 mmHg\n', ...
            P_SAR_ic, P_PAR_ic, P_PVEN_ic);
    fprintf('  Volumes: V_LV=%.1f  V_RV=%.1f  V_PAR=%.1f  V_PVEN=%.1f  V_SVEN=%.1f mL\n', ...
            LVEDV, RVEDV, V_PAR_ic, V_PVEN_ic, V_SVEN_ic);
    fprintf('  E.LV: EB=%.3f EA=%.3f V0.LV=%.2f | E.RV: EB=%.3f EA=%.3f V0.RV=%.2f\n', ...
            params.E.LV.EB, params.E.LV.EA, params.V0.LV, ...
            params.E.RV.EB, params.E.RV.EA, params.V0.RV);
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
