function params = apply_physiological_scaling(params_ref, patient, scaling_mode)
% APPLY_PHYSIOLOGICAL_SCALING
% -----------------------------------------------------------------------
% Applies declared physiological size scaling laws without constructing
% initial conditions or repairing blood volume. Vascular V0 reconciliation
% is handled separately in apply_scaling/params_from_clinical.
%
% INPUTS:
%   params_ref     - reference parameter struct                      [-]
%   patient        - patient demographics struct                     [-]
%   scaling_mode   - 'zhang' (default)                              [-]
%
% OUTPUTS:
%   params         - scaled parameter struct with unchanged params.ic [-]
%
% NOTES:
%   This function may modify only physiological parameters:
%     HR, R, C, E, V0, valve parameters, timing fields
%   It must not construct the state vector or overwrite params.ic.V.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3 || isempty(scaling_mode)
    scaling_mode = 'zhang';
end

scaling_mode = lower(char(scaling_mode));
switch scaling_mode
    case 'zhang'
        params = apply_zhang_scaling(params_ref, patient);
    otherwise
        error('apply_physiological_scaling:unknownMode', ...
            'Unsupported scaling mode: %s', scaling_mode);
end

params.scaling.mode = scaling_mode;

end

function params = apply_zhang_scaling(params_ref, patient)
% APPLY_ZHANG_SCALING -- weight-based physiological scaling only.

W_ref = 70;   % [kg]
w = patient.weight_kg / W_ref;

if isfield(patient, 'BSA') && ~isnan(patient.BSA)
    BSA_patient = patient.BSA;
else
    BSA_patient = sqrt(patient.height_cm * patient.weight_kg / 3600);
end

params = params_ref;
params.scaling.W_ref       = W_ref;
params.scaling.w           = w;
params.scaling.BSA_ref     = 1.73;
params.scaling.BSA_patient = BSA_patient;
params.scaling.patient     = patient;

eHR    = -0.30;
eE_lv  = -0.50;
eE_rv  = -0.75;
eV0    = +0.80;
eR_sys = -0.475;
eR_pul = -0.70;
eC     = +1.00;
eRv_op = -0.90;

params.HR = params_ref.HR * w^eHR;

params.E.LV.EA = params_ref.E.LV.EA * w^eE_lv;
params.E.LV.EB = params_ref.E.LV.EB * w^eE_lv;
params.E.LA.EA = params_ref.E.LA.EA * w^eE_lv;
params.E.LA.EB = params_ref.E.LA.EB * w^eE_lv;
params.E.RV.EA = params_ref.E.RV.EA * w^eE_rv;
params.E.RV.EB = params_ref.E.RV.EB * w^eE_rv;
params.E.RA.EA = params_ref.E.RA.EA * w^eE_rv;
params.E.RA.EB = params_ref.E.RA.EB * w^eE_rv;

params.V0.LV = params_ref.V0.LV * w^eV0;
params.V0.RV = params_ref.V0.RV * w^eV0;
params.V0.LA = params_ref.V0.LA * w^eV0;
params.V0.RA = params_ref.V0.RA * w^eV0;
params.V0.SAR = params_ref.V0.SAR * w^eV0;
params.V0.SC = params_ref.V0.SC * w^eV0;
params.V0.SVEN = params_ref.V0.SVEN * w^eV0;
params.V0.PAR = params_ref.V0.PAR * w^eV0;
params.V0.PVEN = params_ref.V0.PVEN * w^eV0;
if isfield(params_ref.V0, 'PCOX'), params.V0.PCOX = params_ref.V0.PCOX * w^eV0; end
if isfield(params_ref.V0, 'PCNO'), params.V0.PCNO = params_ref.V0.PCNO * w^eV0; end

params.R.SAR  = params_ref.R.SAR  * w^eR_sys;
params.R.SC   = params_ref.R.SC   * w^eR_sys;
params.R.SVEN = params_ref.R.SVEN * w^eR_sys;
params.R.PAR  = params_ref.R.PAR  * w^eR_pul;
params.R.PCOX = params_ref.R.PCOX * w^eR_pul;
params.R.PCNO = params_ref.R.PCNO * w^eR_pul;
params.R.PVEN = params_ref.R.PVEN * w^eR_pul;

params.C.SAR  = params_ref.C.SAR  * w^eC;
params.C.SC   = params_ref.C.SC   * w^eC;
params.C.SVEN = params_ref.C.SVEN * w^eC;
params.C.PAR  = params_ref.C.PAR  * w^eC;
params.C.PCOX = params_ref.C.PCOX * w^eC;
params.C.PCNO = params_ref.C.PCNO * w^eC;
params.C.PVEN = params_ref.C.PVEN * w^eC;
if isfield(params_ref.C, 'RA'), params.C.RA = params_ref.C.RA * w^eC; end
if isfield(params_ref.C, 'LA'), params.C.LA = params_ref.C.LA * w^eC; end

params.Rvalve.open   = params_ref.Rvalve.open * w^eRv_op;
params.Rvalve.closed = params_ref.Rvalve.closed;

T_HB = 60 / params.HR;
params = recompute_timing(params, T_HB);
params.scaling.zhang_exponents = struct( ...
    'HR', eHR, 'E_lv', eE_lv, 'E_rv', eE_rv, 'V0', eV0, ...
    'R_sys', eR_sys, 'R_pul', eR_pul, 'C', eC, 'Rvalve_open', eRv_op);
end

function params = recompute_timing(params, T_HB)
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
