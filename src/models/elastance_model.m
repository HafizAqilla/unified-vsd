function [E_LV, E_RV, E_LA, E_RA] = elastance_model(t, params)
% ELASTANCE_MODEL
% -----------------------------------------------------------------------
% Time-varying elastance for all four cardiac chambers.
%
% Implements the piecewise activation formulation from Valenti's thesis:
%   Eq. (2.2):  E_ch(t) = E_A · e(t) + E_B         [mmHg/mL]
%   Eq. (2.3):  e_v(t)  — ventricular activation
%   Eq. (2.4):  e_a(t)  — atrial activation (time-shifted)
%
% Fully vectorised: t may be scalar or any array.
%
% INPUTS:
%   t       - time, scalar or array                         [s]
%   params  - parameter struct (from apply_scaling.m)
%             Required fields:
%               params.HR              [bpm]
%               params.Tc_LV, Tr_LV    [s]
%               params.Tc_RV, Tr_RV    [s]
%               params.t_ac_LA, Tc_LA, t_ar_LA, Tr_LA  [s]
%               params.t_ac_RA, Tc_RA, t_ar_RA, Tr_RA  [s]
%               params.E.LV.EA, E.LV.EB   [mmHg/mL]
%               params.E.RV.EA, E.RV.EB
%               params.E.LA.EA, E.LA.EB
%               params.E.RA.EA, E.RA.EB
%
% OUTPUTS:
%   E_LV, E_RV, E_LA, E_RA  — instantaneous elastance (same size as t)
%                                                           [mmHg/mL]
%
% ASSUMPTIONS:
%   - Ventricles begin contraction at the start of each beat (phase = 0).
%   - Atrial kick occurs near end-diastole (t_ac > 0.75 T_HB).
%   - All timing parameters in absolute seconds (computed by apply_scaling.m
%     or params_from_clinical.m from fractional values in default_parameters.m).
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eqs. (2.2)–(2.4), Table 2.1.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

T_HB = 60 / params.HR;          % heartbeat period [s]
phi  = mod(t, T_HB);             % phase within current beat [0, T_HB)

%% Ventricular activation  e_v(phi)  [Eq. 2.3]
e_LV = ev_activation(phi, params.Tc_LV, params.Tr_LV, T_HB);
e_RV = ev_activation(phi, params.Tc_RV, params.Tr_RV, T_HB);

%% Atrial activation  e_a(phi)  [Eq. 2.4]
e_LA = ea_activation(phi, params.t_ac_LA, params.Tc_LA, ...
                          params.t_ar_LA, params.Tr_LA, T_HB);
e_RA = ea_activation(phi, params.t_ac_RA, params.Tc_RA, ...
                          params.t_ar_RA, params.Tr_RA, T_HB);

%% Elastance  E_ch = E_A * e + E_B  [Eq. 2.2]
E_LV = params.E.LV.EA * e_LV + params.E.LV.EB;
E_RV = params.E.RV.EA * e_RV + params.E.RV.EB;
E_LA = params.E.LA.EA * e_LA + params.E.LA.EB;
E_RA = params.E.RA.EA * e_RA + params.E.RA.EB;

end  % elastance_model

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function ev = ev_activation(phi, Tvc, Tvr, ~)
% EV_ACTIVATION — Ventricular normalised activation  e_v ∈ [0,1]  [Eq. 2.3]
%   Fully vectorised: phi may be scalar or any array.
%   Phases: contraction [0, Tvc] / relaxation (Tvc, Tvc+Tvr] / diastole (rest=0)
ev = zeros(size(phi));
m1 = phi <= Tvc;
m2 = (phi > Tvc) & (phi <= Tvc + Tvr);
ev(m1) = 0.5 * (1 - cos(pi * phi(m1) / Tvc));
ev(m2) = 0.5 * (1 + cos(pi * (phi(m2) - Tvc) / Tvr));
% diastolic phase: ev = 0 (already initialised by zeros)
end  % ev_activation

function ea = ea_activation(phi, t_ac, Tc, t_ar, Tr, ~)
% EA_ACTIVATION — Atrial normalised activation  e_a ∈ [0,1]  [Eq. 2.4]
%   Fully vectorised: phi may be scalar or any array.
%   t_ac = contraction start, t_ar = relaxation start (= t_ac + Tc)
ea = zeros(size(phi));
m1 = (phi >= t_ac) & (phi < t_ar);
m2 = (phi >= t_ar) & (phi < t_ar + Tr);
ea(m1) = 0.5 * (1 - cos(pi * (phi(m1) - t_ac) / Tc));
ea(m2) = 0.5 * (1 + cos(pi * (phi(m2) - t_ar) / Tr));
% resting phase: ea = 0 (already initialised by zeros)
end  % ea_activation
