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
% VERSION:  1.1  (ea_activation updated 2026-04-07: full Eq. 2.4 wraparound)
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

function ev = ev_activation(phi, Tvc, Tvr, T_HB)
% EV_ACTIVATION — Ventricular normalised activation  e_v ∈ [0,1]  [Eq. 2.3]
%   phi ∈ [0, T_HB),  phases: contraction / relaxation / diastole
ev = zeros(size(phi));
for i = 1:numel(phi)
    p = phi(i);
    if p <= Tvc
        ev(i) = 0.5 * (1 - cos(pi * p / Tvc));
    elseif p <= Tvc + Tvr
        ev(i) = 0.5 * (1 + cos(pi * (p - Tvc) / Tvr));
    else
        ev(i) = 0;
    end
end
end  % ev_activation

function ea = ea_activation(phi, t_ac, Tc, t_ar, Tr, T_HB)
% EA_ACTIVATION — Atrial normalised activation  e_a ∈ [0,1]  [Valenti Eq. 2.4]
%
% Implements the exact 4-part piecewise formula, including the wraparound
% continuation of relaxation from the PREVIOUS heartbeat cycle.
%
% The 4 cases (phi = phase within current beat, phi ∈ [0, T_HB)):
%
%  Case 1  0 ≤ phi ≤ t_ar + Tr − T_HB               [wraparound tail]
%          ea = 0.5·(1 + cos(π·(phi + T_HB − t_ar) / Tr))
%          Continuation of last cycle's relaxation into the new beat.
%          Only active when t_ar + Tr > T_HB.
%          Example LA: t_ac=0.75·T_HB, Tr=0.80·T_HB → wrap spans 0→0.65·T_HB.
%
%  Case 2  t_ar + Tr − T_HB < phi ≤ t_ac             [diastolic baseline]
%          ea = 0   (fully relaxed; quiescent before next kick)
%
%  Case 3  t_ac < phi ≤ t_ar  (= t_ac + Tc)          [contraction]
%          ea = 0.5·(1 − cos(π·(phi − t_ac) / Tc))
%
%  Case 4  t_ar < phi ≤ T_HB                          [relaxation]
%          ea = 0.5·(1 + cos(π·(phi − t_ar) / Tr))
%
% INPUTS:
%   phi   - phase within beat, scalar or array           [s]  phi ∈ [0, T_HB)
%   t_ac  - contraction start time                       [s]
%   Tc    - contraction duration                         [s]
%   t_ar  - relaxation start time (= t_ac + Tc)          [s]
%   Tr    - relaxation duration                          [s]
%   T_HB  - heartbeat period                             [s]
%
% REFERENCE:  Valenti (2023). Thesis. Eq. (2.4).

ea      = zeros(size(phi));
t_wrap  = t_ar + Tr - T_HB;   % upper limit of wraparound region [s]
                               % ≤ 0 means no wraparound this beat

for i = 1:numel(phi)
    p = phi(i);

    if t_wrap > 0 && p <= t_wrap
        % Case 1 — tail of PREVIOUS cycle's relaxation (wraparound)
        ea(i) = 0.5 * (1 + cos(pi * (p + T_HB - t_ar) / Tr));

    elseif p <= t_ac
        % Case 2 — diastolic baseline (fully relaxed)
        ea(i) = 0;

    elseif p <= t_ar
        % Case 3 — atrial contraction
        ea(i) = 0.5 * (1 - cos(pi * (p - t_ac) / Tc));

    else
        % Case 4 — atrial relaxation (current cycle)
        ea(i) = 0.5 * (1 + cos(pi * (p - t_ar) / Tr));
    end
end
end  % ea_activation
