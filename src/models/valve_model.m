function Q = valve_model(P_up, P_down, params)
% VALVE_MODEL
% -----------------------------------------------------------------------
% Non-ideal diode cardiac valve model (Valenti 2023).
%
% Implements the two-resistance diode formulation from Valenti (2023) Eq.
% (2.6) and Table 3.3.  Each valve is characterised by a minimum (open)
% resistance R_min and a maximum (closed) resistance R_max:
%
%   R_i(p1, p2) = { R_min,  p1 > p2   (katup terbuka / valve open)
%                 { R_max,  p1 <= p2   (katup tertutup / valve closed)
%
% Parameter values (Valenti 2023, Table 3.3 — healthy adult):
%   R_min = 6.2872e-3 mmHg·s/mL   (open  — params.Rvalve.open)
%   R_max = 9.4168e+4 mmHg·s/mL   (closed — params.Rvalve.closed)
%
% VALVE LOGIC:
%   dP > 0  ->  valve OPEN   -> R_eff = R_min  (forward flow)
%   dP <= 0 ->  valve CLOSED -> R_eff = R_max  (reverse pressure blocked)
%
% Note: R_max >> R_min so a tiny reverse bleed exists numerically.
% This is intentional — hard zero flow introduces a true discontinuity
% that can degrade ode15s accuracy during isovolumetric phases.
%
% INPUTS:
%   P_up    - upstream pressure                             [mmHg]
%   P_down  - downstream pressure                           [mmHg]
%   params  - parameter struct with:
%               params.Rvalve.open      R_min (open)    [mmHg·s/mL]
%               params.Rvalve.closed    R_max (closed)  [mmHg·s/mL]
%
% OUTPUTS:
%   Q       - volumetric flow through valve                 [mL/s]
%
% SIGN CONVENTIONS:
%   Q > 0 : forward (upstream -> downstream)
%   Q < 0 : tiny reverse bleed through R_max
%
% EXPECTED EFFECT vs SMOOTH MODEL:
%   - Sharper PV loop corners (more rectangular shape)
%   - More vertical isovolumetric contraction / relaxation phases
%   - More abrupt valve opening and closing
%   - Possible increase in solver step rejections — monitor simulation time
%
% SOLVER NOTE:
%   Hard switching may increase ode15s step-size rejections if the valve
%   event occurs mid-step.  If simulations become very slow or unstable,
%   consider reverting to the tanh model (epsilon_valve = 0.5 mmHg) or
%   using odeset('Events', ...) to trigger step resets at dP = 0.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eq. (2.6) — two-resistance diode model.
%   [2] Valenti (2023). Thesis. Table 3.3 — R_min=6.2872e-3, R_max=9.4168e4.
%   [3] Guardrail §8.4 — numerical robustness note on valve switching.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-07
% VERSION:  2.1  (Valenti Table 3.3 R_min/R_max; replaces Bozkurt-averaged open R)
% -----------------------------------------------------------------------

dP = P_up - P_down;   % pressure gradient across valve  [mmHg]

if dP > 0
    R_eff = params.Rvalve.open;     % [mmHg·s/mL]  R_min — valve open
else
    R_eff = params.Rvalve.closed;   % [mmHg·s/mL]  R_max — valve closed
end

Q = dP / R_eff;   % volumetric flow  [mL/s]

end
