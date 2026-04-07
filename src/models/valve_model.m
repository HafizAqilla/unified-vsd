function Q = valve_model(P_up, P_down, params)
% VALVE_MODEL
% -----------------------------------------------------------------------
% Ideal hard-switch (diode) cardiac valve model.
%
% Replaces the earlier smooth tanh-blended model (Valenti Eq. 2.6 variant)
% with a strict pressure-threshold switch to test whether the epsilon_valve
% smoothing was blunting isovolumetric phases and distorting PV loop shape.
%
% VALVE LOGIC:
%   dP > 0  →  valve OPEN  → R_eff = R_open  (forward flow)
%   dP ≤ 0  →  valve CLOSED → R_eff = R_closed (reverse pressure blocked)
%
% Note: R_closed >> R_open (≈ 9.4×10⁴ vs 6.3×10⁻³ mmHg·s/mL) so a tiny
% reverse bleed exists numerically.  This is intentional — hard zero flow
% is not used because it introduces a true discontinuity that can degrade
% solver accuracy on the LV isovolumetric phases.
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
%   Q > 0 : forward (upstream → downstream)
%   Q < 0 : tiny reverse bleed through R_closed
%
% EXPECTED EFFECT vs SMOOTH MODEL:
%   - Sharper PV loop corners (more rectangular shape)
%   - More vertical isovolumetric contraction / relaxation phases
%   - More abrupt valve opening and closing
%   - Possible increase in solver step rejections — monitor simulation time
%
% ⚠ SOLVER NOTE:
%   Hard switching may increase ode15s step-size rejections if the valve
%   event occurs mid-step.  If simulations become very slow or unstable,
%   consider reverting to the tanh model (epsilon_valve = 0.5 mmHg) or
%   using odeset('Events', ...) to trigger step resets at dP = 0.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eq. (2.6) — two-resistance diode model.
%   [2] Guardrail §8.4 — numerical robustness note on valve switching.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-07
% VERSION:  2.0  (hard-switch; replaces tanh-blended v1.1)
% -----------------------------------------------------------------------

dP = P_up - P_down;   % pressure gradient across valve  [mmHg]

if dP > 0
    R_eff = params.Rvalve.open;     % [mmHg·s/mL]  valve open
else
    R_eff = params.Rvalve.closed;   % [mmHg·s/mL]  valve closed
end

Q = dP / R_eff;   % volumetric flow  [mL/s]

end
