function Q = valve_model(P_up, P_down, params)
% VALVE_MODEL
% -----------------------------------------------------------------------
% Smooth-diode cardiac valve model (tanh-blended two-resistance).
%
% Implements Valenti Eq. (2.6) with a continuous tanh gate so the ODE RHS
% remains differentiable everywhere — preventing step-size collapse in
% ode15s around valve opening/closing events.
%
% VALVE LOGIC:
%   gate = 0.5 + 0.5·tanh(dP / ε)          ∈ [0, 1]
%   Q    = dP · (gate/R_open + (1–gate)/R_closed)
%
%   dP >> 0  →  gate → 1  →  Q ≈ dP / R_open   (forward, open)
%   dP << 0  →  gate → 0  →  Q ≈ dP / R_closed  (reverse, closed — tiny bleed)
%
% INPUTS:
%   P_up    - upstream pressure                             [mmHg]
%   P_down  - downstream pressure                           [mmHg]
%   params  - parameter struct with:
%               params.Rvalve.open      R_min (open)    [mmHg·s/mL]
%               params.Rvalve.closed    R_max (closed)  [mmHg·s/mL]
%               params.epsilon_valve    tanh width ε    [mmHg]
%
% OUTPUTS:
%   Q       - volumetric flow through valve                 [mL/s]
%
% SIGN CONVENTIONS:
%   Q > 0 : forward (upstream → downstream)
%   Q < 0 : tiny reverse bleed through R_closed  (numerical guard only)
%
% NOTE — HOT-PATH INLINING:
%   In system_rhs.m the four cardiac valve calls are inlined directly to
%   avoid function-call overhead in the ODE hot path.  This standalone
%   function is retained for post-processing callers (compute_clinical_indices,
%   plotting_tools) and unit tests.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eq. (2.6) — two-resistance diode model.
%   [2] vsd_shunt_model.m — same gate pattern for the VSD diode.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-14
% VERSION:  3.0  (smooth tanh; re-instated from v1.1; replaces hard-switch v2.0)
% -----------------------------------------------------------------------

dP   = P_up - P_down;                               % pressure gradient  [mmHg]
gate = 0.5 + 0.5 * tanh(dP / params.epsilon_valve); % smooth diode gate  [0…1]
Q    = dP .* (gate / params.Rvalve.open + (1 - gate) / params.Rvalve.closed);

end
