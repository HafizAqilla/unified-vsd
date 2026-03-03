function Q = valve_model(P_up, P_down, params)
% VALVE_MODEL
% -----------------------------------------------------------------------
% Non-ideal diode cardiac valve model (Valenti Eq. 2.6).
%
% The valve opens when upstream pressure exceeds downstream pressure and
% closes otherwise.  Two discrete resistance values model open / closed
% states; no dynamical valve mass is included.
%
% INPUTS:
%   P_up    - upstream pressure                             [mmHg]
%   P_down  - downstream pressure                           [mmHg]
%   params  - parameter struct with:
%               params.Rvalve.open      R_min (open)    [mmHg·s/mL]
%               params.Rvalve.closed    R_max (closed)  [mmHg·s/mL]
%               params.epsilon_valve    smooth-switch width  [mmHg]
%
% OUTPUTS:
%   Q       - volumetric flow through valve                 [mL/s]
%
% SIGN CONVENTIONS:
%   Q > 0 : forward (upstream → downstream)
%   Q < 0 : small reverse flow through R_closed (regurgitation prevented
%            by R_closed >> R_open; not set to hard zero to avoid
%            discontinuities in the ODE).
%
% NUMERICAL METHOD (Guardrail §8.4):
%   Hard if/else switching creates near-discontinuities that collapse the
%   ODE solver step size.  This implementation uses a smooth blending:
%
%     weight = 0.5 + 0.5·tanh(dP / epsilon_valve)   ∈ (0,1)
%     R_eff  = weight·R_open + (1-weight)·R_closed
%     Q      = dP / R_eff
%
%   where epsilon_valve [mmHg] controls the width of the transition zone.
%   Smaller epsilon → sharper physiological switch but increased stiffness.
%   Larger epsilon → smoother numerics but blunted valve closure.
%   Default epsilon_valve = 0.5 mmHg (< 1% of typical dP across open valve).
%
%   At dP = 0:  weight = 0.5; R_eff = 0.5·(R_open+R_closed) ≈ R_closed/2
%               → Q ≈ 0 (effectively closed).
%   At dP >> epsilon:  R_eff → R_open  (fully open).
%   At dP << -epsilon: R_eff → R_closed (fully closed).
%
% ASSUMPTIONS:
%   - Incompressible flow; inertia neglected at valve level.
%   - All four cardiac valves share identical open/closed resistances.
%     Per-valve overrides can be added by extending params.Rvalve.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eq. (2.6).
%   [2] Guardrail §8.4 — smooth switching for stiff ODE robustness.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.1
% -----------------------------------------------------------------------

dP     = P_up - P_down;                           % pressure gradient  [mmHg]
eps_v  = params.epsilon_valve;                    % transition width   [mmHg]

% Smooth weight: 1 = fully open, 0 = fully closed
weight = 0.5 + 0.5 * tanh(dP / eps_v);           % [dimensionless]

% Effective resistance: blends open and closed states continuously
R_eff  = weight * params.Rvalve.open + (1 - weight) * params.Rvalve.closed;  % [mmHg·s/mL]

Q = dP / R_eff;                                   % volumetric flow    [mL/s]

end
