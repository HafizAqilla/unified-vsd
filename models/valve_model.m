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
%               params.Rvalve.open    R_min (open)    [mmHg·s/mL]
%               params.Rvalve.closed  R_max (closed)  [mmHg·s/mL]
%
% OUTPUTS:
%   Q       - volumetric flow through valve                 [mL/s]
%
% SIGN CONVENTIONS:
%   Q > 0 : forward (upstream → downstream)
%   Q < 0 : would imply regurgitation; prevented by R_closed >> R_open.
%
% ASSUMPTIONS:
%   - Incompressible flow; inertia neglected at valve level.
%   - Numerical closure by large R (no hard floor on Q).
%   - All four cardiac valves use identical open/closed resistances here;
%     per-valve overrides can be added by extending params.Rvalve.
%
% REFERENCES:
%   [1] Valenti (2023). Thesis. Eq. (2.6).
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

if P_up > P_down
    R = params.Rvalve.open;
else
    R = params.Rvalve.closed;
end
Q = (P_up - P_down) / R;

end
