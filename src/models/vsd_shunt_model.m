function Q_VSD = vsd_shunt_model(P_LV, P_RV, params)
% VSD_SHUNT_MODEL — Smooth-diode VSD shunt flow (L→R only)
%   Q_VSD = gate(ΔP) · ΔP / R_vsd     [mL/s]
%   gate  = 0.5 + 0.5·tanh(ΔP / ε)    sigmoid, ε = params.epsilon_vsd
%
% INPUTS:
%   P_LV   - left ventricular pressure  [mmHg]  (scalar or column vector)
%   P_RV   - right ventricular pressure [mmHg]  (scalar or column vector)
%   params - parameter struct; must contain:
%              params.R.vsd           [mmHg·s/mL]
%              params.epsilon_vsd     [mmHg]
%
% OUTPUTS:
%   Q_VSD  - VSD shunt flow  [mL/s]
%            Positive = left-to-right (physiologically typical).
%            Near-zero for R.vsd >> 1e4 (post-surgery: effectively closed).
%
% DIODE GATE RATIONALE:
%   epsilon_vsd = 0.1 mmHg is chosen to preserve >99% of diastolic
%   L→R flow at ΔP ≈ 2 mmHg (restrictive pediatric VSD).
%   At 0.5 mmHg the gate clips ~10% of diastolic shunt — too aggressive.
%   This VSD gate is intentionally tighter than cardiac valves
%   (epsilon_valve = 0.5 mmHg in valve_model/system_rhs).
%   The tanh gate replaces a hard if-statement so the ODE RHS remains
%   continuously differentiable, preventing step-size collapse in ode15s.
%
% REFERENCES:
%   [1] system_rhs.m — VSD shunt ODE context.
%   [2] Valenti (2023). Thesis. §2.8.3 (VSD shunt model).
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-25
% VERSION:  1.0
% -----------------------------------------------------------------------

dP   = P_LV - P_RV;                                       % [mmHg]
gate = 0.5 + 0.5 * tanh(dP / params.epsilon_vsd);         % smooth L→R diode
Q_VSD = gate .* dP / params.R.vsd;                        % [mL/s]
end
