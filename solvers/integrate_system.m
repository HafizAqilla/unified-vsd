function sim = integrate_system(params)
% INTEGRATE_SYSTEM
% -----------------------------------------------------------------------
% Numerically integrates the 14-state cardiovascular ODE to periodic
% steady-state, then returns only the last N cardiac cycles.
%
% INPUTS:
%   params  - parameter struct (from apply_scaling.m + params_from_clinical.m)
%             Required fields:
%               params.HR                  [bpm]
%               params.sim.nCyclesSteady   number of cycles simulated
%               params.sim.nCyclesKeep     number of last cycles returned
%               params.sim.rtol            ode15s relative tolerance
%               params.sim.atol            ode15s absolute tolerance
%               params.ic.V                14×1 initial condition vector
%
% OUTPUTS:
%   sim     - struct with fields:
%               .t   time vector (last nCyclesKeep cycles)   [s]
%               .V   state matrix, n×14                       (units per state)
%
% ASSUMPTIONS:
%   - ode15s (stiff solver) is appropriate for the RLC ODE; the large
%     closed-valve resistance makes the system mildly stiff.
%   - Periodic steady-state is reached after nCyclesSteady cycles.
%     Increase this value if warnings about periodicity appear.
%
% REFERENCES:
%   [1] MATLAB ode15s documentation.
%   [2] system_rhs.m — ODE right-hand side.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

T_HB  = 60 / params.HR;           % cardiac cycle period  [s]
t_end = params.sim.nCyclesSteady * T_HB;

odefun = @(t_ode, X) system_rhs(t_ode, X, params);
opts   = odeset('RelTol', params.sim.rtol, 'AbsTol', params.sim.atol);

[t_full, V_full] = ode15s(odefun, [0, t_end], params.ic.V(:), opts);

%% Trim to last nCyclesKeep cycles
keep_duration = params.sim.nCyclesKeep * T_HB;
t_start_keep  = t_end - keep_duration;
idx_keep      = t_full >= t_start_keep;

sim.t = t_full(idx_keep)  - t_start_keep;   % re-zero time axis
sim.V = V_full(idx_keep, :);

end
