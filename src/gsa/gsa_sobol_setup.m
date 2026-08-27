function cfg = gsa_sobol_setup(params0, scenario, N_override, registry_context)
% GSA_SOBOL_SETUP
% -----------------------------------------------------------------------
% Configure the direct Saltelli/Jansen Sobol sensitivity analysis.
%
% This direct Sobol path uses the same uncertain-parameter space as the
% PCE-GSA path through gsa_parameter_space().
%
% INPUTS:
%   params0    - scaled parameter struct used as nominal point for bounds
%   scenario   - 'pre_surgery' | 'post_surgery'
%   N_override - optional explicit Sobol base sample count N
%
% ENVIRONMENT:
%   GSA_SOBOL_N - optional positive integer Sobol base sample count
%
% OUTPUTS:
%   cfg        - struct with parameter bounds, sample matrices, metrics,
%                and reduced-fidelity simulation overrides
%
% REFERENCES:
%   Saltelli et al. (2010). Variance based sensitivity analysis of model
%   output. Design and estimator for the total sensitivity index.
%   Jansen (1999). Analysis of variance designs for model output.
% -----------------------------------------------------------------------

cfg = struct();
cfg.scenario = scenario;
if nargin < 4 || isempty(registry_context)
    registry_context = struct();
end

% Base Sobol sample count per Saltelli A/B block.
default_N = 256;
cfg.N = default_N;
N_source = 'default';

if nargin >= 3 && ~isempty(N_override)
    cfg.N = max(1, round(N_override));
    N_source = 'function_override';
else
    env_N = getenv('GSA_SOBOL_N');
    if ~isempty(env_N)
        env_N_num = str2double(env_N);
        if ~isnan(env_N_num) && isfinite(env_N_num) && env_N_num > 0
            cfg.N = round(env_N_num);
            N_source = 'environment';
        end
    end
end

cfg.N_source = N_source;

% GSA evaluations use reduced warmup for screening speed.
cfg.gsa_sim_overrides.nCyclesSteady = 10;
cfg.gsa_sim_overrides.ss_tol_P = 1.0;
cfg.gsa_sim_overrides.ss_tol_V = 1.0;

%% Shared uncertain parameters and bounds
space = gsa_parameter_space(params0, scenario);
cfg.names = space.names;
cfg.x0 = space.x0;
cfg.lb = space.lb;
cfg.ub = space.ub;
gsa_bounds = build_gsa_registry_bounds(params0, scenario, cfg.names, registry_context);
cfg.x0 = gsa_bounds.x0;
cfg.lb = gsa_bounds.lb;
cfg.ub = gsa_bounds.ub;
cfg.bounds_policy = gsa_bounds.policy;
cfg.bounds_table = gsa_bounds.table;
cfg.registry_context = registry_context;

d = numel(cfg.names);
lb = cfg.lb;
ub = cfg.ub;

%% Scenario-specific output metrics
switch scenario
    case 'pre_surgery'
        cfg.primary_metrics = {'RAP_mean', 'PAP_mean', 'SAP_mean', 'QpQs', 'CO_Lmin'};
        cfg.secondary_metrics = {'RAP_mean', 'LVEDV', 'RVEDV', 'LVEF', 'SVR', 'RVEF'};
    case 'post_surgery'
        cfg.primary_metrics = {'RAP_mean', 'PAP_mean', 'SAP_mean', 'QpQs', 'CO_Lmin'};
        cfg.secondary_metrics = {'LVEF', 'RVEF', 'SVR', 'LVEDV', 'RVEDV'};
    otherwise
        error('gsa_sobol_setup:unknownScenario', ...
            'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

cfg.all_metrics = unique([
    cfg.primary_metrics, cfg.secondary_metrics, ...
    {'RAP_mean', 'PAP_mean', 'SAP_mean', 'SVR', 'PVR', 'QpQs', ...
     'CO_Lmin', 'LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVEF', 'RVEF'}
], 'stable');

%% Saltelli sample matrices
rng(42, 'combRecursive');

sob = sobolset(2*d, 'Skip', 1e3, 'Leap', 1e2);
raw = soob_or_rand(sob, cfg.N);

A_01 = raw(:, 1:d);
B_01 = raw(:, d+1:2*d);

A = bsxfun(@plus, lb', bsxfun(@times, A_01, (ub - lb)'));
B = bsxfun(@plus, lb', bsxfun(@times, B_01, (ub - lb)'));

AB = cell(d, 1);
for i = 1:d
    AB_i = A;
    AB_i(:, i) = B(:, i);
    AB{i} = AB_i;
end

cfg.saltelli = struct('A', A, 'B', B, 'AB', {AB});

fprintf('[gsa_sobol_setup] d=%d params | N=%d samples | scenario=%s\n', ...
    d, cfg.N, scenario);
fprintf('[gsa_sobol_setup] N source: %s\n', cfg.N_source);

if bitand(cfg.N, cfg.N - 1) ~= 0
    fprintf(2, '[gsa_sobol_setup] Warning: N=%d is not a power of two; Sobol stability may degrade.\n', cfg.N);
end

end

function M = soob_or_rand(sobolobj, N)
try
    M = net(sobolobj, N);
catch
    d2 = numel(sobolobj);
    M = rand(N, d2);
end
end
