function cfg = gsa_pce_setup(params0, scenario, uqlab_path, sobios_path)
% GSA_PCE_SETUP
% -----------------------------------------------------------------------
% Configure the PCE-based Global Sensitivity Analysis using the SoBioS
% toolbox pattern and UQLab.
%
% INPUTS:
%   params0     - scaled parameter struct (baseline operating point)
%   scenario    - 'pre_surgery' | 'post_surgery'
%   uqlab_path  - optional full path to the UQLab core folder
%   sobios_path - optional full path to the SoBioS folder
%
% ENVIRONMENT:
%   UNIFIED_VSD_GSA_PCE_N - optional positive integer training sample count
%
% OUTPUTS:
%   cfg         - struct with GSA parameters, bounds, metrics, UQLab input,
%                 PCE options, and reduced-fidelity simulation overrides
%
% REFERENCES:
%   Tosin M., Cortes A.M.A., Cunha Jr A. (2020). A Tutorial on Sobol'
%   Global Sensitivity Analysis Applied to Biological Models. Springer.
%   https://doi.org/10.1007/978-3-030-51862-2_6
% -----------------------------------------------------------------------

%% Path and UQLab init
if nargin >= 3 && ~isempty(uqlab_path)
    addpath(uqlab_path);
end
if nargin >= 4 && ~isempty(sobios_path)
    addpath(genpath(sobios_path));
end

uqlab('-nosplash');

cfg = struct();
cfg.scenario = scenario;
cfg.params0 = params0;

%% Shared uncertain parameters and bounds
space = gsa_parameter_space(params0, scenario);
cfg.names = space.names;
cfg.x0 = space.x0;
cfg.lb = space.lb;
cfg.ub = space.ub;

d = numel(cfg.names);
lb = cfg.lb;
ub = cfg.ub;

%% Scenario-specific metrics
switch scenario
    case 'pre_surgery'
        cfg.primary_metrics = {'QpQs', 'PAP_mean', 'PVR', 'SAP_mean', 'CO_Lmin'};
        cfg.secondary_metrics = {'RAP_mean', 'LVEDV', 'LVESV', 'RVEDV', 'RVESV', ...
            'RVP_min', 'RVP_max', 'RVP_mean', ...
            'PAP_min', 'PAP_max', ...
            'SAP_min', 'SAP_max', ...
            'PVP_mean', 'SVR', 'LVEF', 'RVEF', ...
            'LVP_min', 'LVP_max', 'LVP_mean'};
    case 'post_surgery'
        cfg.primary_metrics = {'QpQs', 'PAP_mean', 'PVR', 'SAP_mean', 'CO_Lmin'};
        cfg.secondary_metrics = {'LVEF', 'RVEF', 'SVR', ...
            'RVP_min', 'RVP_max', 'RVP_mean', ...
            'PAP_min', 'PAP_max', ...
            'PVP_mean', ...
            'LVP_min', 'LVP_max', 'LVP_mean'};
    otherwise
        error('gsa_pce_setup:unknownScenario', ...
            'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

cfg.all_metrics = unique([
    cfg.primary_metrics, cfg.secondary_metrics, ...
    {'RAP_mean', 'LAP_mean', ...
     'PAP_min', 'PAP_max', 'PAP_mean', ...
     'PVP_mean', ...
     'RVP_min', 'RVP_max', 'RVP_mean', ...
     'LVP_min', 'LVP_max', 'LVP_mean', ...
     'SAP_min', 'SAP_max', 'SAP_mean', ...
     'SVR', 'PVR', 'QpQs', ...
     'LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVEF', 'RVEF'}
], 'stable');

%% UQLab probabilistic input
InputOpts = struct();
for i = 1:d
    InputOpts.Marginals(i).Name = strrep(cfg.names{i}, '.', '_');
    InputOpts.Marginals(i).Type = 'Uniform';
    InputOpts.Marginals(i).Parameters = [lb(i), ub(i)];
end
cfg.Input = uq_createInput(InputOpts);

%% UQLab PCE metamodel options
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
PCEOpts.Method = 'LARS';
PCEOpts.Degree = 1:3;
PCEOpts.TruncOptions.qNorm = 0.75;

default_N_train = 128;
env_N_train = getenv('UNIFIED_VSD_GSA_PCE_N');
if ~isempty(env_N_train)
    env_N_train_num = str2double(env_N_train);
    if ~isnan(env_N_train_num) && isfinite(env_N_train_num) && env_N_train_num > 0
        default_N_train = max(1, round(env_N_train_num));
    else
        fprintf(2, ['[gsa_pce_setup] Ignoring invalid UNIFIED_VSD_GSA_PCE_N=''%s''; ', ...
            'using N_train=%d.\n'], env_N_train, default_N_train);
    end
end

PCEOpts.ExpDesign.NSamples = default_N_train;
PCEOpts.ExpDesign.Sampling = 'Halton';
PCEOpts.Input = cfg.Input;
cfg.PCEOpts = PCEOpts;

% Reduced-fidelity ODE overrides for GSA sampling.
cfg.gsa_sim_overrides.nCyclesSteady = 10;
cfg.gsa_sim_overrides.ss_tol_P = 1.0;
cfg.gsa_sim_overrides.ss_tol_V = 1.0;

fprintf('[gsa_pce_setup] d=%d params | N_train=%d | scenario=%s\n', ...
    d, PCEOpts.ExpDesign.NSamples, scenario);

end
