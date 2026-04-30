function cfg = gsa_pce_setup(params0, scenario, uqlab_path, sobios_path)
% GSA_PCE_SETUP
% -----------------------------------------------------------------------
% Configure the PCE-based Global Sensitivity Analysis using the
% SoBioS toolbox pattern (Tosin, Côrtes, Cunha Jr, 2020) and UQLab.
%
% Follows the structure of SoBioS_PCE_QoI_spike_nfkb_7vars.m from the
% SoBioS-1.0 package. UQLab is called as the PCE engine.
%
% INPUTS:
%   params0       - scaled parameter struct (baseline operating point)
%   scenario      - 'pre_surgery' | 'post_surgery'
%   uqlab_path    - (optional) full path to the UQLab core folder
%   sobios_path   - (optional) full path to the SoBioS folder
%
% OUTPUTS:
%   cfg           - struct with fields:
%                     .names              parameter names  (d × 1 cell)
%                     .lb / .ub / .x0    bounds / nominal (d × 1)
%                     .primary_metrics / .secondary_metrics / .all_metrics
%                     .params0            baseline param struct (for wrapper)
%                     .scenario           scenario string
%                     .Input              UQLab Input object
%                     .PCEOpts            UQLab PCE metamodel options
%
% REFERENCES:
%   Tosin M., Côrtes A.M.A., Cunha Jr A. (2020). A Tutorial on Sobol'
%   Global Sensitivity Analysis Applied to Biological Models. Springer.
%   https://doi.org/10.1007/978-3-030-51862-2_6
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-03-10
% VERSION:  2.0  (SoBioS-pattern rewrite)
% -----------------------------------------------------------------------

%% --- Path & UQLab init --------------------------------------------------
if nargin >= 3 && ~isempty(uqlab_path)
    addpath(uqlab_path);
end
if nargin >= 4 && ~isempty(sobios_path)
    % Add SoBioS (for reference / future QoI wrappers)
    addpath(genpath(sobios_path));
end

% Call uqlab to initialise the framework.
uqlab('-nosplash');

cfg          = struct();
cfg.scenario = scenario;
cfg.params0  = params0;   % stored so gsa_run_pce wrapper can access it

%% =====================================================================
%  UNCERTAIN PARAMETERS  (same list as gsa_sobol_setup for comparability)
%% =====================================================================
cfg.names = {
    % Resistances
    'R.SAR'
    'R.SC'
    'R.SVEN'
    'R.PAR'
    'R.PCOX'
    'R.PVEN'
    % Compliances
    'C.SAR'
    'C.SVEN'
    'C.PAR'
    'C.PVEN'
    % Ventricular elastances
    'E.LV.EA'
    'E.LV.EB'
    'E.RV.EA'
    'E.RV.EB'
    % Atrial elastances (active calibration params — must appear in GSA for mask)
    'E.LA.EA'
    'E.RA.EA'
    % Unstressed volumes
    'V0.LV'
    'V0.RV'
    'V0.LA'
    'V0.RA'
};

if strcmp(scenario, 'pre_surgery')
    cfg.names{end+1} = 'R.vsd';
end

d = numel(cfg.names);

%% =====================================================================
%  NOMINAL VALUES AND BOUNDS
%% =====================================================================
x0 = zeros(d, 1);
lb = zeros(d, 1);
ub = zeros(d, 1);

for i = 1:d
    nm    = cfg.names{i};
    x0(i) = get_param_by_name(params0, nm);

    if startsWith(nm, 'R.')
        if strcmp(nm, 'R.vsd')
            lb(i) = max(0.001, 0.05 * x0(i));
            ub(i) = min(500,   20.0 * x0(i));
        else
            lb(i) = 0.4 * x0(i);
            ub(i) = 2.5 * x0(i);
        end
    elseif startsWith(nm, 'C.')
        lb(i) = 0.5 * x0(i);  ub(i) = 2.0 * x0(i);
    elseif startsWith(nm, 'E.')
        lb(i) = 0.5 * x0(i);  ub(i) = 2.0 * x0(i);
    elseif startsWith(nm, 'V0.')
        lb(i) = 0.6 * x0(i);  ub(i) = 1.7 * x0(i);
    else
        lb(i) = 0.7 * x0(i);  ub(i) = 1.3 * x0(i);
    end
end
cfg.x0 = x0;  cfg.lb = lb;  cfg.ub = ub;

%% =====================================================================
%  SCENARIO-SPECIFIC METRICS
%% =====================================================================
switch scenario
    case 'pre_surgery'
        % Primary metrics drive the GSA mask — parameters with high ST
        % for ANY primary metric survive the Sobol threshold cut.
        % LVEF and SAP_mean added so E.LV.EA and C.SAR are never masked out.
        cfg.primary_metrics   = {'QpQs', 'PAP_mean', 'LVEF', 'SAP_mean', 'CO_Lmin'};
        cfg.secondary_metrics = {'RAP_mean', 'LVEDV', 'LVESV', 'RVEDV', 'RVESV', ...
                                 'RVP_min', 'RVP_max', 'RVP_mean', ...
                                 'PAP_min', 'PAP_max', ...
                                 'SAP_min', 'SAP_max', ...
                                 'PVP_mean', 'SVR', ...
                                 'LVP_min', 'LVP_max', 'LVP_mean'};
    case 'post_surgery'
        cfg.primary_metrics   = {'LVEF', 'SAP_mean', 'SVR'};
        cfg.secondary_metrics = {'RVEF', 'QpQs', 'PAP_mean', 'PVR', ...
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
     'PAP_min',  'PAP_max',  'PAP_mean', ...
     'PVP_mean', ...
     'RVP_min',  'RVP_max',  'RVP_mean', ...
     'LVP_min',  'LVP_max',  'LVP_mean', ...
     'SAP_min',  'SAP_max',  'SAP_mean', ...
     'SVR', 'PVR', 'QpQs', ...
     'LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVEF', 'RVEF'}
], 'stable');

%% =====================================================================
%  UQLAB PROBABILISTIC INPUT  (Uniform marginals — SoBioS pattern)
%  Matches InputOpts structure from SoBioS_PCE_QoI_spike_nfkb_7vars.m
%% =====================================================================
InputOpts = struct();
for i = 1:d
    InputOpts.Marginals(i).Name       = strrep(cfg.names{i}, '.', '_');
    InputOpts.Marginals(i).Type       = 'Uniform';
    InputOpts.Marginals(i).Parameters = [lb(i), ub(i)];
end
cfg.Input = uq_createInput(InputOpts);

%% =====================================================================
%  UQLAB PCE METAMODEL OPTIONS  (LARS — SoBioS pattern)
%  PCEOpts.FullModel is NOT set here: set in gsa_run_pce per metric.
%% =====================================================================
PCEOpts.Type                   = 'Metamodel';
PCEOpts.MetaType               = 'PCE';
PCEOpts.Method                 = 'LARS';
PCEOpts.Degree                 = 1:3;           % reduced from 1:4 → faster LARS sweep
PCEOpts.TruncOptions.qNorm     = 0.75;          % single value (was 0.5:0.1:1.0)
PCEOpts.ExpDesign.NSamples     = 100;           % reduced from 200 → 2x ODE speedup
PCEOpts.ExpDesign.Sampling     = 'Halton';
PCEOpts.Input                  = cfg.Input;

cfg.PCEOpts = PCEOpts;

% Reduced-fidelity ODE overrides for GSA sampling (same as gsa_sobol_setup).
% Without this, each of the 100 ODE runs uses the default 80 warm-up cycles.
% 10 cycles at relaxed tolerances is sufficient for Sobol index screening.
cfg.gsa_sim_overrides.nCyclesSteady = 10;   % [cycles]  default 80 → 10
cfg.gsa_sim_overrides.ss_tol_P      = 1.0;  % [mmHg]    relaxed from 0.1
cfg.gsa_sim_overrides.ss_tol_V      = 1.0;  % [mL]      relaxed from 0.1

fprintf('[gsa_pce_setup] d=%d params | N_train=%d | scenario=%s\n', ...
        d, PCEOpts.ExpDesign.NSamples, scenario);

end  % gsa_pce_setup


% =========================================================================
%  LOCAL HELPER
% =========================================================================
function v = get_param_by_name(params, name)
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts), v = v.(parts{k}); end
end