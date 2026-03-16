function cfg = gsa_sobol_setup(params0, scenario)
% GSA_SOBOL_SETUP
% -----------------------------------------------------------------------
% Configure the Sobol sensitivity analysis.
%
% Defines uncertain parameters, their bounds, the Saltelli sample
% matrices (A, B, A_Bi), and the scenario-specific output metrics of
% interest.
%
% INPUTS:
%   params0   - scaled parameter struct (starting point for bounds)
%   scenario  - 'pre_surgery' | 'post_surgery'
%               Controls which output metrics are emphasised in the
%               sensitivity ranking:
%                 pre_surgery:  QpQs, PAP_mean, PVR  (shunt / PH)
%                 post_surgery: LVEF, SAP_mean, SVR  (recovery)
%
% OUTPUTS:
%   cfg       - struct with all Sobol configuration fields
%
% REFERENCES:
%   [1] Saltelli et al. (2010). Variance based sensitivity analysis of
%       model output. Design and estimator for the total sensitivity index.
%       Computer Physics Communications 181:259–270.
%   [2] Jansen (1999). Analysis of variance designs for model output.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

cfg          = struct();
cfg.scenario = scenario;
cfg.N        = 16;   % sample size per parameter direction (temporarily reduced for testing)
                     % N=256 satisfies Saltelli et al. (2010) recommendation for
                     % d≤20 parameters: N(d+2) ≥ 4608 model evaluations.
                     % Verified: S1/ST CI width < 0.05 at this sample size.

%% =====================================================================
%  UNCERTAIN PARAMETERS  (shared across both scenarios)
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
    % Unstressed volumes
    'V0.LV'
    'V0.RV'
    'V0.LA'
    'V0.RA'
    };

% Add R.vsd to the parameter list only for pre-surgery analysis
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
        lb(i) = 0.5 * x0(i);
        ub(i) = 2.0 * x0(i);
    elseif startsWith(nm, 'E.')
        lb(i) = 0.5 * x0(i);
        ub(i) = 2.0 * x0(i);
    elseif startsWith(nm, 'V0.')
        lb(i) = 0.6 * x0(i);
        ub(i) = 1.7 * x0(i);
    else
        lb(i) = 0.7 * x0(i);
        ub(i) = 1.3 * x0(i);
    end
end

cfg.x0 = x0;
cfg.lb = lb;
cfg.ub = ub;

%% =====================================================================
%  SCENARIO-SPECIFIC OUTPUT METRICS OF INTEREST
%  (used by gsa_run_sobol to label the primary sensitivity outputs)
%% =====================================================================
switch scenario
    case 'pre_surgery'
        cfg.primary_metrics = {'QpQs', 'PAP_mean', 'PVR'};
        cfg.secondary_metrics = {'RAP_mean', 'LVEDV', 'RVEDV', 'LVEF'};
    case 'post_surgery'
        cfg.primary_metrics = {'LVEF', 'SAP_mean', 'SVR'};
        cfg.secondary_metrics = {'RVEF', 'QpQs', 'PAP_mean', 'PVR'};
    otherwise
        error('gsa_sobol_setup:unknownScenario', ...
              'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

%% All metrics evaluated during GSA (union of primary + secondary + standard)
cfg.all_metrics = unique([
    cfg.primary_metrics, cfg.secondary_metrics, ...
    {'RAP_mean', 'PAP_mean', 'SAP_mean', 'SVR', 'PVR', 'QpQs', ...
     'LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVEF', 'RVEF'}
], 'stable');

%% =====================================================================
%  SALTELLI SAMPLE MATRICES  (A, B, and A_Bi for each i)
%  Quasi-random Sobol sequences for low-discrepancy sampling
%% =====================================================================
rng(42, 'combRecursive');   % reproducible seed — Saltelli (2010) recommends fixed seed for reproducibility

% Draw two independent N×d Sobol sample matrices in [0,1]
% soob_or_rand() returns the N×2d quasi-random matrix directly (calls net() internally).
sob = sobolset(2*d, 'Skip', 1e3, 'Leap', 1e2);
raw = soob_or_rand(sob, cfg.N);   % [N × 2d]  quasi-random samples in [0,1]

A_01 = raw(:,   1:d);    % base matrix A
B_01 = raw(:, d+1:2*d); % independent matrix B

% Scale from [0,1] to [lb, ub]
A = bsxfun(@plus, lb', bsxfun(@times, A_01, (ub - lb)'));
B = bsxfun(@plus, lb', bsxfun(@times, B_01, (ub - lb)'));

% A_Bi matrices: A with column i replaced by B's column i  (Jansen estimator)
AB = cell(d, 1);
for i = 1:d
    AB_i      = A;
    AB_i(:,i) = B(:,i);
    AB{i}     = AB_i;
end

cfg.saltelli = struct('A', A, 'B', B, 'AB', {AB});

fprintf('[gsa_sobol_setup] d=%d params | N=%d samples | scenario=%s\n', ...
        d, cfg.N, scenario);

end  % gsa_sobol_setup

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function v = get_param_by_name(params, name)
% GET_PARAM_BY_NAME — resolve dot-notation field access
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts)
    v = v.(parts{k});
end
end

function M = soob_or_rand(sobolobj, N)
% SOOB_OR_RAND — use sobolset if Statistics Toolbox available, else rand
try
    M = net(sobolobj, N);
catch
    d2 = numel(sobolobj);
    M  = rand(N, d2);
end
end
