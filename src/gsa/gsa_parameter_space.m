function space = gsa_parameter_space(params0, scenario)
% GSA_PARAMETER_SPACE
% -----------------------------------------------------------------------
% Shared uncertain-parameter definitions for PCE and direct Sobol GSA.
%
% Keeping names and bounds in one helper prevents the PCE and direct Sobol
% paths from drifting apart when a new physiology parameter is added.
%
% INPUTS:
%   params0  - scaled parameter struct used as the nominal point
%   scenario - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   space    - struct with .names, .x0, .lb, .ub
% -----------------------------------------------------------------------

names = {
    'R.SAR'
    'R.SC'
    'R.SVEN'
    'R.PAR'
    'R.PCOX'
    'R.PVEN'
    'C.SAR'
    'C.SVEN'
    'C.PAR'
    'C.PVEN'
    'E.LV.EA'
    'E.LV.EB'
    'E.RV.EA'
    'E.RV.EB'
    'E.LA.EA'
    'E.RA.EA'
    'V0.LV'
    'V0.RV'
    'V0.LA'
    'V0.RA'
};

if strcmp(scenario, 'pre_surgery')
    if isfield(params0, 'vsd') && isfield(params0.vsd, 'mode') && ...
            strcmpi(params0.vsd.mode, 'orifice_bidirectional')
        names{end+1} = 'vsd.Cd';
    else
        names{end+1} = 'R.vsd';
    end
end

d = numel(names);
x0 = zeros(d, 1);
lb = zeros(d, 1);
ub = zeros(d, 1);

for i = 1:d
    nm    = names{i};
    x0(i) = get_param_by_name(params0, nm);

    if startsWith(nm, 'R.')
        if strcmp(nm, 'R.vsd')
            lb(i) = max(0.001, 0.05 * x0(i));
            ub(i) = min(500,   20.0 * x0(i));
        else
            lb(i) = 0.4 * x0(i);
            ub(i) = 2.5 * x0(i);
        end
    elseif strcmp(nm, 'vsd.Cd')
        lb(i) = 0.8 * x0(i);
        ub(i) = 1.2 * x0(i);
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

space = struct('names', {names}, 'x0', x0, 'lb', lb, 'ub', ub);
end

function v = get_param_by_name(params, name)
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts)
    v = v.(parts{k});
end
end
