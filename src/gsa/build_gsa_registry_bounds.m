function bounds = build_gsa_registry_bounds(params0, scenario, names)
% BUILD_GSA_REGISTRY_BOUNDS
% -----------------------------------------------------------------------
% Builds GSA sampling bounds from the central parameter registry.
% Falls back per-parameter only when registry metadata is unavailable.
% -----------------------------------------------------------------------

names = names(:);
n = numel(names);
x0 = nan(n, 1);
lb = nan(n, 1);
ub = nan(n, 1);
source = strings(n, 1);
note = strings(n, 1);

registry = table();
try
    registry = build_parameter_registry(params0, params0, params0, ...
        scenario, struct(), names);
catch ME
    registry_warning = ME.message;
end

for idx = 1:n
    name = names{idx};
    x0(idx) = get_param_by_name(params0, name);
    row = [];
    if ~isempty(registry)
        row = find(strcmp(registry.name, name), 1);
    end
    if ~isempty(row) && isfinite(registry.lb(row)) && isfinite(registry.ub(row)) && ...
            registry.lb(row) < registry.ub(row)
        lb(idx) = registry.lb(row);
        ub(idx) = registry.ub(row);
        source(idx) = "registry";
        note(idx) = string(registry.notes{row});
    else
        [lb(idx), ub(idx)] = fallback_bounds(name, x0(idx));
        source(idx) = "fallback_multiplier";
        if exist('registry_warning', 'var')
            note(idx) = string(registry_warning);
        else
            note(idx) = "Registry row unavailable; used legacy multiplier.";
        end
    end
    if ~(isfinite(lb(idx)) && isfinite(ub(idx)) && lb(idx) < ub(idx))
        [lb(idx), ub(idx)] = fallback_bounds(name, x0(idx));
        source(idx) = "fallback_multiplier";
        note(idx) = "Registry bounds invalid; used legacy multiplier.";
    end
end

bounds = struct();
bounds.names = names;
bounds.x0 = x0;
bounds.lb = lb;
bounds.ub = ub;
bounds.policy = 'registry_backed_with_per_parameter_fallback';
bounds.table = table(names, x0, lb, ub, source, note, ...
    'VariableNames', {'name','x0','lb','ub','source','note'});
end

function [lb, ub] = fallback_bounds(name, x0)
if startsWith(name, 'R.')
    if strcmp(name, 'R.vsd')
        lb = max(0.001, 0.05 * x0);
        ub = min(500, 20.0 * x0);
    else
        lb = 0.4 * x0;
        ub = 2.5 * x0;
    end
elseif strcmp(name, 'vsd.Cd')
    lb = 0.8 * x0;
    ub = 1.2 * x0;
elseif startsWith(name, 'C.') || startsWith(name, 'E.')
    lb = 0.5 * x0;
    ub = 2.0 * x0;
elseif startsWith(name, 'V0.')
    lb = 0.6 * x0;
    ub = 1.7 * x0;
else
    lb = 0.7 * x0;
    ub = 1.3 * x0;
end
end

function v = get_param_by_name(params, name)
parts = strsplit(name, '.');
v = params;
for k = 1:numel(parts)
    v = v.(parts{k});
end
end

