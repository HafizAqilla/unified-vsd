function T = make_gsa_summary_table(gsa_out, top_n)
% MAKE_GSA_SUMMARY_TABLE
% -----------------------------------------------------------------------
% Produce a compact summary table showing the top-N most influential
% parameters for each output metric in the GSA result.
%
% INPUTS:
%   gsa_out  - struct from gsa_run_sobol.m
%   top_n    - number of top parameters to show per metric  (default: 5)
%
% OUTPUTS:
%   T        - MATLAB table with columns:
%              Metric | Primary | Parameter | Sobol_S1 | Sobol_ST
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-02-26
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2, top_n = 5; end

metrics = fieldnames(gsa_out);
metrics = metrics(~ismember(metrics, {'scenario', 'cfg'}));

rows = {};
for m = 1:numel(metrics)
    mf  = metrics{m};
    res = gsa_out.(mf);
    if ~isstruct(res) || ~isfield(res, 'table'), continue; end

    tbl   = res.table;
    n_top = min(top_n, height(tbl));
    for r = 1:n_top
        rows(end+1, :) = { ...
            mf, ...
            res.primary, ...
            tbl.Parameter{r}, ...
            round(tbl.Sobol_S1(r), 4), ...
            round(tbl.Sobol_ST(r), 4) ...
            }; %#ok<AGROW>
    end
end

if isempty(rows)
    T = table();
    return;
end

T = cell2table(rows, 'VariableNames', ...
    {'Metric', 'IsPrimary', 'Parameter', 'Sobol_S1', 'Sobol_ST'});

%% Print highlighted primary-metric rows
fprintf('\n=== Sobol GSA Summary — Scenario: %s ===\n', gsa_out.scenario);
fprintf('    (showing top %d parameters per metric)\n\n', top_n);

primary_rows = T(logical([T.IsPrimary{:}]'), :);
fprintf('--- PRIMARY METRICS ---\n');
disp(primary_rows(:, {'Metric','Parameter','Sobol_S1','Sobol_ST'}));

secondary_rows = T(~logical([T.IsPrimary{:}]'), :);
if ~isempty(secondary_rows)
    fprintf('--- SECONDARY METRICS ---\n');
    disp(secondary_rows(:, {'Metric','Parameter','Sobol_S1','Sobol_ST'}));
end

end
