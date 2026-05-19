function validate_bounds(parameter_registry, scenario)
% VALIDATE_BOUNDS
% -----------------------------------------------------------------------
% Sanity-checks registry-driven calibration bounds before optimisation.
%
% INPUTS:
%   parameter_registry - registry table from build_parameter_registry.m   [-]
%   scenario           - 'pre_surgery' | 'post_surgery'                  [-]
%
% OUTPUTS:
%   none; throws if the registry is internally inconsistent.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(scenario)
    scenario = 'pre_surgery';
end
if nargin < 1 || ~istable(parameter_registry)
    error('validate_bounds:invalidRegistry', ...
        'parameter_registry must be a table.');
end

required_vars = {'name','baseline_scaled','seeded_value','bound_anchor', ...
    'is_calibratable','lb','ub'};
missing_vars = required_vars(~ismember(required_vars, parameter_registry.Properties.VariableNames));
if ~isempty(missing_vars)
    error('validate_bounds:missingColumns', ...
        'Registry is missing required columns: %s', strjoin(missing_vars, ', '));
end

for row_idx = 1:height(parameter_registry)
    row = parameter_registry(row_idx, :);
    name = row.name{1};
    lb = row.lb;
    ub = row.ub;
    anchor = row.bound_anchor;
    seeded_value = row.seeded_value;
    baseline_scaled = row.baseline_scaled;

    if ~isfinite(lb) || ~isfinite(ub)
        error('validate_bounds:nonFiniteBounds', ...
            'Non-finite bounds for %s.', name);
    end
    if lb >= ub
        error('validate_bounds:invalidOrder', ...
            'Lower bound must be below upper bound for %s.', name);
    end

    if row.is_calibratable
        if ~isfinite(anchor) || ~(lb < anchor && anchor < ub)
            error('validate_bounds:anchorOutsideBounds', ...
                'Registry anchor must lie inside bounds for calibratable %s.', name);
        end
        if ~isfinite(seeded_value) || ~(lb <= seeded_value && seeded_value <= ub)
            error('validate_bounds:seedOutsideBounds', ...
                'Seeded value must lie inside bounds for calibratable %s.', name);
        end
    end

    if startsWith(name, 'R.') || startsWith(name, 'C.') || startsWith(name, 'L.') || ...
            startsWith(name, 'E.') || strcmp(name, 'group.R_sys_scale') || ...
            strcmp(name, 'group.R_pul_scale') || ...
            strcmp(name, 'vsd.Cd')
        if lb <= 0
            error('validate_bounds:nonPositiveLowerBound', ...
                'Lower bound must stay positive for %s.', name);
        end
    end

    if strcmp(name, 'R.vsd') && strcmp(scenario, 'post_surgery') && row.is_calibratable
        error('validate_bounds:postClosureVsdActive', ...
            'R.vsd must not be calibratable in post_surgery.');
    end

    is_pulse_pressure_anchor = ismember('bound_type', parameter_registry.Properties.VariableNames) && ...
        strcmp(row.bound_type{1}, 'pulse_pressure_anchor');

    is_patient_specific_vsd_anchor = ...
        (strcmp(name, 'R.vsd') && strcmp(scenario, 'pre_surgery')) || ...
        (strcmp(name, 'vsd.Cd') && strcmp(scenario, 'pre_surgery'));

    if row.is_calibratable && ~isnan(baseline_scaled) && ...
            ~is_pulse_pressure_anchor && ...
            ~is_patient_specific_vsd_anchor && ...
            ~(lb <= baseline_scaled && baseline_scaled <= ub)
        error('validate_bounds:scaledBaselineOutsideBounds', ...
            'Scaled baseline must lie inside bounds for %s.', name);
    end
end
end
