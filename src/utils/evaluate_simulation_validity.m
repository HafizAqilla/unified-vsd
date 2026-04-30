function validity = evaluate_simulation_validity(sim, params, metrics, scenario, clinical)
% EVALUATE_SIMULATION_VALIDITY
% -----------------------------------------------------------------------
% Applies hard physiological validity gates to a completed simulation.
%
% INPUTS:
%   sim       - simulation struct from integrate_system                 [-]
%   params    - parameter struct                                        [-]
%   metrics   - clinical metric struct                                  [-]
%   scenario  - scenario string                                         [-]
%   clinical  - optional clinical struct for scenario-aware bounds      [-]
%
% OUTPUTS:
%   validity  - struct with flags, summary, and scalar penalty          [-]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 5
    clinical = [];
end

XV = sim.V;
sidx = params.idx;

vol_idx = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV ...
           sidx.V_SAR sidx.V_SC sidx.V_SVEN sidx.V_PAR sidx.V_PVEN];
chamber_idx = [sidx.V_RA sidx.V_RV sidx.V_LA sidx.V_LV];

flags = struct();
flags.steady_state_not_reached = ~isfield(sim, 'ss_reached') || ~sim.ss_reached;
flags.nonfinite_state = any(~isfinite(XV(:)));
flags.negative_chamber_volume = any(any(XV(:, chamber_idx) < 0));
flags.negative_total_volume_state = any(any(XV(:, vol_idx) < 0));
flags.unrealistic_v0_adjusted_volume = any(any( ...
    XV(:, sidx.V_LV) - params.V0.LV < -1 | ...
    XV(:, sidx.V_RV) - params.V0.RV < -1));

flags.flow_collapse = false;
if isfield(metrics, 'Qs_mean_mLs') && isfield(metrics, 'QpQs')
    flags.flow_collapse = isfinite(metrics.QpQs) && abs(metrics.QpQs) > 0 && abs(metrics.Qs_mean_mLs) < 1e-3;
end

clinical_targets = struct();
if ~isempty(clinical) && isstruct(clinical)
    try
        clinical_targets = targets_to_struct(get_calibration_targets(scenario, clinical));
    catch
        clinical_targets = struct();
    end
end

[pvr_lb, pvr_ub] = metric_bounds('PVR', 0.1, 20, clinical_targets);
[svr_lb, svr_ub] = metric_bounds('SVR', 0.5, 40, clinical_targets);
[lvef_lb, lvef_ub] = metric_bounds('LVEF', 0.05, 0.90, clinical_targets);
[rvef_lb, rvef_ub] = metric_bounds('RVEF', 0.05, 0.90, clinical_targets);

flags.pvr_out_of_bounds = metric_out_of_bounds(metrics, 'PVR', pvr_lb, pvr_ub);
flags.svr_out_of_bounds = metric_out_of_bounds(metrics, 'SVR', svr_lb, svr_ub);
flags.ef_out_of_bounds = any([
    metric_out_of_bounds(metrics, 'LVEF', lvef_lb, lvef_ub), ...
    metric_out_of_bounds(metrics, 'RVEF', rvef_lb, rvef_ub)]);

flags.vsd_geometry_invalid = false;
if isfield(params, 'vsd') && isfield(params.vsd, 'area_mm2')
    flags.vsd_geometry_invalid = params.vsd.area_mm2 < 0 || params.vsd.area_mm2 > 400;
end
if isfield(params, 'R') && isfield(params.R, 'vsd')
    flags.vsd_geometry_invalid = flags.vsd_geometry_invalid || params.R.vsd < 0;
end

flags.reverse_shunt_missing = false;
if strcmp(scenario, 'pre_surgery') && isfield(params, 'vsd') && isfield(params.vsd, 'mode')
    if contains(lower(params.vsd.mode), 'bidirectional')
        [P, Q] = reconstruct_signals_validity(sim, params);
        reverse_possible = any(P.RV > P.LV);
        reverse_present = any(Q.VSD < 0);
        flags.reverse_shunt_missing = reverse_possible && ~reverse_present;
    end
end

flag_names = fieldnames(flags);
flag_values = struct2cell(flags);
flag_values = cellfun(@logical, flag_values);

validity.flags = flags;
validity.failed_flags = flag_names(flag_values);
validity.is_valid = ~any(flag_values);
validity.penalty = 0;
if ~validity.is_valid
    validity.penalty = 1e3 * sum(flag_values);
end
validity.metric_bounds = struct( ...
    'PVR', [pvr_lb, pvr_ub], ...
    'SVR', [svr_lb, svr_ub], ...
    'LVEF', [lvef_lb, lvef_ub], ...
    'RVEF', [rvef_lb, rvef_ub]);
validity.metric_values = struct( ...
    'PVR', field_or_nan(metrics, 'PVR'), ...
    'SVR', field_or_nan(metrics, 'SVR'), ...
    'LVEF', field_or_nan(metrics, 'LVEF'), ...
    'RVEF', field_or_nan(metrics, 'RVEF'));
validity.bound_summary = make_bound_summary(validity.metric_values, validity.metric_bounds);

end

function [P, Q] = reconstruct_signals_validity(sim, params)
[P, Q] = reconstruct_hemodynamic_signals(sim.t(:), sim.V, params);
end

function tf = metric_out_of_bounds(metrics, field_name, lb, ub)
tf = false;
if ~isfield(metrics, field_name) || ~isfinite(metrics.(field_name))
    return;
end
tf = metrics.(field_name) < lb || metrics.(field_name) > ub;
end

function [lb, ub] = metric_bounds(metric_name, default_lb, default_ub, clinical_targets)
lb = default_lb;
ub = default_ub;
if ~isfield(clinical_targets, metric_name)
    return;
end

target = clinical_targets.(metric_name);
if ~isfinite(target)
    return;
end

switch metric_name
    case {'PVR', 'SVR'}
        lb = min(default_lb, 0.25 * target);
        ub = max(default_ub, 2.50 * target);
    case {'LVEF', 'RVEF'}
        lb = max(default_lb, target - 0.35);
        ub = min(default_ub, target + 0.35);
end
end

function targets_struct = targets_to_struct(targets)
targets_struct = struct();
for i = 1:numel(targets)
    targets_struct.(targets(i).Metric) = targets(i).ClinicalValue;
end
end

function value = field_or_nan(data, field_name)
value = NaN;
if isfield(data, field_name)
    value = data.(field_name);
end
end

function summary = make_bound_summary(metric_values, metric_bounds)
metric_names = fieldnames(metric_bounds);
n_metrics = numel(metric_names);
summary = table('Size', [n_metrics, 4], ...
    'VariableTypes', {'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Metric', 'Value', 'LowerBound', 'UpperBound'});
for i = 1:n_metrics
    name = metric_names{i};
    bounds = metric_bounds.(name);
    summary.Metric(i) = string(name);
    summary.Value(i) = metric_values.(name);
    summary.LowerBound(i) = bounds(1);
    summary.UpperBound(i) = bounds(2);
end
end
