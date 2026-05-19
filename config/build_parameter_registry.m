function registry = build_parameter_registry(params_adult, params_scaled, params_seeded, scenario, case_profile, parameter_names)
% BUILD_PARAMETER_REGISTRY
% -----------------------------------------------------------------------
% Builds centralized metadata for calibration variables, including anchor
% values, registry-driven bounds, and evidence notes for plausibility
% checks. The registry is scenario- and case-profile-aware, but it does
% not decide the final active GSA mask.
%
% INPUTS:
%   params_adult     - adult baseline parameter struct                  [-]
%   params_scaled    - scaled pediatric baseline before clinical seeding[-]
%   params_seeded    - clinically seeded parameter struct              [-]
%   scenario         - 'pre_surgery' | 'post_surgery'                  [-]
%   case_profile     - calibration governance profile                  [-]
%   parameter_names  - calibration variable names to include           [cellstr]
%
% OUTPUTS:
%   registry         - table with parameter metadata and bounds         [-]
%
% NOTES:
%   - For grouped calibration variables, baseline_scaled is defined as a
%     dimensionless scale of 1.0 relative to the scaled pediatric anchor.
%   - For pre-surgery R.vsd, the physiologic anchor is patient-specific,
%     so baseline_scaled and seeded_value are both taken from the seeded
%     shunt estimate rather than the healthy-adult closed-shunt default.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-07
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 6 || isempty(parameter_names)
    error('build_parameter_registry:missingNames', ...
        'parameter_names must be supplied explicitly.');
end
if nargin < 5 || isempty(case_profile)
    case_profile = struct();
end
if nargin < 4 || isempty(scenario)
    scenario = 'pre_surgery';
end
if nargin < 3 || isempty(params_seeded)
    params_seeded = params_scaled;
end
if nargin < 2 || isempty(params_scaled)
    params_scaled = params_seeded;
end
if nargin < 1 || isempty(params_adult)
    params_adult = params_scaled;
end

names = parameter_names(:);
n_params = numel(names);

group_col = cell(n_params, 1);
unit_col = cell(n_params, 1);
baseline_adult_col = nan(n_params, 1);
baseline_scaled_col = nan(n_params, 1);
seeded_value_col = nan(n_params, 1);
bound_anchor_col = nan(n_params, 1);
scale_method_col = cell(n_params, 1);
is_calibratable_col = false(n_params, 1);
lb_col = nan(n_params, 1);
ub_col = nan(n_params, 1);
bound_type_col = cell(n_params, 1);
bound_source_col = cell(n_params, 1);
confidence_col = cell(n_params, 1);
notes_col = cell(n_params, 1);

for idx = 1:n_params
    name = names{idx};
    group_info = get_coupled_parameter_group(name, case_profile);

    group_col{idx} = resolve_parameter_group(name, group_info);
    unit_col{idx} = resolve_parameter_unit(name, group_info);
    is_calibratable_col(idx) = resolve_calibratable_flag(name, scenario, case_profile);

    if group_info.found
        baseline_adult = NaN;
        baseline_scaled = 1.0;
        seeded_value = get_calibration_param_value(params_seeded, params_scaled, name, case_profile);
        scale_method = sprintf('group_scale_relative_to_scaled(%s)', strjoin(group_info.memberNames, ','));
        note_parts = {sprintf('Coupled calibration scale for %s.', strjoin(group_info.memberNames, ', '))};
    else
        baseline_adult = get_calibration_param_value(params_adult, params_adult, name, case_profile);
        baseline_scaled = get_calibration_param_value(params_scaled, params_scaled, name, case_profile);
        seeded_value = get_calibration_param_value(params_seeded, params_scaled, name, case_profile);
        scale_method = resolve_scale_method(name, scenario);
        note_parts = {''};
        if strcmp(name, 'R.vsd') && strcmp(scenario, 'pre_surgery')
            baseline_adult = NaN;
            baseline_scaled = resolve_preclosure_vsd_anchor(seeded_value);
            note_parts{end + 1} = 'Healthy-adult closed-shunt baseline is not used as the plausibility anchor for pre-closure R.vsd.';
        end
    end

    [lb, ub, bound_anchor, bound_type, bound_source, confidence_level, bound_note] = ...
        resolve_bounds(name, group_info, baseline_scaled, seeded_value, scenario, case_profile);
    [lb, ub, widened_note] = ensure_seed_inside_bounds(lb, ub, seeded_value, name);

    baseline_adult_col(idx) = baseline_adult;
    baseline_scaled_col(idx) = baseline_scaled;
    seeded_value_col(idx) = seeded_value;
    bound_anchor_col(idx) = bound_anchor;
    scale_method_col{idx} = scale_method;
    lb_col(idx) = lb;
    ub_col(idx) = ub;
    bound_type_col{idx} = bound_type;
    bound_source_col{idx} = bound_source;
    confidence_col{idx} = confidence_level;

    note_parts{end + 1} = bound_note;
    note_parts{end + 1} = widened_note;
    notes_col{idx} = join_nonempty(note_parts, ' ');
end

registry = table( ...
    names, group_col, unit_col, baseline_adult_col, baseline_scaled_col, seeded_value_col, ...
    bound_anchor_col, scale_method_col, is_calibratable_col, lb_col, ub_col, ...
    bound_type_col, bound_source_col, confidence_col, notes_col, ...
    'VariableNames', {'name','group','unit','baseline_adult','baseline_scaled', ...
    'seeded_value','bound_anchor','scale_method','is_calibratable','lb','ub', ...
    'bound_type','bound_source','confidence_level','notes'});
end

function group_name = resolve_parameter_group(name, group_info)
if group_info.found
    switch group_info.groupType
        case 'R'
            group_name = 'coupled_resistance_scale';
        case 'C'
            group_name = 'coupled_compliance_scale';
        otherwise
            group_name = 'coupled_parameter_scale';
    end
    return;
end

if startsWith(name, 'R.')
    if strcmp(name, 'R.vsd')
        group_name = 'vsd_resistance';
    else
        group_name = 'resistance';
    end
elseif startsWith(name, 'C.')
    group_name = 'compliance';
elseif startsWith(name, 'L.')
    group_name = 'inertance';
elseif startsWith(name, 'E.')
    group_name = 'elastance';
elseif startsWith(name, 'V0.')
    group_name = 'unstressed_volume';
elseif strcmp(name, 'vsd.Cd')
    group_name = 'vsd_orifice';
else
    group_name = 'other';
end
end

function unit_name = resolve_parameter_unit(name, group_info)
if group_info.found
    unit_name = 'ratio';
    return;
end
if startsWith(name, 'R.')
    unit_name = 'mmHg*s/mL';
elseif startsWith(name, 'C.')
    unit_name = 'mL/mmHg';
elseif startsWith(name, 'L.')
    unit_name = 'mmHg*s^2/mL';
elseif startsWith(name, 'E.')
    unit_name = 'mmHg/mL';
elseif startsWith(name, 'V0.')
    unit_name = 'mL';
elseif strcmp(name, 'vsd.Cd')
    unit_name = 'ratio';
else
    unit_name = '-';
end
end

function scale_method = resolve_scale_method(name, scenario)
if strcmp(name, 'R.vsd') && strcmp(scenario, 'pre_surgery')
    scale_method = 'patient_specific_shunt_seed';
elseif strcmp(name, 'vsd.Cd')
    scale_method = 'patient_specific_orifice_seed';
elseif startsWith(name, 'R.')
    scale_method = 'allometric_scaled_then_clinically_seeded';
elseif startsWith(name, 'C.')
    scale_method = 'allometric_scaled_then_clinically_seeded';
elseif startsWith(name, 'L.')
    scale_method = 'allometric_scaled_or_fixed_baseline';
elseif startsWith(name, 'E.')
    scale_method = 'allometric_scaled_then_clinically_seeded';
elseif startsWith(name, 'V0.')
    scale_method = 'allometric_scaled_then_clinically_seeded';
else
    scale_method = 'scenario_specific';
end
end

function tf = resolve_calibratable_flag(name, scenario, case_profile)
if strcmp(name, 'R.vsd') && strcmp(scenario, 'post_surgery')
    tf = false;
    return;
end
if isfield(case_profile, 'allowedFreeParameters') && ~isempty(case_profile.allowedFreeParameters)
    tf = ismember(name, case_profile.allowedFreeParameters(:));
else
    tf = true;
end
end

function [lb, ub, bound_anchor, bound_type, bound_source, confidence_level, note_text] = ...
        resolve_bounds(name, group_info, baseline_scaled, seeded_value, scenario, case_profile)

if group_info.found
    switch group_info.groupType
        case 'R'
            [lb, ub] = deal(0.25, 2.80);
            bound_type = 'group_scale_multiplier';
            bound_source = 'kung2013_bsa_rc_coupled_resistance_prior';
            confidence_level = 'medium';
        case 'C'
            [lb, ub] = deal(0.20, 5.00);
            bound_type = 'group_scale_multiplier';
            bound_source = 'provisional_multiplier_scaled_baseline';
            confidence_level = 'medium';
        otherwise
            [lb, ub] = deal(0.30, 4.00);
            bound_type = 'group_scale_multiplier';
            bound_source = 'provisional_multiplier_scaled_baseline';
            confidence_level = 'low';
    end
    bound_anchor = 1.0;
    note_text = 'Grouped scale bound is dimensionless relative to the BSA-scaled baseline. Vascular compliance is derived from resistance using the Kung/Pennati-Fumero resting R-C relation.';
    [lb, ub] = apply_profile_tightening(name, group_info, lb, ub, bound_anchor, case_profile);
    return;
end

bound_anchor = baseline_scaled;
if strcmp(name, 'R.vsd') && strcmp(scenario, 'pre_surgery')
    anchor = resolve_preclosure_vsd_anchor(seeded_value);
    lb = max(0.001, 0.05 * anchor);
    ub = min(500, 20.0 * anchor);
    bound_anchor = anchor;
    bound_type = 'patient_specific_wide_box';
    bound_source = 'geometry_gradient_flow_seed_provisional_box';
    confidence_level = 'low';
    note_text = 'Pre-closure VSD resistance uses a wide patient-specific box around the seeded shunt estimate.';
elseif strcmp(name, 'vsd.Cd')
    lb = 0.20;
    ub = 1.20;
    bound_anchor = seeded_value;
    bound_type = 'physiologic_fraction';
    bound_source = 'orifice_discharge_coefficient_literature_plus_expert_box';
    confidence_level = 'medium';
    note_text = 'Discharge coefficient is constrained to a physiologic orifice-flow interval.';
elseif strcmp(name, 'C.SAR') || strcmp(name, 'C.PAR')
    if isfinite(seeded_value) && seeded_value > 0
        bound_anchor = seeded_value;
    end
    [lower_mult, upper_mult, note_text] = resolve_arterial_compliance_prior(name);
    lb = lower_mult * bound_anchor;
    ub = upper_mult * bound_anchor;
    bound_type = 'pulse_pressure_anchor';
    bound_source = 'windkessel_stroke_volume_over_pulse_pressure_anchor';
    confidence_level = 'high';
elseif startsWith(name, 'R.')
    [lower_mult, upper_mult, note_text] = resolve_resistance_prior(name);
    lb = lower_mult * baseline_scaled;
    ub = upper_mult * baseline_scaled;
    bound_type = 'class_specific_resistance_prior';
    bound_source = 'kung2013_sharp2000_scaled_resistance_prior';
    confidence_level = 'medium';
elseif startsWith(name, 'C.')
    [lower_mult, upper_mult, note_text] = resolve_compliance_prior(name);
    lb = lower_mult * baseline_scaled;
    ub = upper_mult * baseline_scaled;
    bound_type = 'class_specific_compliance_prior';
    bound_source = 'kung2013_rc_coupled_compliance_prior';
    confidence_level = 'medium';
elseif startsWith(name, 'L.')
    lb = 0.25 * baseline_scaled;
    ub = 4.00 * baseline_scaled;
    bound_type = 'baseline_multiplier';
    bound_source = 'literature_informed_provisional_multiplier';
    confidence_level = 'low';
    note_text = 'Inertance usually remains fixed; this bound is a narrow architecture-dependent placeholder.';
elseif startsWith(name, 'E.')
    [lower_mult, upper_mult, note_text] = resolve_elastance_prior(name);
    lb = lower_mult * baseline_scaled;
    ub = upper_mult * baseline_scaled;
    bound_type = 'class_specific_elastance_prior';
    bound_source = 'zhang2019_pediatric_elastance_anchor';
    confidence_level = 'medium';
elseif startsWith(name, 'V0.')
    [lower_mult, upper_mult, note_text] = resolve_v0_prior(name);
    lb = lower_mult * baseline_scaled;
    ub = upper_mult * baseline_scaled;
    bound_type = 'class_specific_v0_prior';
    bound_source = 'blood_volume_preload_consistency_prior';
    confidence_level = 'medium';
else
    lb = 0.30 * baseline_scaled;
    ub = 4.00 * baseline_scaled;
    bound_type = 'baseline_multiplier';
    bound_source = 'provisional_multiplier_scaled_baseline';
    confidence_level = 'low';
    note_text = 'Fallback multiplier bound.';
end

[lb, ub] = apply_profile_tightening(name, group_info, lb, ub, bound_anchor, case_profile);
end

function [lb, ub] = apply_profile_tightening(name, group_info, lb, ub, anchor_value, case_profile)
if ~isfield(case_profile, 'boundScale') || isempty(case_profile.boundScale) || ...
        ~isfield(case_profile.boundScale, 'names')
    return;
end

idx = find(strcmp(case_profile.boundScale.names, name), 1, 'first');
if isempty(idx)
    return;
end

if group_info.found
    lb = max(lb, case_profile.boundScale.lower(idx));
    ub = min(ub, case_profile.boundScale.upper(idx));
else
    lb = max(lb, case_profile.boundScale.lower(idx) * anchor_value);
    ub = min(ub, case_profile.boundScale.upper(idx) * anchor_value);
end
end

function [lb, ub, note_text] = ensure_seed_inside_bounds(lb, ub, seeded_value, name)
note_text = '';
if ~isfinite(seeded_value)
    return;
end

span = max(ub - lb, 1e-9);
margin = max(0.05 * abs(seeded_value), 0.02 * span);

if seeded_value < lb
    lb = max(seeded_value - margin, eps);
    note_text = sprintf(['Registry lower bound widened to include clinically seeded start for %s. ', ...
        'Scaled baseline remains the plausibility anchor.'], name);
elseif seeded_value > ub
    ub = seeded_value + margin;
    note_text = sprintf(['Registry upper bound widened to include clinically seeded start for %s. ', ...
        'Scaled baseline remains the plausibility anchor.'], name);
end
end

function anchor = resolve_preclosure_vsd_anchor(seeded_value)
if ~isfinite(seeded_value) || seeded_value >= 1e5
    anchor = 1.0;
else
    anchor = max(seeded_value, 0.01);
end
end

function [lower_mult, upper_mult, note_text] = resolve_resistance_prior(name)
lower_mult = 0.40;
upper_mult = 2.50;
note_text = ['Resistance prior is anchored to the BSA-scaled baseline and ', ...
    'clinical SVR/PVR seeding.'];

if strcmp(name, 'R.SVEN') || strcmp(name, 'R.PVEN')
    lower_mult = 0.40;
    upper_mult = 3.00;
    note_text = ['Venous resistance prior remains slightly wider while paired ', ...
        'compliance is derived through the resting R-C coupling relation.'];
end
end

function [lower_mult, upper_mult, note_text] = resolve_compliance_prior(name)
lower_mult = 0.60;
upper_mult = 1.60;
note_text = ['Compliance prior is no longer a free wide box; it is anchored to ', ...
    'the BSA-scaled baseline and expected to follow the resting R-C coupling relation.'];

if strcmp(name, 'C.SVEN') || strcmp(name, 'C.PVEN')
    lower_mult = 0.50;
    upper_mult = 1.80;
end
end

function [lower_mult, upper_mult, note_text] = resolve_arterial_compliance_prior(name)
if strcmp(name, 'C.SAR')
    lower_mult = 0.75;
    upper_mult = 1.35;
    note_text = ['Systemic arterial compliance is anchored to stroke volume over ', ...
        'systemic pulse pressure, with only a narrow residual tuning interval.'];
else
    lower_mult = 0.70;
    upper_mult = 1.45;
    note_text = ['Pulmonary arterial compliance is anchored to stroke volume over ', ...
        'pulmonary pulse pressure, with a narrow residual tuning interval.'];
end
end

function [lower_mult, upper_mult, note_text] = resolve_elastance_prior(name)
lower_mult = 0.60;
upper_mult = 2.20;
note_text = 'Elastance prior is tightened around the allometric pediatric anchor.';

if contains(name, '.RV.')
    lower_mult = 0.55;
    upper_mult = 2.60;
elseif contains(name, '.LA.') || contains(name, '.RA.')
    lower_mult = 0.50;
    upper_mult = 2.50;
elseif endsWith(name, '.EB')
    lower_mult = 0.60;
    upper_mult = 2.50;
end
end

function [lower_mult, upper_mult, note_text] = resolve_v0_prior(name)
lower_mult = 0.70;
upper_mult = 1.35;
note_text = ['Unstressed-volume prior is tied to blood-volume and preload ', ...
    'consistency rather than a broad free search box.'];

if strcmp(name, 'V0.LV')
    lower_mult = 0.75;
    upper_mult = 1.40;
elseif strcmp(name, 'V0.RV')
    lower_mult = 0.70;
    upper_mult = 1.35;
end
end

function text = join_nonempty(parts, delimiter)
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    text = '';
else
    text = strjoin(parts, delimiter);
end
end
