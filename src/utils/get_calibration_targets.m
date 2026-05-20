function targets = get_calibration_targets(scenario, clinical)
% GET_CALIBRATION_TARGETS
% -----------------------------------------------------------------------
% Single source of truth for model metric to clinical field mappings.
%
% INPUTS:
%   scenario  - clinical scenario: 'pre_surgery' or 'post_surgery'
%   clinical  - unified clinical struct from patient_template()
%
% OUTPUTS:
%   targets   - struct array with metric metadata and clinical values
%               UncertaintyFraction is relative 1-sigma uncertainty [-];
%               UncertaintyAbs is absolute 1-sigma uncertainty [Unit].
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(clinical)
    clinical = patient_template();
end

switch scenario
    case 'pre_surgery'
        src = clinical.pre_surgery;
        rows = {
            'RAP_min',      'RAP_dia_mmHg',   'mmHg',  'Right atrial diastolic-like minimum', false, false, false, 'Low'
            'RAP_mean',     'RAP_mean_mmHg',  'mmHg',  'Right atrial mean pressure',           true,  false, true,  'High'
            'RAP_max',      'RAP_sys_mmHg',   'mmHg',  'Right atrial systolic-like maximum',   false, false, false, 'Low'
            'LAP_min',      'LAP_dia_mmHg',   'mmHg',  'Left atrial diastolic-like minimum',   false, false, false, 'Low'
            'LAP_mean',     'LAP_mean_mmHg',  'mmHg',  'Left atrial mean pressure',            true,  false, true,  'Moderate'
            'LAP_max',      'LAP_sys_mmHg',   'mmHg',  'Left atrial systolic-like maximum',    false, false, false, 'Low'
            'PAP_min',      'PAP_dia_mmHg',   'mmHg',  'PA diastolic pressure',                true,  false, true,  'High'
            'PAP_max',      'PAP_sys_mmHg',   'mmHg',  'PA systolic pressure',                 true,  false, true,  'High'
            'PAP_mean',     'PAP_mean_mmHg',  'mmHg',  'PA mean pressure',                     true,  true,  false, 'High'
            'SAP_min',      'SAP_dia_mmHg',   'mmHg',  'Systemic arterial diastolic pressure', true,  false, true,  'Moderate'
            'SAP_max',      'SAP_sys_mmHg',   'mmHg',  'Systemic arterial systolic pressure',  true,  false, true,  'Moderate'
            'SAP_mean',     'SAP_mean_mmHg',  'mmHg',  'Mean arterial pressure',               true,  true,  false, 'High'
            'QpQs',         'QpQs',           '-',     'Pulmonary/systemic flow ratio',        true,  true,  false, 'High'
            'Q_shunt_Lmin', 'Q_shunt_Lmin',   'L/min', 'Net shunt flow derived as Qp minus Qs', true,  false, false, 'Derived'
            'PVR',          'PVR_WU',         'WU',    'Pulmonary vascular resistance',        true,  true,  false, 'High'
            'SVR',          'SVR_WU',         'WU',    'Systemic vascular resistance',         true,  false, true,  'High'
            'CO_Lmin',      'CO_Lmin',        'L/min', 'Effective systemic cardiac output (Qs)', true, false, true, 'High'
            'VSD_frac_pct', 'VSD_frac_pct',   '%',     'VSD shunt fraction of Qp',             false, false, false, 'Derived'
            'LVEDV',        'LVEDV_mL',       'mL',    'LV end-diastolic volume',              true,  false, true,  'Moderate'
            'LVESV',        'LVESV_mL',       'mL',    'LV end-systolic volume',               true,  false, true,  'Moderate'
            'RVEDV',        'RVEDV_mL',       'mL',    'RV end-diastolic volume',              true,  false, true,  'Moderate'
            'RVESV',        'RVESV_mL',       'mL',    'RV end-systolic volume',               true,  false, true,  'Moderate'
            'LVEF',         'EF',             '-',     'LV ejection fraction',                 true,  false, true,  'Moderate'
            'RVEF',         'RVEF',           '-',     'RV ejection fraction',                 false, false, true,  'Moderate'
            };

    case 'post_surgery'
        src = clinical.post_surgery;
        rows = {
            'RAP_min',      'RAP_dia_mmHg',  'mmHg',  'Right atrial diastolic-like minimum', false, false, false, 'Low'
            'RAP_mean',     'RAP_mean_mmHg', 'mmHg',  'Right atrial mean pressure',           true,  false, true,  'High'
            'RAP_max',      'RAP_sys_mmHg',  'mmHg',  'Right atrial systolic-like maximum',   false, false, false, 'Low'
            'LAP_min',      'LAP_dia_mmHg',  'mmHg',  'Left atrial diastolic-like minimum',   false, false, false, 'Low'
            'LAP_mean',     'LAP_mean_mmHg', 'mmHg',  'Left atrial mean pressure',            false, false, true,  'Moderate'
            'LAP_max',      'LAP_sys_mmHg',  'mmHg',  'Left atrial systolic-like maximum',    false, false, false, 'Low'
            'PAP_min',      'PAP_dia_mmHg',  'mmHg',  'PA diastolic pressure',                false, false, true,  'High'
            'PAP_max',      'PAP_sys_mmHg',  'mmHg',  'PA systolic pressure',                 false, false, true,  'High'
            'PAP_mean',     'PAP_mean_mmHg', 'mmHg',  'PA mean pressure',                     true,  false, true,  'High'
            'SAP_min',      'SAP_dia_mmHg',  'mmHg',  'Systemic arterial diastolic pressure', true,  false, true,  'Moderate'
            'SAP_max',      'SAP_sys_mmHg',  'mmHg',  'Systemic arterial systolic pressure',  true,  false, true,  'Moderate'
            'SAP_mean',     'MAP_mmHg',      'mmHg',  'Mean arterial pressure',               true,  true,  false, 'High'
            'QpQs',         'QpQs',          '-',     'Qp/Qs ratio after VSD closure',        true,  true,  false, 'High'
            'PVR',          'PVR_WU',        'WU',    'Pulmonary vascular resistance',        true,  true,  false, 'High'
            'SVR',          'SVR_WU',        'WU',    'Systemic vascular resistance',         true,  false, true,  'High'
            'CO_Lmin',      'CO_Lmin',       'L/min', 'Effective systemic cardiac output (Qs)', true, false, true, 'High'
            'VSD_frac_pct', 'VSD_frac_pct',  '%',     'Residual VSD shunt fraction of Qp',    false, false, false, 'Derived'
            'LVEDV',        'LVEDV_mL',      'mL',    'LV end-diastolic volume',              true,  false, true,  'Moderate'
            'RVEDV',        'RVEDV_mL',      'mL',    'RV end-diastolic volume',              true,  false, true,  'Moderate'
            'LVEF',         'EF',            '-',     'LV ejection fraction',                 true,  false, true,  'Moderate'
            'RVEF',         'RVEF',          '-',     'RV ejection fraction',                 true,  false, true,  'Moderate'
            };

    otherwise
        error('get_calibration_targets:unknownScenario', ...
              'scenario must be ''pre_surgery'' or ''post_surgery''.');
end

targets = repmat(empty_target(), size(rows, 1), 1);
for i = 1:size(rows, 1)
    targets(i).Metric            = rows{i, 1};
    targets(i).ClinicalField     = rows{i, 2};
    targets(i).Unit              = rows{i, 3};
    targets(i).Description       = rows{i, 4};
    targets(i).UseForCalibration = rows{i, 5};
    targets(i).MandatoryPrimary  = rows{i, 6};
    targets(i).CandidatePrimary  = rows{i, 7};
    targets(i).Reliability       = rows{i, 8};
    [targets(i).UncertaintyFraction, targets(i).UncertaintyAbs, targets(i).UncertaintyNote] = ...
        reliability_uncertainty(rows{i, 8}, rows{i, 3});
    targets(i).Comparator        = rows{i, 1};
    if isfield(src, targets(i).ClinicalField)
        targets(i).ClinicalValue = src.(targets(i).ClinicalField);
    end
end

targets = apply_patient_specific_uncertainty(targets, scenario, clinical);
targets = apply_source_declared_uncertainty(targets, scenario, clinical);

end


function target = empty_target()
target = struct( ...
    'Metric', '', ...
    'ClinicalField', '', ...
    'Unit', '', ...
    'Description', '', ...
    'UseForCalibration', false, ...
    'MandatoryPrimary', false, ...
    'CandidatePrimary', false, ...
    'Reliability', '', ...
    'UncertaintyFraction', NaN, ...
    'UncertaintyAbs', NaN, ...
    'UncertaintyNote', '', ...
    'Comparator', '', ...
    'ClinicalValue', NaN);
end

function [frac, abs_sigma, note] = reliability_uncertainty(reliability, unit)
% RELIABILITY_UNCERTAINTY - map qualitative reliability to default sigma.
switch lower(strtrim(reliability))
    case 'high'
        frac = 0.05;
    case 'moderate'
        frac = 0.10;
    otherwise
        frac = 0.20;
end
abs_sigma = NaN;
note = sprintf('%s reliability mapped to %.0f%% relative uncertainty.', ...
    reliability, 100 * frac);
if strcmp(unit, '-')
    note = sprintf('%s Dimensionless target.', note);
end
end

function targets = apply_patient_specific_uncertainty(targets, scenario, clinical)
% APPLY_PATIENT_SPECIFIC_UNCERTAINTY - record known target inconsistencies.
if ~strcmp(scenario, 'pre_surgery') || ~isstruct(clinical) || ...
        ~isfield(clinical, 'common') || ~isfield(clinical.common, 'patient_name')
    return;
end

case_id = lower(char(clinical.common.patient_name));
if ~strcmp(case_id, 'reyna')
    return;
end

idx = find(strcmp({targets.Metric}, 'CO_Lmin'), 1, 'first');
if isempty(idx)
    return;
end

targets(idx).Reliability = 'Moderate';
targets(idx).UncertaintyFraction = 0.15;
targets(idx).UncertaintyAbs = 0.50;
targets(idx).Comparator = 'Qs_Lmin';
targets(idx).UncertaintyNote = [ ...
    'Reyna catheter Fick CO is treated as systemic flow Qs with widened ' ...
    'uncertainty because echo-derived LV stroke volume is internally lower.'];
end

function targets = apply_source_declared_uncertainty(targets, scenario, clinical)
% APPLY_SOURCE_DECLARED_UNCERTAINTY - honor uncertainty metadata in profiles.
if ~isstruct(clinical) || ~isfield(clinical, scenario)
    return;
end
src = clinical.(scenario);
idx = find(strcmp({targets.Metric}, 'CO_Lmin'), 1, 'first');
if isempty(idx)
    return;
end
if isfield(src, 'CO_uncertainty_Lmin') && isfinite(src.CO_uncertainty_Lmin) && ...
        src.CO_uncertainty_Lmin > 0
    targets(idx).UncertaintyAbs = src.CO_uncertainty_Lmin;
end
if isfield(src, 'CO_comparator') && ~isempty(src.CO_comparator)
    targets(idx).Comparator = char(src.CO_comparator);
end
end
