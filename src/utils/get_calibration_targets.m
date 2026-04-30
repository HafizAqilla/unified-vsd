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
            'PVR',          'PVR_WU',         'WU',    'Pulmonary vascular resistance',        true,  true,  false, 'High'
            'SVR',          'SVR_WU',         'WU',    'Systemic vascular resistance',         true,  false, true,  'High'
            'CO_Lmin',      'CO_Lmin',        'L/min', 'Effective systemic cardiac output (Qs)', true, false, true, 'High'
            'VSD_frac_pct', 'VSD_frac_pct',   '%',     'VSD shunt fraction of Qp',             false, false, false, 'Derived'
            'LVEDV',        'LVEDV_mL',       'mL',    'LV end-diastolic volume',              true,  false, true,  'Moderate'
            'LVESV',        'LVESV_mL',       'mL',    'LV end-systolic volume',               true,  false, true,  'Moderate'
            'RVEDV',        'RVEDV_mL',       'mL',    'RV end-diastolic volume',              true,  false, true,  'Moderate'
            'RVESV',        'RVESV_mL',       'mL',    'RV end-systolic volume',               true,  false, true,  'Moderate'
            'LVEF',         'LVEF',           '-',     'LV ejection fraction',                 true,  false, true,  'Moderate'
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
            'LVEF',         'LVEF',          '-',     'LV ejection fraction',                 true,  false, true,  'Moderate'
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
    if isfield(src, targets(i).ClinicalField)
        targets(i).ClinicalValue = src.(targets(i).ClinicalField);
    end
end

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
    'ClinicalValue', NaN);
end
