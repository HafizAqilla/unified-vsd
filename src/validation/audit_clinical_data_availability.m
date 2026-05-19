function evidence_table = audit_clinical_data_availability(clinical, scenario)
% AUDIT_CLINICAL_DATA_AVAILABILITY
% -----------------------------------------------------------------------
% Ranks scenario clinical fields from directly available evidence to
% unavailable fields so patient configs can be reviewed against protocol.
%
% INPUTS:
%   clinical - unified clinical struct from config/patient_*.m          [-]
%   scenario - scenario string: 'pre_surgery' | 'post_surgery'          [-]
%
% OUTPUTS:
%   evidence_table - table of field availability, provenance, and use    [-]
%
% ASSUMPTIONS:
%   - Directly reported/measured values are strongest evidence.
%   - Values derived algebraically from directly reported values are valid
%     evidence, but should not be double-counted as independent targets.
%   - Model-control fields are useful for reproducibility but are not
%     clinical evidence.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%   [2] AGENTS.md Sections 9-10.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-18
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(scenario)
    scenario = 'pre_surgery';
end

rows = expected_field_rows(scenario);
rows = append_unlisted_fields(rows, clinical, scenario);
evidence_table = build_evidence_table(rows, clinical, scenario);
end

function evidence_table = build_evidence_table(rows, clinical, scenario)
% BUILD_EVIDENCE_TABLE - resolve registry rows into ranked evidence table.
n_rows = size(rows, 1);
rank_col = nan(n_rows, 1);
available_col = false(n_rows, 1);
section_col = cell(n_rows, 1);
field_col = cell(n_rows, 1);
label_col = cell(n_rows, 1);
unit_col = cell(n_rows, 1);
evidence_class_col = cell(n_rows, 1);
planned_class_col = cell(n_rows, 1);
recommended_use_col = cell(n_rows, 1);
value_col = cell(n_rows, 1);
note_col = cell(n_rows, 1);

for idx = 1:n_rows
    section_name = rows{idx, 1};
    field_name = rows{idx, 2};
    planned_class = rows{idx, 5};
    value = resolve_field_value(clinical, scenario, section_name, field_name);
    is_available = is_available_value(value);
    evidence_class = resolved_evidence_class(planned_class, is_available);

    rank_col(idx) = evidence_rank(evidence_class);
    available_col(idx) = is_available;
    section_col{idx} = section_name;
    field_col{idx} = field_name;
    label_col{idx} = rows{idx, 3};
    unit_col{idx} = rows{idx, 4};
    evidence_class_col{idx} = evidence_class;
    planned_class_col{idx} = planned_class;
    recommended_use_col{idx} = recommended_use(evidence_class, planned_class);
    value_col{idx} = value_to_text(value);
    note_col{idx} = rows{idx, 6};
end

evidence_table = table(rank_col, available_col, section_col, field_col, ...
    label_col, unit_col, value_col, evidence_class_col, planned_class_col, ...
    recommended_use_col, note_col, ...
    'VariableNames', {'AvailabilityRank','Available','Section','Field', ...
    'Label','Unit','Value','EvidenceClass','PlannedClass', ...
    'RecommendedUse','Note'});
evidence_table = sortrows(evidence_table, ...
    {'AvailabilityRank','Section','Field'}, {'ascend','ascend','ascend'});
end

function rows = expected_field_rows(scenario)
% EXPECTED_FIELD_ROWS - protocol-facing and model-control field registry.
rows = common_field_rows();

if strcmp(scenario, 'pre_surgery')
    rows = [rows; pre_surgery_field_rows()];
elseif strcmp(scenario, 'post_surgery')
    rows = [rows; post_surgery_field_rows()];
end
end

function rows = common_field_rows()
% COMMON_FIELD_ROWS - demographics and model-control registry rows.
rows = {
    'common', 'age_years', 'Age', 'years', 'direct_reported', 'Protocol row 1.'
    'common', 'sex', 'Sex', '-', 'direct_reported', 'Protocol row 2; 0=female, 1=male.'
    'common', 'height_cm', 'Height', 'cm', 'direct_reported', 'Protocol row 3.'
    'common', 'weight_kg', 'Weight', 'kg', 'direct_reported', 'Protocol row 4.'
    'common', 'BSA', 'Body surface area', 'm^2', 'directly_derived', 'Protocol row 5, derived from height and weight.'
    'common', 'HR', 'Heart rate', 'bpm', 'direct_reported', 'Protocol row 6.'
    'common', 'patient_name', 'Patient label', '-', 'model_setting', 'Run/report label, not clinical evidence.'
    'common', 'patient_id', 'Patient identifier', '-', 'model_setting', 'Optional safe identifier, not a fitting target.'
    'common', 'maturation_mode', 'Maturation mode', '-', 'model_setting', 'Model-control field for scaling.'
    };
end

function rows = pre_surgery_field_rows()
% PRE_SURGERY_FIELD_ROWS - pre-release occluder clinical evidence registry.
rows = {
    'pre_surgery', 'VSD_diameter_mm', 'Effective VSD diameter', 'mm', 'directly_derived', 'Protocol row 7, effective RV-side midpoint.'
    'pre_surgery', 'VSD_gradient_mmHg', 'VSD pressure gradient', 'mmHg', 'directly_derived', 'Protocol row 9, LV systolic minus RV systolic.'
    'pre_surgery', 'Q_shunt_Lmin', 'Net shunt flow', 'L/min', 'directly_derived', 'Derived as Qp-Qs from rows 21 and 23.'
    'pre_surgery', 'QpQs', 'Qp/Qs ratio', '-', 'direct_reported', 'Protocol row 23.'
    'pre_surgery', 'VSD_mode', 'VSD model mode', '-', 'model_setting', 'Model choice, not protocol evidence.'
    'pre_surgery', 'PAP_sys_mmHg', 'Pulmonary arterial systolic pressure', 'mmHg', 'direct_reported', 'Protocol row 16.'
    'pre_surgery', 'PAP_dia_mmHg', 'Pulmonary arterial diastolic pressure', 'mmHg', 'direct_reported', 'Protocol row 17.'
    'pre_surgery', 'PAP_mean_mmHg', 'Pulmonary arterial mean pressure', 'mmHg', 'direct_reported', 'Protocol row 18.'
    'pre_surgery', 'PVR_WU', 'Pulmonary vascular resistance', 'WU', 'directly_derived', 'Protocol row 24 if calculated.'
    'pre_surgery', 'SAP_sys_mmHg', 'Systemic arterial systolic pressure', 'mmHg', 'direct_reported', 'Protocol row 12, RFA catheter.'
    'pre_surgery', 'SAP_dia_mmHg', 'Systemic arterial diastolic pressure', 'mmHg', 'direct_reported', 'Protocol row 13, RFA catheter.'
    'pre_surgery', 'SAP_mean_mmHg', 'Systemic arterial mean pressure', 'mmHg', 'directly_derived', 'Derived from RFA systolic/diastolic pressure.'
    'pre_surgery', 'SVR_WU', 'Systemic vascular resistance', 'WU', 'directly_derived', 'Protocol row 25 if calculated.'
    'pre_surgery', 'RAP_mean_mmHg', 'Right atrial mean pressure', 'mmHg', 'direct_reported', 'Protocol row 19.'
    'pre_surgery', 'LAP_mean_mmHg', 'Left atrial mean pressure', 'mmHg', 'direct_reported', 'Only if measured/reported.'
    'pre_surgery', 'LVEDP_mmHg', 'LV end-diastolic pressure', 'mmHg', 'direct_reported', 'Only if measured/reported.'
    'pre_surgery', 'LVEDV_mL', 'LV end-diastolic volume', 'mL', 'direct_reported', 'Protocol row 26.'
    'pre_surgery', 'LVESV_mL', 'LV end-systolic volume', 'mL', 'direct_reported', 'Protocol row 27.'
    'pre_surgery', 'RVEDV_mL', 'RV end-diastolic volume', 'mL', 'direct_reported', 'Protocol row 28.'
    'pre_surgery', 'RVESV_mL', 'RV end-systolic volume', 'mL', 'direct_reported', 'Protocol row 29.'
    'pre_surgery', 'EF', 'Ejection fraction', '-', 'directly_derived', 'Derived from LVEDV and LVESV when EF is not separately reported.'
    'pre_surgery', 'override_IC', 'Initial-condition override', '-', 'model_setting', 'Model-control field, not protocol evidence.'
    'pre_surgery', 'CO_comparator', 'CO comparator', '-', 'model_setting', 'Model-control field, not protocol evidence.'
    'pre_surgery', 'CO_uncertainty_Lmin', 'CO uncertainty', 'L/min', 'model_setting', 'Model-control uncertainty for fitting.'
    'pre_surgery', 'CO_Lmin', 'Systemic flow Qs', 'L/min', 'directly_derived', 'Back-calculated from reported Qp and Qp/Qs.'
    };
end

function rows = post_surgery_field_rows()
% POST_SURGERY_FIELD_ROWS - post-closure clinical evidence registry.
rows = {
    'post_surgery', 'QpQs', 'Post-closure Qp/Qs ratio', '-', 'direct_reported', 'Post-op shunt ratio if reported.'
    'post_surgery', 'PAP_sys_mmHg', 'Pulmonary arterial systolic pressure', 'mmHg', 'direct_reported', 'Post-op catheter/echo pressure if reported.'
    'post_surgery', 'PAP_dia_mmHg', 'Pulmonary arterial diastolic pressure', 'mmHg', 'direct_reported', 'Post-op catheter/echo pressure if reported.'
    'post_surgery', 'PAP_mean_mmHg', 'Pulmonary arterial mean pressure', 'mmHg', 'direct_reported', 'Post-op catheter/echo pressure if reported.'
    'post_surgery', 'PVR_WU', 'Pulmonary vascular resistance', 'WU', 'directly_derived', 'Post-op calculated resistance if source pressures/flows are available.'
    'post_surgery', 'SAP_sys_mmHg', 'Systemic arterial systolic pressure', 'mmHg', 'direct_reported', 'Post-op pressure if reported.'
    'post_surgery', 'SAP_dia_mmHg', 'Systemic arterial diastolic pressure', 'mmHg', 'direct_reported', 'Post-op pressure if reported.'
    'post_surgery', 'MAP_mmHg', 'Mean arterial pressure', 'mmHg', 'directly_derived', 'Post-op calculated/reported MAP.'
    'post_surgery', 'SVR_WU', 'Systemic vascular resistance', 'WU', 'directly_derived', 'Post-op calculated resistance if source pressures/flows are available.'
    'post_surgery', 'RAP_mean_mmHg', 'Right atrial mean pressure', 'mmHg', 'direct_reported', 'Post-op catheter pressure if reported.'
    'post_surgery', 'LAP_mean_mmHg', 'Left atrial mean pressure', 'mmHg', 'direct_reported', 'Post-op catheter/PCWP pressure if reported.'
    'post_surgery', 'LVEDV_mL', 'LV end-diastolic volume', 'mL', 'direct_reported', 'Post-op echo/MRI volume if reported.'
    'post_surgery', 'LVESV_mL', 'LV end-systolic volume', 'mL', 'direct_reported', 'Post-op echo/MRI volume if reported.'
    'post_surgery', 'RVEDV_mL', 'RV end-diastolic volume', 'mL', 'direct_reported', 'Post-op echo/MRI volume if reported.'
    'post_surgery', 'RVESV_mL', 'RV end-systolic volume', 'mL', 'direct_reported', 'Post-op echo/MRI volume if reported.'
    'post_surgery', 'EF', 'Ejection fraction', '-', 'directly_derived', 'Use as derived if LVEDV/LVESV are also available.'
    'post_surgery', 'RVEF', 'RV ejection fraction', '-', 'directly_derived', 'Use as derived if RVEDV/RVESV are also available.'
    'post_surgery', 'CO_Lmin', 'Systemic cardiac output', 'L/min', 'direct_reported', 'Post-op Fick/thermodilution output if reported.'
    };
end

function rows = append_unlisted_fields(rows, clinical, scenario)
% APPEND_UNLISTED_FIELDS - keep newly added config fields visible.
known = strcat(rows(:, 1), '.', rows(:, 2));
sections = {'common', scenario};
for section_idx = 1:numel(sections)
    section_name = sections{section_idx};
    if ~isfield(clinical, section_name) || ~isstruct(clinical.(section_name))
        continue;
    end
    names = fieldnames(clinical.(section_name));
    for idx = 1:numel(names)
        key = sprintf('%s.%s', section_name, names{idx});
        if ismember(key, known)
            continue;
        end
        rows(end + 1, :) = {section_name, names{idx}, names{idx}, '-', ...
            'unclassified', 'Field is not in the clinical evidence registry yet.'}; %#ok<AGROW>
        known{end + 1} = key; %#ok<AGROW>
    end
end
end

function value = resolve_field_value(clinical, scenario, section_name, field_name)
value = NaN;
if strcmp(section_name, 'common')
    if isfield(clinical, 'common') && isfield(clinical.common, field_name)
        value = clinical.common.(field_name);
    end
elseif strcmp(section_name, scenario)
    if isfield(clinical, scenario) && isfield(clinical.(scenario), field_name)
        value = clinical.(scenario).(field_name);
    end
end
end

function tf = is_available_value(value)
if isnumeric(value) || islogical(value)
    tf = isscalar(value) && isfinite(double(value));
elseif ischar(value) || isstring(value)
    tf = strlength(string(value)) > 0;
else
    tf = ~isempty(value);
end
end

function evidence_class = resolved_evidence_class(planned_class, is_available)
if ~is_available
    evidence_class = 'unavailable';
else
    evidence_class = planned_class;
end
end

function rank = evidence_rank(evidence_class)
switch evidence_class
    case 'direct_reported'
        rank = 1;
    case 'directly_derived'
        rank = 2;
    case 'model_setting'
        rank = 3;
    case 'unclassified'
        rank = 4;
    otherwise
        rank = 5;
end
end

function use_text = recommended_use(evidence_class, planned_class)
switch evidence_class
    case 'direct_reported'
        use_text = 'candidate_fit_or_validation_target';
    case 'directly_derived'
        use_text = 'use_with_dependency_check';
    case 'model_setting'
        use_text = 'reproducibility_only_not_clinical_target';
    case 'unclassified'
        use_text = 'review_before_use';
    otherwise
        if strcmp(planned_class, 'model_setting')
            use_text = 'leave_unset_or_default_model_control';
        else
            use_text = 'leave_as_NaN';
        end
end
end

function text = value_to_text(value)
if isnumeric(value)
    if isscalar(value)
        if isfinite(value)
            text = sprintf('%.10g', value);
        else
            text = 'NaN';
        end
    else
        text = sprintf('[%s]', strjoin(arrayfun(@(x) sprintf('%.10g', x), ...
            value(:)', 'UniformOutput', false), ','));
    end
elseif islogical(value)
    if value
        text = 'true';
    else
        text = 'false';
    end
elseif isstring(value)
    text = char(value);
elseif ischar(value)
    text = value;
else
    text = class(value);
end
end
