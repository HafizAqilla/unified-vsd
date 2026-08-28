function ranges = clinical_reference_ranges(scenario, clinical, case_profile)
% CLINICAL_REFERENCE_RANGES
% -----------------------------------------------------------------------
% Broad pediatric VSD plausibility ranges used to gate baseline outputs.
% These are screening ranges, not a replacement for target-tier governance.
% -----------------------------------------------------------------------

if nargin < 1 || isempty(scenario)
    scenario = 'pre_surgery';
end
if nargin < 2
    clinical = struct();
end
if nargin < 3
    case_profile = struct(); %#ok<NASGU>
end

common = struct();
if isstruct(clinical) && isfield(clinical, 'common')
    common = clinical.common;
end
age_years = get_field_or(common, 'age_years', NaN);

[hr_soft, hr_hard] = pediatric_hr_ranges(age_years);

rows = {
    'HR',       'bpm',   hr_soft(1),  hr_soft(2),  hr_hard(1),  hr_hard(2),  'Heart rate plausibility for pediatric baseline.'
    'RAP_mean', 'mmHg',  1,           10,          0,           15,          'Right atrial pressure baseline screen.'
    'PAP_mean', 'mmHg',  8,           35,          5,           55,          'Mean pulmonary arterial pressure baseline screen.'
    'SAP_mean', 'mmHg',  50,          95,          35,          120,         'Mean systemic pressure baseline screen.'
    'CO_Lmin',  'L/min', 1.0,         6.5,         0.4,         9.0,         'Systemic cardiac output baseline screen.'
    'QpQs',     '-',     0.7,         4.5,         0.3,         8.0,         'Pulmonary-to-systemic flow ratio baseline screen.'
    'LVEF',     '-',     0.40,        0.85,        0.20,        0.95,        'Left ventricular ejection fraction screen.'
    'RVEF',     '-',     0.30,        0.80,        0.15,        0.95,        'Right ventricular ejection fraction screen.'
    'LVEDV',    'mL',    5,           120,         1,           220,         'LV end-diastolic volume broad pediatric screen.'
    'LVESV',    'mL',    1,           80,          0.2,         180,         'LV end-systolic volume broad pediatric screen.'
    'RVEDV',    'mL',    5,           140,         1,           260,         'RV end-diastolic volume broad pediatric screen.'
    'RVESV',    'mL',    1,           90,          0.2,         200,         'RV end-systolic volume broad pediatric screen.'
    'PVR',      'WU',    0.3,         12,          0.05,        25,          'Pulmonary vascular resistance baseline screen.'
    'SVR',      'WU',    8,           140,         2,           260,         'Systemic vascular resistance baseline screen.'
    };

if strcmpi(scenario, 'post_surgery')
    rows(strcmp(rows(:, 1), 'QpQs'), 3:6) = {0.85, 1.25, 0.5, 2.0};
end

ranges = cell2table(rows, 'VariableNames', ...
    {'Metric','Unit','SoftLow','SoftHigh','HardLow','HardHigh','Rationale'});
end

function [soft, hard] = pediatric_hr_ranges(age_years)
if ~isfinite(age_years)
    soft = [60 140];
    hard = [40 220];
elseif age_years < 1
    soft = [90 170];
    hard = [60 240];
elseif age_years < 5
    soft = [70 150];
    hard = [45 220];
elseif age_years < 12
    soft = [60 130];
    hard = [40 200];
else
    soft = [50 120];
    hard = [35 180];
end
end

function value = get_field_or(s, field_name, fallback)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name)) && ...
        isfinite(s.(field_name))
    value = s.(field_name);
else
    value = fallback;
end
end

