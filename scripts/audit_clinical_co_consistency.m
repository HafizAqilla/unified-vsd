function audit_clinical_co_consistency(clinical, scenario)
% AUDIT_CLINICAL_CO_CONSISTENCY
% -----------------------------------------------------------------------
% Audits whether clinical cardiac-output targets agree with chamber-volume
% targets before calibration. This is an input-data consistency diagnostic;
% it does not run the ODE solver or change model parameters.
%
% INPUTS:
%   clinical  - unified clinical struct from config/ patient profiles    [-]
%   scenario  - scenario string: 'pre_surgery' | 'post_surgery'          [-]
%
% OUTPUTS:
%   results/tables/clinical_co_consistency_<patient>_<scenario>_*.csv
%
% ASSUMPTIONS:
%   - Catheter Fick CO in unrepaired VSD is systemic flow Qs [L/min].
%   - LV stroke output should approximate Qp, not Qs, when left-to-right
%     shunting is present.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-08
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

if nargin < 1 || isempty(clinical)
    clinical = patient_reyna();
end
if nargin < 2 || isempty(scenario)
    scenario = 'pre_surgery';
end
if ~isfield(clinical, scenario)
    error('audit_clinical_co_consistency:unknownScenario', ...
        'Clinical struct does not contain scenario "%s".', scenario);
end

src = clinical.(scenario);
patient_label = resolve_patient_label(clinical);
hr_bpm = field_or_nan(clinical.common, 'HR');                 % [bpm]
co_qs_Lmin = field_or_nan(src, 'CO_Lmin');                    % [L/min]
qpqs = field_or_nan(src, 'QpQs');                             % [-]
lvsv_mL = field_or_nan(src, 'LVEDV_mL') - field_or_nan(src, 'LVESV_mL');
rvsv_mL = field_or_nan(src, 'RVEDV_mL') - field_or_nan(src, 'RVESV_mL');
lvco_from_volumes_Lmin = lvsv_mL * hr_bpm / 1000;             % [L/min]
rvco_from_volumes_Lmin = rvsv_mL * hr_bpm / 1000;             % [L/min]
qp_from_qs_qpqs_Lmin = co_qs_Lmin * qpqs;                     % [L/min]
qshunt_from_qpqs_Lmin = co_qs_Lmin * (qpqs - 1);              % [L/min]

rows = {
    'Fick_Qs_target_Lmin', co_qs_Lmin, 'L/min', ...
        'Systemic flow target used for CO_Lmin comparator'
    'Qp_from_Qs_times_QpQs_Lmin', qp_from_qs_qpqs_Lmin, 'L/min', ...
        'Expected pulmonary/LV output implied by Fick Qs and Qp/Qs'
    'Qshunt_from_Qp_minus_Qs_Lmin', qshunt_from_qpqs_Lmin, 'L/min', ...
        'Expected left-to-right shunt flow implied by Fick Qs and Qp/Qs'
    'LVSV_from_volumes_mL', lvsv_mL, 'mL', ...
        'LV stroke volume implied by clinical LVEDV-LVESV'
    'LVCO_from_volumes_Lmin', lvco_from_volumes_Lmin, 'L/min', ...
        'LVSV times heart rate; should approximate Qp in unrepaired VSD'
    'RVSV_from_volumes_mL', rvsv_mL, 'mL', ...
        'RV stroke volume implied by clinical RVEDV-RVESV'
    'RVCO_from_volumes_Lmin', rvco_from_volumes_Lmin, 'L/min', ...
        'RVSV times heart rate; should approximate Qp at steady state'
    'LVCO_minus_Qp_target_Lmin', lvco_from_volumes_Lmin - qp_from_qs_qpqs_Lmin, 'L/min', ...
        'Negative value indicates volume-derived LV output below Fick-implied Qp'
    'LVCO_to_Qp_target_ratio', lvco_from_volumes_Lmin / max(qp_from_qs_qpqs_Lmin, 1e-9), '-', ...
        'Ratio close to 1 means LV volumes and Fick/QpQs agree'
    };

audit_tbl = cell2table(rows, 'VariableNames', ...
    {'Check','Value','Unit','Interpretation'});

tables_dir = fullfile(root, 'results', 'tables');
if ~exist(tables_dir, 'dir')
    mkdir(tables_dir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
csv_path = fullfile(tables_dir, sprintf( ...
    'clinical_co_consistency_%s_%s_%s.csv', patient_label, scenario, timestamp));
txt_path = strrep(csv_path, '.csv', '.txt');
writetable(audit_tbl, csv_path);
write_text_summary(txt_path, patient_label, scenario, audit_tbl);

fprintf('[audit_clinical_co_consistency] Audit exported:\n  %s\n  %s\n', ...
    csv_path, txt_path);
disp(audit_tbl);
end

function write_text_summary(txt_path, patient_label, scenario, audit_tbl)
% WRITE_TEXT_SUMMARY - plain-language audit summary for thesis notes.
fid = fopen(txt_path, 'w');
if fid < 0
    warning('audit_clinical_co_consistency:openFailed', ...
        'Unable to write text summary.');
    return;
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Clinical CO Consistency Audit\n');
fprintf(fid, '=============================\n');
fprintf(fid, 'Patient: %s\n', patient_label);
fprintf(fid, 'Scenario: %s\n\n', scenario);
for idx = 1:height(audit_tbl)
    fprintf(fid, '%s: %.6g %s\n', ...
        audit_tbl.Check{idx}, audit_tbl.Value(idx), audit_tbl.Unit{idx});
    fprintf(fid, '  %s\n', audit_tbl.Interpretation{idx});
end
end

function label = resolve_patient_label(clinical)
% RESOLVE_PATIENT_LABEL - folder-safe patient label.
label = 'patient_unknown';
if isfield(clinical, 'common') && isfield(clinical.common, 'patient_name') && ...
        ~isempty(clinical.common.patient_name)
    label = lower(regexprep(char(clinical.common.patient_name), '[^a-zA-Z0-9_]+', '_'));
end
end

function value = field_or_nan(s, field_name)
% FIELD_OR_NAN - scalar numeric field lookup with NaN fallback.
if isstruct(s) && isfield(s, field_name) && isnumeric(s.(field_name)) && ...
        isscalar(s.(field_name)) && isfinite(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

function project_path = build_clean_project_path(root)
% BUILD_CLEAN_PROJECT_PATH - omit hidden worktrees that may shadow files.
project_paths = strsplit(genpath(root), pathsep);
project_paths = project_paths(~cellfun('isempty', project_paths));
is_shadow = contains(project_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, project_paths);
project_path = strjoin(project_paths(~is_shadow & is_existing), pathsep);
end
