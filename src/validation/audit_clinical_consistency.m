function audit = audit_clinical_consistency(clinical, scenario, config)
% AUDIT_CLINICAL_CONSISTENCY
% -----------------------------------------------------------------------
% Audits whether flow-derived and volume-derived stroke-volume estimates
% are mutually consistent before patient-specific calibration.
%
% INPUTS:
%   clinical  - unified clinical struct from config/ patient profiles    [-]
%   scenario  - scenario string: 'pre_surgery' | 'post_surgery'          [-]
%   config    - optional threshold/configuration struct                  [-]
%
% OUTPUTS:
%   audit     - struct containing stroke volumes, pairwise differences,
%               severity flags, notes, and target-tier recommendations   [-]
%
% ASSUMPTIONS:
%   - CO_Lmin is interpreted according to the clinical profile comparator;
%     for Reyna pre-surgery it is systemic flow Qs.
%   - In unrepaired left-to-right VSD, Qp = Qs * QpQs and ventricular
%     ejected stroke volumes should be checked against Qp, not Qs alone.
%   - Large disagreement between catheter-derived flow and echo-derived
%     chamber volumes is treated as data-governance uncertainty, not as
%     proof that either modality is wrong.
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%   [2] AGENTS.md Section 9: clinical data reliability and uncertainty.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-14
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(scenario)
    scenario = 'pre_surgery';
end
if nargin < 3 || isempty(config)
    config = default_audit_config();
else
    config = merge_struct(default_audit_config(), config);
end

audit = empty_audit(config);
if nargin < 1 || isempty(clinical) || ~isstruct(clinical) || ...
        ~isfield(clinical, scenario)
    audit.notes = {'Clinical struct or scenario is unavailable; audit not evaluable.'};
    return;
end

src = clinical.(scenario);                          % [-]
common = struct();                                  % [-]
if isfield(clinical, 'common')
    common = clinical.common;
end

HR_bpm = field_or_nan(common, 'HR');                % [bpm]
CO_Lmin = field_or_nan(src, 'CO_Lmin');             % [L/min]
QpQs = field_or_nan(src, 'QpQs');                   % [-]
LVEDV_mL = field_or_nan(src, 'LVEDV_mL');           % [mL]
LVESV_mL = field_or_nan(src, 'LVESV_mL');           % [mL]
RVEDV_mL = field_or_nan(src, 'RVEDV_mL');           % [mL]
RVESV_mL = field_or_nan(src, 'RVESV_mL');           % [mL]

SV_Qs_mL = NaN;                                     % [mL/beat]
SV_Qp_mL = NaN;                                     % [mL/beat]
if isfinite(HR_bpm) && HR_bpm > 0 && isfinite(CO_Lmin)
    SV_Qs_mL = CO_Lmin * 1000 / HR_bpm;
    if isfinite(QpQs)
        SV_Qp_mL = CO_Lmin * QpQs * 1000 / HR_bpm;
    end
end

SV_LV_mL = LVEDV_mL - LVESV_mL;                     % [mL/beat]
SV_RV_mL = RVEDV_mL - RVESV_mL;                     % [mL/beat]
if ~all(isfinite([LVEDV_mL, LVESV_mL]))
    SV_LV_mL = NaN;
end
if ~all(isfinite([RVEDV_mL, RVESV_mL]))
    SV_RV_mL = NaN;
end

audit.HR_bpm = HR_bpm;
audit.CO_Lmin = CO_Lmin;
audit.QpQs = QpQs;
audit.SV_Qs = SV_Qs_mL;
audit.SV_Qp = SV_Qp_mL;
audit.SV_LV = SV_LV_mL;
audit.SV_RV = SV_RV_mL;
audit.stroke_volume_table = table( ...
    {'SV_Qs'; 'SV_Qp'; 'SV_LV'; 'SV_RV'}, ...
    [SV_Qs_mL; SV_Qp_mL; SV_LV_mL; SV_RV_mL], ...
    {'mL/beat'; 'mL/beat'; 'mL/beat'; 'mL/beat'}, ...
    {'CO_Lmin * 1000 / HR'; ...
     'CO_Lmin * QpQs * 1000 / HR'; ...
     'LVEDV - LVESV'; ...
     'RVEDV - RVESV'}, ...
    'VariableNames', {'Estimate', 'Value', 'Unit', 'Formula'});

pair_names = {'SV_Qs_vs_SV_Qp'; 'SV_Qp_vs_SV_LV'; 'SV_Qp_vs_SV_RV'; ...
    'SV_Qs_vs_SV_LV'; 'SV_Qs_vs_SV_RV'; 'SV_LV_vs_SV_RV'};
pair_a = [SV_Qs_mL; SV_Qp_mL; SV_Qp_mL; SV_Qs_mL; SV_Qs_mL; SV_LV_mL];
pair_b = [SV_Qp_mL; SV_LV_mL; SV_RV_mL; SV_LV_mL; SV_RV_mL; SV_RV_mL];
rel_diff = nan(numel(pair_names), 1);
severity = cell(numel(pair_names), 1);
for idx = 1:numel(pair_names)
    rel_diff(idx) = relative_difference(pair_a(idx), pair_b(idx));
    severity{idx} = severity_from_difference(rel_diff(idx), config);
end
audit.pairwise_relative_differences = cell2struct(num2cell(rel_diff), pair_names, 1);
audit.pairwise_table = table(pair_names, pair_a, pair_b, rel_diff, severity, ...
    'VariableNames', {'Comparison', 'ValueA_mL', 'ValueB_mL', ...
    'RelativeDifference', 'Severity'});

valid_rel_diff = rel_diff(isfinite(rel_diff));
if isempty(valid_rel_diff)
    audit.severity = 'not_evaluable';
else
    audit.max_relative_difference = max(valid_rel_diff);
    audit.severity = severity_from_difference(audit.max_relative_difference, config);
end

rv_qp_diff = relative_difference(SV_Qp_mL, SV_RV_mL);
lv_qp_diff = relative_difference(SV_Qp_mL, SV_LV_mL);
audit.flags = {};
audit.notes = {};
audit.recommended_target_tier_changes = struct( ...
    'metric', {}, 'tier', {}, 'reason', {});

if isfinite(rv_qp_diff) && rv_qp_diff > config.strong_threshold
    audit.flags{end + 1} = 'RV_STROKE_VOLUME_CONFLICT';
    audit.notes{end + 1} = sprintf(['RVEDV/RVESV imply SV_RV=%.3f mL, ', ...
        'but catheter Qs*QpQs/HR implies SV_Qp=%.3f mL (relative difference %.1f%%).'], ...
        SV_RV_mL, SV_Qp_mL, 100 * rv_qp_diff);
    audit.recommended_target_tier_changes(end + 1) = struct( ...
        'metric', 'RVEDV', ...
        'tier', 'consistency_check_only', ...
        'reason', ['Echo-derived RVEDV conflicts with catheter-derived ', ...
        'pulmonary stroke volume; keep visible but exclude from hard fitting.']);
end

if isfinite(lv_qp_diff) && lv_qp_diff > config.strong_threshold
    audit.flags{end + 1} = 'LV_STROKE_VOLUME_CONFLICT';
    audit.notes{end + 1} = sprintf(['LVEDV/LVESV imply SV_LV=%.3f mL, ', ...
        'but catheter Qs*QpQs/HR implies SV_Qp=%.3f mL (relative difference %.1f%%).'], ...
        SV_LV_mL, SV_Qp_mL, 100 * lv_qp_diff);
    audit.recommended_target_tier_changes(end + 1) = struct( ...
        'metric', 'LVEDV', ...
        'tier', 'consistency_check_only', ...
        'reason', ['Echo-derived LVEDV conflicts with catheter-derived ', ...
        'pulmonary stroke volume; keep visible but exclude from hard fitting.']);
    audit.recommended_target_tier_changes(end + 1) = struct( ...
        'metric', 'LVESV', ...
        'tier', 'consistency_check_only', ...
        'reason', ['Echo-derived LVESV conflicts with catheter-derived ', ...
        'pulmonary stroke volume; keep visible but exclude from hard fitting.']);
end

if isempty(audit.flags)
    audit.flags = {'NONE'};
    audit.notes = {'No strong stroke-volume inconsistency was detected with the available data.'};
end

audit.summary = sprintf('Clinical consistency audit severity: %s; max relative SV difference = %.1f%%.', ...
    audit.severity, 100 * audit.max_relative_difference);
end

function config = default_audit_config()
% DEFAULT_AUDIT_CONFIG - explicit thresholds for clinical interpretability.
config = struct();
config.mild_threshold = 0.20;       % [-] relative difference above 20%
config.strong_threshold = 0.35;     % [-] relative difference above 35%
config.critical_threshold = 0.50;   % [-] relative difference above 50%
end

function audit = empty_audit(config)
audit = struct();
audit.config = config;
audit.HR_bpm = NaN;
audit.CO_Lmin = NaN;
audit.QpQs = NaN;
audit.SV_Qs = NaN;
audit.SV_Qp = NaN;
audit.SV_LV = NaN;
audit.SV_RV = NaN;
audit.max_relative_difference = NaN;
audit.severity = 'not_evaluable';
audit.flags = {'NOT_EVALUABLE'};
audit.notes = {};
audit.recommended_target_tier_changes = struct( ...
    'metric', {}, 'tier', {}, 'reason', {});
audit.stroke_volume_table = table();
audit.pairwise_table = table();
audit.pairwise_relative_differences = struct();
end

function severity = severity_from_difference(rel_diff, config)
if ~isfinite(rel_diff)
    severity = 'not_evaluable';
elseif rel_diff > config.critical_threshold
    severity = 'critical';
elseif rel_diff > config.strong_threshold
    severity = 'strong';
elseif rel_diff > config.mild_threshold
    severity = 'mild';
else
    severity = 'none';
end
end

function rel_diff = relative_difference(a, b)
if ~isfinite(a) || ~isfinite(b)
    rel_diff = NaN;
    return;
end
denom = max((abs(a) + abs(b)) / 2, 1e-9);
rel_diff = abs(a - b) / denom;
end

function value = field_or_nan(s, field_name)
if isstruct(s) && isfield(s, field_name) && isnumeric(s.(field_name)) && ...
        isscalar(s.(field_name)) && isfinite(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

function merged = merge_struct(defaults, overrides)
merged = defaults;
fields = fieldnames(overrides);
for idx = 1:numel(fields)
    merged.(fields{idx}) = overrides.(fields{idx});
end
end
