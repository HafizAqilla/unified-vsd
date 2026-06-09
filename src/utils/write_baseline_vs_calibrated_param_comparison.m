function csv_path = write_baseline_vs_calibrated_param_comparison( ...
        params_baseline, params_calibrated, calib_out, tables_dir, scenario, warm_start, clinical)
% WRITE_BASELINE_VS_CALIBRATED_PARAM_COMPARISON
% -----------------------------------------------------------------------
% Generates CSV and Excel side-by-side comparison tables of parameters
% across Adult Baseline, Demographic Scaled, PreOp Calibrated,
% PostOp Baseline, and PostOp Calibrated states.
%
% INPUTS:
%   params_baseline    - post-surgery baseline parameter struct           [-]
%   params_calibrated  - post-surgery calibrated parameter struct         [-]
%   calib_out          - calibration diagnostic structure                 [-]
%   tables_dir         - directory to write reporting tables              [-]
%   scenario           - current run scenario ('post_surgery')            [-]
%   warm_start         - warm-start audit structure                       [-]
%   clinical           - clinical targets and pre-op parameter seed       [-]
%
% OUTPUTS:
%   csv_path           - path to the generated comparison CSV file        [-]
%
% REFERENCES:
%   [1] apply_post_surgery_warm_start.m
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-06-06
% VERSION:  2.0
% -----------------------------------------------------------------------

csv_path = '';

if nargin < 6 || isempty(warm_start)
    warm_start = struct('applied', false, 'copiedFields', {{}});
end
if nargin < 7
    clinical = struct();
end

%% ---- Reconstruct Adult and Demographic Scaled parameters ------------
params_ref = default_parameters();

patient = struct();
if isfield(clinical, 'common') && isfield(clinical.common, 'weight_kg')
    patient.age_years = clinical.common.age_years;
    patient.age_days = clinical.common.age_years * 365.25;
    patient.weight_kg = clinical.common.weight_kg;
    patient.height_cm = clinical.common.height_cm;
    patient.sex = clinical.common.sex;
    patient.maturation_mode = 'normal';
    patient.scaling_mode = resolve_reporting_scaling_mode( ...
        clinical, params_baseline, 'lundquist_bsa');
elseif isfield(params_baseline, 'scaling') && isfield(params_baseline.scaling, 'patient')
    patient = params_baseline.scaling.patient;
else
    patient.age_years = 3.0;
    patient.age_days = 3.0 * 365.25;
    patient.weight_kg = 14.0;
    patient.height_cm = 95.0;
    patient.sex = 0;
    patient.maturation_mode = 'normal';
    patient.scaling_mode = 'lundquist_bsa';
end

params_scaled = apply_scaling(params_ref, patient);

%% ---- Retrieve PreOp Calibrated parameters ---------------------------
params_pre_cal = struct();
if isfield(clinical, 'pre_surgery') && isfield(clinical.pre_surgery, 'CalibParams')
    params_pre_cal = clinical.pre_surgery.CalibParams;
end

%% ---- Extract all parameter structures to flat cell lists ------------
[names_ref, vals_ref] = extract_scalars(params_ref, '');
[names_scaled, vals_scaled] = extract_scalars(params_scaled, '');
[names_pre_cal, vals_pre_cal] = extract_scalars(params_pre_cal, '');
[names_post_base, vals_post_base] = extract_scalars(params_baseline, '');
[names_post_cal, vals_post_cal] = extract_scalars(params_calibrated, '');

%% ---- Build Map containers for fast value lookups -------------------
map_ref = containers.Map(names_ref, vals_ref);
map_scaled = containers.Map(names_scaled, vals_scaled);
map_pre_cal = containers.Map(names_pre_cal, vals_pre_cal);
map_post_base = containers.Map(names_post_base, vals_post_base);
map_post_cal = containers.Map(names_post_cal, vals_post_cal);

%% ---- Collect all unique parameter keys ------------------------------
all_keys = unique([keys(map_ref), keys(map_scaled), keys(map_pre_cal), ...
    keys(map_post_base), keys(map_post_cal)]);

%% ---- Filter and populate columns ------------------------------------
Parameter = {};
Adult_Baseline = [];
Demographic_Scaled = [];
PreOp_Calibrated = [];
PostOp_Baseline = [];
PostOp_Calibrated = [];

for k_idx = 1:numel(all_keys)
    key = all_keys{k_idx};
    
    % Skip metadata, internal indices, timing-internal arrays, and overrides
    if startsWith(key, 'idx.') || startsWith(key, 'conv.') || ...
       startsWith(key, 'sim.') || startsWith(key, 'scaling.patient.') || ...
       startsWith(key, 'scaling.zhang_exponents.') || startsWith(key, 'scaling.lundquist_exponents.') || ...
       contains(key, 'clinical_override.') || startsWith(key, 'ic.V')
        continue;
    end
    
    Parameter{end+1} = key; %#ok<AGROW>
    
    if isKey(map_ref, key)
        Adult_Baseline(end+1) = map_ref(key); %#ok<AGROW>
    else
        Adult_Baseline(end+1) = NaN; %#ok<AGROW>
    end
    
    if isKey(map_scaled, key)
        Demographic_Scaled(end+1) = map_scaled(key); %#ok<AGROW>
    else
        Demographic_Scaled(end+1) = NaN; %#ok<AGROW>
    end
    
    if isKey(map_pre_cal, key)
        PreOp_Calibrated(end+1) = map_pre_cal(key); %#ok<AGROW>
    else
        PreOp_Calibrated(end+1) = NaN; %#ok<AGROW>
    end
    
    if isKey(map_post_base, key)
        PostOp_Baseline(end+1) = map_post_base(key); %#ok<AGROW>
    else
        PostOp_Baseline(end+1) = NaN; %#ok<AGROW>
    end
    
    if isKey(map_post_cal, key)
        PostOp_Calibrated(end+1) = map_post_cal(key); %#ok<AGROW>
    else
        PostOp_Calibrated(end+1) = NaN; %#ok<AGROW>
    end
end

% Transpose to column vectors
Parameter = Parameter(:);
Adult_Baseline = Adult_Baseline(:);
Demographic_Scaled = Demographic_Scaled(:);
PreOp_Calibrated = PreOp_Calibrated(:);
PostOp_Baseline = PostOp_Baseline(:);
PostOp_Calibrated = PostOp_Calibrated(:);

% Compute percentage changes of interest
PctChange_Adult_to_Scaled = (Demographic_Scaled - Adult_Baseline) ./ abs(Adult_Baseline) * 100;
PctChange_Scaled_to_PreCal = (PreOp_Calibrated - Demographic_Scaled) ./ abs(Demographic_Scaled) * 100;
PctChange_PostBase_to_PostCal = (PostOp_Calibrated - PostOp_Baseline) ./ abs(PostOp_Baseline) * 100;

%% ---- Create output Table -------------------------------------------
T = table(Parameter, Adult_Baseline, Demographic_Scaled, PreOp_Calibrated, ...
    PostOp_Baseline, PostOp_Calibrated, ...
    PctChange_Adult_to_Scaled, PctChange_Scaled_to_PreCal, PctChange_PostBase_to_PostCal, ...
    'VariableNames', { ...
        'Parameter', ...
        'Adult_Baseline', ...
        'Demographic_Scaled', ...
        'PreOp_Calibrated', ...
        'PostOp_Baseline', ...
        'PostOp_Calibrated', ...
        'PctChange_Adult_to_Scaled', ...
        'PctChange_Scaled_to_PreCal', ...
        'PctChange_PostBase_to_PostCal'});

% Sort alphabetically by parameter path name
T = sortrows(T, 'Parameter');

%% ---- Write output files ---------------------------------------------
if ~exist(tables_dir, 'dir')
    mkdir(tables_dir);
end

xls_file = fullfile(tables_dir, ...
    sprintf('param_comparison_detailed_%s.xlsx', scenario));
csv_path = fullfile(tables_dir, ...
    sprintf('param_comparison_detailed_%s.csv', scenario));

try
    writetable(T, xls_file);
    fprintf('[param_comparison] Excel saved to: %s\n', xls_file);
catch ME
    warning('write_baseline_vs_calibrated_param_comparison:xlsError', ...
        'Failed to save Excel file: %s', ME.message);
end

try
    writetable(T, csv_path);
    fprintf('[param_comparison] CSV saved to: %s\n', csv_path);
catch ME
    warning('write_baseline_vs_calibrated_param_comparison:csvError', ...
        'Failed to save CSV file: %s', ME.message);
end

end

%% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function [names, vals] = extract_scalars(s, prefix)
% EXTRACT_SCALARS - recursive extractor for structure numeric fields.
names = {};
vals = {};
if ~isstruct(s) || isempty(fieldnames(s))
    return;
end

fields = fieldnames(s);
for i = 1:numel(fields)
    name = fields{i};
    val = s.(name);
    
    if isempty(prefix)
        path = name;
    else
        path = sprintf('%s.%s', prefix, name);
    end
    
    if isstruct(val)
        [sub_names, sub_vals] = extract_scalars(val, path);
        names = [names; sub_names(:)]; %#ok<AGROW>
        vals = [vals; sub_vals(:)]; %#ok<AGROW>
    elseif isnumeric(val) && isscalar(val) && isfinite(val)
        names{end+1, 1} = path; %#ok<AGROW>
        vals{end+1, 1} = double(val); %#ok<AGROW>
    end
end
end

function scaling_mode = resolve_reporting_scaling_mode(clinical, params_baseline, default_mode)
% RESOLVE_REPORTING_SCALING_MODE - preserve the run's demographic scaling mode.
scaling_mode = default_mode;
if isfield(params_baseline, 'scaling') && isfield(params_baseline.scaling, 'mode') && ...
        ~isempty(params_baseline.scaling.mode)
    scaling_mode = params_baseline.scaling.mode;
elseif isfield(params_baseline, 'scaling') && ...
        isfield(params_baseline.scaling, 'requested_mode') && ...
        ~isempty(params_baseline.scaling.requested_mode)
    scaling_mode = params_baseline.scaling.requested_mode;
elseif isfield(clinical, 'common') && isfield(clinical.common, 'scaling_mode') && ...
        ~isempty(clinical.common.scaling_mode)
    scaling_mode = clinical.common.scaling_mode;
end

scaling_mode = lower(strtrim(char(scaling_mode)));
if strcmp(scaling_mode, 'lundquist') || strcmp(scaling_mode, 'lundqvist')
    scaling_mode = 'lundquist_bsa';
end
end
