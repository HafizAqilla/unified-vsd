function seed_resolution = resolve_pre_to_post_seed(protocol, root_dir)
% RESOLVE_PRE_TO_POST_SEED  Resolve pinned pre-op seed with safe fallback.

if nargin < 2 || isempty(root_dir)
    root_dir = pwd;
end
protocol = load_run_protocol(protocol, 'RootDir', root_dir);

seed_resolution = struct();
seed_resolution.ProtocolID = protocol.ProtocolID;
seed_resolution.Mode = protocol.Mode;
seed_resolution.AllowSeedFallback = protocol.AllowSeedFallback;
seed_resolution.ResolvedPath = '';
seed_resolution.Source = '';
seed_resolution.Status = 'unresolved';
seed_resolution.Warning = '';
seed_resolution.Metadata = struct();

if ~isempty(protocol.PreSeedPath)
    seed_path = char(protocol.PreSeedPath);
    if ~isfile(seed_path) && ~isabsolute_path(seed_path)
        seed_path = fullfile(root_dir, seed_path);
    end
    if ~isfile(seed_path)
        error('resolve_pre_to_post_seed:missingPinnedSeed', ...
            'Pinned PreSeedPath does not exist: %s', seed_path);
    end
    seed_resolution = fill_seed_resolution(seed_resolution, seed_path, 'pinned');
    return;
end

if strcmpi(protocol.Mode, 'publication') && ~protocol.AllowSeedFallback
    error('resolve_pre_to_post_seed:publicationSeedRequired', ...
        'Publication mode requires an explicit PreSeedPath.');
end

if ~protocol.AllowSeedFallback
    error('resolve_pre_to_post_seed:fallbackDisabled', ...
        'No PreSeedPath supplied and fallback is disabled.');
end

seed_path = find_latest_seed(root_dir, get_field_or(protocol, 'PatientLabel', ''));
if isempty(seed_path)
    error('resolve_pre_to_post_seed:noFallbackSeed', ...
        'No pre-to-post seed found under %s.', root_dir);
end
seed_resolution = fill_seed_resolution(seed_resolution, seed_path, 'latest_fallback');
seed_resolution.Warning = ['Fallback seed discovery was used. Publication ' ...
    'runs should pin PreSeedPath explicitly.'];
end

function seed_resolution = fill_seed_resolution(seed_resolution, seed_path, source)
seed_resolution.ResolvedPath = seed_path;
seed_resolution.Source = source;
seed_resolution.Status = 'resolved';
seed_resolution.Metadata = read_seed_metadata(seed_path);
end

function metadata = read_seed_metadata(seed_path)
metadata = struct('File', seed_path);
try
    vars = whos('-file', seed_path);
    metadata.Variables = {vars.name};
    if any(strcmp({vars.name}, 'pre_to_post_seed'))
        loaded = load(seed_path, 'pre_to_post_seed');
        seed = loaded.pre_to_post_seed;
        metadata.SourceScenario = get_field_or(seed, 'source_scenario', '');
        metadata.Timestamp = get_field_or(seed, 'timestamp', '');
        metadata.CalibrationStatus = get_nested_field_or(seed, ...
            {'calibration_status','label'}, '');
        if isfield(seed, 'scaling_policy')
            metadata.ScalingMode = get_field_or(seed.scaling_policy, 'ScalingMode', '');
            metadata.ScalingRole = get_field_or(seed.scaling_policy, 'ScalingRole', '');
        end
    end
catch ME
    metadata.ReadWarning = ME.message;
end
end

function seed_path = find_latest_seed(root_dir, patient_label)
patterns = {};
if ~isempty(patient_label)
    patterns{end + 1} = fullfile(root_dir, 'results', 'runs', ...
        sprintf('*_%s_pre_surgery', patient_label), 'mat', 'pre_to_post_seed_latest.mat');
    patterns{end + 1} = fullfile(root_dir, 'results', 'runs', ...
        sprintf('*_%s_pre_surgery', patient_label), 'mat', 'pre_to_post_seed_*.mat');
end
patterns{end + 1} = fullfile(root_dir, 'results', 'runs', ...
    '*_pre_surgery', 'mat', 'pre_to_post_seed_latest.mat');
patterns{end + 1} = fullfile(root_dir, 'results', 'runs', ...
    '*_pre_surgery', 'mat', 'pre_to_post_seed_*.mat');

matches = [];
for idx = 1:numel(patterns)
    matches = [matches; dir(patterns{idx})]; %#ok<AGROW>
end
if isempty(matches)
    seed_path = '';
    return;
end
[~, order] = sort([matches.datenum], 'descend');
matches = matches(order);
seed_path = fullfile(matches(1).folder, matches(1).name);
end

function tf = isabsolute_path(path_value)
path_value = char(path_value);
tf = startsWith(path_value, filesep) || ...
    (~isempty(regexp(path_value, '^[A-Za-z]:[\\/]', 'once')));
end

function value = get_field_or(s, field_name, fallback)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    value = s.(field_name);
else
    value = fallback;
end
end

function value = get_nested_field_or(s, fields, fallback)
value = fallback;
cursor = s;
for idx = 1:numel(fields)
    if isstruct(cursor) && isfield(cursor, fields{idx})
        cursor = cursor.(fields{idx});
    else
        return;
    end
end
if ~isempty(cursor)
    value = cursor;
end
end

