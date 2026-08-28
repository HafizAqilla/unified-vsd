function result = run_pre_post_protocol(protocol_ref, varargin)
% RUN_PRE_POST_PROTOCOL  Dry-run friendly pre/post publication orchestrator.

root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root_dir));

opts = parse_options(varargin{:});
protocol = load_run_protocol(protocol_ref, 'RootDir', root_dir);
protocol = apply_runtime_overrides(protocol, opts);

seed_resolution = struct('Status', 'not_requested');
if protocol.RunPostPrediction || protocol.RunPostCalibration
    try
        seed_resolution = resolve_pre_to_post_seed(protocol, root_dir);
    catch ME
        if opts.DryRun
            seed_resolution = struct();
            seed_resolution.ProtocolID = protocol.ProtocolID;
            seed_resolution.Mode = protocol.Mode;
            seed_resolution.AllowSeedFallback = protocol.AllowSeedFallback;
            seed_resolution.ResolvedPath = '';
            seed_resolution.Source = '';
            seed_resolution.Status = 'missing_required';
            seed_resolution.Warning = ME.message;
            seed_resolution.ErrorIdentifier = ME.identifier;
        else
            rethrow(ME);
        end
    end
end

out_dir = fullfile(root_dir, 'results', 'protocols', protocol.ProtocolID);
manifest_file = write_protocol_manifest(out_dir, protocol, seed_resolution);

result = struct();
result.Protocol = protocol;
result.SeedResolution = seed_resolution;
result.ManifestFile = manifest_file;
result.DryRun = opts.DryRun;

fprintf('[run_pre_post_protocol] Protocol: %s\n', protocol.ProtocolID);
fprintf('[run_pre_post_protocol] Mode: %s | Scaling: %s (%s)\n', ...
    protocol.Mode, protocol.ScalingMode, protocol.ScalingRole);
fprintf('[run_pre_post_protocol] Manifest: %s\n', manifest_file);
if isfield(seed_resolution, 'ResolvedPath') && ~isempty(seed_resolution.ResolvedPath)
    fprintf('[run_pre_post_protocol] Seed: %s\n', seed_resolution.ResolvedPath);
elseif isfield(seed_resolution, 'Warning') && ~isempty(seed_resolution.Warning)
    fprintf('[run_pre_post_protocol] Seed warning: %s\n', seed_resolution.Warning);
end

if opts.DryRun
    fprintf('[run_pre_post_protocol] Dry run only; no simulation launched.\n');
    return;
end

error('run_pre_post_protocol:notImplementedForHeavyRun', ...
    ['Heavy execution is intentionally not launched by this scaffold yet. ' ...
     'Use DryRun=true or implement the phase-specific execution calls.']);
end

function opts = parse_options(varargin)
parser = inputParser();
addParameter(parser, 'DryRun', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'RunPreCalibration', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(parser, 'RunPostPrediction', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(parser, 'RunPostCalibration', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.DryRun = logical(opts.DryRun);
end

function protocol = apply_runtime_overrides(protocol, opts)
fields = {'RunPreCalibration','RunPostPrediction','RunPostCalibration'};
for idx = 1:numel(fields)
    field_name = fields{idx};
    if ~isempty(opts.(field_name))
        protocol.(field_name) = logical(opts.(field_name));
    end
end
end
