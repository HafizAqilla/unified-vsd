function protocol = load_run_protocol(protocol_ref, varargin)
% LOAD_RUN_PROTOCOL  Load and normalize a publication/exploratory protocol.

opts = parse_options(varargin{:});
if isstruct(protocol_ref)
    protocol = protocol_ref;
elseif ischar(protocol_ref) || isstring(protocol_ref)
    protocol_ref = char(protocol_ref);
    if isfile(protocol_ref)
        [folder, name] = fileparts(protocol_ref);
        addpath(folder);
        cleanup = onCleanup(@() rmpath(folder)); %#ok<NASGU>
        protocol_fn = str2func(name);
        protocol = protocol_fn();
    else
        protocol_fn = str2func(protocol_ref);
        protocol = protocol_fn();
    end
else
    error('load_run_protocol:badInput', ...
        'Protocol must be a struct, function name, or .m file path.');
end

protocol = apply_protocol_defaults(protocol, opts.RootDir);
end

function protocol = apply_protocol_defaults(protocol, root_dir)
if ~isfield(protocol, 'ProtocolID') || isempty(protocol.ProtocolID)
    protocol.ProtocolID = 'unnamed_protocol';
end
if ~isfield(protocol, 'Mode') || isempty(protocol.Mode)
    protocol.Mode = 'exploratory';
end
if ~isfield(protocol, 'ScalingMode')
    protocol.ScalingMode = '';
end
policy = resolve_scaling_policy(protocol.ScalingMode, protocol.Mode, ...
    'ScalingRole', get_field_or(protocol, 'ScalingRole', ''), ...
    'OverrideRationale', get_field_or(protocol, 'ScalingOverrideRationale', ''));

protocol.ScalingMode = policy.ScalingMode;
protocol.ScalingRole = policy.ScalingRole;
protocol.ScalingCitation = policy.ScalingCitation;
protocol.ImplementationVariant = policy.ImplementationVariant;
protocol.DeviationFromCitation = policy.DeviationFromCitation;

protocol.RootDir = root_dir;
protocol.AllowSeedFallback = logical(get_field_or(protocol, 'AllowSeedFallback', ...
    ~strcmpi(protocol.Mode, 'publication')));
protocol.PreSeedPath = char(string(get_field_or(protocol, 'PreSeedPath', '')));
protocol.RunPreCalibration = logical(get_field_or(protocol, 'RunPreCalibration', true));
protocol.RunPostPrediction = logical(get_field_or(protocol, 'RunPostPrediction', true));
protocol.RunPostCalibration = logical(get_field_or(protocol, 'RunPostCalibration', false));
protocol.PreScenario = char(string(get_field_or(protocol, 'PreScenario', 'pre_surgery')));
protocol.PostScenario = char(string(get_field_or(protocol, 'PostScenario', 'post_surgery')));
end

function opts = parse_options(varargin)
parser = inputParser();
addParameter(parser, 'RootDir', pwd, @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.RootDir = char(opts.RootDir);
end

function value = get_field_or(s, field_name, fallback)
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    value = s.(field_name);
else
    value = fallback;
end
end

