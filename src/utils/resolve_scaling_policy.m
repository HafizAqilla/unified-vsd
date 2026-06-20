function policy = resolve_scaling_policy(scaling_mode, run_mode, varargin)
% RESOLVE_SCALING_POLICY
% -----------------------------------------------------------------------
% Resolves scaling mode, role, citation, and implementation caveats.
%
% Publication mode defaults to Zhang as the primary prior. Lundquist-BSA is
% retained as comparator unless an explicit override-primary rationale is
% supplied.
% -----------------------------------------------------------------------

opts = parse_policy_options(varargin{:});

if nargin < 2 || isempty(run_mode)
    run_mode = getenv('UNIFIED_VSD_RUN_MODE');
end
if isempty(run_mode)
    run_mode = getenv('UNIFIED_VSD_PROTOCOL_MODE');
end
if isempty(run_mode)
    run_mode = 'exploratory';
end
run_mode = lower(strtrim(char(run_mode)));

if nargin < 1 || isempty(scaling_mode)
    scaling_mode = getenv('UNIFIED_VSD_SCALING_MODE');
end
if isempty(scaling_mode)
    if strcmp(run_mode, 'publication')
        scaling_mode = 'zhang';
    else
        scaling_mode = 'lundquist_bsa';
    end
end
scaling_mode = normalize_scaling_mode(scaling_mode);

registry = scaling_method_registry();
if ~isfield(registry, scaling_mode)
    error('resolve_scaling_policy:unknownMode', ...
        'Unsupported scaling mode: %s', scaling_mode);
end
entry = registry.(scaling_mode);

if isempty(opts.ScalingRole)
    if strcmp(run_mode, 'publication')
        scaling_role = entry.DefaultPublicationRole;
    else
        scaling_role = entry.DefaultExploratoryRole;
    end
else
    scaling_role = lower(strtrim(char(opts.ScalingRole)));
end

override_rationale = strtrim(char(opts.OverrideRationale));
if strcmp(run_mode, 'publication') && strcmp(scaling_mode, 'lundquist_bsa') && ...
        any(strcmp(scaling_role, {'primary_prior', 'override_primary'}))
    if isempty(override_rationale)
        error('resolve_scaling_policy:missingOverrideRationale', ...
            ['Lundquist-BSA cannot be primary in publication mode without ' ...
             'an explicit override rationale.']);
    end
    scaling_role = 'override_primary';
end

policy = entry;
policy.RunMode = run_mode;
policy.ScalingMode = scaling_mode;
policy.ScalingRole = scaling_role;
policy.ScalingCitation = build_citation_text(entry);
policy.OverrideRationale = override_rationale;
policy.AllowedDownstream = true;
policy.BaselineGateStatus = 'not_evaluated';

end

function opts = parse_policy_options(varargin)
parser = inputParser();
addParameter(parser, 'ScalingRole', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OverrideRationale', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
end

function mode = normalize_scaling_mode(mode)
mode = lower(strtrim(char(mode)));
switch mode
    case {'lundquist', 'lundqvist'}
        mode = 'lundquist_bsa';
end
end

function text = build_citation_text(entry)
parts = {entry.CitationLabel};
if isfield(entry, 'DOI') && ~isempty(entry.DOI)
    parts{end + 1} = ['DOI:' entry.DOI]; %#ok<AGROW>
end
if isfield(entry, 'PMID') && ~isempty(entry.PMID)
    parts{end + 1} = ['PMID:' entry.PMID]; %#ok<AGROW>
end
text = strjoin(parts, '; ');
end

