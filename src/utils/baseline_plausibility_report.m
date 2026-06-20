function report = baseline_plausibility_report(metrics, clinical, scenario, case_profile, scaling_policy, varargin)
% BASELINE_PLAUSIBILITY_REPORT
% -----------------------------------------------------------------------
% Screens baseline model outputs before GSA/calibration.
% -----------------------------------------------------------------------

opts = parse_options(varargin{:});
if nargin < 5 || isempty(scaling_policy)
    scaling_policy = resolve_scaling_policy('', 'exploratory');
end
if nargin < 4
    case_profile = struct();
end

ranges = clinical_reference_ranges(scenario, clinical, case_profile);
n_rows = height(ranges);
value_col = nan(n_rows, 1);
status_col = strings(n_rows, 1);
message_col = strings(n_rows, 1);

for idx = 1:n_rows
    metric = ranges.Metric{idx};
    if isstruct(metrics) && isfield(metrics, metric) && isfinite(metrics.(metric))
        value = metrics.(metric);
        value_col(idx) = value;
        if value < ranges.HardLow(idx) || value > ranges.HardHigh(idx)
            status_col(idx) = "hard_fail";
            message_col(idx) = sprintf('%s %.4g outside hard range [%.4g, %.4g].', ...
                metric, value, ranges.HardLow(idx), ranges.HardHigh(idx));
        elseif value < ranges.SoftLow(idx) || value > ranges.SoftHigh(idx)
            status_col(idx) = "soft_warning";
            message_col(idx) = sprintf('%s %.4g outside soft range [%.4g, %.4g].', ...
                metric, value, ranges.SoftLow(idx), ranges.SoftHigh(idx));
        else
            status_col(idx) = "pass";
            message_col(idx) = sprintf('%s within soft range.', metric);
        end
    else
        status_col(idx) = "not_available";
        message_col(idx) = sprintf('%s not available in baseline metrics.', metric);
    end
end

table_out = ranges;
table_out.Value = value_col;
table_out.Status = status_col;
table_out.Message = message_col;
table_out.ScalingMode = repmat(string(scaling_policy.ScalingMode), n_rows, 1);
table_out.ScalingRole = repmat(string(scaling_policy.ScalingRole), n_rows, 1);
table_out.ImplementationVariant = repmat( ...
    string(scaling_policy.ImplementationVariant), n_rows, 1);

has_hard = any(status_col == "hard_fail");
has_soft = any(status_col == "soft_warning");
if has_hard && any(strcmp(scaling_policy.ScalingRole, {'comparator', 'exploratory'}))
    gate_status = 'comparator_only';
elseif has_hard
    gate_status = 'hard_fail';
elseif has_soft
    gate_status = 'soft_warning';
else
    gate_status = 'pass';
end

publication_mode = strcmp(get_policy_field(scaling_policy, 'RunMode', ''), 'publication');
allowed_downstream = ~(publication_mode && any(strcmp(gate_status, ...
    {'hard_fail', 'comparator_only'})));

report = struct();
report.status = gate_status;
report.allowedDownstream = allowed_downstream;
report.table = table_out;
report.failedMetrics = cellstr(table_out.Metric(status_col == "hard_fail"));
report.warningMetrics = cellstr(table_out.Metric(status_col == "soft_warning"));
report.scalingPolicy = scaling_policy;
report.summary = sprintf('Baseline plausibility: %s (%d hard fail, %d soft warning).', ...
    gate_status, nnz(status_col == "hard_fail"), nnz(status_col == "soft_warning"));

if ~isempty(opts.ResultsDir)
    if ~exist(opts.ResultsDir, 'dir')
        mkdir(opts.ResultsDir);
    end
    report.tableFile = fullfile(opts.ResultsDir, 'baseline_plausibility.csv');
    writetable(table_out, report.tableFile);
else
    report.tableFile = '';
end
end

function opts = parse_options(varargin)
parser = inputParser();
addParameter(parser, 'ResultsDir', '', @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});
opts = parser.Results;
opts.ResultsDir = char(opts.ResultsDir);
end

function value = get_policy_field(policy, field_name, fallback)
if isstruct(policy) && isfield(policy, field_name) && ~isempty(policy.(field_name))
    value = policy.(field_name);
else
    value = fallback;
end
end

