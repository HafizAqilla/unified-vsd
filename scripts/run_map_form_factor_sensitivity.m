function record = run_map_form_factor_sensitivity(form_factor, varargin)
% RUN_MAP_FORM_FACTOR_SENSITIVITY
% -----------------------------------------------------------------------
% Runs one pre-surgery calibration arm with the systemic mean pressure
% recomputed under an explicit MAP form factor, and reports the resulting
% acceptance gate.
%
% Why this exists
% ---------------
% Reyna's SAP_mean is not an independent measurement. The catheter recorded
% SAP_sys = 100 and SAP_dia = 57; SAP_mean = 71.3 was then COMPUTED as
% dia + (sys - dia)/3, the resting-adult one-third rule.
%
% The model does not use that rule. It reports SAP_mean as a true time
% average of its arterial waveform, whose form factor is a model output, not
% an assumption. Observed on 2026-08-28: the model's own form factor was
% 0.460 while the clinical target assumed 0.333.
%
% The consequence is that sys = 100, dia = 57 and mean = 71.3 cannot all be
% satisfied simultaneously by a waveform whose mean is a true time average.
% Pinned to 71.3, the model lowers the whole systemic waveform, which pushes
% SAP_min far below the measured diastolic pressure. That is a target
% definition conflict, not a model deficiency, and it must be reported as a
% sensitivity rather than tuned away.
%
% The one-third rule assumes roughly 60-70 bpm. As heart rate rises, systole
% occupies a larger fraction of the cycle and the true form factor increases,
% so the rule progressively underestimates MAP. Reyna's heart rate is
% 119 bpm.
%
% INPUTS:
%   form_factor - MAP form factor k in MAP = dia + k*(sys - dia)      [-]
%                 1/3 reproduces the value recorded in the protocol.
%   varargin    - optional name/value pairs:
%       'Scenario'   scenario string, default 'pre_surgery'
%       'ScalingMode' scaling prior, default 'zhang'
%
% OUTPUTS:
%   record - struct with the applied target, run folder, and gate counts
%
% USAGE:
%   addpath(genpath(pwd));
%   record = run_map_form_factor_sensitivity(0.40);
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 6)
%   [2] config/patient_reyna.m lines 64-71 (catheter RFA source values)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------

parser = inputParser;
parser.FunctionName = mfilename;
addRequired(parser, 'form_factor', ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0 && x < 1);
addParameter(parser, 'Scenario', 'pre_surgery', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ScalingMode', 'zhang', @(x) ischar(x) || isstring(x));
parse(parser, form_factor, varargin{:});
opts = parser.Results;
scenario = char(opts.Scenario);

clinical = patient_reyna();
src = clinical.(scenario);

sys_mmHg = src.SAP_sys_mmHg;
dia_mmHg = src.SAP_dia_mmHg;
if ~isfinite(sys_mmHg) || ~isfinite(dia_mmHg)
    error('run_map_form_factor_sensitivity:missingCuffPressures', ...
        'Systolic and diastolic catheter pressures are required.');
end

original_mean = src.SAP_mean_mmHg;
applied_mean = dia_mmHg + form_factor * (sys_mmHg - dia_mmHg);
clinical.(scenario).SAP_mean_mmHg = applied_mean;

fprintf('\n=== MAP FORM FACTOR SENSITIVITY ===\n');
fprintf('  Catheter sys/dia   : %.1f / %.1f mmHg (pulse %.1f)\n', ...
    sys_mmHg, dia_mmHg, sys_mmHg - dia_mmHg);
fprintf('  Form factor applied: %.3f\n', form_factor);
fprintf('  SAP_mean target    : %.2f mmHg (protocol value %.2f)\n', ...
    applied_mean, original_mean);
fprintf('  Heart rate         : %.0f bpm\n', clinical.common.HR);
fprintf(['  NOTE: SAP_mean is derived from the two catheter pressures, ', ...
    'not measured\n        independently. This arm varies the estimator, ', ...
    'not the measurement.\n']);

setenv('UNIFIED_VSD_SCALING_MODE', char(opts.ScalingMode));
setenv('UNIFIED_VSD_DISABLE_HISTORICAL_SEEDS', '1');

main_run(scenario, clinical);

record = collect_gate_record(scenario, form_factor, applied_mean, original_mean);
end

% =========================================================================
function record = collect_gate_record(scenario, form_factor, applied_mean, original_mean)
% COLLECT_GATE_RECORD - read the newest exported gate table for this run.
record = struct( ...
    'form_factor', form_factor, ...
    'sap_mean_target', applied_mean, ...
    'sap_mean_protocol', original_mean, ...
    'run_dir', '', ...
    'governed_within_gate', NaN, ...
    'governed_total', NaN, ...
    'metrics_within_gate', NaN, ...
    'metrics_total', NaN, ...
    'worst_metric', '', ...
    'worst_abs_error_pct', NaN);

pattern = fullfile('results', 'runs', '**', ...
    sprintf('full_metric_gate_%s.csv', scenario));
found = dir(pattern);
if isempty(found)
    warning('run_map_form_factor_sensitivity:noGateTable', ...
        'No exported metric gate table was found.');
    return;
end
[~, newest] = max([found.datenum]);
gate_path = fullfile(found(newest).folder, found(newest).name);
record.run_dir = found(newest).folder;

tbl = readtable(gate_path);
governed = logical(tbl.InPrimaryRMSE);
within = logical(tbl.WithinGate);
record.governed_total = nnz(governed);
record.governed_within_gate = nnz(governed & within);
record.metrics_total = height(tbl);
record.metrics_within_gate = nnz(within);

[worst_pct, worst_ix] = max(tbl.AbsError_pct);
if isfinite(worst_pct)
    record.worst_metric = tbl.Metric{worst_ix};
    record.worst_abs_error_pct = worst_pct;
end

fprintf('\n=== SENSITIVITY RESULT (form factor %.3f) ===\n', form_factor);
fprintf('  Governed gate : %d of %d within band\n', ...
    record.governed_within_gate, record.governed_total);
fprintf('  All targets   : %d of %d within band\n', ...
    record.metrics_within_gate, record.metrics_total);
fprintf('  Worst metric  : %s (%.2f%%)\n', ...
    record.worst_metric, record.worst_abs_error_pct);
end
