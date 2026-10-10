function out_dir = export_run_supplement(run_dir, clinical, varargin)
% EXPORT_RUN_SUPPLEMENT
% -----------------------------------------------------------------------
% Writes every raw result of one pre-closure run as plain CSV/JSON into
% <run_dir>/supplement/, so the run can be audited and reproduced without
% MATLAB. Called immediately after each run finishes.
%
% FILES WRITTEN:
%   protocol_comparison_pre.csv   every pre-form protocol entry with its role
%                                 and the baseline/calibrated model value
%   protocol_comparison_post.csv  every post-form entry vs the zero-shot
%                                 post-closure prediction
%   inputs.csv                    values imposed on the model (not validated)
%   parameters_calibration.csv    all calibration candidates: scaled, seeded,
%                                 bounds, GSA selection, calibrated value
%   parameters_full_baseline.csv  every numeric model parameter, baseline
%   parameters_full_calibrated.csv  ... calibrated (accepted candidate)
%   parameters_full_post.csv      ... post-closure prediction state
%   gsa_indices_initial.csv       S1 and ST per output and input, before calibration
%   gsa_indices_final.csv         ... at the calibrated operating point
%   gsa_surrogate_quality.csv     PCE design size, degree, sparsity, LOO error
%   gsa_mask.csv                  screening result per candidate parameter
%   multistart.csv                every start: RMSE, failures, selection
%   identifiability_columns.csv   sensitivity-matrix column norms/correlations
%   identifiability_pairs.csv     pairwise column correlations
%   identifiability_reduction.csv best-conditioned subset per size
%   timeseries_pre_baseline.csv   last cardiac cycle, all states and signals
%   timeseries_pre_calibrated.csv
%   timeseries_post_prediction.csv
%   run_metadata.json             seed, status, fit summary, provenance
%
% INPUTS:
%   run_dir   - results/runs/<timestamp>_<patient>_pre_surgery            [-]
%   clinical  - clinical struct including clinical.protocol               [-]
%   'OutDir'  - override output folder (default <run_dir>/supplement)     [-]
% -----------------------------------------------------------------------

p = inputParser;
addParameter(p, 'OutDir', fullfile(run_dir, 'supplement'));
parse(p, varargin{:});
out_dir = p.Results.OutDir;
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

pkg = dir(fullfile(run_dir, 'mat', 'run_package_pre_surgery_*.mat'));
assert(~isempty(pkg), 'export_run_supplement:noPackage', 'No run package in %s', run_dir);
rp = load(fullfile(pkg(1).folder, pkg(1).name));
rp = rp.run_package;
ac = rp.accepted_candidate;

% ---- post-closure zero-shot prediction (HR from the post form) --------
pred = evaluate_post_closure_prediction( ...
    fullfile(run_dir, 'mat', 'params_accepted_candidate_pre_surgery.mat'), clinical, false);

% ---- protocol comparisons ---------------------------------------------
tiers = read_tiers(run_dir);
pre_tbl = protocol_table('pre', clinical, rp.metrics_baseline, ac.metrics, tiers);
% The post-closure state has more than one declared dataset (they disagree on
% ventricular volume); every one is evaluated and labelled, none overwrites
% another. The first variant is the primary one.
variants = post_variants(clinical);
post_parts = cell(1, numel(variants));
for v = 1:numel(variants)
    cv = clinical;
    cv.post_surgery = variants(v).data;
    t = protocol_table('post', cv, [], pred.metrics, table());
    t = addvars(t, repmat({variants(v).name}, height(t), 1), ...
        repmat({variants(v).provenance}, height(t), 1), 'Before', 1, ...
        'NewVariableNames', {'Dataset', 'DatasetProvenance'});
    post_parts{v} = t;
    if v == 1
        writetable(t, fullfile(out_dir, 'protocol_comparison_post.csv'));
    end
end
post_tbl = vertcat(post_parts{:});
writetable(pre_tbl, fullfile(out_dir, 'protocol_comparison_pre.csv'));
writetable(post_tbl, fullfile(out_dir, 'protocol_comparison_post_all_datasets.csv'));
writetable(inputs_table(clinical, ac.params, pred.params), fullfile(out_dir, 'inputs.csv'));

% ---- parameters ---------------------------------------------------------
reference = rp.params_baseline;
if isfield(rp, 'gsa_init_cfg') && isfield(rp.gsa_init_cfg, 'registry_context') && ...
        isfield(rp.gsa_init_cfg.registry_context, 'params_scaled')
    reference = rp.gsa_init_cfg.registry_context.params_scaled;
end
writetable(parameters_table(rp, ac.params, reference), ...
    fullfile(out_dir, 'parameters_calibration.csv'));
writetable(flatten_params(rp.params_baseline), fullfile(out_dir, 'parameters_full_baseline.csv'));
writetable(flatten_params(ac.params), fullfile(out_dir, 'parameters_full_calibrated.csv'));
writetable(flatten_params(pred.params), fullfile(out_dir, 'parameters_full_post.csv'));

% ---- global sensitivity analysis ---------------------------------------
[gi, qi] = gsa_tables(rp, 'gsa_init_out', 'gsa_init_cfg', 'before_calibration');
[gf, qf] = gsa_tables(rp, 'gsa_final_out', 'gsa_final_cfg', 'after_calibration');
writetable(gi, fullfile(out_dir, 'gsa_indices_initial.csv'));
writetable(gf, fullfile(out_dir, 'gsa_indices_final.csv'));
writetable([qi; qf], fullfile(out_dir, 'gsa_surrogate_quality.csv'));
writetable(mask_table(rp, gi), fullfile(out_dir, 'gsa_mask.csv'));

% ---- multi-start ---------------------------------------------------------
ms = rp.calib_out.multistart.table;
ms.Selected = ms.Start == rp.calib_out.multistart.best_start_index;
ms.MultistartSeed = repmat(rp.calib_out.multistart.seed, height(ms), 1);
writetable(ms, fullfile(out_dir, 'multistart.csv'));

% ---- identifiability -----------------------------------------------------
copy_if_exists(fullfile(run_dir, 'tables', 'parameter_identifiability_pre_surgery.csv'), ...
    fullfile(out_dir, 'identifiability_columns.csv'));
copy_if_exists(fullfile(run_dir, 'tables', 'parameter_identifiability_pairs_pre_surgery.csv'), ...
    fullfile(out_dir, 'identifiability_pairs.csv'));
red = analyse_parameter_reduction(run_dir, 'pre_surgery', 1);
writetable(reduction_table(red), fullfile(out_dir, 'identifiability_reduction.csv'));

% ---- time series (last cardiac cycle) -------------------------------------
writetable(last_cycle_table(rp.sim_baseline, rp.params_baseline), ...
    fullfile(out_dir, 'timeseries_pre_baseline.csv'));
writetable(last_cycle_table(ac.sim, ac.params), ...
    fullfile(out_dir, 'timeseries_pre_calibrated.csv'));
writetable(last_cycle_table(pred.sim, pred.params), ...
    fullfile(out_dir, 'timeseries_post_prediction.csv'));

% ---- metadata ---------------------------------------------------------------
meta = struct();
meta.run_folder = rp.run_folder;
meta.timestamp = rp.timestamp;
meta.multistart_seed = int64(rp.calib_out.multistart.seed);
meta.num_starts = rp.calib_out.multistart.num_starts;
meta.calibration_status = ac.calibration_status.label;
meta.n_calibrated_parameters = numel(rp.calib_out.names);
meta.calibrated_parameters = rp.calib_out.names;
meta.gsa_threshold = rp.sobol_ST_threshold;
meta.heart_rate_pre_bpm = ac.params.HR;
meta.heart_rate_post_bpm = pred.params.HR;
meta.scaling_mode = ac.params.scaling.mode;
meta.post_within_10pct = pred.within10;
meta.post_n_compared = pred.n;
meta.post_chi2_per_n = pred.chi2 / max(pred.n, 1);
meta.steady_state_pre_calibrated = logical(ac.sim.ss_reached);
meta.steady_state_post_prediction = logical(pred.sim.ss_reached);
meta.exported_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
write_json(fullfile(out_dir, 'run_metadata.json'), meta);

fprintf('[export_run_supplement] %d files written to %s\n', ...
    numel(dir(fullfile(out_dir, '*.*'))) - 2, out_dir);
end

% =======================================================================
function T = protocol_table(form, clinical, metrics_base, metrics_model, tiers)
% One row per protocol entry. Role says how the entry relates to the model.
% 'config' rows compare against the exact value the calibration/evaluation
% used (get_calibration_targets); ReadingsMean shows the raw form readings.
rec = clinical.protocol.(form);
spec = protocol_roles(form);
scenario = ternary(strcmp(form, 'pre'), 'pre_surgery', 'post_surgery');
targets = get_calibration_targets(scenario, clinical);
rows = {};
for i = 1:numel(rec)
    r = rec(i);
    key = sprintf('%d|%s', r.row, r.name);
    if ~isKey(spec, key)
        error('export_run_supplement:unmappedEntry', ...
            'Protocol entry "%s" (%s row %d) has no declared role.', r.name, form, r.row);
    end
    s = spec(key);
    rmean = NaN;
    if ~isempty(r.readings), rmean = mean(r.readings); end
    value = rmean;
    if strcmp(s.rule, 'config')
        value = target_value(targets, s.metric);
    end
    rows(end+1, :) = comparison_row(form, r.row, r.name, r.unit, r.method, s.role, ...
        tier_of(tiers, s.metric, s.role), readings_text(r.readings), rmean, value, s.rule, ...
        s.metric, metrics_base, metrics_model, r.note, s.note); %#ok<AGROW>
end
dd = derived_entries(form, clinical);
for k = 1:numel(dd)
    d = dd(k);
    rows(end+1, :) = comparison_row(form, NaN, d.name, d.unit, d.method, d.role, ...
        tier_of(tiers, d.metric, d.role), '', NaN, d.value, d.rule, d.metric, ...
        metrics_base, metrics_model, '', d.note); %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Form','ProtocolRow','Variable','Unit', ...
    'Method','Role','TierInRun','Readings','ReadingsMean','ComparedValue','ValueRule', ...
    'ModelMetric','ModelBaseline','ModelValue','ErrorPctBaseline','ErrorPct', ...
    'AbsErrorPct','Within10Pct','FormNote','RoleNote'});
end

function row = comparison_row(form, prow, name, unit, method, role, tier, readings, ...
    rmean, value, rule, metric, metrics_base, metrics_model, form_note, role_note)
y_base = metric_value(metrics_base, metric);
y_model = metric_value(metrics_model, metric);
err_base = 100 * (y_base - value) / abs(value);
err = 100 * (y_model - value) / abs(value);
within = double(abs(err) <= 10);
within(isnan(err)) = NaN;
row = {form, prow, name, unit, method, role, tier, readings, rmean, value, rule, ...
    metric, y_base, y_model, err_base, err, abs(err), within, form_note, role_note};
end

function v = target_value(targets, metric)
v = NaN;
idx = find(strcmp({targets.Metric}, metric), 1);
if ~isempty(idx), v = targets(idx).ClinicalValue; end
end

function tier = tier_of(tiers, metric, role)
% Tier describes how a target was fitted, so it only applies to fitted rows;
% evaluation-only rows sharing the same model output must not inherit it.
tier = '';
if ~strcmp(role, 'fitted'), return; end
if ~isempty(tiers) && ~isempty(metric) && any(strcmp(tiers.Metric, metric))
    tier = tiers.Tier{find(strcmp(tiers.Metric, metric), 1)};
end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function spec = protocol_roles(form)
% Declared role of every protocol entry. Roles:
%   input                 imposed on the model; not a validation comparison
%   fitted                calibration target (tier in TierInRun)
%   evaluation_only       compared, never fitted
%   predicted_evaluated   post-closure measurement vs zero-shot prediction
%   excluded_echo         echocardiographic entry excluded by study decision
%   excluded_timepoint    measured, but belongs to the other state
%   missing               blank on the form
spec = containers.Map();
add = @(row, name, role, metric, rule, note) put_role(spec, row, name, role, metric, rule, note);
if strcmp(form, 'pre')
    add(1, 'Age', 'input', '', 'readings', 'scaling input');
    add(3, 'Height', 'input', '', 'readings', 'scaling input');
    add(4, 'Weight', 'input', '', 'readings', 'scaling input');
    add(5, 'Body surface area', 'input', '', 'readings', 'reported; scaling is weight-based');
    add(6, 'Heart rate', 'input', '', 'readings', 'cycle length input');
    add(7, 'VSD diameter, LV side', 'input', '', 'readings', 'not used; RV side sets the orifice');
    add(7, 'VSD diameter, RV side', 'input', '', 'readings', 'orifice diameter = mean of range');
    add(9, 'LV systolic pressure', 'evaluation_only', 'LVP_max', 'mean_of_readings', '');
    add(9, 'RV systolic pressure', 'evaluation_only', 'RVP_max', 'mean_of_readings', '');
    add(10, 'Systolic pressure, NIBP', 'evaluation_only', 'SAP_max', 'mean_of_readings', 'non-invasive cuff; model has one systemic artery');
    add(11, 'Diastolic pressure, NIBP', 'evaluation_only', 'SAP_min', 'mean_of_readings', 'non-invasive cuff');
    add(12, 'Systolic pressure, RFA', 'fitted', 'SAP_max', 'config', '');
    add(13, 'Diastolic pressure, RFA', 'fitted', 'SAP_min', 'config', '');
    add(14, 'Systolic pressure, DAo', 'evaluation_only', 'SAP_max', 'mean_of_readings', 'second catheter site; same model output as RFA');
    add(15, 'Diastolic pressure, DAo', 'evaluation_only', 'SAP_min', 'mean_of_readings', 'second catheter site');
    add(16, 'PA systolic pressure', 'fitted', 'PAP_max', 'config', '');
    add(17, 'PA diastolic pressure', 'fitted', 'PAP_min', 'config', '');
    add(18, 'PA mean pressure', 'fitted', 'PAP_mean', 'config', 'config uses the modal reading');
    add(19, 'RA mean pressure', 'fitted', 'RAP_mean', 'config', '');
    add(20, 'Mean arterial pressure, form', 'evaluation_only', 'SAP_mean', 'mean_of_readings', 'form MAP from NIBP; fitted MAP is the derived RFA row');
    add(21, 'Pulmonary flow Qp', 'evaluation_only', 'Qp_Lmin', 'mean_of_readings', 'enters the fit only through derived Qs and shunt flow');
    add(22, 'Systemic flow Qs', 'missing', '', 'none', 'derived value listed separately');
    add(23, 'Qp/Qs', 'fitted', 'QpQs', 'config', '');
    add(24, 'PVR', 'missing', '', 'none', '');
    add(25, 'SVR', 'missing', '', 'none', 'derived value listed separately');
    add(26, 'LVEDV', 'excluded_timepoint', '', 'none', 'evaluated in post-closure table');
    add(27, 'LVESV', 'excluded_timepoint', '', 'none', 'evaluated in post-closure table');
    add(28, 'RVEDV', 'excluded_timepoint', '', 'none', 'evaluated in post-closure table');
    add(29, 'RVESV', 'excluded_timepoint', '', 'none', 'evaluated in post-closure table');
    add(30, 'Stroke volume', 'missing', '', 'none', '');
    add(31, 'Cardiac output, echo', 'missing', '', 'none', '');
    add(32, 'Ejection fraction', 'missing', '', 'none', '');
    add(33, 'LVOT diameter', 'excluded_echo', '', 'none', '');
    add(34, 'LVOT VTI', 'missing', '', 'none', '');
    add(35, 'RVOT diameter', 'excluded_echo', '', 'none', '');
    add(36, 'RVOT VTI', 'missing', '', 'none', '');
else
    add(6, 'Heart rate', 'input', '', 'readings', 'cycle length input of the prediction');
    add(7, 'Systolic pressure, NIBP', 'evaluation_only', 'SAP_max', 'mean_of_readings', 'non-invasive cuff');
    add(8, 'Diastolic pressure, NIBP', 'evaluation_only', 'SAP_min', 'mean_of_readings', 'non-invasive cuff');
    add(9, 'Systolic pressure, RFA', 'predicted_evaluated', 'SAP_max', 'config', '');
    add(10, 'Diastolic pressure, RFA', 'predicted_evaluated', 'SAP_min', 'config', '');
    add(11, 'Systolic pressure, DAo', 'evaluation_only', 'SAP_max', 'mean_of_readings', 'second catheter site');
    add(12, 'Diastolic pressure, DAo', 'evaluation_only', 'SAP_min', 'mean_of_readings', 'second catheter site');
    add(13, 'PA systolic pressure', 'predicted_evaluated', 'PAP_max', 'config', '');
    add(14, 'PA diastolic pressure', 'predicted_evaluated', 'PAP_min', 'config', '');
    add(15, 'PA mean pressure', 'predicted_evaluated', 'PAP_mean', 'config', '');
    add(16, 'RA mean pressure', 'predicted_evaluated', 'RAP_mean', 'config', '');
    add(17, 'Mean arterial pressure, form', 'evaluation_only', 'SAP_mean', 'mean_of_readings', 'form MAP from NIBP');
    for rr = 18:29
        names = {'Pulmonary flow Qp','Systemic flow Qs','Qp/Qs','PVR','SVR', ...
            'LVEDV','LVESV','RVEDV','RVESV','Stroke volume','Cardiac output','Ejection fraction'};
        add(rr, names{rr - 17}, 'missing', '', 'none', '');
    end
end
end

function put_role(spec, row, name, role, metric, rule, note)
% containers.Map is a handle, so this mutates the caller's map.
spec(sprintf('%d|%s', row, name)) = struct('role', role, 'metric', metric, ...
    'rule', rule, 'note', note);
end

function d = derived_entries(form, clinical)
% Values derived from protocol rows by stated formulas.
e = @(name, value, unit, metric, role, rule, method, note) struct('name', name, ...
    'value', value, 'unit', unit, 'metric', metric, 'role', role, 'rule', rule, ...
    'method', method, 'note', note);
if strcmp(form, 'pre')
    pre = clinical.pre_surgery;
    d = [ ...
        e('Mean arterial pressure (RFA, protocol formula)', pre.SAP_mean_mmHg, 'mmHg', 'SAP_mean', ...
          'fitted', 'config', 'derived', '(sys + 2 dia)/3, protocol row 20 formula on rows 12-13')
        e('Systemic flow Qs (derived)', pre.CO_Lmin, 'L/min', 'CO_Lmin', ...
          'fitted', 'config', 'derived', 'Qp / (Qp/Qs), rows 21 and 23')
        e('Shunt flow (derived)', pre.Q_shunt_Lmin, 'L/min', 'Q_shunt_Lmin', ...
          'fitted', 'config', 'derived', 'Qp - Qs; fitted, excluded from the acceptance gate')
        e('SVR (derived)', (pre.SAP_mean_mmHg - pre.RAP_mean_mmHg) / pre.CO_Lmin, 'WU', 'SVR', ...
          'evaluation_only', 'formula', 'derived', '(MAP - RAP) / Qs')
        e('LV-RV systolic gradient (derived)', 94 - 25, 'mmHg', 'LV_RV_systolic_gradient', ...
          'evaluation_only', 'formula', 'derived', 'row 9 LV minus RV systolic pressure')];
else
    post = clinical.post_surgery;
    d = [ ...
        e('Mean arterial pressure (RFA, protocol formula)', post.MAP_mmHg, 'mmHg', 'SAP_mean', ...
          'predicted_evaluated', 'config', 'derived', '(sys + 2 dia)/3 on rows 9-10')
        e('LVEDV (medical record)', post.LVEDV_mL, 'mL', 'LVEDV', ...
          'predicted_evaluated', 'config', 'echocardiography', 'post-closure echocardiography; form cell blank')
        e('LVESV (medical record)', post.LVESV_mL, 'mL', 'LVESV', ...
          'predicted_evaluated', 'config', 'echocardiography', 'post-closure echocardiography; form cell blank')
        e('RVEDV (medical record)', post.RVEDV_mL, 'mL', 'RVEDV', ...
          'predicted_evaluated', 'config', 'echocardiography', 'post-closure echocardiography; form cell blank')
        e('RVESV (medical record)', post.RVESV_mL, 'mL', 'RVESV', ...
          'predicted_evaluated', 'config', 'echocardiography', 'post-closure echocardiography; form cell blank')
        e('LVEF (derived from the volumes)', post.EF, '-', 'LVEF', ...
          'predicted_evaluated', 'config', 'derived', '(LVEDV - LVESV)/LVEDV')
        e('RVEF (derived from the volumes)', post.RVEF, '-', 'RVEF', ...
          'predicted_evaluated', 'config', 'derived', '(RVEDV - RVESV)/RVEDV')
        e('Cardiac output (medical record)', post.CO_Lmin, 'L/min', 'CO_Lmin', ...
          'predicted_evaluated', 'config', 'echocardiography estimate', 'form cell blank')
        e('Qp/Qs (medical record)', post.QpQs, '-', 'QpQs', ...
          'predicted_evaluated', 'config', 'echocardiography estimate', 'form cell blank')
        e('SVR (medical record)', post.SVR_WU, 'WU', 'SVR', ...
          'predicted_evaluated', 'config', 'derived', '(MAP - RAP)/cardiac output')];
end
end

function y = metric_value(metrics, name)
y = NaN;
if isempty(metrics) || isempty(name), return; end
if strcmp(name, 'LV_RV_systolic_gradient')
    if isfield(metrics, 'LVP_max') && isfield(metrics, 'RVP_max')
        y = metrics.LVP_max - metrics.RVP_max;
    end
    return;
end
if isfield(metrics, name), y = metrics.(name); end
end

function t = readings_text(x)
if isempty(x), t = ''; else, t = strjoin(compose('%g', x), '; '); end
end

function tiers = read_tiers(run_dir)
tiers = table();
f = fullfile(run_dir, 'tables', 'validation_target_tiers_pre_surgery.csv');
if exist(f, 'file')
    tiers = readtable(f, 'TextType', 'char');
    if ~all(ismember({'Metric','Tier'}, tiers.Properties.VariableNames))
        tiers = table();
    end
end
end

% =======================================================================
function T = inputs_table(clinical, params_pre, params_post)
c = clinical.common;
pre = clinical.pre_surgery;
names = {'Age';'Weight';'Height';'Body surface area';'Heart rate, pre-closure'; ...
    'Heart rate, post-closure';'VSD diameter (orifice)'};
units = {'yr';'kg';'cm';'m2';'bpm';'bpm';'mm'};
clin = [c.age_years; c.weight_kg; c.height_cm; c.BSA; pre.HR; clinical.post_surgery.HR; ...
    pre.VSD_diameter_mm];
used = [c.age_years; c.weight_kg; c.height_cm; c.BSA; params_pre.HR; params_post.HR; ...
    params_pre.vsd.diameter_mm];
src = {'pre row 1';'pre row 4';'pre row 3';'pre row 5';'pre row 6';'post row 6'; ...
    'pre row 7, mean of RV-side range'};
T = table(names, units, clin, used, src, 'VariableNames', ...
    {'Input','Unit','ProtocolValue','ValueUsedByModel','Source'});
end

function T = parameters_table(rp, params_cal, reference)
co = rp.calib_out;
reg = co.parameterRegistry;
names = co.names_all(:);
n = numel(names);
cal = nan(n, 1);
for i = 1:n
    cal(i) = get_calibration_param_value(params_cal, reference, names{i}, co.caseProfile);
end
T = table(names, 'VariableNames', {'Parameter'});
cols = {'unit','baseline_adult','baseline_scaled','seeded_value','lb','ub'};
for k = 1:numel(cols)
    if ismember(cols{k}, reg.Properties.VariableNames)
        T.(cols{k}) = reg.(cols{k})(match_rows(reg.name, names));
    end
end
T.GSASelected = logical(rp.optMask(:));
T.Calibrated = ismember(names, co.names(:));
T.ValueUsed = cal;
T.RatioToScaled = cal ./ T.baseline_scaled;
end

function idx = match_rows(reg_names, names)
idx = zeros(numel(names), 1);
for i = 1:numel(names)
    idx(i) = find(strcmp(reg_names, names{i}), 1);
end
end

function T = flatten_params(params)
paths = {}; vals = [];
walk(params, '');
T = table(paths(:), vals(:), 'VariableNames', {'Parameter','Value'});
    function walk(s, prefix)
        f = fieldnames(s);
        for k = 1:numel(f)
            v = s.(f{k});
            name = f{k};
            if ~isempty(prefix), name = [prefix '.' f{k}]; end
            if any(strcmp(name, {'idx','ic','scaling','clinical_override','maturation'}))
                continue;
            end
            if isstruct(v) && isscalar(v)
                walk(v, name);
            elseif (isnumeric(v) || islogical(v)) && isscalar(v)
                paths{end+1} = name; vals(end+1) = double(v); %#ok<AGROW>
            end
        end
    end
end

% =======================================================================
function [T, Q] = gsa_tables(rp, out_field, cfg_field, stage)
T = table(); Q = table();
if ~isfield(rp, out_field) || isempty(rp.(out_field)), return; end
g = rp.(out_field);
names = rp.(cfg_field).names(:);
N = NaN;
if isfield(g, 'checkpoint_signature') && isfield(g.checkpoint_signature, 'N_train')
    N = g.checkpoint_signature.N_train;
end
metrics = fieldnames(g);
rows = {}; qrows = {};
for m = 1:numel(metrics)
    v = g.(metrics{m});
    if ~isstruct(v) || ~isfield(v, 'ST'), continue; end
    for i = 1:numel(names)
        rows(end+1, :) = {stage, metrics{m}, names{i}, v.S1(i), v.ST(i)}; %#ok<AGROW>
    end
    [loo, deg, nnzc] = deal(NaN);
    try, loo = v.surrogate.Error.ModifiedLOO; catch, end
    try, deg = v.surrogate.PCE.Basis.Degree; catch, end
    try, nnzc = nnz(v.surrogate.PCE.Coefficients); catch, end
    qrows(end+1, :) = {stage, metrics{m}, N, deg, nnzc, loo}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Stage','Output','Parameter','S1','ST'});
Q = cell2table(qrows, 'VariableNames', {'Stage','Output','DesignSize', ...
    'PolynomialDegree','NonzeroCoefficients','ModifiedLOO'});
end

function T = mask_table(rp, gsa_init)
co = rp.calib_out;
names = co.names_all(:);
metrics = {};
if isfield(rp, 'mask_metrics'), metrics = rp.mask_metrics; end
maxST = nan(numel(names), 1);
if ~isempty(gsa_init)
    for i = 1:numel(names)
        members = expand_calibration_parameter_names(names{i}, co.caseProfile);
        sel = ismember(gsa_init.Parameter, members) & ismember(gsa_init.Output, metrics);
        if any(sel), maxST(i) = max(gsa_init.ST(sel)); end
    end
end
T = table(names, maxST, logical(rp.optMask(:)), logical(co.mask(:)), repmat(rp.sobol_ST_threshold, numel(names), 1), ...
    repmat({strjoin(cellstr(metrics), ';')}, numel(names), 1), ...
    'VariableNames', {'Parameter','MaxST_MaskMetrics','SelectedByGSA','Calibrated','Threshold','MaskMetrics'});
end

function T = reduction_table(red)
rows = cell(numel(red.by_size), 4);
for i = 1:numel(red.by_size)
    r = red.by_size{i};
    rows(i, :) = {r.p, r.dof, r.cond, strjoin(r.names(:)', ';')};
end
T = cell2table(rows, 'VariableNames', {'SubsetSize','DegreesOfFreedom','ConditionNumber','Parameters'});
T = [T; {numel(red.column_norms.Parameter), red.n_metrics - numel(red.column_norms.Parameter), ...
    red.full_cond, strjoin(red.column_norms.Parameter(:)', ';')}];
T = sortrows(T, 'SubsetSize', 'descend');
end

% =======================================================================
function T = last_cycle_table(sim, params)
t = sim.t(:);
T_HB = 60 / params.HR;
mask = t >= t(end) - T_HB;
tc = t(mask);
XV = sim.V(mask, :);
[P, Q] = reconstruct_hemodynamic_signals(tc, XV, params);
T = table(tc - tc(1), (tc - tc(1)) / T_HB, 'VariableNames', {'time_s','phase'});
for f = fieldnames(P)'
    T.(['P_' f{1} '_mmHg']) = P.(f{1})(:);
end
for f = fieldnames(Q)'
    T.(['Q_' f{1} '_mLs']) = Q.(f{1})(:);
end
for f = fieldnames(params.idx)'
    T.(['state_' f{1}]) = XV(:, params.idx.(f{1}));
end
end

function v = post_variants(clinical)
% Declared post-closure datasets; falls back to the single struct when a
% patient has only one.
if isfield(clinical, 'post_surgery_variants') && ~isempty(clinical.post_surgery_variants)
    v = clinical.post_surgery_variants;
else
    v = struct('name', 'post_surgery', 'source', 'clinical.post_surgery', ...
        'provenance', 'single declared dataset', 'data', clinical.post_surgery);
end
end

function copy_if_exists(src, dst)
if exist(src, 'file'), copyfile(src, dst); end
end

function write_json(path, s)
fid = fopen(path, 'w');
fprintf(fid, '%s', jsonencode(s, 'PrettyPrint', true));
fclose(fid);
end
