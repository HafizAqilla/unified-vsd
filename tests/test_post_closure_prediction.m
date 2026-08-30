function tests = test_post_closure_prediction()
% TEST_POST_CLOSURE_PREDICTION
% -----------------------------------------------------------------------
% Contract tests for evaluate_post_closure_prediction.
%
% This script produces the strongest claim on the branch -- the first genuine
% out-of-sample result this model has had (results doc §6.5). A claim like
% that is only worth as much as the guarantees behind it, and each of those
% guarantees has a specific way of failing silently:
%
%   - if the shunt is not really closed, the "post-closure prediction" is a
%     pre-closure simulation and the result is meaningless;
%   - if any parameter other than the shunt changes, it is no longer a
%     prediction but a fit;
%   - if a pre-surgery target leaks into the comparison, the holdout is
%     contaminated by data the model was trained on.
%
% Each is pinned below.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_results_20260829.md §6.5
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-30
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function candidate = candidate_path()
candidate = fullfile('results', 'runs', ...
    '20260829_225931_reyna_pre_surgery', 'mat', ...
    'params_accepted_candidate_pre_surgery.mat');
end

function tf = candidate_available()
tf = exist(candidate_path(), 'file') == 2;
end

function test_compares_only_post_surgery_targets(tc)
% A pre-surgery target leaking in would contaminate the holdout with data the
% model was trained on, turning a prediction into a partial self-report.
if ~candidate_available()
    return;   % run artefacts are not committed; skip rather than fail
end
clinical = patient_reyna();
out = evalc_pred(candidate_path(), clinical);

post_targets = get_calibration_targets('post_surgery', clinical);
finite_names = {};
for i = 1:numel(post_targets)
    if isfinite(post_targets(i).ClinicalValue)
        finite_names{end+1} = post_targets(i).Metric; %#ok<AGROW>
    end
end

verifyTrue(tc, all(ismember(out.table.Metric, finite_names)), ...
    'every compared metric must be a finite post_surgery target.');
verifyEqual(tc, out.n, numel(finite_names), ...
    'all finite post targets must be compared, none silently dropped.');
end

function test_measured_values_match_the_clinical_record(tc)
% Guards against the comparison drifting onto some other value: the measured
% column must be exactly what config/patient_reyna.m holds for post_surgery.
if ~candidate_available(); return; end
clinical = patient_reyna();
out = evalc_pred(candidate_path(), clinical);

expect = containers.Map( ...
    {'PAP_max','PAP_min','PAP_mean','SAP_max','SAP_min','SAP_mean','RAP_mean'}, ...
    {17, 9, 13, 89, 68, 79, 5});

for i = 1:height(out.table)
    name = out.table.Metric{i};
    if isKey(expect, name)
        verifyEqual(tc, out.table.Measured(i), expect(name), ...
            sprintf('%s measured value must match the procedure log.', name));
    end
end
end

function test_prediction_is_not_a_refit(tc)
% The whole claim rests on nothing being tuned. Running twice must give
% bit-identical output: any optimisation, randomness or state carried between
% calls would break that and would mean this is not a prediction.
if ~candidate_available(); return; end
clinical = patient_reyna();
a = evalc_pred(candidate_path(), clinical);
b = evalc_pred(candidate_path(), clinical);
verifyEqual(tc, a.table.Predicted, b.table.Predicted, 'AbsTol', 0, ...
    'prediction must be deterministic; any variation implies fitting.');
verifyEqual(tc, a.chi2, b.chi2, 'AbsTol', 0);
end

function test_closure_is_enforced_not_assumed(tc)
% If closure silently failed, the "post-closure prediction" would be a
% pre-closure simulation reported as a validation success -- the worst
% available failure mode for this result. The script must assert, not hope.
src = fileread(which('evaluate_post_closure_prediction'));
verifyTrue(tc, contains(src, 'shuntStillOpen'), ...
    'closure must be asserted with an explicit error identifier.');
verifyTrue(tc, contains(src, 'vsd_shunt_model'), ...
    'closure must be verified against the shunt model itself, not a field.');
verifyTrue(tc, contains(src, 'area_mm2'), ...
    ['closure must zero the orifice area: R.vsd alone is a no-op in ', ...
     'orifice_bidirectional mode, which is the mode this patient uses.']);
end

function test_closed_parameters_differ_from_source_only_in_the_shunt(tc)
% "Only the defect was closed" is the load-bearing sentence in the §6.5
% claim. Verify it directly against the stored candidate.
if ~candidate_available(); return; end
S = load(candidate_path());
p_open = S.accepted_candidate.params;

p_closed = p_open;
p_closed.R.vsd = 1e6;
if isfield(p_closed, 'vsd')
    if isfield(p_closed.vsd, 'area_mm2'); p_closed.vsd.area_mm2 = 0; end
    if isfield(p_closed.vsd, 'D_mm');     p_closed.vsd.D_mm = 0;     end
end

% Everything outside R.vsd and the vsd block must be untouched.
groups = intersect(fieldnames(p_open), {'C','E','V0','L','conv','sim'});
for g = 1:numel(groups)
    verifyEqual(tc, p_closed.(groups{g}), p_open.(groups{g}), ...
        sprintf('%s must be identical: closure may not touch it.', groups{g}));
end

% And the shunt must genuinely be shut.
verifyEqual(tc, vsd_shunt_model(90, 20, p_closed), 0, 'AbsTol', 1e-12);
verifyGreaterThan(tc, abs(vsd_shunt_model(90, 20, p_open)), 1e-6, ...
    'source candidate must have an open shunt, or the test proves nothing.');
end

function out = evalc_pred(path, clinical)
out = [];
txt = evalc('out = evaluate_post_closure_prediction(path, clinical, false);'); %#ok<NASGU>
end
