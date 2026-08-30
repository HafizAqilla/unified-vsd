function out = evaluate_post_closure_prediction(candidate_mat, clinical, verbose)
% EVALUATE_POST_CLOSURE_PREDICTION
% -----------------------------------------------------------------------
% Genuine out-of-sample validation: does a model calibrated on PRE-closure
% data alone predict the POST-closure state?
%
% WHY THIS IS A REAL HOLDOUT
% --------------------------
% docs/reyna_statistical_calibration_results_20260829.md §5 records that
% nothing in the Reyna recipe is presently a genuine validation holdout:
% every finite target is either fitted, or is algebra over quantities that
% are fitted (SVR, Q_shunt), so "predicting" it demonstrates nothing the
% model was not already told.
%
% The post-closure pressures obtained from the 06/04/2026 procedure log are
% different in kind. They are:
%   - independent measurements, not derived from the pre-closure rows;
%   - taken in a genuinely different haemodynamic state;
%   - never seen by the pre-only calibration that produced these parameters.
%
% So this is the first true prediction test available for this model. The
% only intervention applied is closing the defect -- no refitting, no
% parameter adjustment of any kind.
%
% WHAT A FAILURE WOULD MEAN
% -------------------------
% A large error here is informative rather than merely disappointing: it
% would say the shared-parameter assumption underlying Phase 4 (one patient,
% one parameter set, only the defect changes) does not hold, which is a
% finding that must be reported rather than tuned away.
%
% INPUTS:
%   candidate_mat - path to a params_*_candidate_pre_surgery.mat            [-]
%   clinical      - unified clinical struct (default: patient_reyna())      [-]
%   verbose       - print the comparison table (default true)               [-]
%
% OUTPUTS:
%   out.table     - per-metric predicted vs measured, with z-scores         [-]
%   out.chi2      - sum of z^2 over post targets                            [-]
%   out.n         - number of post targets compared                    [count]
%   out.within10  - count within 10%                                   [count]
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-30
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 2 || isempty(clinical); clinical = patient_reyna(); end
if nargin < 3; verbose = true; end

S = load(candidate_mat);
params = extract_params(S);

% The ONLY change: close the defect. Mode-aware, because vsd_shunt_model
% ignores R.vsd entirely in orifice_bidirectional mode -- see
% objective_joint_pre_post.m close_vsd for the full rationale.
params_post = params;
if isfield(params_post, 'R') && isfield(params_post.R, 'vsd')
    params_post.R.vsd = 1e6;
end
if isfield(params_post, 'vsd')
    if isfield(params_post.vsd, 'area_mm2'); params_post.vsd.area_mm2 = 0; end
    if isfield(params_post.vsd, 'D_mm');     params_post.vsd.D_mm = 0;     end
end

% Sanity: the shunt must actually be shut, or this measures nothing.
q_res = vsd_shunt_model(90, 20, params_post);
assert(abs(q_res) < 1e-9, ...
    'evaluate_post_closure_prediction:shuntStillOpen', ...
    'VSD closure failed; the prediction would be meaningless.');

sim = integrate_system(params_post);
metrics = compute_clinical_indices(sim, params_post);

targets = get_calibration_targets('post_surgery', clinical);

names = {}; measured = []; predicted = []; err_pct = []; zs = [];
for i = 1:numel(targets)
    t = targets(i);
    if ~isfinite(t.ClinicalValue) || ~isfield(metrics, t.Metric); continue; end
    y = metrics.(t.Metric);
    if ~isfinite(y); continue; end
    sigma = resolve_sigma(t);
    names{end+1} = t.Metric;                                     %#ok<AGROW>
    measured(end+1) = t.ClinicalValue;                           %#ok<AGROW>
    predicted(end+1) = y;                                        %#ok<AGROW>
    err_pct(end+1) = 100 * (y - t.ClinicalValue) / abs(t.ClinicalValue); %#ok<AGROW>
    zs(end+1) = (y - t.ClinicalValue) / sigma;                   %#ok<AGROW>
end

out = struct();
out.table = table(names(:), measured(:), predicted(:), err_pct(:), zs(:), ...
    'VariableNames', {'Metric','Measured','Predicted','Error_pct','ZScore'});
out.chi2 = sum(zs.^2);
out.n = numel(zs);
out.within10 = sum(abs(err_pct) <= 10);

if verbose
    fprintf('\n=== POST-CLOSURE PREDICTION (OUT-OF-SAMPLE) ===\n');
    fprintf('  Source: %s\n', candidate_mat);
    fprintf('  These post-closure measurements were NEVER used to fit these\n');
    fprintf('  parameters. Only the defect was closed. No refitting.\n\n');
    disp(out.table);
    fprintf('  Within 10%%: %d of %d\n', out.within10, out.n);
    fprintf('  chi2 = %.3f over %d targets  (chi2/N = %.3f)\n', ...
        out.chi2, out.n, out.chi2 / max(out.n, 1));
end
end

% =======================================================================
function params = extract_params(S)
fn = fieldnames(S);
for i = 1:numel(fn)
    v = S.(fn{i});
    if isstruct(v) && isfield(v, 'R') && isfield(v, 'vsd')
        params = v;
        return;
    end
    if isstruct(v) && isfield(v, 'params') && isstruct(v.params)
        params = v.params;
        return;
    end
end
error('evaluate_post_closure_prediction:noParams', ...
    'No parameter struct found in %s', strjoin(fn, ', '));
end

% =======================================================================
function sigma = resolve_sigma(t)
sigma = NaN;
if isfield(t, 'UncertaintyAbs') && isfinite(t.UncertaintyAbs) && t.UncertaintyAbs > 0
    sigma = t.UncertaintyAbs;
elseif isfield(t, 'UncertaintyFraction') && isfinite(t.UncertaintyFraction) && ...
        t.UncertaintyFraction > 0
    sigma = t.UncertaintyFraction * abs(t.ClinicalValue);
end
if ~isfinite(sigma) || sigma <= 0
    sigma = max(abs(t.ClinicalValue), 1e-6) * 0.10;
end
sigma = max(sigma, 1e-9);
end
