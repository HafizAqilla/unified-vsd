function [J, info] = objective_joint_pre_post(x, params0_pre, params0_post, ...
    clinical, calib, pce_surrogate)
% OBJECTIVE_JOINT_PRE_POST
% -----------------------------------------------------------------------
% Joint pre/post-closure calibration objective (PRD
% reyna_statistical_calibration_v1 Phase 4).
%
% WHY
% ---
% One patient has one set of structural parameters. Between pre-closure and
% post-closure essentially only the defect changes. Fitting both haemodynamic
% states from a single shared parameter vector therefore roughly doubles the
% observation count without adding free parameters, which is the only
% available route to positive degrees of freedom for this case:
%
%   pre-only : N = 9  governed observations, p = 12  ->  dof = -3
%   joint    : N = 16 governed observations, p = 12  ->  dof = +4
%
% See docs/reyna_statistical_calibration_results_20260829.md §6.0 for the
% clinical record that makes this possible: the post-closure pressures come
% from the SAME catheterisation session as the pre-closure ones (device
% released 12.06.23, readings 12.15-12.32, unchanged anaesthesia), so the two
% states are directly comparable rather than being two separate studies.
%
% WHAT DIFFERS BETWEEN THE TWO SIMULATIONS
% ----------------------------------------
% Shared  : chamber elastances, V0, vascular R and C, group scales.
% Differs : the shunt only. Pre keeps its geometry-derived orifice; post is
%           closed via R.vsd = VSD_CLOSED_RESISTANCE, matching the convention
%           already used for the pre-to-post seed package in main_run.m:839.
%
% INPUTS:
%   x             - shared calibration vector over calib.names             [-]
%   params0_pre   - baseline pre-closure parameter struct                  [-]
%   params0_post  - baseline post-closure parameter struct                 [-]
%   clinical      - unified clinical struct with both scenarios            [-]
%   calib         - calibration config (calibration_param_sets)            [-]
%   pce_surrogate - optional PCE surrogate, passed through                 [-]
%
% OUTPUTS:
%   J    - chi2_pre + chi2_post, plus shared-vector penalties counted ONCE  [-]
%   info - breakdown struct:
%          .chi2_pre, .chi2_post, .chi2_total
%          .n_pre, .n_post, .n_total   governed observation counts    [count]
%          .penalty                    shared regularisation/plausibility [-]
%          .post_available             false when post targets are absent
%          .valid                      false if either simulation failed
%
% ASSUMPTIONS:
%   - Regularisation and plausibility terms describe the SHARED vector and
%     are therefore applied once, not once per scenario. Applying them twice
%     would silently double the regularisation weight relative to the
%     single-scenario objective and make the two incomparable.
%   - Degrades to pre-only (with a warning) when the post scenario carries no
%     finite targets, so this objective is safe to call for any patient.
%
% REFERENCES:
%   [1] docs/reyna_statistical_calibration_prd.md §7
%   [2] docs/reyna_statistical_calibration_results_20260829.md §0.2, §6.0
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-30
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 6
    pce_surrogate = [];
end

info = init_info();

% ---- Apply the SHARED vector to both scenario parameter structs ---------
params_pre  = apply_shared_vector(params0_pre,  x, calib);
params_post = apply_shared_vector(params0_post, x, calib);

% The shunt is the only thing allowed to differ.
params_post = close_vsd(params_post);

% ---- Pre-closure data term ---------------------------------------------
[chi2_pre, n_pre, ok_pre] = scenario_chi2(params_pre, clinical, 'pre_surgery', ...
    calib, pce_surrogate);
info.chi2_pre = chi2_pre;
info.n_pre = n_pre;

% ---- Post-closure data term --------------------------------------------
[chi2_post, n_post, ok_post] = scenario_chi2(params_post, clinical, 'post_surgery', ...
    calib, pce_surrogate);
info.chi2_post = chi2_post;
info.n_post = n_post;
info.post_available = n_post > 0;

if ~info.post_available
    % Not an error: a patient may genuinely have no post-closure record.
    % Reducing to pre-only keeps this objective callable everywhere, but the
    % caller must know the DOF argument above no longer applies.
    warning('objective_joint_pre_post:noPostTargets', ...
        ['No finite post_surgery targets; joint objective has reduced to ', ...
         'pre-only. The degrees-of-freedom benefit of joint inversion does ', ...
         'NOT apply to this result.']);
    chi2_post = 0;
    ok_post = true;
end

info.valid = ok_pre && ok_post;
if ~info.valid
    J = calib.invalidPenaltyScale;
    info.chi2_total = Inf;
    return;
end

info.chi2_total = chi2_pre + chi2_post;
info.n_total = n_pre + n_post;

% ---- Shared-vector penalties, applied ONCE -----------------------------
info.penalty = shared_vector_penalty(x, calib);

J = info.chi2_total + info.penalty;
end

% =======================================================================
function info = init_info()
info = struct('chi2_pre', Inf, 'chi2_post', Inf, 'chi2_total', Inf, ...
    'n_pre', 0, 'n_post', 0, 'n_total', 0, 'penalty', 0, ...
    'post_available', false, 'valid', false);
end

% =======================================================================
function params = apply_shared_vector(params, x, calib)
% APPLY_SHARED_VECTOR - write the shared calibration vector into a scenario
% parameter struct. Identical values go into both structs by construction,
% which is what makes the parameters "shared" rather than merely similar.
for idx = 1:numel(calib.names)
    params = set_calibration_param_value(params, calib.referenceParams, ...
        calib.names{idx}, x(idx), calib.caseProfile);
end
end

% =======================================================================
function params = close_vsd(params)
% CLOSE_VSD - post-closure shunt state, correct for EVERY vsd mode.
%
% This must be mode-aware, and getting it wrong is silent rather than loud.
% vsd_shunt_model dispatches on params.vsd.mode:
%
%   resistive modes ('*_diode', plain resistive)
%       Q = dP / R.vsd            -> closed by making R.vsd very large
%   'orifice_bidirectional'       (the mode Reyna actually uses)
%       Q = Cd * A * sqrt(...)    -> R.vsd is NEVER READ; closure requires
%                                    vsd.area_mm2 = 0, which is the explicit
%                                    early-return branch in vsd_shunt_model
%
% Setting only R.vsd would therefore leave an orifice-mode patient shunting
% at full strength in the "post-closure" simulation, and the joint fit would
% silently be fitting two open-VSD states. Both channels are closed here.
%
% NOTE: main_run.m's pre-to-post seed package (`post_seed_params.R.vsd = 1e6`)
% closes only the resistive channel, so for an orifice-mode patient that seed
% is not actually a closed-VSD model. Flagged in
% docs/reyna_statistical_calibration_results_20260829.md §6.4.
if isfield(params, 'R') && isfield(params.R, 'vsd')
    params.R.vsd = vsd_closed_resistance();
end
if isfield(params, 'vsd')
    if isfield(params.vsd, 'area_mm2')
        params.vsd.area_mm2 = 0;
    end
    if isfield(params.vsd, 'D_mm')
        params.vsd.D_mm = 0;
    end
end
end

% =======================================================================
function r = vsd_closed_resistance()
% VSD_CLOSED_RESISTANCE - single definition of "closed" for the resistive
% channel, mirroring main_run.m:839 so the two cannot drift apart silently.
r = 1e6;
end

% =======================================================================
function [chi2, n_obs, ok] = scenario_chi2(params, clinical, scenario, calib, pce)
% SCENARIO_CHI2 - sum of squared sigma-normalised residuals over the governed
% targets of one scenario.
%
% Uses the same sigma resolution as Phase 1's sigma-weighted objective so the
% joint statistic is expressed in the same units as the single-scenario one
% and the two remain comparable.
chi2 = Inf;
n_obs = 0;
ok = false;

try
    sim = integrate_system(params);
    metrics = compute_clinical_indices(sim, params);
catch
    return;
end

targets = get_calibration_targets(scenario, clinical);

total = 0;
count = 0;
for idx = 1:numel(targets)
    name = targets(idx).Metric;
    y_clin = targets(idx).ClinicalValue;
    if ~isfinite(y_clin) || ~isfield(metrics, name)
        continue;
    end
    if ~is_governed_target(name, scenario, calib)
        continue;
    end
    y_model = metrics.(name);
    if ~isfinite(y_model)
        continue;
    end
    sigma = resolve_sigma(targets(idx), y_clin);
    z = (y_model - y_clin) / sigma;
    total = total + z^2;
    count = count + 1;
end

chi2 = total;
n_obs = count;
ok = true;
end

% =======================================================================
function tf = is_governed_target(name, scenario, calib)
% IS_GOVERNED_TARGET - only rows inside the governed primary RMSE set count,
% matching compute_chi_squared_report so the joint chi2 and the reported one
% describe the same observation set.
tf = false;
tiers = [];
if isfield(calib, 'targetTiersByScenario') && ...
        isfield(calib.targetTiersByScenario, scenario)
    tiers = calib.targetTiersByScenario.(scenario);
elseif isfield(calib, 'targetTiers')
    tiers = calib.targetTiers;
end
if ~isstruct(tiers) || ~isfield(tiers, 'table') || isempty(tiers.table)
    return;
end
row = find(strcmp(tiers.table.Metric, name), 1, 'first');
if isempty(row)
    return;
end
tf = logical(tiers.table.IncludedInPrimaryRMSE(row));
end

% =======================================================================
function sigma = resolve_sigma(target, y_clin)
% RESOLVE_SIGMA - declared measurement uncertainty for one target.
%
% Resolution order matches Phase 1's sigma_weighted_residual exactly
% (objective_calibration.m): UncertaintyAbs when finite, else
% UncertaintyFraction * |ClinicalValue|, else a 10% fallback. Keeping the two
% identical is what makes chi2_pre from this objective comparable with the
% chi2 reported by compute_chi_squared_report.
sigma = NaN;
if isfield(target, 'UncertaintyAbs') && isfinite(target.UncertaintyAbs) && ...
        target.UncertaintyAbs > 0
    sigma = target.UncertaintyAbs;
elseif isfield(target, 'UncertaintyFraction') && ...
        isfinite(target.UncertaintyFraction) && target.UncertaintyFraction > 0
    sigma = target.UncertaintyFraction * abs(y_clin);
end
if ~isfinite(sigma) || sigma <= 0
    sigma = max(abs(y_clin), 1e-6) * 0.10;
end
sigma = max(sigma, 1e-9);
end

% =======================================================================
function penalty = shared_vector_penalty(x, calib)
% SHARED_VECTOR_PENALTY - regularisation on the shared vector.
%
% Deliberately computed once. The vector is shared between the two
% simulations, so charging it per scenario would double its weight relative
% to the single-scenario objective and make joint and pre-only results
% incomparable -- see PRD §7.3.
penalty = 0;
if ~isfield(calib, 'regLambda') || ~isfinite(calib.regLambda) || calib.regLambda <= 0
    return;
end
if ~isfield(calib, 'x0') || numel(calib.x0) ~= numel(x)
    return;
end
x0 = calib.x0(:);
xv = x(:);
valid = isfinite(xv) & isfinite(x0) & x0 > 0 & xv > 0;
if any(valid)
    penalty = calib.regLambda * sum(log(xv(valid) ./ x0(valid)).^2);
end
end
