function J = objective_calibration(x, params0, clinical, calib, scenario, pce_surrogate)
if nargin < 6; pce_surrogate = []; end
% OBJECTIVE_CALIBRATION
% -----------------------------------------------------------------------
% Weighted normalised least-squares objective for fmincon.
%
%   J = Σ_k  w_k(err_k) · ((y_k(x) − y_k^{clin}) / y_k^{clin})²
%     + J_CO_guard                    (soft CO range penalty)
%     + J_vol_balance                 (volume under-prediction guard)
%     + λ · Σ_i  ((x_i − x0_i) / x0_i)²
%
% Primary metric penalty (QpQs, SAP_mean, LVEF):
%   w_k is amplified smoothly as error grows above 5%:
%     scale = 1 + 19 * t / (1 + t),  t = max(0, err/0.05 − 1)
%   This replaces a hard 100x cliff.  The smooth curve keeps a finite
%   gradient everywhere so L-BFGS can still improve volume/PAP metrics
%   without hitting an objective wall.
%     t=0 (≤5%) → 1x,  t=1 (10%) → ~10.5x,  t→∞ → 20x
%
% CO guard (§6b):
%   Adds 10·((violation)/limit)² when CO < 1.5 or > 10.0 L/min.
%
% Volume balance guard (§6c):
%   Adds 2·mean(rel_under²) when LVEDV or RVEDV is under-predicted.
%
% INPUTS:
%   x         - parameter vector (free variables, ordered per calib.names)
%   params0   - baseline parameter struct (NOT modified in-place)
%   clinical  - unified clinical struct
%   calib     - configuration struct from calibration_param_sets.m
%   scenario  - 'pre_surgery' | 'post_surgery'
%
% OUTPUTS:
%   J         - scalar objective value  (dimensionless)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-04-08
% VERSION:  2.0
% -----------------------------------------------------------------------

%% 1. Update parameter struct with trial x
params = params0;
for i = 1:numel(calib.names)
    params = set_param_by_name(params, calib.names{i}, x(i));
end

%% 2 & 3. Simulate or Predict
if ~isempty(pce_surrogate)
    x_full = pce_surrogate.cfg.x0;
    for idx_full = 1:numel(pce_surrogate.cfg.names)
        idx = find(strcmp(calib.names, pce_surrogate.cfg.names{idx_full}));
        if ~isempty(idx)
            x_full(idx_full) = x(idx);
        end
    end
    % Clamp x_full to the PCE surrogate training bounds before evaluation.
    % UQLab normalises inputs using Internal.Input.Marginals, NOT cfg.lb/ub.
    % cfg.lb/ub can be stale (e.g. checkpoint from a prior patient run), so
    % we extract the ground-truth bounds from the first available surrogate.
    % If extraction fails, fall back to cfg bounds.
    pce_lb = pce_surrogate.cfg.lb(:);
    pce_ub = pce_surrogate.cfg.ub(:);
    surr_fns = fieldnames(pce_surrogate);
    for kf = 1:numel(surr_fns)
        fn = surr_fns{kf};
        if isstruct(pce_surrogate.(fn)) && isfield(pce_surrogate.(fn), 'surrogate')
            s = pce_surrogate.(fn).surrogate;
            if isfield(s, 'Internal') && isfield(s.Internal, 'Input') && ...
                    isfield(s.Internal.Input, 'Marginals')
                margs = s.Internal.Input.Marginals;
                if numel(margs) == numel(pce_lb)
                    pce_lb = arrayfun(@(m) m.Parameters(1), margs(:));
                    pce_ub = arrayfun(@(m) m.Parameters(2), margs(:));
                end
            end
            break;
        end
    end
    x_full = min(max(x_full(:), pce_lb), pce_ub);
    metrics = struct();
    % Evaluate polynomial chaos surrogate avoiding numerical integration
    for k = 1:numel(calib.metricFields)
        mf = calib.metricFields{k};
        if isfield(pce_surrogate, mf) && isfield(pce_surrogate.(mf), 'surrogate')
            metrics.(mf) = uq_evalModel(pce_surrogate.(mf).surrogate, x_full(:)');
        end
    end
    % Compute derived constraints if not present
    if ~isfield(metrics, 'CO_Lmin') && isfield(metrics, 'QpQs')
        % Fallback for guards if not modelled
        metrics.CO_Lmin = 5.0;
    end
else
    params.sim.nCyclesSteady = 40;
    params.sim.ss_tol_P      = 0.5;
    params.sim.ss_tol_V      = 0.5;
    params.sim.rtol          = 1e-8;
    params.sim.atol          = 1e-10;
    try
        sim = integrate_system(params);
    catch
        x_range = max(calib.ub(:) - calib.lb(:), 1e-6);
        dist_sq = sum(((x(:) - calib.x0(:)) ./ x_range).^2);
        J = 1e4 * (1 + dist_sq);
        return;
    end
    metrics = compute_clinical_indices(sim, params);
end


%% 4. Select clinical target sub-section
switch scenario
    case 'pre_surgery',  src = clinical.pre_surgery;
    case 'post_surgery', src = clinical.post_surgery;
end

%% 5. Build clinical lookup table (field mapping same as validation_report.m)
% Map model metric name → clinical struct field name
field_map = build_field_map(scenario);

%% 6. Weighted squared errors
J = 0;
for k = 1:numel(calib.metricFields)
    mf  = calib.metricFields{k};
    w   = calib.weights.(mf);

    % Retrieve clinical value via field map
    if isfield(field_map, mf) && isfield(src, field_map.(mf))
        y_clin = src.(field_map.(mf));
    elseif isfield(src, mf)
        y_clin = src.(mf);
    else
        continue;
    end

    if isnan(y_clin), continue; end

    if isfield(metrics, mf)
        y_model = metrics.(mf);
        denom   = max(abs(y_clin), 1e-6);
        err_rel = abs(y_model - y_clin) / denom;
        
        % Smooth amplification for primary metrics.
        % Replaces the hard 100x cliff: that jump created a steep objective wall
        % that caused L-BFGS to ignore volume/PAP directions.
        % New rule: scale rises continuously from 1x (at <=5% error) toward
        % ~20x asymptotically, keeping a finite gradient everywhere.
        %   t = 0          → scale = 1x  (error at or below threshold)
        %   t = 1 (10% err)→ scale ≈ 10.5x
        %   t → ∞          → scale → 20x
        if ismember(mf, {'QpQs', 'SAP_mean', 'SAP_max', 'PAP_max', 'PAP_mean', 'SAP_min', 'PAP_min', 'CO_Lmin', 'LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVEF', 'SVR', 'LAP_mean'}) 
            t = max(0, err_rel / 0.05 - 1);   % 0 below 5%, positive above
            w = w * (1 + 19 * t / (1 + t));   % smooth, asymptotes to 20x
        end

        J = J + w * err_rel^2;
    end
end

%% 6b. CO physiological guard
% Soft penalty when cardiac output leaves a viable physiological range.
% Does not require a patient-specific CO target — it guards against
% parameter excursions that produce non-physiological haemodynamics.
%   Paediatric VSD range: 1.5 – 10.0 L/min (covers high-shunt states).
if isfield(metrics, 'CO_Lmin')
    co_val = metrics.CO_Lmin;
    co_lo  = 1.5;   % [L/min]
    co_hi  = 10.0;  % [L/min]
    if co_val < co_lo
        J = J + 10 * ((co_lo - co_val) / co_lo)^2;
    elseif co_val > co_hi
        J = J + 10 * ((co_val - co_hi) / co_hi)^2;
    end
end

%% 6c. Ventricular volume balance guard
% Extra penalty when both LVEDV and RVEDV are simultaneously under-predicted.
% The weighted metric loop already includes these terms; this guard adds a
% small asymmetric push to prevent the optimizer treating large volume under-
% predictions as equally cheap as over-predictions of the same magnitude.
vol_pairs = {{'LVEDV','LVEDV_mL'}, {'RVEDV','RVEDV_mL'}};
vol_under  = 0;
n_vol_valid = 0;
for kv = 1:numel(vol_pairs)
    vf   = vol_pairs{kv}{1};
    cf   = vol_pairs{kv}{2};
    if isfield(metrics, vf) && isfield(src, cf)
        y_c = src.(cf);
        y_m = metrics.(vf);
        if ~isnan(y_c) && ~isnan(y_m) && y_m < y_c
            rel_under = (y_c - y_m) / max(abs(y_c), 1e-6);
            vol_under  = vol_under + rel_under^2;
            n_vol_valid = n_vol_valid + 1;
        end
    end
end
if n_vol_valid > 0
    J = J + 2.0 * vol_under / n_vol_valid;
end

%% 7. Regularisation: pull parameters toward initial guess
if isfield(calib, 'regLambda') && calib.regLambda > 0
    x0  = calib.x0(:);
    J   = J + calib.regLambda * sum(((x(:) - x0) ./ max(abs(x0), 1e-6)).^2);
end

end  % objective_calibration

% =========================================================================
%  LOCAL HELPERS
% =========================================================================

function params = set_param_by_name(params, name, value)
% SET_PARAM_BY_NAME — assign value to nested struct field (dot notation)
parts = strsplit(name, '.');
switch numel(parts)
    case 2, params.(parts{1}).(parts{2}) = value;
    case 3, params.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_param_by_name:unsupportedDepth', ...
              'Only 2- or 3-level notation supported. Got: %s', name);
end
end

function fmap = build_field_map(scenario)
% BUILD_FIELD_MAP — maps model metric name → clinical sub-struct field name
%   This must stay in sync with validation_report.m metric_defs.
fmap = struct();
switch scenario
    case 'pre_surgery'
        fmap.RAP_mean  = 'RAP_mean_mmHg';
        fmap.LAP_mean  = 'LAP_mean_mmHg';   % estimated pulmonary back-pressure
        fmap.PAP_min   = 'PAP_dia_mmHg';    % diastolic PA pressure
        fmap.PAP_max   = 'PAP_sys_mmHg';    % systolic PA pressure
        fmap.PAP_mean  = 'PAP_mean_mmHg';
        fmap.SAP_min   = 'SAP_dia_mmHg';    % diastolic systemic pressure
        fmap.SAP_max   = 'SAP_sys_mmHg';    % systolic systemic pressure
        fmap.SAP_mean  = 'SAP_mean_mmHg';
        fmap.QpQs      = 'QpQs';
        fmap.PVR       = 'PVR_WU';
        fmap.SVR       = 'SVR_WU';
        fmap.LVEDV     = 'LVEDV_mL';
        fmap.LVESV     = 'LVESV_mL';
        fmap.RVEDV     = 'RVEDV_mL';
        fmap.RVESV     = 'RVESV_mL';
        fmap.LVEF      = 'LVEF';
    case 'post_surgery'
        fmap.SAP_min   = 'SAP_sys_mmHg';
        fmap.SAP_max   = 'SAP_sys_mmHg';
        fmap.SAP_mean  = 'MAP_mmHg';
        fmap.SVR       = 'SVR_WU';
        fmap.PVR       = 'PVR_WU';
        fmap.LVEF      = 'LVEF';
        fmap.RVEF      = 'RVEF';
        fmap.QpQs      = 'QpQs';
        fmap.LVEDV     = 'LVEDV_mL';
        fmap.RVEDV     = 'RVEDV_mL';
        fmap.RAP_mean  = 'RAP_mean_mmHg';
        fmap.PAP_mean  = 'PAP_mean_mmHg';
end
end
