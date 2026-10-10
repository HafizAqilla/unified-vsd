function [records, draws, spec] = draw_measurement_noise(clinical, n_replicates, noise_seed)
% DRAW_MEASUREMENT_NOISE
% -----------------------------------------------------------------------
% Builds n_replicates perturbed copies of the pre-closure clinical record
% for the measurement-noise (Monte Carlo) simulation.
%
% Every directly measured pre-closure target is perturbed independently by
% zero-mean Gaussian noise whose standard deviation is that target's
% declared 1-sigma measurement uncertainty, resolved exactly as the
% chi-squared report resolves it (UncertaintyAbs, else UncertaintyFraction
% times the measured value). Quantities that the clinical record itself
% derives from other measurements are recomputed from the perturbed values
% with the same formulas the record uses, so each replicate is an
% internally consistent record:
%   mean arterial pressure  MAP     = (SAP_sys + 2 SAP_dia) / 3  (protocol row 20)
%   shunt flow              Q_shunt = Qs (Qp/Qs - 1)             (= Qp - Qs)
% SVR is left NaN in the record, as in the original, and is derived
% downstream by get_calibration_targets.
%
% Model inputs (age, weight, height, BSA, heart rate, defect diameter) are
% NOT perturbed; neither is any post-closure value (the held-out set).
%
% A draw is rejected and redrawn only if it is physically impossible:
% a non-positive value, diastolic >= systolic, a PA mean outside its own
% diastolic-systolic range, or Qp/Qs <= 1 (a left-to-right shunt was
% measured, and the discharge-coefficient seed needs positive shunt flow).
% The number of rejections is returned so it can be reported.
%
% INPUTS:
%   clinical      - unified clinical struct (patient_reyna())            [-]
%   n_replicates  - number of perturbed records                      [count]
%   noise_seed    - rng seed for the draws                               [-]
%
% OUTPUTS:
%   records - 1 x n cell array of clinical structs (pre_surgery perturbed)
%   draws   - table, one row per replicate: perturbed and derived values,
%             plus the number of rejected draws before acceptance
%   spec    - table, one row per perturbed target: field, measured value,
%             sigma, sigma source
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-26
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 3 || isempty(noise_seed)
    noise_seed = 20260926;
end
validateattributes(n_replicates, {'numeric'}, {'scalar', 'integer', 'positive'});

% Directly measured pre-closure targets that are perturbed.
% {metric in get_calibration_targets, field in clinical.pre_surgery}
perturbed = {
    'RAP_mean', 'RAP_mean_mmHg'
    'PAP_mean', 'PAP_mean_mmHg'
    'PAP_min',  'PAP_dia_mmHg'
    'PAP_max',  'PAP_sys_mmHg'
    'SAP_min',  'SAP_dia_mmHg'
    'SAP_max',  'SAP_sys_mmHg'
    'CO_Lmin',  'CO_Lmin'
    'QpQs',     'QpQs'
    };

targets = get_calibration_targets('pre_surgery', clinical);
target_metrics = {targets.Metric};
n_p = size(perturbed, 1);
mu = nan(n_p, 1);
sigma = nan(n_p, 1);
sigma_source = cell(n_p, 1);
for i = 1:n_p
    t = targets(strcmp(target_metrics, perturbed{i, 1}));
    assert(numel(t) == 1, 'draw_measurement_noise:missingTarget', ...
        'Target %s not found.', perturbed{i, 1});
    mu(i) = clinical.pre_surgery.(perturbed{i, 2});
    assert(isfinite(mu(i)), 'draw_measurement_noise:missingValue', ...
        'Measured value for %s is missing; it cannot be perturbed.', perturbed{i, 2});
    [sigma(i), sigma_source{i}] = resolve_sigma(t, mu(i));
end
spec = table(perturbed(:, 1), perturbed(:, 2), mu, sigma, 100 * sigma ./ abs(mu), sigma_source, ...
    'VariableNames', {'Metric', 'ClinicalField', 'Measured', 'Sigma', 'SigmaPct', 'SigmaSource'});

stream = RandStream('mt19937ar', 'Seed', noise_seed);
values = nan(n_replicates, n_p);
rejected = zeros(n_replicates, 1);
for k = 1:n_replicates
    while true
        x = mu + sigma .* randn(stream, n_p, 1);
        if is_physical(x, perturbed(:, 2))
            break;
        end
        rejected(k) = rejected(k) + 1;
        assert(rejected(k) < 1000, 'draw_measurement_noise:rejectionLoop', ...
            'Could not draw a physical record after 1000 attempts.');
    end
    values(k, :) = x';
end

records = cell(1, n_replicates);
map_col = nan(n_replicates, 1);
shunt_col = nan(n_replicates, 1);
for k = 1:n_replicates
    c = clinical;
    for i = 1:n_p
        c.pre_surgery.(perturbed{i, 2}) = values(k, i);
    end
    c.pre_surgery.SAP_mean_mmHg = (c.pre_surgery.SAP_sys_mmHg + 2 * c.pre_surgery.SAP_dia_mmHg) / 3;
    c.pre_surgery.Q_shunt_Lmin = c.pre_surgery.CO_Lmin * (c.pre_surgery.QpQs - 1);
    c.pre_surgery.noise_replicate = k;
    c.pre_surgery.noise_seed = noise_seed;
    map_col(k) = c.pre_surgery.SAP_mean_mmHg;
    shunt_col(k) = c.pre_surgery.Q_shunt_Lmin;
    records{k} = c;
end

draws = array2table(values, 'VariableNames', perturbed(:, 2)');
draws = addvars(draws, (1:n_replicates)', 'Before', 1, 'NewVariableNames', 'Replicate');
draws.SAP_mean_mmHg = map_col;
draws.Q_shunt_Lmin = shunt_col;
draws.RejectedDraws = rejected;
end

% =======================================================================
function [sigma, source] = resolve_sigma(t, value)
if isfinite(t.UncertaintyAbs) && t.UncertaintyAbs > 0
    sigma = t.UncertaintyAbs;
    source = 'declared absolute';
elseif isfinite(t.UncertaintyFraction) && t.UncertaintyFraction > 0
    sigma = t.UncertaintyFraction * abs(value);
    source = sprintf('declared %.0f%% (%s reliability)', 100 * t.UncertaintyFraction, t.Reliability);
else
    error('draw_measurement_noise:noSigma', 'No declared uncertainty for %s.', t.Metric);
end
end

function ok = is_physical(x, fields)
v = @(name) x(strcmp(fields, name));
ok = all(x > 0) ...
    && v('PAP_dia_mmHg') < v('PAP_mean_mmHg') ...
    && v('PAP_mean_mmHg') < v('PAP_sys_mmHg') ...
    && v('SAP_dia_mmHg') < v('SAP_sys_mmHg') ...
    && v('QpQs') > 1;
end
