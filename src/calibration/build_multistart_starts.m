function plan = build_multistart_starts(x_seed, lb, ub, num_starts, seed)
% BUILD_MULTISTART_STARTS
% -----------------------------------------------------------------------
% Deterministic multi-start design over the calibration parameter box.
%
% A 14-parameter calibration run from a single local start cannot separate
% "this prior fits worse" from "the solver never left its start point".
% Sampling several starts makes that difference measurable and turns the
% reported RMSE into a distribution rather than a single unqualified number.
%
% Start 1 is always the incoming seed, so enabling multi-start can never lose
% the previously reported operating point. Starts 2..N come from a scrambled
% Sobol sequence, which covers the box more evenly than independent uniform
% draws at the small sample counts used here.
%
% INPUTS:
%   x_seed     - seed parameter vector                            [n x 1]
%   lb         - lower bounds                                     [n x 1]
%   ub         - upper bounds                                     [n x 1]
%   num_starts - number of starts, >= 1                         [count]
%   seed       - RNG seed for the fallback sampler                   [-]
%
% OUTPUTS:
%   plan       - struct with:
%       .starts      [n x num_starts] start vectors, column 1 = seed
%       .num_starts  realised number of starts                  [count]
%       .seed        RNG seed actually used                        [-]
%       .sampler     'seed_only' | 'sobol_scrambled' | 'uniform_seeded'
%
% ASSUMPTIONS:
%   - Bounds are finite and ub >= lb; violating entries are clamped.
%   - Sobol sampling needs Statistics and Machine Learning Toolbox; without
%     it the function falls back to a seeded uniform design and says so.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 4 || isempty(num_starts) || ~isfinite(num_starts)
    num_starts = 1;
end
if nargin < 5 || isempty(seed) || ~isfinite(seed)
    seed = 20260828;
end

x_seed = x_seed(:);
lb = lb(:);
ub = ub(:);
n_dim = numel(x_seed);
num_starts = max(1, round(num_starts));

plan = struct();
plan.num_starts = num_starts;
plan.seed = round(seed);
plan.sampler = 'seed_only';
plan.starts = repmat(x_seed, 1, num_starts);

if num_starts < 2 || n_dim == 0
    return;
end

if numel(lb) ~= n_dim || numel(ub) ~= n_dim
    error('build_multistart_starts:boundsMismatch', ...
        'lb and ub must match the seed vector length (%d).', n_dim);
end

span = ub - lb;

% scramble() draws from the global stream, so the design is only
% reproducible if that stream is seeded first. Save and restore the caller's
% RNG state so seeding here cannot perturb anything downstream.
rng_state = rng();
restore_rng = onCleanup(@() rng(rng_state));
rng(plan.seed, 'twister');

try
    stream = sobolset(n_dim, 'Skip', 1);
    stream = scramble(stream, 'MatousekAffineOwen');
    unit_pts = net(stream, num_starts - 1);
    plan.sampler = 'sobol_scrambled';
catch
    unit_pts = rand(num_starts - 1, n_dim);
    plan.sampler = 'uniform_seeded';
end

for idx = 1:(num_starts - 1)
    candidate = lb + unit_pts(idx, :)' .* span;
    plan.starts(:, idx + 1) = min(max(candidate, lb), ub);
end
end
