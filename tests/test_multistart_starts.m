function tests = test_multistart_starts()
% TEST_MULTISTART_STARTS
% -----------------------------------------------------------------------
% Contract tests for build_multistart_starts: the seed is preserved, every
% start is inside bounds, the design is deterministic, and single-start
% behaviour is unchanged so enabling multi-start is opt-in.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 3)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function [x0, lb, ub] = box()
lb = [0.25; 0.25; 0.45; 0.75; 0.70];
ub = [2.80; 2.80; 2.80; 1.35; 1.45];
x0 = [1.00; 0.90; 1.10; 1.00; 1.05];
end

function test_single_start_returns_only_the_seed(tc)
[x0, lb, ub] = box();
plan = build_multistart_starts(x0, lb, ub, 1, 20260828);
verifyEqual(tc, plan.num_starts, 1);
verifyEqual(tc, plan.sampler, 'seed_only');
verifyEqual(tc, plan.starts, x0, 'AbsTol', 1e-12);
end

function test_first_start_is_always_the_seed(tc)
% Enabling multi-start must never lose the previously reported operating
% point, so column one stays exactly the incoming seed.
[x0, lb, ub] = box();
plan = build_multistart_starts(x0, lb, ub, 8, 20260828);
verifyEqual(tc, plan.starts(:, 1), x0, 'AbsTol', 1e-12);
end

function test_all_starts_are_inside_bounds(tc)
[x0, lb, ub] = box();
plan = build_multistart_starts(x0, lb, ub, 16, 20260828);
verifyEqual(tc, size(plan.starts), [numel(x0), 16]);
verifyTrue(tc, all(plan.starts >= lb - 1e-12, 'all'), ...
    'No start may sit below its lower bound.');
verifyTrue(tc, all(plan.starts <= ub + 1e-12, 'all'), ...
    'No start may sit above its upper bound.');
end

function test_design_is_deterministic(tc)
[x0, lb, ub] = box();
a = build_multistart_starts(x0, lb, ub, 8, 20260828);
b = build_multistart_starts(x0, lb, ub, 8, 20260828);
verifyEqual(tc, a.starts, b.starts, 'AbsTol', 1e-12, ...
    'Same seed and count must reproduce the same design.');
end

function test_starts_are_distinct(tc)
% A design that collapses onto the seed would silently reduce to one start.
[x0, lb, ub] = box();
plan = build_multistart_starts(x0, lb, ub, 8, 20260828);
verifyEqual(tc, size(unique(plan.starts', 'rows'), 1), 8, ...
    'Every start must be distinct.');
end

function test_design_grows_by_extension_not_reshuffle(tc)
% Sobol is a sequence: asking for more starts must extend the design rather
% than replace it, so a longer run stays comparable with a shorter one.
[x0, lb, ub] = box();
short = build_multistart_starts(x0, lb, ub, 4, 20260828);
long = build_multistart_starts(x0, lb, ub, 8, 20260828);
if strcmp(short.sampler, 'sobol_scrambled')
    verifyEqual(tc, long.starts(:, 1:4), short.starts, 'AbsTol', 1e-12);
end
end

function test_bounds_mismatch_is_rejected(tc)
[x0, lb, ~] = box();
verifyError(tc, @() build_multistart_starts(x0, lb, [1; 2], 4, 1), ...
    'build_multistart_starts:boundsMismatch');
end

function test_reports_which_sampler_was_used(tc)
[x0, lb, ub] = box();
plan = build_multistart_starts(x0, lb, ub, 4, 20260828);
verifyTrue(tc, ismember(plan.sampler, ...
    {'sobol_scrambled', 'uniform_seeded'}), ...
    'The sampler actually used must be reported for provenance.');
verifyEqual(tc, plan.seed, 20260828);
end
