function tests = test_evidence_timing_governance()
% TEST_EVIDENCE_TIMING_GOVERNANCE
% -----------------------------------------------------------------------
% Evidence recorded at one surgical timepoint must not silently fit a
% scenario at another. For Reyna the LV/RV chamber block is H+1
% post-operative echo, which patient_reyna() marks unavailable for
% pre-surgery fitting.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 5)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function [clinical, recipe] = reyna_fixture()
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
end

function test_reyna_recipe_declares_chamber_block_as_post_operative(tc)
[~, recipe] = reyna_fixture();
verifyTrue(tc, isfield(recipe, 'evidence_timing'), ...
    'Recipe must declare evidence timing for its overrides.');
for field = {'LVEDV_mL','LVESV_mL','RVEDV_mL','RVESV_mL','LVEF'}
    name = field{1};
    verifyTrue(tc, isfield(recipe.evidence_timing, name), ...
        sprintf('%s must declare an evidence timing.', name));
    verifyEqual(tc, recipe.evidence_timing.(name), 'post_operative_H1', ...
        sprintf('%s is H+1 post-operative echo.', name));
end
end

function test_chamber_volumes_are_not_pre_surgery_targets_at_all(tc)
% The block is removed, not merely demoted: no finite pre-surgery comparator
% should exist for any chamber row.
[clinical, ~] = reyna_fixture();
targets = get_calibration_targets('pre_surgery', clinical);

for metric = {'LVEDV','LVESV','RVEDV','RVESV','LVEF','RVEF'}
    name = metric{1};
    idx = find(strcmp({targets.Metric}, name), 1, 'first');
    verifyNotEmpty(tc, idx, sprintf('%s should still be a known metric.', name));
    verifyFalse(tc, isfinite(targets(idx).ClinicalValue), sprintf( ...
        ['%s must have no pre-surgery clinical comparator: it is H+1 ', ...
         'post-operative echo measuring a different loading state.'], name));
end
end

function test_no_post_operative_row_is_fitted_pre_surgery(tc)
[clinical, ~] = reyna_fixture();
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

for metric = {'LVEDV','LVESV','RVEDV','RVESV','LVEF'}
    name = metric{1};
    row = find(strcmp(tbl.Metric, name), 1, 'first');
    verifyNotEmpty(tc, row, sprintf('%s must appear in the tier table.', name));
    verifyFalse(tc, tbl.IncludedInCalibration(row), sprintf( ...
        '%s is post-operative evidence and must not be fitted pre-surgery.', name));
    verifyFalse(tc, tbl.IncludedInPrimaryRMSE(row), sprintf( ...
        '%s must not count toward the governed pre-surgery RMSE.', name));
end
end

function test_excluded_evidence_is_retained_for_provenance(tc)
% Removing the rows as targets must not delete the record of what was
% measured; the values stay documented with their timing and destination.
[~, recipe] = reyna_fixture();
verifyTrue(tc, isfield(recipe, 'excluded_evidence'), ...
    'Excluded chamber evidence must remain documented.');
ev = recipe.excluded_evidence;
verifyEqual(tc, ev.timing, 'post_operative_H1');
verifyEqual(tc, ev.recommended_scenario, 'post_surgery');
verifyEqual(tc, ev.LVEDV_mL, 41.0, 'AbsTol', 1e-12);
verifyEqual(tc, ev.LVESV_mL, 19.3, 'AbsTol', 1e-12);
verifyEqual(tc, ev.RVEDV_mL, 30.5, 'AbsTol', 1e-12);
verifyEqual(tc, ev.RVESV_mL, 12.0, 'AbsTol', 1e-12);
verifyEqual(tc, ev.LVEF, 0.528, 'AbsTol', 1e-12);
verifyNotEmpty(tc, ev.reason);
end

function test_chamber_ic_override_is_disabled_pre_surgery(tc)
% override_IC seeds chamber V0/E from the same post-operative block, a
% second channel independent of target tiers.
[~, recipe] = reyna_fixture();
verifyFalse(tc, recipe.pre_surgery_overrides.override_IC, ...
    'override_IC must stay off so post-operative echo does not seed the fit.');
verifyFalse(tc, recipe.allow_cross_timing_evidence, ...
    'Reyna must not opt in to cross-timing evidence.');
end

function test_governance_reports_no_violations_for_reyna(tc)
[clinical, recipe] = reyna_fixture();
tiers = build_target_tiers(clinical, 'pre_surgery', [], ...
    recipe_tier_config(recipe));
findings = assert_evidence_timing_governance(clinical, 'pre_surgery', ...
    tiers, recipe);

verifyEmpty(tc, findings.violations, sprintf( ...
    'Unexpected cross-timing fits: %s', strjoin(findings.violations, ', ')));
verifyEqual(tc, height(findings.table), 5, ...
    'All five chamber rows should be tracked.');
verifyTrue(tc, all(findings.table.CrossTiming), ...
    'Every chamber row is cross-timing for a pre-surgery run.');
verifyFalse(tc, any(findings.table.Fitted), ...
    'No cross-timing row may carry a fitted tier.');
end

function [clinical, recipe] = reinjected_fixture()
% REINJECTED_FIXTURE - simulate someone putting the chamber block back into
% pre_surgery and tiering RVESV as fitted. The guard must catch it.
[clinical, recipe] = reyna_fixture();
clinical.pre_surgery.RVESV_mL = 12.0;
recipe.soft_metrics = [recipe.soft_metrics, {'RVESV'}];
recipe.consistency_only = setdiff(recipe.consistency_only, {'RVESV'}, 'stable');
end

function test_governance_rejects_a_reinjected_post_operative_row(tc)
% Guard the guard: re-adding RVESV as a fitted target must be refused.
[clinical, bad_recipe] = reinjected_fixture();
tiers = build_target_tiers(clinical, 'pre_surgery', [], ...
    recipe_tier_config(bad_recipe));

verifyError(tc, ...
    @() assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, bad_recipe), ...
    'assert_evidence_timing_governance:crossTimingEvidence');
end

function test_explicit_opt_in_allows_cross_timing(tc)
% The escape hatch must exist, but only when set deliberately.
[clinical, opt_in] = reinjected_fixture();
opt_in.allow_cross_timing_evidence = true;

tiers = build_target_tiers(clinical, 'pre_surgery', [], recipe_tier_config(opt_in));
findings = assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, opt_in);

verifyTrue(tc, findings.allowed);
verifyTrue(tc, ismember('RVESV', findings.violations), ...
    'The opt-in must still report what it permitted.');
end

function test_recipe_without_timing_metadata_is_inert(tc)
[clinical, recipe] = reyna_fixture();
legacy = rmfield(recipe, 'evidence_timing');
tiers = build_target_tiers(clinical, 'pre_surgery', [], recipe_tier_config(legacy));
findings = assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, legacy);
verifyEmpty(tc, findings.violations);
verifyEmpty(tc, findings.table);
end

function config = recipe_tier_config(recipe)
config = struct();
config.policy_name = 'test_recipe_tiers';
config.hard = recipe.primary_metrics;
config.soft = recipe.soft_metrics;
config.consistency_only = recipe.consistency_only;
config.derived_validation = recipe.derived_validation;
config.validation_holdout = recipe.validation_holdout;
config.primary_rmse_holdout = recipe.primary_rmse_holdout;
end
