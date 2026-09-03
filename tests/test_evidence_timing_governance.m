function tests = test_evidence_timing_governance()
% TEST_EVIDENCE_TIMING_GOVERNANCE
% -----------------------------------------------------------------------
% Two things are tested here:
%
%   (A) Reyna's pre-surgery chamber volumes (LVEDV/LVESV/RVEDV/RVESV/LVEF).
%       CORRECTED (publication-readiness reconciliation, 2026-09):
%       these were previously believed to be H+1 POST-operative echo and
%       were excluded from pre-surgery entirely. The IRB-governed protocol
%       form contradicts that: rows 26-29 report them as PRE-release
%       (same pre-surgery catheterisation session), so config/patient_reyna.m
%       now carries them as finite values. They still must never be FITTED
%       — the LV pair is internally implausible (SV_LV = 8.4 mL vs ~34 mL
%       implied by protocol Qp) — so they are governed as consistency-only
%       via recipe.consistency_only, not via a (now factually wrong)
%       cross-timing exclusion.
%
%   (B) The generic assert_evidence_timing_governance mechanism itself,
%       using a synthetic fixture decoupled from Reyna's corrected data, so
%       the cross-timing guard is still exercised even though it is no
%       longer the reason Reyna's chamber block is unfitted.
%
% REFERENCES:
%   [1] config/patient_reyna.m
%   [2] config/calibration_recipes/reyna_pre_surgery.m
%   [3] src/calibration/assert_evidence_timing_governance.m
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-03
% VERSION:  2.0 (rewritten for the corrected pre-release chamber-volume
%           timing; see git history for the pre-rewrite version)
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function [clinical, recipe] = reyna_fixture()
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
end

%% ---- (A) Reyna: consistency-only, not cross-timing ----------------------

function test_reyna_chamber_block_is_finite_pre_surgery(tc)
% The protocol form's pre-release volumes must now be visible, not NaN'd
% out by either patient_reyna() or the recipe overrides.
[clinical, ~] = reyna_fixture();
targets = get_calibration_targets('pre_surgery', clinical);

expected = struct('LVEDV', 32.0, 'LVESV', 23.6, 'RVEDV', 30.5, 'RVESV', 12.0);
for metric = fieldnames(expected)'
    name = metric{1};
    idx = find(strcmp({targets.Metric}, name), 1, 'first');
    verifyNotEmpty(tc, idx, sprintf('%s should be a known metric.', name));
    verifyEqual(tc, targets(idx).ClinicalValue, expected.(name), 'AbsTol', 1e-9, ...
        sprintf('%s should carry the protocol pre-release value.', name));
end
end

function test_reyna_chamber_rows_are_consistency_only_not_fitted(tc)
[clinical, ~] = reyna_fixture();
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

for metric = {'LVEDV','LVESV','RVEDV','RVESV','LVEF'}
    name = metric{1};
    row = find(strcmp(tbl.Metric, name), 1, 'first');
    verifyNotEmpty(tc, row, sprintf('%s must appear in the tier table.', name));
    verifyEqual(tc, tbl.Tier{row}, 'consistency_check_only', sprintf( ...
        '%s must be governed as consistency-only, not silently fitted.', name));
    verifyFalse(tc, tbl.IncludedInCalibration(row), sprintf( ...
        '%s must not be fitted pre-surgery (internally implausible LV pair).', name));
    verifyFalse(tc, tbl.IncludedInPrimaryRMSE(row), sprintf( ...
        '%s must not count toward the governed pre-surgery primary RMSE.', name));
end
end

function test_reyna_superseded_evidence_is_retained_for_provenance(tc)
% The previous (unconfirmed) H+1-labelled LV pair must not be silently
% deleted just because it was superseded by the protocol form. It is kept
% as recipe.excluded_evidence purely for the direction-consistency check in
% predicted_chamber_state_report.m; it is no longer an exclusion record.
[~, recipe] = reyna_fixture();
verifyTrue(tc, isfield(recipe, 'excluded_evidence'), ...
    'Superseded chamber evidence must remain documented, not deleted.');
ev = recipe.excluded_evidence;
verifyEqual(tc, ev.LVEDV_mL, 41.0, 'AbsTol', 1e-12);
verifyEqual(tc, ev.LVESV_mL, 19.3, 'AbsTol', 1e-12);
verifyEqual(tc, ev.RVEDV_mL, 30.5, 'AbsTol', 1e-12);
verifyEqual(tc, ev.RVESV_mL, 12.0, 'AbsTol', 1e-12);
verifyEqual(tc, ev.LVEF, 0.528, 'AbsTol', 1e-12);
verifyNotEmpty(tc, ev.reason);
verifyEqual(tc, ev.recommended_scenario, 'pre_surgery', ...
    'Unlike the old H+1 framing, this evidence now belongs to pre_surgery, not post_surgery.');
end

function test_reyna_chamber_ic_override_is_disabled_pre_surgery(tc)
% override_IC seeds chamber V0/E from the chamber block, a channel
% independent of target tiers; it must stay off regardless of tiering.
[~, recipe] = reyna_fixture();
verifyFalse(tc, recipe.pre_surgery_overrides.override_IC, ...
    'override_IC must stay off so the consistency-only volumes do not seed the fit.');
end

function test_reyna_vsd_diameter_is_not_recipe_overridden(tc)
% D1 regression: the recipe used to force VSD_diameter_mm to 3.025 mm over
% patient_reyna()'s protocol-sourced 3.665 mm. That override must be gone.
[~, recipe] = reyna_fixture();
verifyFalse(tc, isfield(recipe.pre_surgery_overrides, 'VSD_diameter_mm'), ...
    'The recipe must not override VSD_diameter_mm; patient_reyna() (3.665 mm) must win.');

clinical = patient_reyna();
verifyEqual(tc, clinical.pre_surgery.VSD_diameter_mm, 3.665, 'AbsTol', 1e-9);
end

%% ---- (B) Generic assert_evidence_timing_governance mechanism -------------

function [clinical, recipe] = synthetic_cross_timing_fixture()
% A minimal, Reyna-independent fixture that genuinely mixes evidence
% timepoints, to keep the cross-timing guard itself under test even though
% Reyna's own data no longer exercises it.
clinical = patient_template();
clinical.common.patient_name = 'synthetic_cross_timing_case';
clinical.pre_surgery.RAP_mean_mmHg = 5;
clinical.pre_surgery.PAP_mean_mmHg = 15;
clinical.pre_surgery.SAP_mean_mmHg = 70;
clinical.pre_surgery.QpQs = 1.2;
clinical.pre_surgery.CO_Lmin = 3.5;
% This RVESV is deliberately evidence from a different timepoint than the
% pre_surgery scenario being fitted.
clinical.pre_surgery.RVESV_mL = 12.0;

recipe = struct();
recipe.primary_metrics = {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'};
recipe.soft_metrics = {'RVESV'};
recipe.consistency_only = {};
recipe.derived_validation = {};
recipe.validation_holdout = {};
recipe.primary_rmse_holdout = {};
recipe.evidence_timing = struct('RVESV_mL', 'post_operative_H1');
recipe.allow_cross_timing_evidence = false;
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

function test_governance_rejects_a_cross_timing_fitted_row(tc)
[clinical, recipe] = synthetic_cross_timing_fixture();
tiers = build_target_tiers(clinical, 'pre_surgery', [], recipe_tier_config(recipe));

verifyError(tc, ...
    @() assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, recipe), ...
    'assert_evidence_timing_governance:crossTimingEvidence');
end

function test_explicit_opt_in_allows_cross_timing(tc)
[clinical, recipe] = synthetic_cross_timing_fixture();
recipe.allow_cross_timing_evidence = true;

tiers = build_target_tiers(clinical, 'pre_surgery', [], recipe_tier_config(recipe));
findings = assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, recipe);

verifyTrue(tc, findings.allowed);
verifyTrue(tc, ismember('RVESV', findings.violations), ...
    'The opt-in must still report what it permitted.');
end

function test_recipe_without_timing_metadata_is_inert(tc)
[clinical, recipe] = synthetic_cross_timing_fixture();
legacy = rmfield(recipe, 'evidence_timing');
tiers = build_target_tiers(clinical, 'pre_surgery', [], recipe_tier_config(legacy));
findings = assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, legacy);
verifyEmpty(tc, findings.violations);
verifyEmpty(tc, findings.table);
end

function test_governance_is_inert_for_reyna_since_no_cross_timing_remains(tc)
% Reyna's recipe no longer declares evidence_timing at all (nothing is
% cross-timing anymore), so the generic mechanism must be a silent no-op
% for her, not a false positive or false negative.
[clinical, recipe] = reyna_fixture();
verifyFalse(tc, isfield(recipe, 'evidence_timing'), ...
    'Reyna''s chamber block is same-session evidence; it should not declare a cross-timing tag.');

tiers = build_target_tiers(clinical, 'pre_surgery', [], recipe_tier_config(recipe));
findings = assert_evidence_timing_governance(clinical, 'pre_surgery', tiers, recipe);
verifyEmpty(tc, findings.violations);
verifyEmpty(tc, findings.table);
end
