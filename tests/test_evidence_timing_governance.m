function tests = test_evidence_timing_governance()
% TEST_EVIDENCE_TIMING_GOVERNANCE
% -----------------------------------------------------------------------
% Two things are tested here:
%
%   (A) Reyna's chamber volumes (LVEDV/LVESV/RVEDV/RVESV/LVEF/RVEF).
%       CORRECTED TWICE (publication-readiness reconciliation, 2026-09-05
%       then 2026-09-06). Originally believed to be H+1 post-operative
%       echo, excluded from pre-surgery entirely. The 2026-09-05 pass
%       read protocol form rows 26-29 ("PRE RELEASE OCCLUDER") as
%       same-session PRE-surgery evidence instead. That reading was
%       itself wrong, per the study owner: "pre-release occluder" means
%       the closure device is deployed and occluding the defect, just not
%       yet mechanically detached -- i.e. the VSD is already functionally
%       CLOSED at that measurement. So this data belongs to
%       clinical.post_surgery, and there is currently NO confirmed
%       pre-surgery chamber-volume measurement for Reyna at all.
%
%   (B) The generic assert_evidence_timing_governance mechanism itself,
%       using a synthetic fixture decoupled from Reyna's data, so the
%       cross-timing guard stays under test regardless of which way
%       Reyna's own chamber-volume timing gets read.
%
% REFERENCES:
%   [1] config/patient_reyna.m
%   [2] config/calibration_recipes/reyna_pre_surgery.m
%   [3] src/calibration/assert_evidence_timing_governance.m
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-03
% VERSION:  3.0 (rewritten again for the post-closure, not pre-surgery,
%           chamber-volume timing; see git history for prior versions)
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function [clinical, recipe] = reyna_fixture()
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
end

%% ---- (A) Reyna: chamber volumes belong to post_surgery ------------------

function test_reyna_chamber_block_is_unavailable_pre_surgery(tc)
% There is currently no confirmed pre-surgery chamber-volume measurement
% for Reyna at all; patient_reyna() must not fabricate one.
[clinical, ~] = reyna_fixture();
targets = get_calibration_targets('pre_surgery', clinical);

for metric = {'LVEDV','LVESV','RVEDV','RVESV','LVEF'}
    name = metric{1};
    idx = find(strcmp({targets.Metric}, name), 1, 'first');
    verifyNotEmpty(tc, idx, sprintf('%s should still be a known metric.', name));
    verifyFalse(tc, isfinite(targets(idx).ClinicalValue), sprintf( ...
        ['%s must have no pre-surgery clinical comparator: the only ', ...
         'measured chamber-volume evidence is post-closure ("pre-release ', ...
         'occluder" = device deployed and occluding, not yet detached).'], name));
end
end

function test_reyna_chamber_block_is_finite_post_surgery(tc)
% The protocol form's "pre-release occluder" volumes are post-closure
% measurements and belong in clinical.post_surgery.
clinical = patient_reyna();
targets = get_calibration_targets('post_surgery', clinical);

expected = struct('LVEDV', 32.0, 'LVESV', 23.6, 'RVEDV', 30.5, 'RVESV', 12.0);
for metric = fieldnames(expected)'
    name = metric{1};
    idx = find(strcmp({targets.Metric}, name), 1, 'first');
    verifyNotEmpty(tc, idx, sprintf('%s should be a known metric.', name));
    verifyEqual(tc, targets(idx).ClinicalValue, expected.(name), 'AbsTol', 1e-9, ...
        sprintf('%s should carry the protocol post-closure value.', name));
end
end

function test_reyna_chamber_rows_are_unavailable_not_fitted_pre_surgery(tc)
[clinical, ~] = reyna_fixture();
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

for metric = {'LVEDV','LVESV','RVEDV','RVESV','LVEF'}
    name = metric{1};
    row = find(strcmp(tbl.Metric, name), 1, 'first');
    verifyNotEmpty(tc, row, sprintf('%s must appear in the tier table.', name));
    verifyEqual(tc, tbl.Tier{row}, 'unavailable', sprintf( ...
        '%s has no pre-surgery clinical value, so its tier must be unavailable.', name));
    verifyFalse(tc, tbl.IncludedInCalibration(row), sprintf( ...
        '%s must not be fitted pre-surgery.', name));
    verifyFalse(tc, tbl.IncludedInPrimaryRMSE(row), sprintf( ...
        '%s must not count toward the governed pre-surgery primary RMSE.', name));
end
end

function test_reyna_superseded_evidence_is_retained_for_provenance(tc)
% An old, unconfirmed-provenance LV/RV pair (previously mislabelled "H+1
% post-operative echo") must not be silently deleted; it is kept purely as
% a direction-check reference in predicted_chamber_state_report.m, since
% there is no real pre-surgery chamber comparator to use instead.
[~, recipe] = reyna_fixture();
verifyTrue(tc, isfield(recipe, 'excluded_evidence'), ...
    'The historical reference figure must remain documented, not deleted.');
ev = recipe.excluded_evidence;
verifyEqual(tc, ev.LVEDV_mL, 41.0, 'AbsTol', 1e-12);
verifyEqual(tc, ev.LVESV_mL, 19.3, 'AbsTol', 1e-12);
verifyEqual(tc, ev.RVEDV_mL, 30.5, 'AbsTol', 1e-12);
verifyEqual(tc, ev.RVESV_mL, 12.0, 'AbsTol', 1e-12);
verifyEqual(tc, ev.LVEF, 0.528, 'AbsTol', 1e-12);
verifyNotEmpty(tc, ev.reason);
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
