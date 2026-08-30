function tests = test_no_ungoverned_calibration_targets()
% TEST_NO_UNGOVERNED_CALIBRATION_TARGETS
% -----------------------------------------------------------------------
% A clinical target declared UseForCalibration with a finite value must be
% fitted, or excluded on purpose. Silently landing in 'validation_only'
% removes it from the objective while leaving it inside the governed primary
% RMSE, which scores the optimiser on a metric it never fits.
%
% REFERENCES:
%   [1] docs/reyna_zhang_full_metric_prd.md (Phase 1, Defect A)
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-08-28
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function profiles = patient_profile_handles()
profiles = {@patient_reyna, @patient_profile_Razka, @patient_profile_A, ...
    @patient_azzam, @patient_fathan, @patient_jericho, @patient_syabil, ...
    @patient_ali_zhafran, @patient_hasna_azizah, @patient_ibnu_sina, ...
    @patient_salman, @patient_zoya, @patient_healthy_adult};
end

function test_every_patient_profile_governs_its_calibration_targets(tc)
scenarios = {'pre_surgery', 'post_surgery'};
handles = patient_profile_handles();

for h = 1:numel(handles)
    clinical = handles{h}();
    label = func2str(handles{h});
    for s = 1:numel(scenarios)
        scenario = scenarios{s};
        % build_target_tiers asserts internally; a throw is the failure mode.
        tiers = build_target_tiers(clinical, scenario);
        targets = get_calibration_targets(scenario, clinical);

        for idx = 1:numel(targets)
            if ~targets(idx).UseForCalibration || ...
                    ~isfinite(targets(idx).ClinicalValue)
                continue;
            end
            row = find(strcmp(tiers.table.Metric, targets(idx).Metric), 1, 'first');
            verifyNotEmpty(tc, row, sprintf( ...
                '%s/%s: %s missing from the tier table.', ...
                label, scenario, targets(idx).Metric));
            verifyNotEqual(tc, tiers.table.Tier{row}, 'validation_only', ...
                sprintf('%s/%s: %s is UseForCalibration but ungoverned.', ...
                label, scenario, targets(idx).Metric));
        end
    end
end
end

function test_reyna_pulmonary_pressures_are_fitted(tc)
% The specific defect this phase closes: PAP_max/PAP_min were graded inside
% the governed primary RMSE while sitting outside the objective.
clinical = patient_reyna();
[recipe, found] = load_calibration_recipe(clinical, 'pre_surgery');
verifyTrue(tc, found, 'Reyna pre_surgery recipe must exist.');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);

profile = build_case_calibration_profile(clinical, 'pre_surgery');
tiers = profile.targetTiers;

for metric = {'PAP_max', 'PAP_min'}
    name = metric{1};
    row = find(strcmp(tiers.table.Metric, name), 1, 'first');
    verifyNotEmpty(tc, row, sprintf('%s must appear in the tier table.', name));
    verifyEqual(tc, tiers.table.Tier{row}, 'soft', ...
        sprintf('%s must be a soft calibration target.', name));
    verifyTrue(tc, tiers.table.IncludedInCalibration(row), ...
        sprintf('%s must enter the objective.', name));
    verifyTrue(tc, tiers.table.IncludedInPrimaryRMSE(row), ...
        sprintf('%s must remain inside the governed primary RMSE.', name));
end
end

function test_no_metric_is_graded_without_being_fitted(tc)
% Every row inside the governed primary RMSE must also be inside the
% objective, unless it is an explicitly declared report-only tier.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

graded_not_fitted = tbl.IncludedInPrimaryRMSE & ~tbl.IncludedInCalibration;
verifyFalse(tc, any(graded_not_fitted), sprintf( ...
    'Metrics graded but never fitted: %s', ...
    strjoin(tbl.Metric(graded_not_fitted), ', ')));
end

function test_recipe_primary_rmse_holdout_is_honoured(tc)
% Q_shunt_Lmin is algebraically CO_Lmin * (QpQs - 1); it is fitted as a
% consistency term but must not be counted as an independent RMSE row.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

row = find(strcmp(tbl.Metric, 'Q_shunt_Lmin'), 1, 'first');
verifyNotEmpty(tc, row);
verifyTrue(tc, tbl.IncludedInCalibration(row), ...
    'Q_shunt_Lmin should still inform the objective.');
verifyFalse(tc, tbl.IncludedInPrimaryRMSE(row), ...
    'Q_shunt_Lmin is a derived identity and must be held out of primary RMSE.');
end

function test_allowed_metric_fields_match_the_fitted_tiers(tc)
% The tier table says what should be fitted; allowedMetricFields decides what
% the objective actually sees. A disagreement means the governance table is
% describing a fit that is not happening.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

should_fit = tbl.Metric(logical(tbl.IncludedInCalibration))';
missing = setdiff(should_fit, profile.allowedMetricFields, 'stable');
verifyEmpty(tc, missing, sprintf( ...
    'Tiered as fitted but absent from allowedMetricFields: %s', ...
    strjoin(missing, ', ')));
end

function test_rmse_holdout_still_reaches_the_objective(tc)
% primary_rmse_holdout governs the reported RMSE, not the objective.
% Q_shunt_Lmin must stay fitted as a consistency term while being excluded
% from the headline count.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');

verifyTrue(tc, ismember('Q_shunt_Lmin', profile.allowedMetricFields), ...
    'An RMSE holdout must still be available to the objective.');

tbl = profile.targetTiers.table;
row = find(strcmp(tbl.Metric, 'Q_shunt_Lmin'), 1, 'first');
verifyFalse(tc, tbl.IncludedInPrimaryRMSE(row), ...
    'The holdout must still be excluded from the governed RMSE.');
end

function test_assertion_fires_on_an_ungoverned_target(tc)
% Guard the guard: an empty soft list must be rejected, not silently
% downgraded to validation_only.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);

bad_config = struct();
bad_config.policy_name = 'deliberately_ungoverned';
bad_config.hard = {'CO_Lmin'};
bad_config.soft = {};
bad_config.consistency_only = {};
bad_config.derived_validation = {};
bad_config.validation_holdout = {};
bad_config.primary_rmse_holdout = {};

verifyError(tc, ...
    @() build_target_tiers(clinical, 'pre_surgery', [], bad_config), ...
    'build_target_tiers:ungovernedCalibrationTarget');
end

function test_svr_is_excluded_from_calibration_and_primary_rmse(tc)
% PRD reyna_statistical_calibration_v1 Phase 5 asked to relabel SVR as
% validation_holdout, on the premise that it is "a genuine prediction test
% rather than a fitted result". That premise does not survive checking the
% derivation: SVR_target = (SAP_mean - RAP_mean) / CO_Lmin
% (objective_calibration.m build_systemic_bundle, line ~395), and all three
% of SAP_mean, RAP_mean, CO_Lmin are hard-tier targets already fitted. SVR
% carries no information independent of what is already in the objective,
% so relabelling it validation_holdout would misrepresent algebraic closure
% as generalisation -- the same error already correctly avoided for
% Q_shunt_Lmin. SVR therefore stays derived_validation; this test locks in
% the exclusion outcome both tiers share (excluded from fitting and from
% the governed RMSE) without asserting the (deviated-from) holdout label.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
clinical = apply_calibration_recipe_to_clinical(clinical, 'pre_surgery', recipe);
profile = build_case_calibration_profile(clinical, 'pre_surgery');
tbl = profile.targetTiers.table;

row = find(strcmp(tbl.Metric, 'SVR'), 1, 'first');
verifyNotEmpty(tc, row, 'SVR must appear in the tier table.');
verifyFalse(tc, tbl.IncludedInCalibration(row), ...
    'SVR must not be fitted: it is algebra over three already-fitted targets.');
verifyFalse(tc, tbl.IncludedInPrimaryRMSE(row), ...
    'SVR must not count toward the governed primary RMSE.');
verifyEqual(tc, tbl.Tier{row}, 'derived_validation', ...
    'SVR is a derived quantity, not an independent holdout candidate.');
end

function test_no_current_metric_is_designated_a_genuine_validation_holdout(tc)
% Locks in the honest state of Phase 5: nothing in the current Reyna
% pre-surgery recipe is a genuinely independent held-out prediction test.
% This is expected to change only when a metric is identified that is (a)
% an independent measurement, not an algebraic function of other fitted
% targets, and (b) deliberately excluded from the objective for that
% reason. If this test starts failing because recipe.validation_holdout is
% no longer empty, verify the new entry is actually independent before
% updating this assertion.
clinical = patient_reyna();
[recipe, ~] = load_calibration_recipe(clinical, 'pre_surgery');
verifyEmpty(tc, recipe.validation_holdout);
end
