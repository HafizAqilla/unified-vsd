function tests = test_measurement_noise_draws()
% TEST_MEASUREMENT_NOISE_DRAWS
% -----------------------------------------------------------------------
% Contract tests for draw_measurement_noise: the draws are reproducible,
% only measured pre-closure targets move, derived quantities follow the
% record's own formulas, inputs and the held-out post-closure record are
% untouched, and the spread matches the declared uncertainty.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-26
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function test_same_seed_same_draws(tc)
c = patient_reyna();
[~, a] = draw_measurement_noise(c, 5, 7);
[~, b] = draw_measurement_noise(c, 5, 7);
verifyEqual(tc, a, b);
end

function test_inputs_and_post_record_unchanged(tc)
c = patient_reyna();
recs = draw_measurement_noise(c, 3, 11);
for k = 1:3
    r = recs{k};
    verifyEqual(tc, r.common, c.common);
    verifyEqual(tc, r.post_surgery, c.post_surgery);
    verifyEqual(tc, r.pre_surgery.HR, c.pre_surgery.HR);
    verifyEqual(tc, r.pre_surgery.VSD_diameter_mm, c.pre_surgery.VSD_diameter_mm);
end
end

function test_derived_quantities_recomputed(tc)
c = patient_reyna();
recs = draw_measurement_noise(c, 4, 13);
for k = 1:4
    s = recs{k}.pre_surgery;
    verifyEqual(tc, s.SAP_mean_mmHg, (s.SAP_sys_mmHg + 2 * s.SAP_dia_mmHg) / 3, 'AbsTol', 1e-12);
    verifyEqual(tc, s.Q_shunt_Lmin, s.CO_Lmin * (s.QpQs - 1), 'AbsTol', 1e-12);
    verifyGreaterThan(tc, s.Q_shunt_Lmin, 0);
    verifyLessThan(tc, s.PAP_dia_mmHg, s.PAP_mean_mmHg);
    verifyLessThan(tc, s.PAP_mean_mmHg, s.PAP_sys_mmHg);
end
end

function test_spread_matches_declared_sigma(tc)
c = patient_reyna();
[~, draws, spec] = draw_measurement_noise(c, 4000, 17);
for i = 1:height(spec)
    x = draws.(spec.ClinicalField{i});
    verifyEqual(tc, mean(x), spec.Measured(i), 'AbsTol', 4 * spec.Sigma(i) / sqrt(4000));
    verifyEqual(tc, std(x), spec.Sigma(i), 'RelTol', 0.06);
end
end
