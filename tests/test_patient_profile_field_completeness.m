function tests = test_patient_profile_field_completeness()
% TEST_PATIENT_PROFILE_FIELD_COMPLETENESS
% -----------------------------------------------------------------------
% Regression test for defect D4 (publication-readiness cleanup, Phase 2):
% config/patient_reyna.m was written from a bare struct() rather than
% patient_template(), so it silently lacked common.patient_id,
% common.maturation_mode, pre_surgery.RVEDP_mmHg, post_surgery.LVEDP_mmHg
% and post_surgery.RVEDP_mmHg. Most call sites guard field access with
% isfield() so this was inert rather than crashing, but a missing field
% is a silent trap for the next person extending the file.
%
% This test asserts every clinical.common / pre_surgery / post_surgery
% field defined by patient_template() also exists on patient_reyna().
% It does not assert values, only structural completeness, so Phase 3's
% protocol-data reconciliation (HR, VSD diameter, etc.) is free to change
% values without breaking this test.
%
% REFERENCES:
%   [1] config/patient_template.m
%   [2] config/patient_reyna.m
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-03
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function test_reyna_has_every_template_field(tc)
template = patient_template();
reyna = patient_reyna();

check_section_fields(tc, template.common, reyna.common, 'common');
check_section_fields(tc, template.pre_surgery, reyna.pre_surgery, 'pre_surgery');
check_section_fields(tc, template.post_surgery, reyna.post_surgery, 'post_surgery');
end

function check_section_fields(tc, template_section, actual_section, section_name)
expected_fields = fieldnames(template_section);
for i = 1:numel(expected_fields)
    name = expected_fields{i};
    verifyTrue(tc, isfield(actual_section, name), sprintf( ...
        'patient_reyna() is missing clinical.%s.%s, present in patient_template().', ...
        section_name, name));
end
end
