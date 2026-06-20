function cases = patient_cohort_cases()
% PATIENT_COHORT_CASES
% -----------------------------------------------------------------------
% Returns the available VSD patient profiles in a reproducible order.
%
% INPUTS:
%   none
%
% OUTPUTS:
%   cases - struct array with labels, clinical structs, and run flags    [-]
%
% REFERENCES:
%   [1] docs/clinical_data_dictionary.md
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-15
% VERSION:  1.0
% -----------------------------------------------------------------------

profiles = {
    'reyna',       @patient_reyna
    'razka',       @patient_profile_Razka
    'fathan',      @patient_fathan
    'azzam',       @patient_azzam
    'ali_zhafran', @patient_ali_zhafran
    'salman',      @patient_salman
    'syabil',      @patient_syabil
    'jericho',     @patient_jericho
    'hasna_azizah', @patient_hasna_azizah
    };

cases = repmat(struct('label', '', 'clinical', [], 'scenario', 'pre_surgery', ...
    'can_simulate', false, 'skip_reason', ''), size(profiles, 1), 1);
for idx = 1:size(profiles, 1)
    clinical = profiles{idx, 2}();
    cases(idx).label = profiles{idx, 1};
    cases(idx).clinical = clinical;
    cases(idx).scenario = 'pre_surgery';
    [cases(idx).can_simulate, cases(idx).skip_reason] = can_simulate_case(clinical);
end
end

function [tf, reason] = can_simulate_case(clinical)
% CAN_SIMULATE_CASE - check minimum anthropometry required for scaling.
tf = isfinite(clinical.common.weight_kg) && clinical.common.weight_kg > 0 && ...
    isfinite(clinical.common.height_cm) && clinical.common.height_cm > 0;
if tf
    reason = '';
else
    reason = 'Missing weight and/or height; cannot apply pediatric scaling.';
end
end
