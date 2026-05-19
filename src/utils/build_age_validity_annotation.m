function age_validity = build_age_validity_annotation(age_years, scaling_mode, maturation_mode)
% BUILD_AGE_VALIDITY_ANNOTATION
% -----------------------------------------------------------------------
% Tags a run with the age-validity regime implied by the current pediatric
% scaling evidence. This does not alter parameters; it documents how much
% literature support exists for the age range being simulated.
%
% INPUTS:
%   age_years       - patient age                                      [years]
%   scaling_mode    - scaling mode string                              [-]
%   maturation_mode - maturation mode string                           [-]
%
% OUTPUTS:
%   age_validity    - struct with regime/evidence/summary fields       [-]
%
% REFERENCES:
%   [1] Bozkurt S. PLoS ONE 2019;14(10):e0224663.
%   [2] Hiebing AA et al. J Biomech Eng. 2023;145(10).
%   [3] Colebank MJ et al. arXiv:2602.04008.
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-08
% VERSION:  1.0
% -----------------------------------------------------------------------

if nargin < 1 || isempty(age_years) || ~isfinite(age_years)
    age_years = NaN;
end
if nargin < 2 || isempty(scaling_mode)
    scaling_mode = 'unknown';
end
if nargin < 3 || isempty(maturation_mode)
    maturation_mode = 'unknown';
end

age_validity = struct();
age_validity.age_years = age_years;
age_validity.scaling_mode = char(string(scaling_mode));
age_validity.maturation_mode = char(string(maturation_mode));
age_validity.reference_windows = { ...
    'Infant maturation evidence: birth to ~2 years (Hiebing/Colebank lineage).', ...
    'Generic child geometry/hemodynamics anchor: BSA~1.0 m^2, ~8-12 years (Bozkurt 2019).'};

if isnan(age_years)
    age_validity.regime = 'unknown_age';
    age_validity.evidence_level = 'unknown';
    age_validity.summary = 'Age not available; cannot annotate scaling validity regime.';
    age_validity.action = 'report_only';
    return;
end

if age_years <= 2.0
    age_validity.regime = 'infant_supported';
    age_validity.evidence_level = 'moderate';
    age_validity.summary = ['Age lies within the strongest support window for infant maturation-based ', ...
        'scaling and VSD adaptation studies.'];
    age_validity.action = 'use_infant_scaling_with_standard_caution';
elseif age_years < 8.0
    age_validity.regime = 'preschool_extrapolated';
    age_validity.evidence_level = 'low';
    age_validity.summary = ['Age lies between the infant maturation literature and the generic ', ...
        'school-age child anchor, so scaling is an extrapolated prior rather than a directly validated regime.'];
    age_validity.action = 'treat_as_extrapolated_and_report_uncertainty';
elseif age_years <= 12.0
    age_validity.regime = 'generic_child_supported';
    age_validity.evidence_level = 'moderate';
    age_validity.summary = ['Age overlaps the generic child regime used in Bozkurt 2019, but still ', ...
        'requires caution because that study was not a patient-specific congenital calibration framework.'];
    age_validity.action = 'use_generic_child_anchor_plus_case_specific_uncertainty';
else
    age_validity.regime = 'adolescent_extrapolated';
    age_validity.evidence_level = 'low';
    age_validity.summary = ['Age exceeds the main infant and generic-child support windows used by the ', ...
        'current scaling literature chain, so the baseline should be treated as an extrapolated prior.'];
    age_validity.action = 'report_as_extrapolated_and_prior_sensitive';
end
end
