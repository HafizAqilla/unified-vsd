%% test_baseline_provenance_and_age_validity.m
% Smoke tests for provenance export helpers and age-regime annotation.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
restoredefaultpath();
addpath(genpath(project_root));

fprintf('==========================================\n');
fprintf('  UNIFIED VSD MODEL - Provenance/Age Test\n');
fprintf('==========================================\n\n');

params_ref = default_parameters();
provenance = build_baseline_provenance(params_ref);

required_cols = {'Parameter','Group','Unit','CurrentValue', ...
    'FoundationSource','EvidenceScope','ConfidenceLevel','Notes'};
assert(all(ismember(required_cols, provenance.Properties.VariableNames)), ...
    'Baseline provenance table is missing required columns.');
assert(any(strcmp(provenance.Parameter, 'R.SAR')), ...
    'Baseline provenance table should include R.SAR.');
assert(any(strcmp(provenance.Parameter, 'R.vsd')), ...
    'Baseline provenance table should include R.vsd.');

age_reyna = build_age_validity_annotation(3.17, 'lundquist_bsa', 'normal');
assert(strcmp(age_reyna.regime, 'preschool_extrapolated'), ...
    'Reyna-age example should map to preschool_extrapolated.');
assert(strcmp(age_reyna.evidence_level, 'low'), ...
    'Preschool extrapolated regime should carry low evidence level.');

age_infant = build_age_validity_annotation(0.75, 'lundquist_bsa', 'normal');
assert(strcmp(age_infant.regime, 'infant_supported'), ...
    'Infant example should map to infant_supported.');

age_child = build_age_validity_annotation(9.0, 'lundquist_bsa', 'normal');
assert(strcmp(age_child.regime, 'generic_child_supported'), ...
    'School-age child example should map to generic_child_supported.');

fprintf('[PASS] Baseline provenance and age-validity helpers are structurally consistent.\n');
