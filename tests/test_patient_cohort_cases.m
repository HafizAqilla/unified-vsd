%% test_patient_cohort_cases.m
% Confirms official cohort includes P1-P9 including Hasna.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

cases = patient_cohort_cases();
labels = {cases.label};

assert(numel(cases) == 9, 'Expected nine official cohort cases.');
assert(strcmp(labels{1}, 'reyna'), 'P1 should be Reyna.');
assert(strcmp(labels{9}, 'hasna_azizah'), 'P9 should be Hasna Azizah.');
assert(all(arrayfun(@(c) isfield(c.clinical, 'common'), cases)), ...
    'Each cohort case should expose clinical.common.');

fprintf('  [PASS] patient_cohort_cases returns P1-P9 including Hasna.\n');

