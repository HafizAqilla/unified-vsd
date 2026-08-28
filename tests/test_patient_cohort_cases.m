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
assert(isfield(cases(1).clinical, 'common'), 'Reyna should expose clinical.common.');
for idx = 1:numel(cases)
    if isempty(cases(idx).clinical)
        assert(~isempty(cases(idx).skip_reason), ...
            'Unavailable cohort profiles must have an explicit skip reason.');
    else
        assert(isfield(cases(idx).clinical, 'common'), ...
            'Available cohort profiles should expose clinical.common.');
    end
end

fprintf('  [PASS] patient_cohort_cases returns P1-P9 including Hasna.\n');
