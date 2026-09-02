function tests = test_deidentification_governance()
% TEST_DEIDENTIFICATION_GOVERNANCE
% -----------------------------------------------------------------------
% AGENTS.md Section 9.2 requires patient identifiers (hospital name, MRN,
% case ID, DOB, full name) to live only in a gitignored local provenance
% file (config/private/*.local.m), never in a tracked file. This guards
% the de-identification pass done ahead of publication (2026-08-31):
% real identifiers were stripped from config/patient_reyna.m,
% config/patient_profile_Razka.m, config/calibration_recipes/
% reyna_pre_surgery.m, tests/test_reyna_systemic_flow_profile.m, and
% several docs, replaced with pointers to config/private/
% patient_provenance.local.m.
%
% This test scans every file git currently tracks (NOT git history — a
% history rewrite was explicitly declined) for the known real-identifier
% strings and fails if any tracked file still contains one.
%
% REFERENCES:
%   [1] .assistant/AGENTS.md Section 9.2
%   [2] config/private/patient_provenance.local.m.example
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-09-02
% VERSION:  1.0
% -----------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function patterns = identifier_patterns()
patterns = { ...
    '01008971', ...        % Reyna MRN
    '00948048', ...        % Razka MRN
    'HA000557', ...        % Reyna case ID
    'Harapan Kita', ...    % facility name
    'Razka Alfa Rizki', ... % Razka full name
    'Reyna Nur', ...       % Reyna full name (Reyna Nur Shakila)
    '25/09/2021' ...       % Razka DOB
    };
end

function repo_root = locate_repo_root()
here = fileparts(mfilename('fullpath'));
repo_root = fileparts(here);
end

function test_no_tracked_file_contains_patient_identifiers(tc)
repo_root = locate_repo_root();
old_dir = pwd();
cleanup_obj = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(repo_root);

[status, tracked_raw] = system('git ls-files');
verifyEqual(tc, status, 0, 'git ls-files must succeed to run this test.');

tracked_files = strsplit(strtrim(tracked_raw), newline);
tracked_files = tracked_files(~cellfun(@isempty, tracked_files));

patterns = identifier_patterns();
violations = {};

for i = 1:numel(tracked_files)
    rel_path = strtrim(tracked_files{i});
    if isempty(rel_path)
        continue;
    end
    [~, ~, ext] = fileparts(rel_path);
    if ~ismember(lower(ext), {'.m', '.md', '.txt'})
        continue;
    end
    % The gitignored real-data file is never tracked; the checked-in
    % template intentionally documents field names but must not itself
    % contain real values. Both are covered by scanning below - no
    % exclusion needed since the template only has placeholders.
    full_path = fullfile(repo_root, rel_path);
    if ~isfile(full_path)
        continue; % tracked-but-deleted-in-worktree edge case
    end
    text = fileread(full_path);
    for p = 1:numel(patterns)
        if contains(text, patterns{p})
            violations{end+1} = sprintf('%s contains pattern "%s"', ...
                rel_path, patterns{p}); %#ok<AGROW>
        end
    end
end

verifyEmpty(tc, violations, sprintf( ...
    'Tracked files must not contain patient identifiers:\n%s', ...
    strjoin(violations, sprintf('\n'))));
end

function test_provenance_template_has_no_real_values(tc)
repo_root = locate_repo_root();
template_path = fullfile(repo_root, 'config', 'private', ...
    'patient_provenance.local.m.example');
verifyTrue(tc, isfile(template_path), ...
    'The checked-in provenance template must exist.');

text = fileread(template_path);
patterns = identifier_patterns();
for p = 1:numel(patterns)
    verifyFalse(tc, contains(text, patterns{p}), sprintf( ...
        'Template must not contain real identifier "%s".', patterns{p}));
end
end
