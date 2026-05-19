function audit_baseline_provenance()
% AUDIT_BASELINE_PROVENANCE
% -----------------------------------------------------------------------
% Exports a source-of-truth audit table for the current default baseline.
% The audit is meant to document which parameters come from Bozkurt-style
% adult physiology, which arise from the Valenti architecture, and which
% reflect later consistency corrections in this repo.
%
% OUTPUTS:
%   results/tables/baseline_provenance_audit.csv
%   results/tables/baseline_provenance_summary.txt
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-08
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
results_dir = fullfile(root, 'results', 'tables');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

params_ref = default_parameters();
provenance = build_baseline_provenance(params_ref);

csv_path = fullfile(results_dir, 'baseline_provenance_audit.csv');
txt_path = fullfile(results_dir, 'baseline_provenance_summary.txt');
writetable(provenance, csv_path);

group_summary = groupsummary(provenance, 'FoundationSource');
fid = fopen(txt_path, 'w');
if fid < 0
    error('audit_baseline_provenance:openFailed', ...
        'Unable to write provenance summary file.');
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Unified VSD Baseline Provenance Audit\n');
fprintf(fid, '===================================\n');
fprintf(fid, 'CSV: %s\n\n', csv_path);
fprintf(fid, 'Total rows: %d\n\n', height(provenance));
fprintf(fid, 'Counts by foundation source:\n');
for idx = 1:height(group_summary)
    fprintf(fid, '  %s : %d\n', ...
        group_summary.FoundationSource{idx}, group_summary.GroupCount(idx));
end

fprintf('\n[audit_baseline_provenance] Audit exported:\n  %s\n  %s\n', ...
    csv_path, txt_path);
disp(provenance(:, {'Parameter','CurrentValue','FoundationSource','EvidenceScope','ConfidenceLevel'}));
end
