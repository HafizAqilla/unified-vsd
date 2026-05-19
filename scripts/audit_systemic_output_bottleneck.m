function audit_systemic_output_bottleneck(run_package_path)
% AUDIT_SYSTEMIC_OUTPUT_BOTTLENECK
% -----------------------------------------------------------------------
% Diagnoses why effective systemic output (CO_Lmin = Qs) remains low in a
% saved Reyna run. The audit compares systemic flow definitions, target
% consistency, pressure gradients, and shunt/load metrics for the accepted
% and best candidates saved by main_run.
%
% INPUTS:
%   run_package_path - optional path to a run_package_*.mat file          [-]
%
% OUTPUTS:
%   <run>/tables/systemic_output_bottleneck_audit.csv
%   <run>/tables/systemic_output_bottleneck_audit.txt
%   <run>/tables/co_definition_audit_<scenario>.csv
%
% AUTHOR:   Unified VSD Model
% DATE:     2026-05-08
% VERSION:  1.0
% -----------------------------------------------------------------------

root = fileparts(fileparts(mfilename('fullpath')));
restoredefaultpath();
addpath(build_clean_project_path(root));

if nargin < 1 || isempty(run_package_path)
    run_package_path = find_latest_reyna_run_package(root);
end

loaded = load(run_package_path, 'run_package');
if ~isfield(loaded, 'run_package')
    error('audit_systemic_output_bottleneck:missingRunPackage', ...
        'Selected MAT file does not contain run_package.');
end
run_package = loaded.run_package;
clinical = run_package.clinical;
scenario = run_package.scenario;
src = clinical.(scenario);

candidate_names = {'baseline','accepted_candidate','best_candidate'};
metrics_list = cell(size(candidate_names));
metrics_list{1} = run_package.metrics_baseline;
metrics_list{2} = run_package.accepted_candidate.metrics;
metrics_list{3} = run_package.best_candidate.metrics;

target_co = src.CO_Lmin;
target_sap_mean = src.SAP_mean_mmHg;
target_rap_mean = src.RAP_mean_mmHg;
target_svr = src.SVR_WU;
derived_svr = (target_sap_mean - target_rap_mean) / max(target_co, 1e-9);
target_qp = NaN;
if isfield(src, 'QpQs') && isfinite(src.QpQs)
    target_qp = target_co * src.QpQs;
end

rows = cell(numel(candidate_names), 21);
for idx = 1:numel(candidate_names)
    metrics = metrics_list{idx};
    q_sys = field_or_nan(metrics, 'CO_Lmin');
    q_pul = field_or_nan(metrics, 'Qp_Lmin');
    q_ao = field_or_nan(metrics, 'Qao_Lmin');
    lv_co = field_or_nan(metrics, 'LVCO_Lmin');
    q_vsd = field_or_nan(metrics, 'Qvsd_Lmin');
    sap_mean = field_or_nan(metrics, 'SAP_mean');
    rap_mean = field_or_nan(metrics, 'RAP_mean');
    svr_model = field_or_nan(metrics, 'SVR');
    expected_q_from_model_svr = (sap_mean - rap_mean) / max(svr_model, 1e-9);
    expected_q_from_target_svr = (sap_mean - rap_mean) / max(target_svr, 1e-9);
    svr_from_qao = (sap_mean - rap_mean) / max(q_ao, 1e-9);
    svr_from_lvco = (sap_mean - rap_mean) / max(lv_co, 1e-9);
    lvco_qp_residual = lv_co - q_pul;
    lvco_target_qp_residual = lv_co - target_qp;

    rows(idx, :) = { ...
        candidate_names{idx}, ...
        q_sys, ...
        q_pul, ...
        q_ao, ...
        lv_co, ...
        q_vsd, ...
        field_or_nan(metrics, 'QpQs'), ...
        field_or_nan(metrics, 'Qs_ao_gap_Lmin'), ...
        sap_mean, ...
        rap_mean, ...
        sap_mean - rap_mean, ...
        svr_model, ...
        svr_from_qao, ...
        svr_from_lvco, ...
        expected_q_from_model_svr, ...
        expected_q_from_target_svr, ...
        target_qp, ...
        lvco_qp_residual, ...
        lvco_target_qp_residual, ...
        100 * (q_sys - target_co) / max(abs(target_co), 1e-9)};
end

audit_tbl = cell2table(rows, 'VariableNames', { ...
    'Candidate','CO_Qs_Lmin','Qp_Lmin','Qao_Lmin','LVCO_Lmin','Qvsd_Lmin', ...
    'QpQs','QaoMinusQs_Lmin','SAP_mean','RAP_mean','MAPMinusRAP', ...
    'SVR_from_Qs','SVR_from_Qao','SVR_from_LVCO','Q_from_model_SVR', ...
    'Q_if_target_SVR','Clinical_Qp_Target_Lmin','LVCO_minus_Qp_Lmin', ...
    'LVCO_minus_clinical_Qp_Lmin','CO_ErrorPct'});

tables_dir = run_package.paths.tables_dir;
csv_path = fullfile(tables_dir, 'systemic_output_bottleneck_audit.csv');
txt_path = fullfile(tables_dir, 'systemic_output_bottleneck_audit.txt');
writetable(audit_tbl, csv_path);
co_definition_csv_path = fullfile(tables_dir, sprintf('co_definition_audit_%s.csv', scenario));
writetable(audit_tbl, co_definition_csv_path);

fid = fopen(txt_path, 'w');
if fid < 0
    error('audit_systemic_output_bottleneck:openFailed', ...
        'Unable to write systemic-output audit text file.');
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Systemic Output Bottleneck Audit\n');
fprintf(fid, '================================\n');
fprintf(fid, 'Run package: %s\n\n', run_package_path);
fprintf(fid, 'Clinical CO_Lmin target: %.4f L/min\n', target_co);
fprintf(fid, 'Clinical Qp target from Qs*QpQs: %.4f L/min\n', target_qp);
fprintf(fid, 'Clinical SAP_mean/RAP_mean gradient: %.4f mmHg\n', target_sap_mean - target_rap_mean);
fprintf(fid, 'Documented SVR target: %.4f WU\n', target_svr);
fprintf(fid, 'Derived SVR from MAP-RAP-CO: %.4f WU\n\n', derived_svr);
for idx = 1:height(audit_tbl)
    fprintf(fid, '%s\n', audit_tbl.Candidate{idx});
    fprintf(fid, '  CO_Qs_Lmin: %.4f L/min (error %.2f%%)\n', ...
        audit_tbl.CO_Qs_Lmin(idx), audit_tbl.CO_ErrorPct(idx));
    fprintf(fid, '  Qao_Lmin / LVCO_Lmin / Qvsd_Lmin: %.4f / %.4f / %.4f\n', ...
        audit_tbl.Qao_Lmin(idx), audit_tbl.LVCO_Lmin(idx), audit_tbl.Qvsd_Lmin(idx));
    fprintf(fid, '  Qp_Lmin / LVCO-Qp: %.4f / %.4f L/min\n', ...
        audit_tbl.Qp_Lmin(idx), audit_tbl.LVCO_minus_Qp_Lmin(idx));
    fprintf(fid, '  Qao-Qs gap: %.4f L/min\n', audit_tbl.QaoMinusQs_Lmin(idx));
    fprintf(fid, '  MAP-RAP: %.4f mmHg, SVR(Qs/Qao/LVCO): %.4f / %.4f / %.4f WU\n', ...
        audit_tbl.MAPMinusRAP(idx), audit_tbl.SVR_from_Qs(idx), ...
        audit_tbl.SVR_from_Qao(idx), audit_tbl.SVR_from_LVCO(idx));
    fprintf(fid, '  Q implied by target SVR at model pressure: %.4f L/min\n\n', ...
        audit_tbl.Q_if_target_SVR(idx));
end

fprintf('[audit_systemic_output_bottleneck] Audit exported:\n  %s\n  %s\n', ...
    csv_path, txt_path);
fprintf('  %s\n', co_definition_csv_path);
disp(audit_tbl);
end

function run_package_path = find_latest_reyna_run_package(root)
listing = dir(fullfile(root, 'results', 'runs', '*reyna_pre_surgery', 'mat', 'run_package_*.mat'));
if isempty(listing)
    error('audit_systemic_output_bottleneck:noRunPackage', ...
        'No Reyna pre-surgery run package found.');
end
[~, order] = sort([listing.datenum], 'descend');
run_package_path = fullfile(listing(order(1)).folder, listing(order(1)).name);
end

function value = field_or_nan(s, field_name)
if isstruct(s) && isfield(s, field_name) && isfinite(s.(field_name))
    value = s.(field_name);
else
    value = NaN;
end
end

function project_path = build_clean_project_path(root)
project_paths = strsplit(genpath(root), pathsep);
project_paths = project_paths(~cellfun('isempty', project_paths));
is_shadow = contains(project_paths, [filesep '.claude' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.clone' filesep], 'IgnoreCase', true) | ...
            contains(project_paths, [filesep '.git' filesep], 'IgnoreCase', true);
is_existing = cellfun(@isfolder, project_paths);
project_path = strjoin(project_paths(~is_shadow & is_existing), pathsep);
end
