%% test_reyna_scaling_experiment_contract.m
% Verify the locked Luna experiment contract and dry-run output.

clear; clc;
root = fileparts(mfilename('fullpath'));
project_root = fullfile(root, '..');
addpath(genpath(project_root));

contract = reyna_p1_scaling_v1();
assert(strcmp(contract.scenario, 'pre_surgery'));
assert(contract.gsa_training_samples == 128);
assert(isequal(contract.primary_metrics, ...
    {'RAP_mean','PAP_mean','SAP_mean','QpQs','CO_Lmin'}));
assert(numel(contract.arms) == 4);
assert(sum([contract.arms.historical_seeds_disabled]) == 2);
assert(strcmp(contract.arms(1).scaling_mode, 'zhang'));
assert(strcmp(contract.arms(2).scaling_mode, 'lundquist_bsa'));

result = run_reyna_scaling_experiment( ...
    'config/experiments/reyna_p1_scaling_v1.m', 'DryRun', true);
assert(strcmp(result.decision.label, 'dry_run'));
assert(height(result.summary_table) == 0);
assert(exist(result.contract_file, 'file') == 2);

fprintf('  [PASS] Reyna experiment contract locks arms, targets, and GSA N=128.\n');
