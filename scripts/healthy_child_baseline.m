function T = healthy_child_baseline(out_csv)
% HEALTHY_CHILD_BASELINE
% -----------------------------------------------------------------------
% Simulates the nominal parameter set scaled to a healthy 3-year-old
% (14 kg, 98 cm, WHO median; no defect) with the weight-based scaling of
% the paper, and writes the outputs used in Supplementary Section S6.
% -----------------------------------------------------------------------
if nargin < 1
    out_csv = 'results/runs/20260920_005419_reyna_pre_surgery/supplement/healthy_child_baseline.csv';
end
setenv('UNIFIED_VSD_SCALING_MODE', 'zhang');
patient = struct('weight_kg', 14, 'height_cm', 98, 'age_years', 3, ...
    'BSA', sqrt(14 * 98 / 3600), 'sex', 1, 'maturation_mode', 'normal', ...
    'scaling_mode', 'zhang', 'run_mode', 'exploratory', 'scaling_role', 'exploratory');
p = apply_scaling(default_parameters(), patient);
p.R.vsd = 1e6; p.vsd.area_mm2 = 0;
sim = integrate_system(p);
m = compute_clinical_indices(sim, p);
names = {'HR','CO_Lmin','LVSV','LVEF','RVEF','LVEDV','LVESV','RVEDV','RVESV','SVR','PVR', ...
    'SAP_max','SAP_min','SAP_mean','PAP_max','PAP_min','PAP_mean','LAP_mean','RAP_mean','QpQs'};
val = nan(numel(names), 1);
for i = 1:numel(names)
    if strcmp(names{i}, 'HR'), val(i) = p.HR;
    elseif isfield(m, names{i}), val(i) = m.(names{i});
    end
end
T = table(names(:), val, 'VariableNames', {'Output', 'Simulated'});
T.SteadyState = repmat(logical(sim.ss_reached), height(T), 1);
writetable(T, out_csv);
disp(T);
end
